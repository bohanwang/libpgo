#include "cudaNewtonSolver.h"
#include "cudaBlas.h"

#include <libuipc/linear_system/pcg_kernels.cuh>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace pgo::NonlinearOptimization
{
namespace
{
using CudaUtilities::DeviceBuffer;
using CudaUtilities::checkCuda;
constexpr int blockSize = 256;

__global__ void negateKernel(int n, const double *x, double *out)
{
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < n)
    out[i] = -x[i];
}

__global__ void jacobiKernel(int n, const double *diagonal, double shift,
  double floor, const double *r, double *z)
{
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < n)
    z[i] = r[i] / fmax(diagonal[i] + shift, floor);
}

__global__ void addShiftKernel(int n, const double *v, double shift, double *out)
{
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < n)
    out[i] += shift * v[i];
}

__global__ void trialKernel(int n, const double *x, const double *dx, double alpha, double *trial)
{
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < n)
    trial[i] = x[i] + alpha * dx[i];
}

__global__ void differenceKernel(int n, const double *base, const double *trial, double *difference)
{
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < n)
    difference[i] = trial[i] - base[i];
}

__global__ void trueResidualKernel(int n, const double *gradient, const double *Adx, double *r)
{
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < n)
    r[i] = -gradient[i] - Adx[i];
}

__global__ void residualScoreKernel(int n, const double *Hg, const double *diagonal, double floor, double *scores)
{
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < n)
    scores[i] = fabs(Hg[i]) / fmax(fabs(diagonal[i]), floor);
}

__global__ void selectCoordinateKernel(int index, double *basis, double *scores)
{
  basis[index] = 1;
  scores[index] = 0;
}

__global__ void coordinateTrialKernel(int index, const double *x, double step, double *trial)
{
  trial[index] = x[index] + step;
}

enum class LinearStatus
{
  Converged,
  IterationLimit,
  NonPositiveCurvature,
  NonFinite,
  Running,
};

struct DevicePCGState
{
  LinearStatus status = LinearStatus::Running;
  int iterations = 0;
  double rz = 0, alpha = 0, beta = 0, residualNorm = 0;
};

__global__ void pcgAlphaKernel(DevicePCGState *state, const double *curvature)
{
  if (state->status != LinearStatus::Running)
    return;
  if (!isfinite(*curvature))
    state->status = LinearStatus::NonFinite;
  else if (*curvature <= 0)
    state->status = LinearStatus::NonPositiveCurvature;
  else {
    state->alpha = state->rz / *curvature;
    if (!isfinite(state->alpha))
      state->status = LinearStatus::NonFinite;
  }
}

__global__ void pcgDeviceUpdateXR(int n, const DevicePCGState *state,
  double *x, const double *p, double *r, const double *Ap)
{
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < n && state->status == LinearStatus::Running) {
    x[i] += state->alpha * p[i];
    r[i] -= state->alpha * Ap[i];
  }
}

__global__ void pcgResidualKernel(DevicePCGState *state, const double *norm, double tolerance)
{
  if (state->status != LinearStatus::Running)
    return;
  ++state->iterations;
  state->residualNorm = *norm;
  if (!isfinite(*norm))
    state->status = LinearStatus::NonFinite;
  else if (*norm <= tolerance)
    state->status = LinearStatus::Converged;
}

__global__ void pcgDeviceJacobi(int n, const DevicePCGState *state,
  const double *diagonal, double shift, double floor, const double *r, double *z)
{
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < n && state->status == LinearStatus::Running)
    z[i] = r[i] / fmax(diagonal[i] + shift, floor);
}

__global__ void pcgBetaKernel(DevicePCGState *state, const double *rzNew)
{
  if (state->status != LinearStatus::Running)
    return;
  if (!isfinite(*rzNew) || *rzNew <= 0)
    state->status = LinearStatus::NonFinite;
  else {
    state->beta = *rzNew / state->rz;
    if (!isfinite(state->beta))
      state->status = LinearStatus::NonFinite;
    else
      state->rz = *rzNew;
  }
}

__global__ void pcgDeviceUpdateP(int n, const DevicePCGState *state, double *p, const double *z)
{
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < n && state->status == LinearStatus::Running)
    p[i] = z[i] + state->beta * p[i];
}
}  // namespace

class CudaNewtonSolver::Impl
{
public:
  Impl(std::vector<std::shared_ptr<CudaPotentialEnergy>> models, const SolverParam &params):
    energies(std::move(models)), parameters(params)
  {
    if (energies.empty() || !energies.front())
      throw std::invalid_argument("CUDA Newton requires at least one energy.");
    n = energies.front()->getNumDOFs();
    if (n <= 0)
      throw std::invalid_argument("CUDA Newton requires a positive number of DOFs.");
    for (const auto &model : energies)
      if (!model || model->getNumDOFs() != n)
        throw std::invalid_argument("CUDA Newton energies must have matching DOFs.");
    if (parameters.maxPCGIterations < 0 || parameters.pcgReductionBatchSize < 0 || parameters.pcgReductionBatchSize > 64 ||
      parameters.stagnationWindow < 0 || parameters.stagnationWindow > 64 ||
      !std::isfinite(parameters.stagnationRelativeTolerance) || parameters.stagnationRelativeTolerance <= 0 || parameters.stagnationRelativeTolerance >= 1 ||
      parameters.maxLineSearchIterations < 1 ||
      parameters.maxRegularizationAttempts < 0 || parameters.maxResidualCorrectionCandidates < 0 ||
      !std::isfinite(parameters.pcgRelativeTolerance) || parameters.pcgRelativeTolerance <= 0 || parameters.pcgRelativeTolerance >= 1 ||
      !std::isfinite(parameters.pcgAbsoluteTolerance) || parameters.pcgAbsoluteTolerance < 0 ||
      !std::isfinite(parameters.armijoCoeff) || parameters.armijoCoeff <= 0 || parameters.armijoCoeff >= 1 ||
      !std::isfinite(parameters.initialRegularization) || parameters.initialRegularization <= 0)
      throw std::invalid_argument("CUDA Newton parameters are outside their valid ranges.");
    x.resize(n);
    gradient.resize(n);
    trial.resize(n);
    trialGradient.resize(n);
    dx.resize(n);
    diagonal.resize(n);
    r.resize(n);
    z.resize(n);
    p.resize(n);
    Ap.resize(n);
    if (parameters.checkInputChord)
      initialX.resize(n);
    if (parameters.checkInputChord || parameters.checkTrialChord)
      inputChord.resize(n);
    if (parameters.pcgReductionBatchSize) {
      deviceBlas = std::make_unique<detail::DeviceBlasHandle>();
      pcgState.resize(1);
      pcgReductionScalars.resize(3);
    }
  }

  int blocks() const { return (n + blockSize - 1) / blockSize; }
  void copy(double *out, const double *in)
  {
    checkCuda(cudaMemcpyAsync(out, in, sizeof(double) * n, cudaMemcpyDeviceToDevice), "CUDA vector copy");
  }

  double energy(const double *at, std::vector<double> *values = nullptr)
  {
    double sum = 0;
    if (values)
      values->clear();
    for (const auto &model : energies) {
      const double value = model->func(at);
      if (values)
        values->push_back(value);
      sum += value;
    }
    return sum;
  }

  CudaPotentialEnergy::EnergyChange energyChange(const double *base, const double *next)
  {
    CudaPotentialEnergy::EnergyChange sum{ 0, 0 };
    for (const auto &model : energies) {
      const auto change = model->energyChange(base, next);
      if (!std::isfinite(change.value) || !std::isfinite(change.roundoffScale) || change.roundoffScale < 0)
        return { std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN() };
      sum.value += change.value;
      sum.roundoffScale += std::max(std::abs(change.value), change.roundoffScale);
    }
    return sum;
  }

  void evaluateGradient(const double *at, DeviceBuffer<double> &out)
  {
    out.zero();
    for (const auto &model : energies)
      model->addGradient(at, out.data());
    checkCuda(cudaGetLastError(), "CUDA total gradient");
  }

  void applyHessian(const double *v, double shift, DeviceBuffer<double> &out)
  {
    out.zero();
    for (const auto &model : energies)
      model->addHessianVector(v, out.data());
    if (shift != 0)
      addShiftKernel<<<blocks(), blockSize>>>(n, v, shift, out.data());
    checkCuda(cudaGetLastError(), "CUDA Hessian-vector product");
  }

  bool prepareHessian(bool exact, double &diagonalScale)
  {
    diagonal.zero();
    for (const auto &model : energies) {
      if (exact)
        model->prepareExactHessian(x.data());
      else
        model->prepareHessian(x.data());
      model->addDiagonal(diagonal.data());
    }
    if (!std::isfinite(blas.norm(n, diagonal.data())))
      return false;
    diagonalScale = std::max(1e-12, blas.maxAbs(n, diagonal.data()));
    return true;
  }

  bool certifyChord(const double *from, const double *to, double &lastBound)
  {
    differenceKernel<<<blocks(), blockSize>>>(n, from, to, inputChord.data());
    checkCuda(cudaGetLastError(), "CUDA input-chord displacement");
    lastBound = 1;
    for (const auto &model : energies) {
      const double bound = model->computeMaxStepSize(from, inputChord.data());
      if (!std::isfinite(bound) || bound < 0 || bound > 1) {
        lastBound = bound;
        return false;
      }
      lastBound = std::min(lastBound, bound);
    }
    // Only the energy's certificate for the complete original chord is valid.
    // Rebuilding suffixes or rounded midpoints can change the path and miss a
    // narrow collision. CCD handles conservative advancement internally while
    // retaining its original endpoints and arithmetic error bounds.
    return lastBound == 1;
  }

  bool recoverRoundedDirection(double tolerance, double gradientNorm2, double gradientInfinity,
    double &trialEnergy, double &trialEnergyChange, double &trialGradientNorm2, double &trialGradientInfinity,
    std::vector<double> *trialEnergyValues, int verbose)
  {
    if (parameters.maxResidualCorrectionCandidates == 0)
      return false;
    ++diagnostics.residualCorrectionAttempts;
    double diagonalScale = 0;
    if (!prepareHessian(true, diagonalScale))
      return false;
    // For Phi(x)=.5||g(x)||^2, its gradient is H*g. Rank coordinates by
    // the diagonally scaled residual gradient, without transferring vectors.
    applyHessian(gradient.data(), 0, r);
    if (!std::isfinite(blas.norm(n, r.data())))
      return false;
    residualScoreKernel<<<blocks(), blockSize>>>(n, r.data(), diagonal.data(),
      std::max(1e-30, diagonalScale * 1e-12), p.data());
    checkCuda(cudaGetLastError(), "CUDA residual coordinate scores");
    for (int candidate = 0; candidate < std::min(n, parameters.maxResidualCorrectionCandidates); ++candidate) {
      const int coordinate = blas.maxAbsIndex(n, p.data());
      if (coordinate < 0)
        break;
      double score = 0;
      checkCuda(cudaMemcpy(&score, p.data() + coordinate, sizeof(double), cudaMemcpyDeviceToHost), "CUDA residual coordinate score");
      if (!(score > 0) || !std::isfinite(score))
        break;
      z.zero();
      selectCoordinateKernel<<<1, 1>>>(coordinate, z.data(), p.data());
      applyHessian(z.data(), 0, Ap);
      const double norm = blas.norm(n, Ap.data());
      const double numerator = blas.dot(n, gradient.data(), Ap.data());
      if (!(norm > 0) || !std::isfinite(norm) || !std::isfinite(numerator))
        continue;
      // This minimizes ||g+s*H[:,j]||^2 for the linearized residual. The
      // measured residual below decides whether the represented step works.
      const double step = -(numerator / norm) / norm;
      if (!std::isfinite(step))
        continue;
      copy(trial.data(), x.data());
      coordinateTrialKernel<<<1, 1>>>(coordinate, x.data(), step, trial.data());
      differenceKernel<<<blocks(), blockSize>>>(n, x.data(), trial.data(), dx.data());
      checkCuda(cudaGetLastError(), "CUDA residual coordinate trial");
      if (!(blas.maxAbs(n, dx.data()) > 0))
        continue;
      ++diagnostics.residualCorrectionCandidates;
      bool feasible = true;
      for (const auto &model : energies) {
        if (model->computeMaxStepSize(x.data(), dx.data()) != 1.0) {
          feasible = false;
          break;
        }
      }
      if (!feasible) {
        if (parameters.checkTrialChord)
          ++diagnostics.trialChordRejections;
        continue;
      }
      if (parameters.checkInputChord) {
        double bound = 1;
        if (!certifyChord(initialX.data(), trial.data(), bound)) {
          ++diagnostics.inputChordRejections;
          continue;
        }
      }
      trialEnergy = energy(trial.data(), trialEnergyValues);
      if (!std::isfinite(trialEnergy))
        continue;
      const auto change = energyChange(x.data(), trial.data());
      const double roundoffBound = 16 * std::numeric_limits<double>::epsilon() * change.roundoffScale;
      const double work = blas.dot(n, gradient.data(), dx.data());
      if (!std::isfinite(roundoffBound) || !std::isfinite(change.value) || !std::isfinite(work) ||
        std::abs(change.value) > roundoffBound || std::abs(work) > roundoffBound)
        continue;
      evaluateGradient(trial.data(), trialGradient);
      trialGradientNorm2 = blas.norm(n, trialGradient.data());
      trialGradientInfinity = blas.maxAbs(n, trialGradient.data());
      if (!std::isfinite(trialGradientNorm2) || !(trialGradientNorm2 < gradientNorm2) ||
        !(trialGradientInfinity < tolerance || trialGradientInfinity <= 0.9 * gradientInfinity))
        continue;
      trialEnergyChange = change.value;
      ++diagnostics.residualCorrectionSteps;
      if (verbose >= 2)
        std::cout << "        CUDA residual correction coordinate=" << coordinate
                  << "; step=" << step << "; ||g||_inf=" << trialGradientInfinity
                  << "; energyChange=" << change.value << std::endl;
      return true;
    }
    return false;
  }

  LinearStatus solvePCGBatched(double shift, double diagonalFloor, double tolerance,
    double rz, double rhsNorm, int maxIterations)
  {
    DevicePCGState state;
    state.rz = rz;
    state.residualNorm = rhsNorm;
    pcgState.upload(&state, 1);
    auto *deviceState = pcgState.data();
    auto *scalars = pcgReductionScalars.data();
    for (int queued = 0; queued < maxIterations;) {
      const int count = std::min(parameters.pcgReductionBatchSize, maxIterations - queued);
      for (int iteration = 0; iteration < count; ++iteration) {
        // These operators only write scratch. They may still execute after a
        // terminal flag, while all solution/residual/direction state is frozen.
        applyHessian(p.data(), shift, Ap);
        deviceBlas->dot(n, p.data(), Ap.data(), scalars);
        pcgAlphaKernel<<<1, 1>>>(deviceState, scalars);
        pcgDeviceUpdateXR<<<blocks(), blockSize>>>(n, deviceState, dx.data(), p.data(), r.data(), Ap.data());
        deviceBlas->norm(n, r.data(), scalars + 1);
        pcgResidualKernel<<<1, 1>>>(deviceState, scalars + 1, tolerance);
        pcgDeviceJacobi<<<blocks(), blockSize>>>(n, deviceState, diagonal.data(), shift, diagonalFloor, r.data(), z.data());
        deviceBlas->dot(n, r.data(), z.data(), scalars + 2);
        pcgBetaKernel<<<1, 1>>>(deviceState, scalars + 2);
        pcgDeviceUpdateP<<<blocks(), blockSize>>>(n, deviceState, p.data(), z.data());
        checkCuda(cudaGetLastError(), "CUDA batched PCG iteration");
      }
      queued += count;
      pcgState.download(&state, 1);
      if (state.status != LinearStatus::Running)
        break;
    }
    diagnostics.pcgIterations += state.iterations;
    diagnostics.lastLinearResidualNorm = state.residualNorm;
    return state.status == LinearStatus::Running ? LinearStatus::IterationLimit : state.status;
  }

  // Adapted from libuipc LinearPCG::pcg; see the imported kernels' provenance.
  // Matrix and preconditioner are fixed throughout this solve. Only scalars
  // leave the device, including the explicit curvature and residual checks.
  LinearStatus solvePCG(double shift, double diagonalScale, double newtonTolerance)
  {
    dx.zero();
    negateKernel<<<blocks(), blockSize>>>(n, gradient.data(), r.data());
    const double rhsNorm = blas.norm(n, r.data());
    if (!std::isfinite(rhsNorm))
      return LinearStatus::NonFinite;
    // A stricter requested Newton tolerance must not be blocked by the default
    // absolute PCG tolerance returning a zero direction for a nonzero gradient.
    const double absoluteTolerance = std::min(parameters.pcgAbsoluteTolerance, 0.1 * newtonTolerance);
    const double tolerance = std::max(absoluteTolerance, parameters.pcgRelativeTolerance * rhsNorm);
    diagnostics.lastLinearResidualNorm = rhsNorm;
    if (rhsNorm <= tolerance)
      return LinearStatus::Converged;

    const double diagonalFloor = std::max(1e-30, diagonalScale * 1e-12);
    auto precondition = [&]() {
      jacobiKernel<<<blocks(), blockSize>>>(n, diagonal.data(), shift, diagonalFloor, r.data(), z.data());
      checkCuda(cudaGetLastError(), "CUDA Jacobi preconditioner");
    };
    precondition();
    copy(p.data(), z.data());
    double rz = blas.dot(n, r.data(), z.data());
    if (!std::isfinite(rz) || rz <= 0)
      return LinearStatus::NonFinite;

    const int maxIterations = parameters.maxPCGIterations > 0 ? parameters.maxPCGIterations :
                                                                static_cast<int>(std::min(10000LL, 2LL * n));
    LinearStatus result = LinearStatus::IterationLimit;
    if (parameters.pcgReductionBatchSize) {
      result = solvePCGBatched(shift, diagonalFloor, tolerance, rz, rhsNorm, maxIterations);
      if (result == LinearStatus::NonFinite || result == LinearStatus::NonPositiveCurvature)
        return result;
    }
    else
      for (int iteration = 0; iteration < maxIterations; ++iteration) {
        applyHessian(p.data(), shift, Ap);
        const double curvature = blas.dot(n, p.data(), Ap.data());
        if (!std::isfinite(curvature))
          return LinearStatus::NonFinite;
        if (curvature <= 0)
          return LinearStatus::NonPositiveCurvature;
        const double alpha = rz / curvature;
        if (!std::isfinite(alpha))
          return LinearStatus::NonFinite;
        detail::pcgUpdateXR<<<blocks(), blockSize>>>(alpha, dx.data(), p.data(), r.data(), Ap.data(), n);
        checkCuda(cudaGetLastError(), "CUDA PCG x/residual update");
        ++diagnostics.pcgIterations;
        const double residualNorm = blas.norm(n, r.data());
        diagnostics.lastLinearResidualNorm = residualNorm;
        if (!std::isfinite(residualNorm))
          return LinearStatus::NonFinite;
        if (residualNorm <= tolerance) {
          result = LinearStatus::Converged;
          break;
        }
        precondition();
        const double rzNew = blas.dot(n, r.data(), z.data());
        if (!std::isfinite(rzNew) || rzNew <= 0)
          return LinearStatus::NonFinite;
        const double beta = rzNew / rz;
        if (!std::isfinite(beta))
          return LinearStatus::NonFinite;
        detail::pcgUpdateP<<<blocks(), blockSize>>>(p.data(), z.data(), beta, n);
        checkCuda(cudaGetLastError(), "CUDA PCG direction update");
        rz = rzNew;
      }

    // Recurrence residuals can drift. Report the actual residual and never
    // label a linear solve converged solely from the recurrence estimate.
    applyHessian(dx.data(), shift, Ap);
    trueResidualKernel<<<blocks(), blockSize>>>(n, gradient.data(), Ap.data(), r.data());
    diagnostics.lastLinearResidualNorm = blas.norm(n, r.data());
    if (!std::isfinite(diagnostics.lastLinearResidualNorm) || !std::isfinite(blas.norm(n, dx.data())))
      return LinearStatus::NonFinite;
    if (diagnostics.lastLinearResidualNorm > tolerance)
      result = LinearStatus::IterationLimit;
    return result;
  }

  int solve(double *hostX, int maxIterations, double tolerance, int verbose)
  {
    if (!hostX || maxIterations < 0 || !std::isfinite(tolerance) || tolerance <= 0)
      throw std::invalid_argument("CUDA Newton requires valid state, nonnegative iterations and positive finite tolerance.");
    diagnostics = {};
    diagnostics.stopReason = "iteration_limit";
    diagnostics.finalEnergy = std::numeric_limits<double>::quiet_NaN();
    diagnostics.finalGradientInfinityNorm = std::numeric_limits<double>::quiet_NaN();
    double initialCoordinateScale = 0;
    const bool watchStagnation = parameters.stagnationWindow > 0 && parameters.checkInputChord;
    for (int i = 0; i < n; ++i) {
      if (!std::isfinite(hostX[i])) {
        diagnostics.stopReason = "non_finite_initial_state";
        return SOLVE_NUMERICAL_FAILURE;
      }
      if (watchStagnation)
        initialCoordinateScale = std::max(initialCoordinateScale, std::abs(hostX[i]));
    }
    x.upload(hostX, n);
    // A stationary but infeasible input is not a converged solution. Query
    // enabled model guards before evaluating energy or observing the input;
    // a disabled guard certifies the zero-length chord with a bound of one.
    dx.zero();
    for (const auto &model : energies) {
      const double bound = model->computeMaxStepSize(x.data(), dx.data());
      if (bound != 1.0) {
        diagnostics.lastFeasibleAlpha = bound;
        diagnostics.stopReason = "infeasible_initial_state";
        return SOLVE_NUMERICAL_FAILURE;
      }
    }
    if (parameters.checkInputChord)
      copy(initialX.data(), x.data());
    std::vector<double> currentEnergyValues, trialEnergyValues;
    double currentEnergy = energy(x.data(), verbose >= 3 ? &currentEnergyValues : nullptr);
    evaluateGradient(x.data(), gradient);
    double gradientNorm2 = blas.norm(n, gradient.data());
    if (!std::isfinite(currentEnergy) || !std::isfinite(gradientNorm2)) {
      diagnostics.stopReason = "non_finite_initial_energy_or_gradient";
      return SOLVE_NUMERICAL_FAILURE;
    }
    double gradientInfinity = blas.maxAbs(n, gradient.data());
    if (acceptedIterateObserver)
      acceptedIterateObserver(0, x.data(), currentEnergy, gradientInfinity);
    int status = SOLVE_NOT_CONVERGED;
    int stagnationCount = 0;
    double stagnationAnchorNorm2 = 0, stagnationAnchorInfinity = 0;

    for (int iteration = 0; iteration < maxIterations; ++iteration) {
      if (gradientInfinity < tolerance) {
        diagnostics.stopReason = "gradient_tolerance";
        break;
      }
      if (verbose >= 2)
        std::cout << "    CUDA Newton iter=" << iteration << "; E=" << currentEnergy
                  << "; ||g||_inf=" << gradientInfinity << std::endl;
      bool directionReady = false;
      double directionDerivative = 0;
      double diagonalScale = 0;
      // Local PSD projection can greatly overestimate the curvature of a
      // positive global system. After prolonged convergence, try the exact
      // tangent without changing the objective or line-search safeguards.
      // Any unusable exact solve rebuilds the original projected operator.
      if (parameters.exactHessianStartIteration >= 0 && iteration >= parameters.exactHessianStartIteration) {
        ++diagnostics.exactHessianAttempts;
        if (prepareHessian(true, diagonalScale)) {
          const LinearStatus exactStatus = solvePCG(0, diagonalScale, tolerance);
          if (exactStatus == LinearStatus::Converged || exactStatus == LinearStatus::IterationLimit) {
            directionDerivative = blas.dot(n, gradient.data(), dx.data());
            directionReady = std::isfinite(directionDerivative) && directionDerivative < 0;
          }
        }
        if (directionReady)
          ++diagnostics.exactHessianDirections;
        else
          ++diagnostics.exactHessianFallbacks;
        if (verbose >= 2)
          std::cout << "        CUDA exact Hessian=" << (directionReady ? "descent_direction" : "projected_fallback") << std::endl;
      }
      if (!directionReady && !prepareHessian(false, diagonalScale)) {
        diagnostics.stopReason = "non_finite_hessian_diagonal";
        status = SOLVE_NUMERICAL_FAILURE;
        break;
      }
      double shift = 0;
      for (int attempt = 0; !directionReady && attempt <= parameters.maxRegularizationAttempts; ++attempt) {
        if (attempt > 0) {
          shift = attempt == 1 ? parameters.initialRegularization * diagonalScale : shift * 10.0;
          ++diagnostics.regularizationAttempts;
          diagnostics.maxRegularization = std::max(diagnostics.maxRegularization, shift);
          if (verbose >= 2)
            std::cout << "        CUDA PCG regularization=" << shift << std::endl;
        }
        const LinearStatus linearStatus = solvePCG(shift, diagonalScale, tolerance);
        if (linearStatus == LinearStatus::NonFinite) {
          diagnostics.stopReason = "non_finite_linear_solve";
          status = SOLVE_NUMERICAL_FAILURE;
          break;
        }
        if (linearStatus == LinearStatus::NonPositiveCurvature)
          continue;
        directionDerivative = blas.dot(n, gradient.data(), dx.data());
        if (!std::isfinite(directionDerivative)) {
          diagnostics.stopReason = "non_finite_direction";
          status = SOLVE_NUMERICAL_FAILURE;
          break;
        }
        if (directionDerivative < 0) {
          directionReady = true;
          if (linearStatus == LinearStatus::IterationLimit && verbose >= 2)
            std::cout << "        CUDA PCG returned an inexact descent direction; residual="
                      << diagnostics.lastLinearResidualNorm << std::endl;
          break;
        }
      }
      if (!directionReady) {
        if (status != SOLVE_NUMERICAL_FAILURE)
          diagnostics.stopReason = "linear_curvature_or_descent_failure";
        break;
      }
      const double rawStepInfinity = verbose >= 2 || watchStagnation ? blas.maxAbs(n, dx.data()) : 0;
      if (verbose >= 2)
        std::cout << "        rawStepMaxNorm=" << rawStepInfinity
                  << "; gradientDotStep=" << directionDerivative << std::endl;
      if (verbose >= 3) {
        std::ostringstream line;
        line << std::setprecision(17) << "        componentGradientDotStep=";
        for (const auto &model : energies) {
          z.zero();
          model->addGradient(x.data(), z.data());
          line << '[' << blas.dot(n, z.data(), dx.data()) << ']';
        }
        std::cout << line.str() << std::endl;
      }

      double feasibleAlpha = 1;
      bool validBound = true;
      for (const auto &model : energies) {
        const double bound = model->computeMaxStepSize(x.data(), dx.data());
        if (!std::isfinite(bound) || bound < 0) {
          validBound = false;
          break;
        }
        feasibleAlpha = std::min(feasibleAlpha, bound);
      }
      diagnostics.lastFeasibleAlpha = feasibleAlpha;
      if (!validBound) {
        diagnostics.stopReason = "invalid_feasible_step";
        status = SOLVE_NUMERICAL_FAILURE;
        break;
      }
      if (feasibleAlpha == 0) {
        diagnostics.stopReason = "zero_feasible_step";
        break;
      }
      double alpha = feasibleAlpha;
      bool accepted = false;
      double trialEnergy = currentEnergy;
      double trialGradientNorm2 = gradientNorm2;
      double trialGradientInfinity = gradientInfinity;
      double trialEnergyChange = 0;
      bool roundoffAccepted = false;
      bool roundoffLimitedTrial = false;
      bool allTrialsChordBlocked = parameters.checkInputChord || parameters.checkTrialChord;
      std::string chordFailureReason;
      const int inputChordRejectionsBefore = diagnostics.inputChordRejections;
      for (int backtrack = 0; backtrack < parameters.maxLineSearchIterations; ++backtrack) {
        trialKernel<<<blocks(), blockSize>>>(n, x.data(), dx.data(), alpha, trial.data());
        checkCuda(cudaGetLastError(), "CUDA line search trial");
        if (parameters.checkTrialChord) {
          double chordBound = 1;
          if (!certifyChord(x.data(), trial.data(), chordBound)) {
            ++diagnostics.trialChordRejections;
            chordFailureReason = "trial_chord_blocked";
            if (verbose >= 2)
              std::cout << "        CUDA trial chord rejected; alpha=" << alpha
                        << "; bound=" << chordBound << std::endl;
            alpha *= 0.5;
            continue;
          }
        }
        if (parameters.checkInputChord) {
          double chordBound = 1;
          if (!certifyChord(initialX.data(), trial.data(), chordBound)) {
            ++diagnostics.inputChordRejections;
            chordFailureReason = "input_chord_blocked";
            if (verbose >= 2)
              std::cout << "        CUDA trial input chord rejected; alpha=" << alpha
                        << "; bound=" << chordBound << std::endl;
            alpha *= 0.5;
            continue;
          }
        }
        allTrialsChordBlocked = false;
        trialEnergy = energy(trial.data(), verbose >= 3 ? &trialEnergyValues : nullptr);
        const auto change = std::isfinite(trialEnergy) ? energyChange(x.data(), trial.data()) :
                                                         CudaPotentialEnergy::EnergyChange{ trialEnergy, std::numeric_limits<double>::quiet_NaN() };
        trialEnergyChange = change.value;
        const double roundoffBound = 16 * std::numeric_limits<double>::epsilon() * change.roundoffScale;
        // Near a stationary point, the expected energy decrease is quadratic
        // in the residual and can be smaller than one energy ULP. Comparing
        // rounded differences alone then rejects useful Newton corrections,
        // or accepts unchanged positions until the solver reports a tiny step.
        // Only within this roundoff interval, use strict measured residual
        // improvement as the merit criterion. The collision/material bound
        // and the requested gradient convergence tolerance remain unchanged.
        const bool roundoffLimited = std::isfinite(roundoffBound) &&
          -alpha * directionDerivative <= roundoffBound &&
          std::abs(trialEnergyChange) <= roundoffBound;
        roundoffLimitedTrial |= roundoffLimited;
        const bool armijo = trialEnergyChange <= parameters.armijoCoeff * alpha * directionDerivative;
        if (verbose >= 3 && std::isfinite(trialEnergy)) {
          evaluateGradient(trial.data(), trialGradient);
          trialGradientNorm2 = blas.norm(n, trialGradient.data());
          trialGradientInfinity = blas.maxAbs(n, trialGradient.data());
          std::ostringstream line;
          line << std::setprecision(17) << "        CUDA trial backtrack=" << backtrack
               << "; alpha=" << alpha << "; E0=" << currentEnergy << "; Etrial=" << trialEnergy
               << "; dE=" << trialEnergyChange << "; absoluteEnergyDifference=" << trialEnergy - currentEnergy
               << "; predictedDecrease=" << -alpha * directionDerivative
               << "; roundoffBound=" << roundoffBound << "; armijo=" << armijo
               << "; roundoffLimited=" << roundoffLimited << "; g0=" << gradientInfinity
               << "; gtrial=" << trialGradientInfinity << "; components=";
          for (std::size_t j = 0; j < currentEnergyValues.size(); ++j)
            line << '[' << currentEnergyValues[j] << ',' << trialEnergyValues[j] << ']';
          std::cout << line.str() << std::endl;
        }
        if (std::isfinite(trialEnergy) && std::isfinite(trialEnergyChange) && (armijo || roundoffLimited) &&
          std::isfinite(blas.norm(n, trial.data()))) {
          if (verbose < 3) {
            evaluateGradient(trial.data(), trialGradient);
            trialGradientNorm2 = blas.norm(n, trialGradient.data());
          }
          if (std::isfinite(trialGradientNorm2)) {
            trialGradientInfinity = blas.maxAbs(n, trialGradient.data());
            if (!roundoffLimited || trialGradientInfinity < tolerance || trialGradientInfinity <= 0.9 * gradientInfinity) {
              accepted = true;
              roundoffAccepted = roundoffLimited;
              break;
            }
          }
        }
        alpha *= 0.5;
      }
      if (!accepted && roundoffLimitedTrial && !allTrialsChordBlocked) {
        accepted = recoverRoundedDirection(tolerance, gradientNorm2, gradientInfinity,
          trialEnergy, trialEnergyChange, trialGradientNorm2, trialGradientInfinity,
          verbose >= 3 ? &trialEnergyValues : nullptr, verbose);
        if (accepted) {
          alpha = 1;
          feasibleAlpha = diagnostics.lastFeasibleAlpha = 1;
          roundoffAccepted = true;
        }
      }
      if (!accepted) {
        diagnostics.stopReason = allTrialsChordBlocked ? chordFailureReason : "line_search_failed";
        break;
      }
      // A small displacement is meaningful when coordinates and stiffness
      // have different scales. Stop only when the represented state did not
      // move; an absolute length cutoff can reject resolvable corrections.
      differenceKernel<<<blocks(), blockSize>>>(n, x.data(), trial.data(), z.data());
      checkCuda(cudaGetLastError(), "CUDA accepted-step difference");
      double acceptedStepInfinity = blas.maxAbs(n, z.data());
      bool tinyClippedCandidate = false;
      if (watchStagnation && diagnostics.inputChordRejections > inputChordRejectionsBefore &&
        acceptedStepInfinity > 0 && trialGradientInfinity >= tolerance) {
        const double scale = std::max({ initialCoordinateScale, blas.maxAbs(n, x.data()), rawStepInfinity });
        const double stepRatio = scale > 0 ? acceptedStepInfinity / scale : std::numeric_limits<double>::infinity();
        if (stagnationCount == 0) {
          stagnationAnchorNorm2 = gradientNorm2;
          stagnationAnchorInfinity = gradientInfinity;
        }
        diagnostics.lastStagnationStepRatio = stepRatio;
        diagnostics.stagnationAnchorNorm2 = stagnationAnchorNorm2;
        diagnostics.stagnationAnchorInfinity = stagnationAnchorInfinity;
        diagnostics.lastStagnationNorm2RelativeChange = std::abs(trialGradientNorm2 - stagnationAnchorNorm2) / stagnationAnchorNorm2;
        diagnostics.lastStagnationInfinityRelativeChange = std::abs(trialGradientInfinity - stagnationAnchorInfinity) / stagnationAnchorInfinity;
        // No world-unit length floor: both sides scale with a change of units.
        // Compare against the window's initial residual, so small cumulative
        // progress cannot disappear through comparisons with adjacent iterates.
        constexpr double epsilon = std::numeric_limits<double>::epsilon();
        tinyClippedCandidate = std::isfinite(scale) && std::isfinite(stepRatio) && stepRatio <= 16 * epsilon &&
          diagnostics.lastStagnationNorm2RelativeChange <= parameters.stagnationRelativeTolerance &&
          diagnostics.lastStagnationInfinityRelativeChange <= parameters.stagnationRelativeTolerance;
        if (tinyClippedCandidate)
          ++stagnationCount;
        if (verbose >= 3) {
          std::ostringstream line;
          line << std::setprecision(17) << "        CUDA stagnation candidate: stepRatio=" << stepRatio
               << "; norm2RelativeChange=" << diagnostics.lastStagnationNorm2RelativeChange
               << "; infinityRelativeChange=" << diagnostics.lastStagnationInfinityRelativeChange
               << "; anchorNorm2=" << stagnationAnchorNorm2 << "; anchorInfinity=" << stagnationAnchorInfinity
               << "; streak=" << (tinyClippedCandidate ? stagnationCount : 0)
               << "; relativeTolerance=" << parameters.stagnationRelativeTolerance;
          std::cout << line.str() << std::endl;
        }
      }
      if (!tinyClippedCandidate)
        stagnationCount = 0;
      if (watchStagnation && stagnationCount >= parameters.stagnationWindow) {
        ++diagnostics.stagnationTriggers;
        // A genuine local energy decrease may coexist with an unchanged full
        // residual. Do not relabel it as energy roundoff or as convergence.
        // Recovery still uses its original energy, residual and CCD safeguards.
        if (!recoverRoundedDirection(tolerance, gradientNorm2, gradientInfinity,
              trialEnergy, trialEnergyChange, trialGradientNorm2, trialGradientInfinity,
              verbose >= 3 ? &trialEnergyValues : nullptr, verbose)) {
          diagnostics.stopReason = "representational_stagnation";
          break;
        }
        stagnationCount = 0;
        alpha = 1;
        feasibleAlpha = diagnostics.lastFeasibleAlpha = 1;
        roundoffAccepted = true;
        differenceKernel<<<blocks(), blockSize>>>(n, x.data(), trial.data(), z.data());
        checkCuda(cudaGetLastError(), "CUDA stagnation recovery difference");
        acceptedStepInfinity = blas.maxAbs(n, z.data());
      }
      const double acceptedEnergyChange = trialEnergyChange;
      x.swap(trial);
      gradient.swap(trialGradient);
      currentEnergy = trialEnergy;
      currentEnergyValues.swap(trialEnergyValues);
      gradientNorm2 = trialGradientNorm2;
      gradientInfinity = trialGradientInfinity;
      ++diagnostics.newtonIterations;
      diagnostics.roundoffAcceptedSteps += roundoffAccepted;
      diagnostics.lastAcceptedAlpha = alpha;
      if (acceptedIterateObserver)
        acceptedIterateObserver(diagnostics.newtonIterations, x.data(), currentEnergy, gradientInfinity);
      if (verbose >= 2)
        std::cout << "        feasibleAlpha=" << feasibleAlpha << "; acceptedAlpha=" << alpha
                  << "; acceptedStepMaxNorm=" << acceptedStepInfinity
                  << "; energyChange=" << acceptedEnergyChange
                  << "; roundoffAcceptance=" << roundoffAccepted << std::endl;
      if (acceptedStepInfinity == 0) {
        diagnostics.stopReason = "tiny_step";
        break;
      }
    }

    diagnostics.finalEnergy = currentEnergy;
    diagnostics.finalGradientInfinityNorm = gradientInfinity;
    if (status != SOLVE_NUMERICAL_FAILURE && gradientInfinity < tolerance) {
      status = SOLVE_CONVERGED;
      diagnostics.stopReason = "gradient_tolerance";
    }
    x.download(hostX, n);
    if (verbose >= 1)
      std::cout << "    CUDA Newton status=" << status << "; stop=" << diagnostics.stopReason
                << "; ||g||_inf=" << gradientInfinity << "; Newton iterations=" << diagnostics.newtonIterations
                << "; PCG iterations=" << diagnostics.pcgIterations
                << "; regularizations=" << diagnostics.regularizationAttempts
                << "; exactHessianDirections=" << diagnostics.exactHessianDirections
                << "; exactHessianFallbacks=" << diagnostics.exactHessianFallbacks
                << "; inputChordRejections=" << diagnostics.inputChordRejections
                << "; trialChordRejections=" << diagnostics.trialChordRejections
                << "; residualCorrections=" << diagnostics.residualCorrectionSteps
                << "; stagnationTriggers=" << diagnostics.stagnationTriggers
                << "; roundoffAcceptances=" << diagnostics.roundoffAcceptedSteps << std::endl;
    return status;
  }

  int n = 0;
  std::vector<std::shared_ptr<CudaPotentialEnergy>> energies;
  SolverParam parameters;
  Diagnostics diagnostics;
  std::function<void(int, const double *, double, double)> acceptedIterateObserver;
  DeviceBuffer<double> x, gradient, trial, trialGradient, dx, diagonal, r, z, p, Ap, initialX, inputChord;
  DeviceBuffer<DevicePCGState> pcgState;
  DeviceBuffer<double> pcgReductionScalars;
  detail::BlasHandle blas;
  std::unique_ptr<detail::DeviceBlasHandle> deviceBlas;
};

CudaNewtonSolver::CudaNewtonSolver(std::vector<std::shared_ptr<CudaPotentialEnergy>> energies):
  CudaNewtonSolver(std::move(energies), SolverParam{})
{
}
CudaNewtonSolver::CudaNewtonSolver(std::vector<std::shared_ptr<CudaPotentialEnergy>> energies, const SolverParam &parameters):
  impl_(std::make_unique<Impl>(std::move(energies), parameters))
{
}
CudaNewtonSolver::~CudaNewtonSolver() = default;
int CudaNewtonSolver::solve(double *hostX, int maxIterations, double gradientTolerance, int verbose)
{
  return impl_->solve(hostX, maxIterations, gradientTolerance, verbose);
}
const CudaNewtonSolver::Diagnostics &CudaNewtonSolver::getDiagnostics() const
{
  return impl_->diagnostics;
}
void CudaNewtonSolver::setAcceptedIterateObserver(std::function<void(int, const double *, double, double)> observer)
{
  impl_->acceptedIterateObserver = std::move(observer);
}

}  // namespace pgo::NonlinearOptimization
