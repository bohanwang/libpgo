#include "NewtonSolver.h"

#include "EigenSupport.h"
#include "lineSearch.h"
#include "pgoLogging.h"
#include "scopedProfileSection.h"

#include <cmath>
#include <iostream>
#include <numeric>
#include <chrono>

using namespace pgo;
using namespace pgo::NonlinearOptimization;
namespace ES = pgo::EigenSupport;
using hclock = std::chrono::high_resolution_clock;

namespace
{
constexpr double kRelTolFactor = 1e-5;       // absolute->relative gradient tolerance scale
constexpr double kLooseRelFactor = 1e-4;     // FP-limit "good enough" relative reduction
constexpr double kGradSmallThreshold = 1e-4; // below this gradient, drop damping entirely
constexpr double kLambdaScaleFloor = 1e-8;   // below this damping scale, snap to zero
constexpr double kDampingDecay = 0.9;        // per-iteration damping decay when gradient is not increasing
constexpr double kStepTooSmallEps = 1e-15;   // accepted step max-norm below this == stalled
constexpr double kHistoryGradNormInit = 1e100;
constexpr double kBacktrackArmijo = 0.0001;
constexpr double kBacktrackShrink = 0.5;
constexpr double kBacktrackInitAlpha = 1.0;
constexpr int kLineSearchMaxIter = 50;
constexpr int kLineSearchMaxIterDescent = 3; // fewer iters when the full step already decreases energy
constexpr int kSimpleLineSearchMaxIter = 100;

class LineSearchScope
{
public:
  LineSearchScope(const PotentialEnergy_const_p &energy, EigenSupport::ConstRefVecXd x,
    EigenSupport::ConstRefVecXd dx, bool active):
    energy_(energy), active_(active)
  {
    if (active_)
      energy_->beginLineSearch(x, dx);
  }

  ~LineSearchScope()
  {
    if (active_)
      energy_->endLineSearch();
  }

  LineSearchScope(const LineSearchScope &) = delete;
  LineSearchScope &operator=(const LineSearchScope &) = delete;

private:
  const PotentialEnergy_const_p &energy_;
  bool active_ = false;
};
}  // namespace

namespace pgo::NonlinearOptimization
{
class LineSearchHandle
{
public:
  std::shared_ptr<LineSearch> nativeLineSearch;
};
}  // namespace pgo::NonlinearOptimization

inline double dura(const hclock::time_point &t1, const hclock::time_point &t2)
{
  return std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1e6;
}

NewtonSolver::NewtonSolver(const double *x_, SolverParam sp, PotentialEnergy_const_p energy_, const std::vector<int> &fixedDOFs_, const double *fixedValues_):
  energy(energy_), solverParam(sp)
{
  n3 = (int)energy->getNumDOFs();
  allDOFs.resize(energy->getNumDOFs());
  std::iota(allDOFs.begin(), allDOFs.end(), 0);

  grad.resize(energy->getNumDOFs());
  x = Eigen::Map<const ES::VXd>(x_, energy->getNumDOFs());
  deltax.resize(energy->getNumDOFs());
  lineSearchx.resize(energy->getNumDOFs());
  historyx.resize(energy->getNumDOFs());
  historyGradNormMin = kHistoryGradNormInit;

  setFixedDOFs(fixedDOFs_, fixedValues_);

  if (solverParam.sst == SST_SUBITERATION_LINE_SEARCH) {
    lineSearchHandle = std::make_shared<LineSearchHandle>();

    LineSearch::EvaluateFunction evalFunc = [this](const double *x, double *f, double *grad) -> int {
      if (f)
        *f = energy->func(Eigen::Map<const ES::VXd>(x, n3));

      if (grad) {
        memset(grad, 0, sizeof(double) * n3);
        energy->gradient(Eigen::Map<const ES::VXd>(x, n3), Eigen::Map<ES::VXd>(grad, n3));
      }

      return 0;
    };

    lineSearchHandle->nativeLineSearch = std::make_shared<LineSearch>(n3, evalFunc);
  }
}

void NewtonSolver::setFixedDOFs(const std::vector<int> &fixedDOFs_, const double *fixedValues_)
{
  const bool sameAsBefore = fixedDOFs_.size() != 0 && fixedDOFs.size() == fixedDOFs_.size() &&
    std::memcmp(fixedDOFs.data(), fixedDOFs_.data(), sizeof(int) * fixedDOFs.size()) == 0;

  if (!sameAsBefore) {
    fixedDOFs = fixedDOFs_;

    // dofs
    rhsb2s.clear();
    rhss2b.clear();
    ES::removeRows(n3, fixedDOFs, rhsb2s, rhss2b);

    rhs.resize(n3 - (int)fixedDOFs.size());
    deltaxSmall.resize(n3 - (int)fixedDOFs.size());

    if (energy->isHessianTopologyFixed()) {
      // sparse matrix
      energy->createHessian(sysFull);
      energy->hessian(x, sysFull);

      ES::removeRowsCols(sysFull, fixedDOFs, A11);
      ES::removeRowsCols(sysFull, A11, fixedDOFs, A11Mapping);

      makeLinearSolver(A11);
    }
  }

  fixedValues = ES::Mp<const ES::VXd>(fixedValues_, fixedDOFs_.size());
}

void NewtonSolver::makeLinearSolver(const ES::SpMatD &A)
{
#if defined(PGO_HAS_MKL) && !defined(PGO_HAS_ORIG_PARDISO)
  solver = std::make_shared<ES::EigenMKLPardisoSupport>(A, ES::EigenMKLPardisoSupport::MatrixType::REAL_SYM_INDEFINITE,
    ES::EigenMKLPardisoSupport::ReorderingType::NESTED_DISSECTION, 0, 0, 0, 0, 0, 0);
  solver->analyze(A);
#elif defined(PGO_HAS_ORIG_PARDISO)
  solver = std::make_shared<ES::EigenOrigPardisoSupport>(A, ES::EigenOrigPardisoSupport::MatrixType::REAL_SYM_INDEFINITE,
    ES::EigenOrigPardisoSupport::ReorderingType::NESTED_DISSECTION_4, 0, 0, 0, 0, 0, 0);
  solver->analyze(A);
#else
  solver = std::make_shared<EigenSupport::SymSolver>();
  solver->analyzePattern(A);
#endif
}

SolverResult NewtonSolver::solve(double *x_, int numIter, double epsilon, int verbose)
{
  hclock::time_point t1 = hclock::now();
  SolveStatus status = SolveStatus::MaxIterations;
  int completedIterations = 0;
  solveDiagnostics.reset();

  x.noalias() = Eigen::Map<ES::VXd>(x_, energy->getNumDOFs());

  if (solverParam.sst == SST_SUBITERATION_ONE || solverParam.sst == SST_SUBITERATION_STATIC_DAMPING) {
    historyGradNormMin = kHistoryGradNormInit;
  }

  applyFixedValues();

  int printGap = 10;
  if (verbose == 2) {
    printGap = 10;
  }
  else if (verbose == 3) {
    printGap = 1;
  }

  int iter = 0;
  double lambdaScale = 1.0;
  double lambda0 = 1.0;
  int lineSearchFailedTimes = 0;
  bool hasInitialGradNorm = false;
  double gradMaxNormLast = 0.0;
  for (; iter < numIter; iter++) {
    if (verbose >= 2 && iter % printGap == 0)
      std::cout << "    Iter=" << iter << std::endl;

    // we solve f(x_i) + K(x_i) deltax = 0
    const IterationState state = evaluateCurrentState(iter, epsilon, lambda0, hasInitialGradNorm);
    if (state.nonFiniteEnergy) {
      status = SolveStatus::NonFinite;
      if (verbose >= 1)
        std::cout << "    Iter=" << iter << "; energy is non-finite; status=" << solveStatusToString(status) << std::endl;
      break;
    }
    if (state.nonFiniteGradient) {
      status = SolveStatus::NonFinite;
      if (verbose >= 1)
        std::cout << "    Iter=" << iter << "; gradient is non-finite; status=" << solveStatusToString(status) << std::endl;
      break;
    }

    const double eng = state.energy;
    const double gradMaxNorm = state.gradMaxNorm;
    if (!hasInitialGradNorm) {
      lambda0 = state.lambda0;
      gradMaxNormLast = gradMaxNorm;
      hasInitialGradNorm = true;
    }

    if (verbose >= 2 && iter % printGap == 0)
      std::cout << "        E= " << eng << "; ||grad||_max=" << gradMaxNorm << "; ||x||=" << x.norm() << std::endl;

    // Convergence test: accept either absolute (||grad||_max < eps) or relative
    // reduction from initial gradient (||grad||_max < lambda0 * relTolFactor).
    // Relative branch handles systems where |E| is large enough that the absolute
    // eps becomes unreachable in double precision (line search saturates at FP floor).
    if (isConverged(state)) {
      status = SolveStatus::Converged;
      completedIterations = iter;
      solveDiagnostics.recordFinalGradientStats(state.gradNorm, gradMaxNorm);
      if (verbose >= 1) {
        std::cout << "    Iter=" << iter << "; ||grad||_max=" << gradMaxNorm
                  << (state.absConverged ? " < eps" : " < lambda0*relTol")
                  << " (eps=" << epsilon << ", relThreshold=" << state.relThreshold
                  << "). Done.; status=" << solveStatusToString(status) << std::endl;
      }

      break;
    }

    if (solverParam.sst == SST_SUBITERATION_ONE || solverParam.sst == SST_SUBITERATION_STATIC_DAMPING) {
      if (gradMaxNorm < historyGradNormMin) {
        historyx.noalias() = x;
        historyGradNormMin = gradMaxNorm;
      }
      else {
        if (solverParam.stopAfterIncrease)
          break;
      }
    }

    lambdaScale = updateDampingScale(lambdaScale, gradMaxNorm, gradMaxNormLast);
    gradMaxNormLast = gradMaxNorm;

    const bool fixedHessianTopology = prepareReducedSystem(lambdaScale, lambda0);

    ensureLinearSolver(fixedHessianTopology);
    if (!solveReducedNewtonDirection(fixedHessianTopology)) {
      status = SolveStatus::NonFinite;
      completedIterations = iter + 1;
      if (verbose >= 1)
        std::cout << "    Iter=" << iter << "; dx is non-finite; status=" << solveStatusToString(status) << std::endl;
      break;
    }

    if (verbose >= 3 && iter % printGap == 0)
      std::cout << (A11 * deltaxSmall - rhs).norm() << ' ' << rhs.norm() << std::endl;

    if (!expandReducedStep()) {
      status = SolveStatus::NonFinite;
      completedIterations = iter + 1;
      if (verbose >= 1)
        std::cout << "    Iter=" << iter << "; dx is non-finite; status=" << solveStatusToString(status) << std::endl;
      break;
    }

    const double rawStepMaxNorm = deltax.cwiseAbs().maxCoeff();
    if (verbose >= 2 && iter % printGap == 0)
      std::cout << "        rawStepMaxNorm=" << rawStepMaxNorm << std::endl;

    // x += alpha delta x ?
    if (solverParam.sst == SST_SUBITERATION_LINE_SEARCH) {
      const StepAcceptance accepted = runLineSearchStep(eng, verbose, printGap, iter);
      if (accepted.nonFinite()) {
        status = SolveStatus::NonFinite;
        completedIterations = iter + 1;
        break;
      }

      if (accepted.acceptedEnergy > eng) {
        // Loose relative fallback: if line search can't find descent but Newton already
        // reduced the gradient by 4+ orders of magnitude from the initial state, treat
        // this as converged-at-FP-limit rather than failure.
        status = resolveFpLimitFallback(SolveStatus::LineSearchFailed, gradMaxNorm, lambda0, epsilon);
        completedIterations = iter + 1;
        if (verbose >= 1) {
          std::cout << "    Iter=" << iter << "; line search failed; ||grad||_max=" << gradMaxNorm
                    << " (lambda0=" << lambda0 << ", looseRel=" << lambda0 * kLooseRelFactor << ")"
                    << "; status=" << solveStatusToString(status) << ". Times: " << lineSearchFailedTimes << std::endl;
        }
        break;
      }

      x += deltax * accepted.lineSearchAlpha;

      const double stepSize = accepted.acceptedStepMaxNorm;
      if (stepSize < kStepTooSmallEps) {
        // Same loose relative fallback as the line-search-failed branch above.
        status = resolveFpLimitFallback(SolveStatus::StepTooSmall, gradMaxNorm, lambda0, epsilon);
        completedIterations = iter + 1;
        if (verbose >= 1) {
          std::cout << "    Iter=" << iter << "; dx = " << stepSize
                    << "; dx too small; ||grad||_max=" << gradMaxNorm
                    << " (lambda0=" << lambda0 << ", looseRel=" << lambda0 * kLooseRelFactor << ")"
                    << "; status=" << solveStatusToString(status) << ". Times: " << lineSearchFailedTimes << std::endl;
        }
        break;
      }
    }
    else if (solverParam.sst == SST_SUBITERATION_ONE) {
      x += deltax;
      historyx.noalias() = x;

      memset(grad.data(), 0, sizeof(double) * grad.size());
      energy->gradient(x, grad);
      filterVector(grad);

      historyGradNormMin = grad.norm();

      if (verbose >= 2 && iter % printGap == 0)
        std::cout << "    f=" << energy->func(x) << std::endl;
    }
    else if (solverParam.sst == SST_SUBITERATION_STATIC_DAMPING) {
      if (verbose >= 2 && iter % printGap == 0)
        std::cout << "        E= " << eng << "; ||grad||_max=" << grad.cwiseAbs().maxCoeff() << "; ||grad||=" << grad.norm() << std::endl;

      x += deltax * solverParam.alpha;
    }

    if (stepFunc) {
      stepFunc(x, iter);
    }

    completedIterations = iter + 1;
  }

  if (solverParam.sst == SST_SUBITERATION_ONE || solverParam.sst == SST_SUBITERATION_STATIC_DAMPING) {
    if (historyGradNormMin < epsilon)
      status = SolveStatus::Converged;

    if (verbose >= 1)
      std::cout << "        Final ||grad||=" << historyGradNormMin << std::endl;

    Eigen::Map<ES::VXd>(x_, energy->getNumDOFs()) = historyx;
  }
  else
    Eigen::Map<ES::VXd>(x_, energy->getNumDOFs()) = x;

  hclock::time_point t2 = hclock::now();

  double timeCost = dura(t1, t2);

  SPDLOG_LOGGER_INFO(Logging::lgr(), "Newton solve time: {}", timeCost);
  SolverResult result;
  result.status = status;
  result.iterations = completedIterations;
  result.rawStatusCode = static_cast<int>(status);
  result.diagnostics = solveDiagnostics;
  if (solveDiagnostics.hasFinalGradientStats) {
    result.hasFinalGradientStats = true;
    result.finalGradientNorm = solveDiagnostics.finalGradientNorm;
    result.finalGradientMaxNorm = solveDiagnostics.finalGradientMaxNorm;
  }
  return result;
}

void NewtonSolver::filterVector(ES::VXd &v)
{
  for (int dof : fixedDOFs)
    v[dof] = 0;
}

void NewtonSolver::applyFixedValues()
{
  for (size_t i = 0; i < fixedDOFs.size(); i++) {
    x[fixedDOFs[i]] = fixedValues[i];
  }
}

NewtonSolver::IterationState NewtonSolver::evaluateCurrentState(int iter, double epsilon, double lambda0, bool hasInitialGradNorm)
{
  IterationState state;
  state.iter = iter;

  memset(grad.data(), 0, sizeof(double) * grad.size());
  state.energy = energy->func_grad_hessian(x, grad, sysFull);
  if (!std::isfinite(state.energy)) {
    state.nonFiniteEnergy = true;
    return state;
  }

  sysFull.makeCompressed();
  filterVector(grad);

  state.gradMaxNorm = grad.cwiseAbs().maxCoeff();
  state.gradNorm = grad.norm();
  if (!grad.allFinite() || !std::isfinite(state.gradMaxNorm)) {
    state.nonFiniteGradient = true;
    return state;
  }

  state.lambda0 = hasInitialGradNorm ? lambda0 : state.gradMaxNorm;
  state.relThreshold = state.lambda0 * kRelTolFactor;
  state.absConverged = state.gradMaxNorm < epsilon;
  state.relConverged = state.gradMaxNorm < state.relThreshold;
  return state;
}

bool NewtonSolver::isConverged(const IterationState &state) const
{
  return state.absConverged || state.relConverged;
}

bool NewtonSolver::looseRelativeConverged(double gradMaxNorm, double lambda0) const
{
  return gradMaxNorm < lambda0 * kLooseRelFactor;
}

double NewtonSolver::updateDampingScale(double lambdaScale, double gradMaxNorm, double gradMaxNormLast) const
{
  if (gradMaxNorm < kGradSmallThreshold)  // grad too small: drop damping entirely
    return 0.0;
  if (lambdaScale < kLambdaScaleFloor)  // damping already negligible: snap to zero
    return 0.0;
  // When the gradient is increasing we deliberately leave lambdaScale unchanged;
  // only decay it once the gradient stops increasing.
  if (gradMaxNorm <= gradMaxNormLast)
    lambdaScale *= kDampingDecay;
  return std::min(lambdaScale, 1.0);
}

SolveStatus NewtonSolver::resolveFpLimitFallback(SolveStatus failStatus, double gradMaxNorm, double lambda0, double epsilon)
{
  const bool converged = gradMaxNorm < epsilon || looseRelativeConverged(gradMaxNorm, lambda0);
  if (converged) {
    solveDiagnostics.recordFinalGradientStats(grad.norm(), gradMaxNorm);
    return SolveStatus::Converged;
  }
  return failStatus;
}

bool NewtonSolver::prepareReducedSystem(double lambdaScale, double lambda0)
{
  const bool fixedHessianTopology = energy->isHessianTopologyFixed();
  if (fixedHessianTopology) {
    ES::transferBigToSmall(sysFull, A11, A11Mapping, 1);
  }
  else {
    ES::removeRowsCols(sysFull, fixedDOFs, A11);
  }

  if (solverParam.addDamping) {
    for (int i = 0; i < A11.rows(); i++) {
      A11.coeffRef(i, i) += lambdaScale * lambda0;
    }
  }

  ES::transferBigToSmall(grad, rhs, rhsb2s, 1);
  rhs *= -1.0;

  return fixedHessianTopology;
}

void NewtonSolver::ensureLinearSolver(bool fixedHessianTopology)
{
  if (!fixedHessianTopology || solver == nullptr) {
    makeLinearSolver(A11);
  }
}

bool NewtonSolver::solveReducedNewtonDirection(bool fixedHessianTopology)
{
#if defined(PGO_HAS_MKL) && !defined(PGO_HAS_ORIG_PARDISO)
  {
    Profiling::ScopedProfileSection scopedProfile("solver.linear_solve");
    solver->factorize(A11);
    solver->solve(A11, deltaxSmall.data(), rhs.data(), 1);
  }
#elif defined(PGO_HAS_ORIG_PARDISO)
  {
    Profiling::ScopedProfileSection scopedProfile("solver.linear_solve");
    solver->factorize(A11);
    solver->solve(A11, deltaxSmall.data(), rhs.data(), 1);
  }
#else
  {
    Profiling::ScopedProfileSection scopedProfile("solver.linear_solve");
    solver->factorize(A11);
    deltaxSmall.noalias() = solver->solve(rhs);
  }
#endif

  if (fixedHessianTopology) {
    solver.reset();  // free symbolic factorization memory since we won't reuse it anymore
  }

  return deltaxSmall.allFinite();
}

bool NewtonSolver::expandReducedStep()
{
  memset(deltax.data(), 0, sizeof(double) * n3);
  ES::transferSmallToBig(deltaxSmall, deltax, rhss2b);
  return deltax.allFinite();
}

NewtonSolver::StepAcceptance NewtonSolver::runLineSearchStep(double currentEnergy, int verbose, int printGap, int iter)
{
  StepAcceptance accepted;

  const MaxStepResult maxStep = energy->computeMaxStepLimit(x, deltax);
  solveDiagnostics.recordMaxStep(maxStep);
  accepted.feasibleAlpha = maxStep.alpha;
  if (!std::isfinite(accepted.feasibleAlpha)) {
    accepted.nonFiniteReason = StepAcceptance::NonFiniteReason::FeasibleAlpha;
    if (verbose >= 1)
      std::cout << "    Iter=" << iter << "; feasible alpha is non-finite; status=" << solveStatusToString(SolveStatus::NonFinite) << std::endl;
    return accepted;
  }
  deltax *= accepted.feasibleAlpha;

  {
    // Backtracking and simple line search only evaluate alpha in [0, 1].
    // Golden/Brent may expand the bracket past 1, which is outside the IPC
    // swept active-set superset built for this Newton step.
    const bool useLineSearchActiveSet = solverParam.lsm == LSM_BACKTRACK || solverParam.lsm == LSM_SIMPLE;
    LineSearchScope lineSearchScope(energy, x, deltax, useLineSearchActiveSet);

    lineSearchx.noalias() = x + deltax;
    accepted.acceptedEnergy = energy->func(lineSearchx);
    if (!std::isfinite(accepted.acceptedEnergy)) {
      accepted.nonFiniteReason = StepAcceptance::NonFiniteReason::TrialEnergy;
      if (verbose >= 1)
        std::cout << "    Iter=" << iter << "; trial energy is non-finite; status=" << solveStatusToString(SolveStatus::NonFinite) << std::endl;
    }
    else {
      int maxIter = kLineSearchMaxIter;
      if (accepted.acceptedEnergy < currentEnergy) {
        maxIter = kLineSearchMaxIterDescent;
      }

      if (solverParam.lsm == LSM_GOLDEN) {
        lineSearchHandle->nativeLineSearch->setMaxIterations(maxIter);
        LineSearch::Result ret = lineSearchHandle->nativeLineSearch->golden(x.data(), deltax.data(), currentEnergy);
        accepted.lineSearchAlpha = ret.alpha;
        accepted.acceptedEnergy = ret.f;
      }
      else if (solverParam.lsm == LSM_BRENTS) {
        lineSearchHandle->nativeLineSearch->setMaxIterations(maxIter);
        LineSearch::Result ret = lineSearchHandle->nativeLineSearch->BrentsMethod(x.data(), deltax.data(), currentEnergy);
        accepted.lineSearchAlpha = ret.alpha;
        accepted.acceptedEnergy = ret.f;
      }
      else if (solverParam.lsm == LSM_BACKTRACK) {
        lineSearchHandle->nativeLineSearch->setMaxIterations(maxIter);
        LineSearch::Result ret = lineSearchHandle->nativeLineSearch->backtrackingWithInitialValue(
          x.data(), deltax.data(), currentEnergy, grad.data(), kBacktrackArmijo, kBacktrackShrink, kBacktrackInitAlpha, accepted.acceptedEnergy);
        accepted.lineSearchAlpha = ret.alpha;
        accepted.acceptedEnergy = ret.f;
      }
      else if (solverParam.lsm == LSM_SIMPLE) {
        accepted.acceptedEnergy = currentEnergy;
        for (int i = 0; i < kSimpleLineSearchMaxIter; i++) {
          lineSearchx.noalias() = x + deltax * accepted.lineSearchAlpha;
          accepted.acceptedEnergy = energy->func(lineSearchx);
          if (!std::isfinite(accepted.acceptedEnergy)) {
            accepted.nonFiniteReason = StepAcceptance::NonFiniteReason::LineSearchResult;
            break;
          }

          if (accepted.acceptedEnergy < currentEnergy) {
            break;
          }

          accepted.lineSearchAlpha *= 0.5;
        }
      }

      if (!std::isfinite(accepted.lineSearchAlpha) || !std::isfinite(accepted.acceptedEnergy))
        accepted.nonFiniteReason = StepAcceptance::NonFiniteReason::LineSearchResult;
    }
  }

  if (accepted.nonFinite()) {
    if (verbose >= 1 && accepted.nonFiniteReason == StepAcceptance::NonFiniteReason::LineSearchResult)
      std::cout << "    Iter=" << iter << "; line search energy is non-finite; status=" << solveStatusToString(SolveStatus::NonFinite) << std::endl;
    return accepted;
  }

  accepted.effectiveAlpha = accepted.feasibleAlpha * accepted.lineSearchAlpha;
  accepted.acceptedStepMaxNorm = std::abs(accepted.lineSearchAlpha) * deltax.cwiseAbs().maxCoeff();
  solveDiagnostics.recordLineSearch(accepted.feasibleAlpha, accepted.lineSearchAlpha, accepted.effectiveAlpha);

  if (verbose >= 2 && iter % printGap == 0) {
    std::cout << "        feasibleAlpha=" << accepted.feasibleAlpha << std::endl;
    if (accepted.feasibleAlpha < 1.0) {
      std::cout << "        feasible alpha clamped: material:" << solveDiagnostics.currentMaterialAlpha
                << " contact:" << solveDiagnostics.currentContactAlpha << std::endl;
    }
    std::cout << "        lineSearchAlpha=" << accepted.lineSearchAlpha << std::endl;
    std::cout << "        effectiveAlpha=" << accepted.effectiveAlpha << std::endl;
    std::cout << "        acceptedStepMaxNorm=" << accepted.acceptedStepMaxNorm << std::endl;
    std::cout << "        acceptedEnergy=" << accepted.acceptedEnergy << std::endl;
  }

  return accepted;
}
