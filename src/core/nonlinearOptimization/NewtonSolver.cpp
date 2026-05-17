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

const char *NewtonSolver::solveStatusToString(int status)
{
  switch (static_cast<SolveStatus>(status)) {
    case SolveStatus::Converged:
      return "Converged";
    case SolveStatus::MaxIterations:
      return "MaxIterations";
    case SolveStatus::LineSearchFailed:
      return "LineSearchFailed";
    case SolveStatus::StepTooSmall:
      return "StepTooSmall";
    case SolveStatus::NonFinite:
      return "NonFinite";
    case SolveStatus::LinearSolveFailed:
      return "LinearSolveFailed";
  }

  return "Unknown";
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
  historyGradNormMin = 1e100;

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
  if (fixedDOFs_.size() != 0 && fixedDOFs.size() == fixedDOFs_.size() &&
    std::memcmp(fixedDOFs.data(), fixedDOFs_.data(), sizeof(int) * fixedDOFs.size()) == 0) {
  }
  else {
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

#if defined(PGO_HAS_MKL) && !defined(PGO_HAS_ORIG_PARDISO)
      solver = std::make_shared<ES::EigenMKLPardisoSupport>(A11, ES::EigenMKLPardisoSupport::MatrixType::REAL_SYM_INDEFINITE,
        ES::EigenMKLPardisoSupport::ReorderingType::NESTED_DISSECTION, 0, 0, 0, 0, 0, 0);
      solver->analyze(A11);
#elif defined(PGO_HAS_ORIG_PARDISO)
      solver = std::make_shared<ES::EigenOrigPardisoSupport>(A11, ES::EigenOrigPardisoSupport::MatrixType::REAL_SYM_INDEFINITE,
        ES::EigenOrigPardisoSupport::ReorderingType::NESTED_DISSECTION_4, 0, 0, 0, 0, 0, 0);
      solver->analyze(A11);
#else
      solver = std::make_shared<EigenSupport::SymSolver>();
      solver->analyzePattern(A11);
#endif
    }
  }

  fixedValues = ES::Mp<const ES::VXd>(fixedValues_, fixedDOFs_.size());
}

int NewtonSolver::solve(double *x_, int numIter, double epsilon, int verbose)
{
  hclock::time_point t1 = hclock::now();
  int status = static_cast<int>(SolveStatus::MaxIterations);
  solveDiagnostics.reset();

  x.noalias() = Eigen::Map<ES::VXd>(x_, energy->getNumDOFs());

  if (solverParam.sst == SST_SUBITERATION_ONE || solverParam.sst == SST_SUBITERATION_STATIC_DAMPING) {
    historyGradNormMin = 1e100;
  }

  for (size_t i = 0; i < fixedDOFs.size(); i++) {
    x[fixedDOFs[i]] = fixedValues[i];
  }

  int printGap = 10;
  if (verbose == 2) {
    printGap = 10;
  }
  else if (verbose == 3) {
    printGap = 1;
  }

  // double error0 = 0;
  int iter = 0;
  double lambdaScale = 1.0;
  double lambda0 = 1.0;
  int lineSearchFailedTimes = 0;
  // compute lambda initial
  memset(grad.data(), 0, sizeof(double) * grad.size());
  energy->gradient(x, grad);
  filterVector(grad);
  if (!grad.allFinite()) {
    if (verbose >= 1)
      std::cout << "    Newton solve failed before iteration; status=" << solveStatusToString(static_cast<int>(SolveStatus::NonFinite)) << std::endl;
    return static_cast<int>(SolveStatus::NonFinite);
  }
  lambda0 = grad.cwiseAbs().maxCoeff();

  double gradMaxNormLast = lambda0;
  for (; iter < numIter; iter++) {
    if (verbose >= 2 && iter % printGap == 0)
      std::cout << "    Iter=" << iter << std::endl;

    // we solve f(x_i) + K(x_i) deltax = 0
    memset(grad.data(), 0, sizeof(double) * grad.size());
    double eng = energy->func_grad_hessian(x, grad, sysFull);
    if (!std::isfinite(eng)) {
      status = static_cast<int>(SolveStatus::NonFinite);
      if (verbose >= 1)
        std::cout << "    Iter=" << iter << "; energy is non-finite; status=" << solveStatusToString(status) << std::endl;
      break;
    }
    sysFull.makeCompressed();
    filterVector(grad);

    double gradMaxNorm = grad.cwiseAbs().maxCoeff();
    if (!grad.allFinite() || !std::isfinite(gradMaxNorm)) {
      status = static_cast<int>(SolveStatus::NonFinite);
      if (verbose >= 1)
        std::cout << "    Iter=" << iter << "; gradient is non-finite; status=" << solveStatusToString(status) << std::endl;
      break;
    }

    if (verbose >= 2 && iter % printGap == 0)
      std::cout << "        E= " << eng << "; ||grad||_max=" << gradMaxNorm << "; ||x||=" << x.norm() << std::endl;

    // Convergence test: accept either absolute (||grad||_max < eps) or relative
    // reduction from initial gradient (||grad||_max < lambda0 * relTolFactor).
    // Relative branch handles systems where |E| is large enough that the absolute
    // eps becomes unreachable in double precision (line search saturates at FP floor).
    constexpr double relTolFactor = 1e-5;
    const double relThreshold = lambda0 * relTolFactor;
    const bool absConverged = gradMaxNorm < epsilon;
    const bool relConverged = gradMaxNorm < relThreshold;
    if (absConverged || relConverged) {
      status = static_cast<int>(SolveStatus::Converged);
      if (verbose >= 1) {
        std::cout << "    Iter=" << iter << "; ||grad||_max=" << gradMaxNorm
                  << (absConverged ? " < eps" : " < lambda0*relTol")
                  << " (eps=" << epsilon << ", relThreshold=" << relThreshold
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

    // grad too small, we don't need damping
    if (gradMaxNorm < 1e-4) {
      lambdaScale = 0.0;
    }
    else {
      // if lambda too small, we set it to zero
      if (lambdaScale < 1e-8) {
        lambdaScale = 0;
      }
      else {
        // if gradient increasing, we increase lambda
        if (gradMaxNorm > gradMaxNormLast) {
        }
        else {
          // otherwise we decrease lambda
          lambdaScale *= 0.9;
        }

        // clamp
        lambdaScale = std::min(lambdaScale, 1.0);
      }
    }
    gradMaxNormLast = gradMaxNorm;

    // std::cout << "        Damping lambda=" << lambda << std::endl;

    // remove column rows
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

    // std::cout << "      rhs: ";
    // for (int kk = 0; kk < 10; kk++) {
    //   std::cout << rhs[kk] << ' ';
    // }
    // std::cout << std::endl;

    // std::cout << "      A11: ";
    // for (int kk = 0; kk < 10; kk++) {
    //   std::cout << A11.valuePtr()[kk] << ' ';
    // }
    // std::cout << std::endl;

    if (!fixedHessianTopology || solver == nullptr) {
#if defined(PGO_HAS_MKL) && !defined(PGO_HAS_ORIG_PARDISO)
      solver = std::make_shared<ES::EigenMKLPardisoSupport>(A11, ES::EigenMKLPardisoSupport::MatrixType::REAL_SYM_INDEFINITE,
        ES::EigenMKLPardisoSupport::ReorderingType::NESTED_DISSECTION, 0, 0, 0, 0, 0, 0);
      solver->analyze(A11);
#elif defined(PGO_HAS_ORIG_PARDISO)
      solver = std::make_shared<ES::EigenOrigPardisoSupport>(A11, ES::EigenOrigPardisoSupport::MatrixType::REAL_SYM_INDEFINITE,
        ES::EigenOrigPardisoSupport::ReorderingType::NESTED_DISSECTION_4, 0, 0, 0, 0, 0, 0);
      solver->analyze(A11);
#else
      solver = std::make_shared<EigenSupport::SymSolver>();
      solver->analyzePattern(A11);
#endif
    }

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

    if (verbose >= 3 && iter % printGap == 0)
      std::cout << (A11 * deltaxSmall - rhs).norm() << ' ' << rhs.norm() << std::endl;

    memset(deltax.data(), 0, sizeof(double) * n3);
    ES::transferSmallToBig(deltaxSmall, deltax, rhss2b);

    if (!deltax.allFinite()) {
      status = static_cast<int>(SolveStatus::NonFinite);
      if (verbose >= 1)
        std::cout << "    Iter=" << iter << "; dx is non-finite; status=" << solveStatusToString(status) << std::endl;
      break;
    }

    // for (int kk = 0; kk < 10; kk++) {
    //   std::cout << deltax[kk] << ' ';
    // }
    // std::cout << std::endl;

    // if (iter == 0)
    //   error0 = deltax.norm();
    // else {
    //   if (deltax.norm() < epsilon) {
    //     if (verbose >= 1)
    //       std::cout << "    dir too small." << std::endl;

    //     break;
    //   }
    // }

    const double rawStepMaxNorm = deltax.cwiseAbs().maxCoeff();
    if (verbose >= 2 && iter % printGap == 0)
      std::cout << "        rawStepMaxNorm=" << rawStepMaxNorm << std::endl;

    // x += alpha delta x ?
    if (solverParam.sst == SST_SUBITERATION_LINE_SEARCH) {
      double alpha = 1;
      double stepSize = 0;

      const MaxStepResult maxStep = energy->computeMaxStepLimit(x, deltax);
      solveDiagnostics.recordMaxStep(maxStep);
      const double feasibleAlpha = maxStep.alpha;
      if (!std::isfinite(feasibleAlpha)) {
        status = static_cast<int>(SolveStatus::NonFinite);
        if (verbose >= 1)
          std::cout << "    Iter=" << iter << "; feasible alpha is non-finite; status=" << solveStatusToString(status) << std::endl;
        break;
      }
      deltax *= feasibleAlpha;

      lineSearchx.noalias() = x + deltax;
      double eng1 = energy->func(lineSearchx);
      if (!std::isfinite(eng1)) {
        status = static_cast<int>(SolveStatus::NonFinite);
        if (verbose >= 1)
          std::cout << "    Iter=" << iter << "; trial energy is non-finite; status=" << solveStatusToString(status) << std::endl;
        break;
      }
      int maxIter = 50;
      if (eng1 < eng) {
        maxIter = 3;
      }

      if (solverParam.lsm == LSM_GOLDEN) {
        lineSearchHandle->nativeLineSearch->setMaxIterations(maxIter);
        LineSearch::Result ret = lineSearchHandle->nativeLineSearch->golden(x.data(), deltax.data(), eng);
        alpha = ret.alpha;
        eng1 = ret.f;
      }
      else if (solverParam.lsm == LSM_BRENTS) {
        lineSearchHandle->nativeLineSearch->setMaxIterations(maxIter);
        LineSearch::Result ret = lineSearchHandle->nativeLineSearch->BrentsMethod(x.data(), deltax.data(), eng);
        alpha = ret.alpha;
        eng1 = ret.f;
      }
      else if (solverParam.lsm == LSM_BACKTRACK) {
        lineSearchHandle->nativeLineSearch->setMaxIterations(maxIter);
        LineSearch::Result ret = lineSearchHandle->nativeLineSearch->backtracking(x.data(), deltax.data(), eng, grad.data(), 0.0001, 0.5);
        alpha = ret.alpha;
        eng1 = ret.f;
      }
      else if (solverParam.lsm == LSM_SIMPLE) {
        eng1 = eng;
        for (int i = 0; i < 100; i++) {
          lineSearchx.noalias() = x + deltax * alpha;
          eng1 = energy->func(lineSearchx);
          if (!std::isfinite(eng1)) {
            status = static_cast<int>(SolveStatus::NonFinite);
            break;
          }

          if (eng1 < eng) {
            break;
          }

          alpha *= 0.5;
        }
      }

      if (!std::isfinite(alpha) || !std::isfinite(eng1))
        status = static_cast<int>(SolveStatus::NonFinite);

      if (status == static_cast<int>(SolveStatus::NonFinite)) {
        if (verbose >= 1)
          std::cout << "    Iter=" << iter << "; line search energy is non-finite; status=" << solveStatusToString(status) << std::endl;
        break;
      }

      if (eng1 > eng) {
        // Loose relative fallback: if line search can't find descent but Newton already
        // reduced the gradient by 4+ orders of magnitude from the initial state, treat
        // this as converged-at-FP-limit rather than failure.
        constexpr double looseRelFactor = 1e-4;
        const bool looseRelConverged = gradMaxNorm < lambda0 * looseRelFactor;
        status = (gradMaxNorm < epsilon || looseRelConverged) ? static_cast<int>(SolveStatus::Converged) : static_cast<int>(SolveStatus::LineSearchFailed);
        if (verbose >= 1) {
          std::cout << "    Iter=" << iter << "; line search failed; ||grad||_max=" << gradMaxNorm
                    << " (lambda0=" << lambda0 << ", looseRel=" << lambda0 * looseRelFactor << ")"
                    << "; status=" << solveStatusToString(status) << ". Times: " << lineSearchFailedTimes << std::endl;
        }
        break;

        // solverParam.addDamping = 1;
        // if (lineSearchFailedTimes == 0) {
        //   lambdaScale = 1.0;
        //   lambda0 = grad.cwiseAbs().maxCoeff();
        // }
        // else {
        //   lambdaScale *= 2.0;
        // }

        // if (lineSearchFailedTimes++ >= 3) {
        //   if (verbose >= 1) {
        //     std::cout << "    Iter=" << iter << "; line search failed too many times. Stop." << std::endl;
        //   }

        //   break;
        // }
        // else {
        //   continue;
        // }
      }

      const double effectiveAlpha = feasibleAlpha * alpha;
      const double acceptedStepMaxNorm = std::abs(alpha) * deltax.cwiseAbs().maxCoeff();
      solveDiagnostics.recordLineSearch(feasibleAlpha, alpha, effectiveAlpha);

      if (verbose >= 2 && iter % printGap == 0) {
        std::cout << "        feasibleAlpha=" << feasibleAlpha << std::endl;
        if (feasibleAlpha < 1.0) {
          std::cout << "        feasible alpha clamped: material:" << solveDiagnostics.currentMaterialAlpha
                    << " contact:" << solveDiagnostics.currentContactAlpha << std::endl;
        }
        std::cout << "        lineSearchAlpha=" << alpha << std::endl;
        std::cout << "        effectiveAlpha=" << effectiveAlpha << std::endl;
        std::cout << "        acceptedStepMaxNorm=" << acceptedStepMaxNorm << std::endl;
        std::cout << "        acceptedEnergy=" << eng1 << std::endl;
      }

      x += deltax * alpha;

      stepSize = acceptedStepMaxNorm;
      if (stepSize < 1e-15) {
        // Same loose relative fallback as the line-search-failed branch above.
        constexpr double looseRelFactor = 1e-4;
        const bool looseRelConverged = gradMaxNorm < lambda0 * looseRelFactor;
        status = (gradMaxNorm < epsilon || looseRelConverged) ? static_cast<int>(SolveStatus::Converged) : static_cast<int>(SolveStatus::StepTooSmall);
        if (verbose >= 1) {
          std::cout << "    Iter=" << iter << "; dx = " << stepSize
                    << "; dx too small; ||grad||_max=" << gradMaxNorm
                    << " (lambda0=" << lambda0 << ", looseRel=" << lambda0 * looseRelFactor << ")"
                    << "; status=" << solveStatusToString(status) << ". Times: " << lineSearchFailedTimes << std::endl;
        }
        break;

        // solverParam.addDamping = 1;
        // if (lineSearchFailedTimes == 0) {
        //   lambda0 = grad.cwiseAbs().maxCoeff();
        //   lambdaScale = 1.0;
        // }
        // else {
        //   lambdaScale *= 2.0;
        // }

        // if (lineSearchFailedTimes++ >= 3) {
        //   if (verbose >= 1) {
        //     std::cout << "    Iter=" << iter << "; line search failed too many times. Stop." << std::endl;
        //   }

        //   break;
        // }
        // else {
        //   continue;
        // }
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
  }

  if (solverParam.sst == SST_SUBITERATION_ONE || solverParam.sst == SST_SUBITERATION_STATIC_DAMPING) {
    if (historyGradNormMin < epsilon)
      status = static_cast<int>(SolveStatus::Converged);

    if (verbose >= 1)
      std::cout << "        Final ||grad||=" << historyGradNormMin << std::endl;

    Eigen::Map<ES::VXd>(x_, energy->getNumDOFs()) = historyx;
  }
  else
    Eigen::Map<ES::VXd>(x_, energy->getNumDOFs()) = x;

  hclock::time_point t2 = hclock::now();

  double timeCost = dura(t1, t2);

  SPDLOG_LOGGER_INFO(Logging::lgr(), "Newton solve time: {}", timeCost);
  return status;
}

void NewtonSolver::filterVector(ES::VXd &v)
{
  for (int dof : fixedDOFs)
    v[dof] = 0;
}

//
// OPK_DROP(lineSearchHandle);
