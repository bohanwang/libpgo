#include "NewtonSolver.h"

#include "EigenSupport.h"
#include "lineSearch.h"
#include "pgoLogging.h"

#include <cmath>
#include <iostream>
#include <limits>
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

const char *NewtonSolver::solveStatusName(int status)
{
  switch (status) {
  case SOLVE_CONVERGED:
    return "converged";
  case SOLVE_NOT_CONVERGED:
    return "not_converged";
  case SOLVE_NUMERICAL_FAILURE:
    return "numerical_failure";
  default:
    return "unknown";
  }
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

  x.noalias() = Eigen::Map<ES::VXd>(x_, energy->getNumDOFs());

  if (solverParam.sst == SST_SUBITERATION_ONE || solverParam.sst == SST_SUBITERATION_STATIC_DAMPING) {
    historyGradNormMin = 1e100;
  }

  double eng0 = energy->func(x);
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
  const char *stopReason = "iteration_limit";
  bool numericalFailure = false;
  double lambdaScale = 1.0;
  double lambda0 = 1.0;
  int lineSearchFailedTimes = 0;
  // compute lambda initial
  memset(grad.data(), 0, sizeof(double) * grad.size());
  energy->gradient(x, grad);
  filterVector(grad);
  lambda0 = grad.cwiseAbs().maxCoeff();

  double gradNormLast = lambda0;
  for (; iter < numIter; iter++) {
    if (verbose >= 2 && iter % printGap == 0)
      std::cout << "    Iter=" << iter << std::endl;

    double eng = energy->func(x);

    // we solve f(x_i) + K(x_i) deltax = 0
    memset(grad.data(), 0, sizeof(double) * grad.size());
    energy->gradient(x, grad);
    filterVector(grad);

    if (!std::isfinite(eng) || !grad.allFinite()) {
      if (verbose >= 1)
        std::cout << "    Iter=" << iter << "; non-finite energy or gradient. Stop." << std::endl;
      stopReason = "non_finite_energy";
      numericalFailure = true;
      break;
    }

    double gradNorm = grad.cwiseAbs().maxCoeff();
    if (verbose >= 2 && iter % printGap == 0)
      std::cout << "        E= " << eng << "; ||grad||_max=" << grad.cwiseAbs().maxCoeff() << "; ||x||=" << x.norm() << "; ||grad||=" << gradNorm << std::endl;

    if (solverParam.sst == SST_SUBITERATION_ONE || solverParam.sst == SST_SUBITERATION_STATIC_DAMPING) {
      if (gradNorm < historyGradNormMin) {
        historyx.noalias() = x;
        historyGradNormMin = gradNorm;
      }
      else {
        if (solverParam.stopAfterIncrease) {
          stopReason = "gradient_increase";
          break;
        }
      }
    }

    if (iter) {
      if (gradNorm < epsilon) {
        stopReason = "gradient_tolerance";
        if (verbose >= 1) {
          std::cout << "    Iter=" << iter << "; ||grad||=" << gradNorm << " < eps. Done." << std::endl;
        }

        break;
      }
    }

    if (energy->isHessianTopologyFixed()) {
      memset(sysFull.valuePtr(), 0, sizeof(double) * sysFull.nonZeros());
      energy->hessian(x, sysFull);
    }
    else {
      energy->hessianDirect(x, sysFull);
      sysFull.makeCompressed();
    }

    // grad too small, we don't need damping
    if (gradNorm < 1e-4) {
      lambdaScale = 0.0;
    }
    else {
      // if lambda too small, we set it to zero
      if (lambdaScale < 1e-8) {
        lambdaScale = 0;
      }
      else {
        // if gradient increasing, we increase lambda
        if (gradNorm > gradNormLast) {
        }
        else {
          // otherwise we decrease lambda
          lambdaScale *= 0.9;
        }

        // clamp
        lambdaScale = std::min(lambdaScale, 1.0);
      }
    }
    gradNormLast = gradNorm;

    // std::cout << "        Damping lambda=" << lambda << std::endl;

    // remove column rows
    if (energy->isHessianTopologyFixed()) {
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

    if (!energy->isHessianTopologyFixed() || solver == nullptr) {
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
    solver->factorize(A11);
    solver->solve(A11, deltaxSmall.data(), rhs.data(), 1);
#elif defined(PGO_HAS_ORIG_PARDISO)
    solver->factorize(A11);
    solver->solve(A11, deltaxSmall.data(), rhs.data(), 1);
#else
    solver->factorize(A11);
    deltaxSmall.noalias() = solver->solve(rhs);
#endif

    if (verbose >= 3 && iter % printGap == 0)
      std::cout << (A11 * deltaxSmall - rhs).norm() << ' ' << rhs.norm() << std::endl;

    memset(deltax.data(), 0, sizeof(double) * n3);
    ES::transferSmallToBig(deltaxSmall, deltax, rhss2b);

    if (!deltax.allFinite()) {
      if (verbose >= 1)
        std::cout << "    Iter=" << iter << "; non-finite Newton step. Stop." << std::endl;
      stopReason = "non_finite_step";
      numericalFailure = true;
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

      const double feasibleAlpha = energy->computeMaxStepSize(x, deltax);
      deltax *= feasibleAlpha;

      lineSearchx.noalias() = x + deltax;
      double eng1 = energy->func(lineSearchx);
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

          if (eng1 < eng) {
            break;
          }

          alpha *= 0.5;
        }
      }

      // The stopping decision must not depend on the verbosity level: a failed
      // line search leaves x at the last accepted iterate and ends the solve.
      if (eng1 > eng) {
        if (verbose >= 1)
          std::cout << "    Iter=" << iter << "; line search failed. Times: " << lineSearchFailedTimes << std::endl;
        stopReason = "line_search_failed";
        break;
      }

      const double effectiveAlpha = feasibleAlpha * alpha;
      const double acceptedStepMaxNorm = std::abs(alpha) * deltax.cwiseAbs().maxCoeff();

      if (verbose >= 2 && iter % printGap == 0) {
        std::cout << "        feasibleAlpha=" << feasibleAlpha << std::endl;
        std::cout << "        lineSearchAlpha=" << alpha << std::endl;
        std::cout << "        effectiveAlpha=" << effectiveAlpha << std::endl;
        std::cout << "        acceptedStepMaxNorm=" << acceptedStepMaxNorm << std::endl;
        std::cout << "        acceptedEnergy=" << eng1 << std::endl;
      }

      x += deltax * alpha;

      stepSize = acceptedStepMaxNorm;
      if (stepSize < 1e-15) {
        if (verbose >= 1)
          std::cout << "    Iter=" << iter << "; dx = " << stepSize << "; dx too small. Times: " << lineSearchFailedTimes << std::endl;
        stopReason = "tiny_step";
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
  }

  if (solverParam.sst == SST_SUBITERATION_ONE || solverParam.sst == SST_SUBITERATION_STATIC_DAMPING) {
    if (verbose >= 1)
      std::cout << "        Final ||grad||=" << historyGradNormMin << std::endl;

    Eigen::Map<ES::VXd>(x_, energy->getNumDOFs()) = historyx;
  }
  else
    Eigen::Map<ES::VXd>(x_, energy->getNumDOFs()) = x;

  // The status describes the returned state, including exits after an accepted
  // tiny step or the final allowed iteration, so the gradient is re-evaluated
  // there independent of the verbosity level.
  const Eigen::Map<const ES::VXd> xFinal(x_, energy->getNumDOFs());
  grad.setZero();
  energy->gradient(xFinal, grad);
  filterVector(grad);
  const bool finalStateFinite = xFinal.allFinite() && grad.allFinite();
  const double finalGradientNorm = finalStateFinite ? grad.cwiseAbs().maxCoeff() : std::numeric_limits<double>::quiet_NaN();
  const bool converged = finalStateFinite && finalGradientNorm < epsilon;

  int status = SOLVE_NOT_CONVERGED;
  if (numericalFailure || !finalStateFinite)
    status = SOLVE_NUMERICAL_FAILURE;
  else if (converged)
    status = SOLVE_CONVERGED;

  lastStopReason = stopReason;
  lastGradientNorm = finalGradientNorm;

  hclock::time_point t2 = hclock::now();
  double timeCost = dura(t1, t2);

  if (auto logger = Logging::lgr(); logger) {
    if (verbose >= 1) {
      SPDLOG_LOGGER_INFO(logger,
        "Newton result: stop={}, iteration={}, gradient_inf={:.17g}, tolerance={:.17g}, converged={}, status={}",
        stopReason, iter, finalGradientNorm, epsilon, converged, solveStatusName(status));
    }
    SPDLOG_LOGGER_INFO(logger, "Newton solve time: {}", timeCost);
  }
  return status;
}

void NewtonSolver::filterVector(ES::VXd &v)
{
  for (int dof : fixedDOFs)
    v[dof] = 0;
}

//
// OPK_DROP(lineSearchHandle);
