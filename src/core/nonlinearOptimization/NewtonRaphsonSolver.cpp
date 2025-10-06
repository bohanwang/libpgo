#include "NewtonRaphsonSolver.h"

#include "EigenSupport.h"
#include "lineSearch.h"
#include "pgoLogging.h"

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

NewtonRaphsonSolver::NewtonRaphsonSolver(const double *x_, SolverParam sp, PotentialEnergy_const_p energy_, const std::vector<int> &fixedDOFs_, const double *fixedValues_):
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

void NewtonRaphsonSolver::setFixedDOFs(const std::vector<int> &fixedDOFs_, const double *fixedValues_)
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

    // sparse matrix
    energy->createHessian(sysFull);
    energy->hessian(x, sysFull);

    ES::removeRowsCols(sysFull, fixedDOFs, A11);
    ES::removeRowsCols(sysFull, A11, fixedDOFs, A11Mapping);

    // std::vector<int> flexibleDOFs;
    // std::set_difference(allDOFs.begin(), allDOFs.end(), fixedDOFs.begin(), fixedDOFs.end(), std::back_inserter(flexibleDOFs));
    // ES::SelectRowsCols(sysFull, flexibleDOFs, fixedDOFs, A12);
    // ES::Big2Small(sysFull, A12, flexibleDOFs, fixedDOFs, A12Mapping, 0);

#if defined(PGO_HAS_MKL)
    solver = std::make_shared<ES::EigenMKLPardisoSupport>(A11, ES::EigenMKLPardisoSupport::MatrixType::REAL_SYM_INDEFINITE,
      ES::EigenMKLPardisoSupport::ReorderingType::NESTED_DISSECTION, 0, 0, 2, 0, 0, 0);
    solver->analyze(A11);
#else
    solver = std::make_shared<EigenSupport::SymSolver>();
    solver->analyzePattern(A11);
#endif
  }

  fixedValues = ES::Mp<const ES::VXd>(fixedValues_, fixedDOFs_.size());
}

int NewtonRaphsonSolver::solve(double *x_, int numIter, double epsilon, int verbose)
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
  for (; iter < numIter; iter++) {
    if (verbose >= 2 && iter % printGap == 0)
      std::cout << "    Iter=" << iter << std::endl;

    double eng = energy->func(x);

    // we solve f(x_i) + K(x_i) deltax = 0
    memset(grad.data(), 0, sizeof(double) * grad.size());
    energy->gradient(x, grad);
    filterVector(grad);

    double gradNorm = grad.cwiseAbs().maxCoeff();
    if (verbose >= 2 && iter % printGap == 0)
      std::cout << "        E= " << eng << "; ||grad||_max=" << grad.cwiseAbs().maxCoeff() << "; ||x||=" << x.norm() << "; ||grad||=" << gradNorm << std::endl;

    if (solverParam.sst == SST_SUBITERATION_ONE || solverParam.sst == SST_SUBITERATION_STATIC_DAMPING) {
      if (gradNorm < historyGradNormMin) {
        historyx.noalias() = x;
        historyGradNormMin = gradNorm;
      }
      else {
        if (solverParam.stopAfterIncrease)
          break;
      }
    }

    if (iter) {
      if (gradNorm < epsilon) {
        if (verbose >= 1) {
          std::cout << "    Iter=" << iter << "; ||grad|| < eps. Done." << std::endl;
        }

        break;
      }
    }

    memset(sysFull.valuePtr(), 0, sizeof(double) * sysFull.nonZeros());
    energy->hessian(x, sysFull);

    // remove column rows
    ES::transferBigToSmall(sysFull, A11, A11Mapping, 1);
    // ES::transferBigToSmall(sysFull, A12, A12Mapping, 1);
    ES::transferBigToSmall(grad, rhs, rhsb2s, 1);

    // std::cout << "@@:\n"
    //   << x.norm() << '\n'
    //   << sysFull.coeff(0, 0) << '\n'
    //   << sysFull.coeff(1, 0) << '\n'
    //   << sysFull.coeff(1, 1) << '\n'
    //   << sysFull.coeff(2, 0) << '\n'
    //   << sysFull.coeff(2, 1) << '\n'
    //   << sysFull.coeff(2, 2) << std::endl;

    // rhs = -grad - A12 * fixedvalue
    // ES::mv(A12, fixedValues, rhs, -1.0, -1.0, 0);
    rhs *= -1.0;
#if defined(PGO_HAS_MKL)
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

    if (deltax.hasNaN()) {
      std::cerr << "        Encounter weird number.\n"
                << std::endl;
      abort();
    }

    // if (iter == 0)
    //   error0 = deltax.norm();
    // else {
    //   if (deltax.norm() < epsilon) {
    //     if (verbose >= 1)
    //       std::cout << "    dir too small." << std::endl;

    //     break;
    //   }
    // }

    if (verbose >= 2 && iter % printGap == 0)
      std::cout << "        ||deltax||_max=" << deltax.cwiseAbs().maxCoeff() << std::endl;

    // x += alpha delta x ?
    if (solverParam.sst == SST_SUBITERATION_LINE_SEARCH) {
      double alpha = 1;
      double stepSize = 0;

      alpha = alphaTestFunc ? alphaTestFunc(x, deltax) : 1.0;

      lineSearchx.noalias() = x + deltax;
      double eng1 = energy->func(lineSearchx);
      int maxIter = 10;
      if (eng1 < eng) {
        maxIter = 3;
      }

      if (solverParam.lsm == LSM_GOLDEN) {
        lineSearchHandle->nativeLineSearch->setMaxIterations(maxIter);
        LineSearch::Result ret = lineSearchHandle->nativeLineSearch->golden(x.data(), deltax.data(), eng);
        alpha = ret.alpha;

        if (verbose >= 2 && iter % printGap == 0) {
          std::cout << "        f=" << ret.f << ";alpha=" << alpha << std::endl;
        }
      }
      else if (solverParam.lsm == LSM_BRENTS) {
        lineSearchHandle->nativeLineSearch->setMaxIterations(maxIter);
        LineSearch::Result ret = lineSearchHandle->nativeLineSearch->BrentsMethod(x.data(), deltax.data(), eng);
        alpha = ret.alpha;

        if (verbose >= 2 && iter % printGap == 0) {
          std::cout << "        f=" << ret.f << ";alpha=" << alpha << std::endl;
        }
      }
      else if (solverParam.lsm == LSM_BACKTRACK) {
        lineSearchHandle->nativeLineSearch->setMaxIterations(maxIter);
        LineSearch::Result ret = lineSearchHandle->nativeLineSearch->backtracking(x.data(), deltax.data(), eng, grad.data(), 0.0001, 0.5);
        alpha = ret.alpha;

        if (verbose >= 2 && iter % printGap == 0) {
          std::cout << "        f=" << ret.f << ";alpha=" << alpha << std::endl;
        }
      }
      else if (solverParam.lsm == LSM_SIMPLE) {
        double eng1 = eng;
        for (int i = 0; i < 50; i++) {
          lineSearchx.noalias() = x + deltax * alpha;
          eng1 = energy->func(lineSearchx);

          if (eng1 < eng) {
            break;
          }

          alpha *= 0.75;
        }

        if (verbose >= 2 && iter % printGap == 0) {
          std::cout << "        f=" << eng1 << ";alpha=" << alpha << std::endl;
        }
      }

      x += deltax * alpha;

      stepSize = std::abs(alpha * deltax.cwiseAbs().maxCoeff());
      if (stepSize < 1e-12) {
        if (verbose >= 1) {
          std::cout << "    Iter=" << iter << "; dx = " << stepSize << "; no delta x improvement." << std::endl;
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
  }

  if (solverParam.sst == SST_SUBITERATION_ONE || solverParam.sst == SST_SUBITERATION_STATIC_DAMPING) {
    if (verbose >= 1)
      std::cout << "        Final ||grad||=" << historyGradNormMin << std::endl;

    Eigen::Map<ES::VXd>(x_, energy->getNumDOFs()) = historyx;
  }
  else
    Eigen::Map<ES::VXd>(x_, energy->getNumDOFs()) = x;

  hclock::time_point t2 = hclock::now();

  double timeCost = dura(t1, t2);

  SPDLOG_LOGGER_INFO(Logging::lgr(), "Newton solve time: {}", timeCost);
  return 0;
}

void NewtonRaphsonSolver::filterVector(ES::VXd &v)
{
  for (int dof : fixedDOFs)
    v[dof] = 0;
}

//
// OPK_DROP(lineSearchHandle);
