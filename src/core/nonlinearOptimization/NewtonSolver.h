#pragma once

#include "potentialEnergy.h"

#if defined(PGO_HAS_MKL) && !defined(PGO_HAS_ORIG_PARDISO)
#  include "EigenMKLPardisoSupport.h"
#elif defined(PGO_HAS_ORIG_PARDISO)
#  include "EigenOrigPardisoSupport.h"
#endif

#include <cfloat>
#include <memory>

namespace pgo
{
namespace NonlinearOptimization
{

class LineSearchHandle;

class NewtonSolver
{
public:
  enum SolverSubiterationType
  {
    SST_SUBITERATION_LINE_SEARCH,
    SST_SUBITERATION_STATIC_DAMPING,
    SST_SUBITERATION_ONE
  };

  enum LineSearchMethod
  {
    LSM_GOLDEN,
    LSM_BRENTS,
    LSM_BACKTRACK,
    LSM_SIMPLE,
  };

  struct SolverParam
  {
    double alpha = 0.5;
    SolverSubiterationType sst = SST_SUBITERATION_LINE_SEARCH;
    LineSearchMethod lsm = LSM_SIMPLE;
    int stopAfterIncrease = 1;
    int addDamping = 0;
  };

  // Return codes of solve(). They describe the returned state, not the stop
  // reason: a solve that exhausts its iterations but ends below the gradient
  // tolerance is converged, and a solve that stops after a failed line search
  // above the tolerance is not.
  enum SolveStatus
  {
    SOLVE_CONVERGED = 0,          // gradient infinity norm below epsilon at the returned state
    SOLVE_NOT_CONVERGED = 1,      // stopped early (iteration limit, tiny step, failed line search) with a finite state
    SOLVE_NUMERICAL_FAILURE = 2,  // non-finite energy, gradient, or step; the last finite iterate is returned
  };
  static const char *solveStatusName(int status);

  NewtonSolver(const double *x, SolverParam sp, PotentialEnergy_const_p energy_,
    const std::vector<int> &fixedDOFs, const double *fixedValues_ = nullptr);

  void setFixedDOFs(const std::vector<int> &fixedDOFs, const double *fixedValues);
  int solve(double *x, int numIter, double epsilon, int verbose);

  using StepFunc = std::function<void(const EigenSupport::VXd &, int)>;
  void setStepFunc(StepFunc func) { stepFunc = func; }

  const EigenSupport::VXd &getx() const { return x; }
  // Why the last solve() left its iteration loop, e.g. "gradient_tolerance".
  const char *getLastStopReason() const { return lastStopReason; }
  double getLastGradientNorm() const { return lastGradientNorm; }

protected:
  void filterVector(EigenSupport::VXd &v);

  PotentialEnergy_const_p energy;
  SolverParam solverParam;
  std::shared_ptr<LineSearchHandle> lineSearchHandle;

  EigenSupport::VXd x, grad, deltax, deltaxSmall, lineSearchx;
  EigenSupport::SpMatD sysFull, A11, A12;
  EigenSupport::SpMatI A11Mapping, A12Mapping;

#if defined(PGO_HAS_MKL) && !defined(PGO_HAS_ORIG_PARDISO)
  std::shared_ptr<EigenSupport::EigenMKLPardisoSupport> solver;
#elif defined(PGO_HAS_ORIG_PARDISO)
  std::shared_ptr<EigenSupport::EigenOrigPardisoSupport> solver;
#else
  std::shared_ptr<EigenSupport::SymSolver> solver;
#endif

  std::vector<int> allDOFs, fixedDOFs;
  std::vector<int> rhss2b, rhsb2s;
  EigenSupport::VXd rhs;
  EigenSupport::VXd fixedValues;
  int n3;

  EigenSupport::VXd historyx;
  double historyGradNormMin;

  StepFunc stepFunc;
  const char *lastStopReason = "not_run";
  double lastGradientNorm = -1.0;
};
}  // namespace NonlinearOptimization
}  // namespace pgo
