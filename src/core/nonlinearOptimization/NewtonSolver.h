#pragma once

#include "potentialEnergy.h"

#include "EigenDef.h"

// Panua PARDISO and Intel MKL PARDISO export the same `pardiso` entry point,
// so a build links one of them: the MKL wrapper is compiled out whenever
// PGO_HAS_ORIG_PARDISO is defined.
#if defined(PGO_HAS_MKL) && !defined(PGO_HAS_ORIG_PARDISO)
#  include "EigenMKLPardisoSupport.h"
#endif
#if defined(PGO_HAS_ORIG_PARDISO)
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

  // Direct solver for the Newton system. LS_DEFAULT resolves to the process
  // defaults (see setLinearSolverDefaults), which prefer Panua PARDISO when
  // the build includes it, then Intel MKL PARDISO, then Eigen's LDLT. A build
  // with Panua PARDISO has no MKL PARDISO option (shared entry point).
  enum LinearSolverType
  {
    LS_DEFAULT,
    LS_PARDISO,      // Panua PARDISO (PGO_HAS_ORIG_PARDISO)
    LS_MKL_PARDISO,  // Intel MKL PARDISO (PGO_HAS_MKL without PGO_HAS_ORIG_PARDISO)
    LS_EIGEN_LDLT,   // EigenSupport::SymSolver
  };
  static const char *linearSolverName(LinearSolverType type);
  static bool isLinearSolverAvailable(LinearSolverType type);
  // Maps LS_DEFAULT to the configured default and throws std::invalid_argument
  // when the requested solver is not part of the build.
  static LinearSolverType resolveLinearSolver(LinearSolverType requested);

  struct LinearSolverDefaults
  {
    LinearSolverType type = LS_DEFAULT;  // LS_DEFAULT: best available
    bool pardisoDeterministic = true;    // Panua PARDISO IPARM(34) mode
    int pardisoThreadLimit = 64;         // Panua PARDISO IPARM(3) for identical results
  };
  static void setLinearSolverDefaults(const LinearSolverDefaults &defaults);
  static LinearSolverDefaults getLinearSolverDefaults();

  struct SolverParam
  {
    double alpha = 0.5;
    SolverSubiterationType sst = SST_SUBITERATION_LINE_SEARCH;
    LineSearchMethod lsm = LSM_SIMPLE;
    int stopAfterIncrease = 1;
    int addDamping = 0;
    LinearSolverType linearSolver = LS_DEFAULT;
    int pardisoDeterministic = -1;  // -1: process default; 0/1 override
    int pardisoThreadLimit = 0;     // 0: process default
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
  LinearSolverType getLinearSolverType() const { return linearSolverType; }
  // Why the last solve() left its iteration loop, e.g. "gradient_tolerance".
  const char *getLastStopReason() const { return lastStopReason; }
  double getLastGradientNorm() const { return lastGradientNorm; }

protected:
  void filterVector(EigenSupport::VXd &v);
  bool hasLinearSolver() const;
  // Identity of the active backend object (tests check that a fixed Hessian
  // topology keeps one symbolic factorization across iterations).
  const void *linearSolverAddress() const;
  void createLinearSolver();
  void solveLinearSystem();

  PotentialEnergy_const_p energy;
  SolverParam solverParam;
  std::shared_ptr<LineSearchHandle> lineSearchHandle;

  EigenSupport::VXd x, grad, deltax, deltaxSmall, lineSearchx;
  EigenSupport::SpMatD sysFull, A11, A12;
  EigenSupport::SpMatI A11Mapping, A12Mapping;

  LinearSolverType linearSolverType = LS_EIGEN_LDLT;
  bool pardisoDeterministic = true;
  int pardisoThreadLimit = 64;
#if defined(PGO_HAS_MKL) && !defined(PGO_HAS_ORIG_PARDISO)
  std::shared_ptr<EigenSupport::EigenMKLPardisoSupport> mklSolver;
#endif
#if defined(PGO_HAS_ORIG_PARDISO)
  std::shared_ptr<EigenSupport::EigenOrigPardisoSupport> pardisoSolver;
#endif
  std::shared_ptr<EigenSupport::SymSolver> eigenSolver;

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
