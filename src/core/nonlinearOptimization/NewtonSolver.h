#pragma once

#include "potentialEnergy.h"
#include "solverResult.h"

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
    LineSearchMethod lsm = LSM_BACKTRACK;
    int stopAfterIncrease = 1;
    int addDamping = 0;
  };

  NewtonSolver(const double *x, SolverParam sp, PotentialEnergy_const_p energy_,
    const std::vector<int> &fixedDOFs, const double *fixedValues_ = nullptr);

  void setFixedDOFs(const std::vector<int> &fixedDOFs, const double *fixedValues);
  SolverResult solve(double *x, int numIter, double epsilon, int verbose);

  using StepFunc = std::function<void(const EigenSupport::VXd &, int)>;
  void setStepFunc(StepFunc func) { stepFunc = func; }

  const EigenSupport::VXd &getx() const { return x; }
  const SolveDiagnostics &getSolveDiagnostics() const { return solveDiagnostics; }

protected:
  struct IterationState
  {
    int iter = 0;
    double energy = 0.0;
    double gradMaxNorm = 0.0;
    double gradNorm = 0.0;
    double lambda0 = 1.0;
    double relThreshold = 0.0;
    bool absConverged = false;
    bool relConverged = false;
    bool nonFiniteEnergy = false;
    bool nonFiniteGradient = false;
  };

  struct StepAcceptance
  {
    double feasibleAlpha = 1.0;
    double lineSearchAlpha = 1.0;
    double effectiveAlpha = 1.0;
    double acceptedEnergy = 0.0;
    double acceptedStepMaxNorm = 0.0;

    enum class NonFiniteReason
    {
      None,
      FeasibleAlpha,
      TrialEnergy,
      LineSearchResult
    };

    NonFiniteReason nonFiniteReason = NonFiniteReason::None;

    bool nonFinite() const { return nonFiniteReason != NonFiniteReason::None; }
  };

  // Mutable per-solve state shared between solve() and the step strategy.
  struct SolveContext
  {
    IterationState state;
    double *xOut = nullptr;
    double lambda0 = 1.0;
    double lambdaScale = 1.0;
    double gradMaxNormLast = 0.0;
    bool hasInitialGradNorm = false;
    double epsilon = 0.0;
    int verbose = 0;
    int printGap = 10;
    int iter = 0;
    SolveStatus status = SolveStatus::MaxIterations;
    int completedIterations = 0;
  };

  // Per-subiteration step strategy. Each concrete strategy owns one mode's full
  // lifecycle (begin -> beforeStep -> takeStep -> finalize). begin/beforeStep
  // default to no-ops; the default finalize writes the working x to xOut.
  // beforeStep/takeStep return true when the iteration loop should stop.
  class StepStrategy
  {
  public:
    explicit StepStrategy(NewtonSolver &solver): s(solver) {}
    virtual ~StepStrategy() = default;
    virtual void begin(SolveContext &ctx) {}
    virtual bool beforeStep(SolveContext &ctx) { return false; }
    virtual bool takeStep(SolveContext &ctx) = 0;
    virtual void finalize(SolveContext &ctx);

  protected:
    NewtonSolver &s;
  };
  class LineSearchStrategy;
  class HistoryTrackingStrategy;
  class SubiterationOneStrategy;
  class StaticDampingStrategy;

  void filterVector(EigenSupport::VXd &v);
  void applyFixedValues();
  IterationState evaluateCurrentState(int iter, double epsilon, double lambda0, bool hasInitialGradNorm);
  bool isConverged(const IterationState &state) const;
  bool prepareReducedSystem(double lambdaScale, double lambda0);
  void ensureLinearSolver(bool fixedHessianTopology);
  bool solveReducedNewtonDirection(bool fixedHessianTopology);
  bool expandReducedStep();
  StepAcceptance runLineSearchStep(double currentEnergy, int verbose, int printGap, int iter);
  bool looseRelativeConverged(double gradMaxNorm, double lambda0) const;
  double updateDampingScale(double lambdaScale, double gradMaxNorm, double gradMaxNormLast) const;
  SolveStatus resolveFpLimitFallback(SolveStatus failStatus, double gradMaxNorm, double lambda0, double epsilon);
  void makeLinearSolver(const EigenSupport::SpMatD &A);

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
  SolveDiagnostics solveDiagnostics;

  std::unique_ptr<StepStrategy> stepStrategy;
  StepFunc stepFunc;
};
}  // namespace NonlinearOptimization
}  // namespace pgo
