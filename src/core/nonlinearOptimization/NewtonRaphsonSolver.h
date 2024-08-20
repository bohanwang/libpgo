#pragma once

#include "potentialEnergy.h"

#if defined(PGO_HAS_MKL)
#include "EigenMKLPardisoSupport.h"
#endif

#include <cfloat>
#include <memory>

namespace pgo
{
namespace NonlinearOptimization
{

class LineSearchHandle;

class NewtonRaphsonSolver
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
  };

  NewtonRaphsonSolver(const double *x, SolverParam sp, PotentialEnergy_const_p energy_,
    const std::vector<int> &fixedDOFs, const double *fixedValues_ = nullptr);

  void setFixedDOFs(const std::vector<int> &fixedDOFs, const double *fixedValues);
  int solve(double *x, int numIter, double epsilon, int verbose);

protected:
  void filterVector(EigenSupport::VXd &v);

  PotentialEnergy_const_p energy;
  SolverParam solverParam;
  std::shared_ptr<LineSearchHandle> lineSearchHandle;

  EigenSupport::VXd x, grad, deltax, deltaxSmall, lineSearchx;
  EigenSupport::SpMatD sysFull, A11, A12;
  EigenSupport::SpMatI A11Mapping, A12Mapping;

  #if defined(PGO_HAS_MKL)
  std::shared_ptr<EigenSupport::EigenMKLPardisoSupport> solver;
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
};
}  // namespace NonlinearOptimization
}  // namespace pgo