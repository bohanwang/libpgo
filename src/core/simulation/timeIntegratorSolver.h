#pragma once

#include "timeIntegratorOptions.h"

#include "EigenSupport.h"

#include <memory>

namespace pgo
{
namespace NonlinearOptimization
{
class PotentialEnergy;
class ConstraintFunctions;
}  // namespace NonlinearOptimization

namespace Simulation
{
class TimeIntegratorSolverData;

class TimeIntegratorSolver
{
public:
  TimeIntegratorSolver();

  int solve(bool needRenew, EigenSupport::VXd &x, EigenSupport::VXd &g, EigenSupport::VXd &lambda,
    const EigenSupport::VXd &xlow, const EigenSupport::VXd &xhi,
    const EigenSupport::VXd &clow, const EigenSupport::VXd &chi,
    std::shared_ptr<const NonlinearOptimization::PotentialEnergy> energy,
    std::shared_ptr<const NonlinearOptimization::ConstraintFunctions> constraints,
    int niter, double eps, int verbose, const char *solverConfigFilename, TimeIntegratorSolverOption op);

  static int solveDirect(EigenSupport::VXd &x, EigenSupport::VXd &g, EigenSupport::VXd &lambda,
    const EigenSupport::VXd &xlow, const EigenSupport::VXd &xhi,
    const EigenSupport::VXd &clow, const EigenSupport::VXd &chi,
    std::shared_ptr<const NonlinearOptimization::PotentialEnergy> energy,
    std::shared_ptr<const NonlinearOptimization::ConstraintFunctions> constraints,
    int niter, double eps, int verbose, const char *solverConfigFilename,
    TimeIntegratorSolverOption solverOption);

protected:
  std::shared_ptr<TimeIntegratorSolverData> da;
};

}  // namespace OptimizationBasedIntegrator
}  // namespace VegaFEM
