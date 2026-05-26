#include "solver/solverRunner.h"

#include "NewtonSolver.h"

#include <stdexcept>
#include <string>
#include <vector>

namespace pgo::RunIPCSim
{

StaticSolveResult solveStaticEnergyStrict(
  NonlinearOptimization::PotentialEnergy_const_p energy,
  int numDofs,
  int maxIterations,
  double epsilon,
  int verbose)
{
  NonlinearOptimization::NewtonSolver::SolverParam solverParam;
  solverParam.addDamping = 1;

  StaticSolveResult result;
  result.u = EigenSupport::VXd::Zero(numDofs);

  NonlinearOptimization::NewtonSolver solver(
    result.u.data(), solverParam, energy, std::vector<int>(), nullptr);
  result.solver = solver.solve(result.u.data(), maxIterations, epsilon, verbose);
  if (result.solver.status != NonlinearOptimization::SolveStatus::Converged) {
    throw std::runtime_error("runIPCSim static solve failed to converge; " +
      NonlinearOptimization::formatSolverResultSummary(result.solver));
  }

  return result;
}

}  // namespace pgo::RunIPCSim
