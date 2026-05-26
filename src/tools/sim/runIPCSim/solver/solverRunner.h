#pragma once

#include "EigenSupport.h"
#include "potentialEnergy.h"
#include "solverResult.h"

namespace pgo::RunIPCSim
{

struct StaticSolveResult
{
  NonlinearOptimization::SolverResult solver;
  EigenSupport::VXd u;
};

StaticSolveResult solveStaticEnergyStrict(
  NonlinearOptimization::PotentialEnergy_const_p energy,
  int numDofs,
  int maxIterations,
  double epsilon,
  int verbose);

}  // namespace pgo::RunIPCSim
