#include "solver/staticSolve.h"

#include "EigenSupport.h"
#include "deformationModelEnergy.h"
#include "linearPotentialEnergy.h"
#include "multiVertexPullingSoftConstraints.h"
#include "potentialEnergies.h"
#include "app/config.h"
#include "app/output.h"
#include "setup/setup.h"
#include "solver/solverRunner.h"

#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace pgo::RunIPCSim
{
namespace ES = pgo::EigenSupport;

namespace
{
void setStaticPullingTargets(IpcSimulationContext &context)
{
  for (std::size_t pi = 0; pi < context.pullingEnergies.size(); ++pi) {
    context.pullingEnergies[pi]->setTargetPos(context.pullingTargets[pi].data());
    std::cout << "Static attachment " << pi << " target: "
              << context.pullingTargets[pi].transpose().head(3) << std::endl;
  }
}

ES::VXd buildStaticExternalForce(const RunIPCSimRuntimeConfig &runtimeConfig,
  const IpcSimulationContext &context)
{
  const int n3 = static_cast<int>(context.simulationRestPosition.size());
  ES::VXd g(n3);
  for (int vi = 0; vi < n3 / 3; ++vi)
    g.segment<3>(vi * 3) = runtimeConfig.gravity;

  ES::VXd fext(n3);
  ES::mv(context.M, g, fext);
  if (context.surfacePressureForceEnabled)
    fext.noalias() += context.surfacePressureSimulationForce;
  return fext;
}
}  // namespace

void runIPCSimStaticSolve(
  const RunIPCSimRuntimeConfig &runtimeConfig,
  IpcSimulationContext &context,
  const RunIPCSimOutput &output)
{
  const int n3 = static_cast<int>(context.simulationRestPosition.size());
  if (n3 <= 0)
    throw std::runtime_error("runIPCSim static solve received an empty simulation state.");

  setStaticPullingTargets(context);

  const ES::VXd staticForce = buildStaticExternalForce(runtimeConfig, context);
  auto externalForcesEnergy =
    std::make_shared<PredefinedPotentialEnergies::LinearPotentialEnergy>(staticForce);

  auto energyAll = std::make_shared<NonlinearOptimization::PotentialEnergies>(n3);
  energyAll->addPotentialEnergy(context.elasticEnergy, 1.0);
  for (const auto &pullingEnergy : context.pullingEnergies)
    energyAll->addPotentialEnergy(pullingEnergy, 1.0);
  energyAll->addPotentialEnergy(externalForcesEnergy, -1.0);
  context.contactBackend->addStaticEnergies(runtimeConfig, context, *energyAll);
  energyAll->init();

  ES::VXd u = ES::VXd::Zero(n3);
  energyAll->printEnergy(u);

  const StaticSolveResult staticResult = solveStaticEnergyStrict(
    energyAll, n3, runtimeConfig.solverMaxIter, runtimeConfig.solverEps, 2);
  u = staticResult.u;

  const ES::VXd uvel = ES::VXd::Zero(n3);
  const ES::VXd uacc = ES::VXd::Zero(n3);
  output.writeStateAndSurfaceFrame(
    0, 0, context, u, uvel, uacc, runtimeConfig.scale, true, true);

  if (runtimeConfig.outputVonMises)
    output.writeVonMisesStressJson(0, runtimeConfig.timestep, context, u);
}
}  // namespace pgo::RunIPCSim
