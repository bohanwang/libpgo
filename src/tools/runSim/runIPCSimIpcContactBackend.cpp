#include "runIPCSimContactBackend.h"

#include "embeddedSurfaceFloorPotentialEnergy.h"
#include "embeddedSurfaceIPCPotentialEnergy.h"
#include "implicitBackwardEulerTimeIntegrator.h"
#include "potentialEnergies.h"
#include "runIPCSimConfig.h"
#include "runIPCSimLogging.h"
#include "runIPCSimSession.h"
#include "runIPCSimSetup.h"

namespace pgo::RunIPCSim
{
class IpcContactBackend final : public RunIPCSimContactBackend
{
public:
  ContactBackendKind kind() const override { return ContactBackendKind::Ipc; }
  std::string description() const override { return "ipc"; }

  void initializeAfterRestart(const RunIPCSimRuntimeConfig &, IpcSimulationContext &, RunIPCSimSession &) override {}

  void beginFrame(int frame, const RunIPCSimRuntimeConfig &, IpcSimulationContext &context, RunIPCSimSession &) override
  {
    for (std::size_t fi = 0; fi < context.floorPotentialEnergies.size(); ++fi)
      context.floorPotentialEnergies[fi]->setFloorHeight(floorHeightAtFrame(context.floorMotionStates[fi], frame));
  }

  void addForces(int frame, const RunIPCSimRuntimeConfig &runtimeConfig,
    IpcSimulationContext &context, RunIPCSimSession &session) override
  {
    const double tCurr = static_cast<double>(frame) * runtimeConfig.timestep;
    context.collisionHandler->setObstacleTime(tCurr + runtimeConfig.timestep);
    session.integrator->addGeneralImplicitForceModel(context.collisionHandler, 0, 0);
    for (const auto &forceModel : context.extraGeneralImplicitForceModels)
      session.integrator->addGeneralImplicitForceModel(forceModel, 0, 0);
  }

  void afterStep(int, const RunIPCSimRuntimeConfig &, IpcSimulationContext &, RunIPCSimSession &) override {}

  void addStaticEnergies(const RunIPCSimRuntimeConfig &runtimeConfig,
    IpcSimulationContext &context, NonlinearOptimization::PotentialEnergies &energyAll) override
  {
    const int finalFrame = runtimeConfig.numSimSteps > 0 ? runtimeConfig.numSimSteps - 1 : 0;
    for (std::size_t fi = 0; fi < context.floorPotentialEnergies.size(); ++fi)
      context.floorPotentialEnergies[fi]->setFloorHeight(floorHeightAtFrame(context.floorMotionStates[fi], finalFrame));

    const double staticObstacleTime = runtimeConfig.timestep * static_cast<double>(finalFrame);
    context.collisionHandler->setObstacleTime(staticObstacleTime);
    energyAll.addPotentialEnergy(context.collisionHandler, 1.0);

    for (const auto &forceModel : context.extraGeneralImplicitForceModels)
      energyAll.addPotentialEnergy(forceModel, 1.0);
  }

  void logSummary(const IpcSimulationContext &context, const RunIPCSimSession &session) const override
  {
    logRunIPCSimMaxStepSummary(context.elasticEnergy, context.collisionHandler, session.integrator);
  }
};

std::shared_ptr<RunIPCSimContactBackend> makeIpcContactBackend()
{
  return std::make_shared<IpcContactBackend>();
}
}  // namespace pgo::RunIPCSim
