#include "runIPCSimLoop.h"

#include "embeddedSurfaceIPCPotentialEnergy.h"
#include "embeddedSurfaceFloorPotentialEnergy.h"
#include "implicitBackwardEulerTimeIntegrator.h"
#include "multiVertexPullingSoftConstraints.h"
#include "runIPCSimLogging.h"

#include <algorithm>
#include <iostream>

namespace pgo::RunIPCSim
{
namespace ES = pgo::EigenSupport;

void runIPCSimLoop(const RunIPCSimRuntimeConfig &runtimeConfig,
  IpcSimulationContext &context,
  RunIPCSimSession &session,
  const RunIPCSimOutput &output)
{
  const double ratioDenom = runtimeConfig.numSimSteps > 1 ? static_cast<double>(runtimeConfig.numSimSteps - 1) : 1.0;

  bool executedStep = false;
  for (int framei = session.frameStart + 1; framei < runtimeConfig.numSimSteps; ++framei) {
    session.integrator->clearGeneralImplicitForceModel();

    const double ratio = runtimeConfig.numSimSteps > 1 ? static_cast<double>(framei) / ratioDenom : 1.0;
    for (std::size_t pi = 0; pi < context.pullingEnergies.size(); ++pi) {
      const ES::VXd curTgt = context.pullingTargetRests[pi] * (1.0 - ratio) + context.pullingTargets[pi] * ratio;
      context.pullingEnergies[pi]->setTargetPos(curTgt.data());
      std::cout << "Frame " << framei << ", attachment " << pi << " target: " << curTgt.transpose().head(3) << std::endl;
    }

    session.integrator->addGeneralImplicitForceModel(context.collisionHandler, 0, 0);
    for (std::size_t fi = 0; fi < context.floorPotentialEnergies.size(); ++fi) {
      context.floorPotentialEnergies[fi]->setFloorHeight(floorHeightAtFrame(context.floorMotionStates[fi], framei));
    }
    for (const auto &forceModel : context.extraGeneralImplicitForceModels)
      session.integrator->addGeneralImplicitForceModel(forceModel, 0, 0);
    if (context.surfacePressureForceEnabled) {
      const double ramp = std::min(1.0, static_cast<double>(framei + 1) / static_cast<double>(context.surfacePressureRampSteps));
      session.fext.noalias() = session.gravityForce + ramp * context.surfacePressureSimulationForce;
      session.integrator->setExternalForce(session.fext.data());
    }
    session.integrator->setqState(session.u, session.uvel, session.uacc);
    const double tCurr = static_cast<double>(framei) * runtimeConfig.timestep;
    context.collisionHandler->setObstacleTime(tCurr + runtimeConfig.timestep);
    session.integrator->doTimestep(1, 3, 1);
    executedStep = true;
    session.integrator->getq(session.u);
    session.integrator->getqvel(session.uvel);
    session.integrator->getqacc(session.uacc);
    logRunIPCSimMaxStepSummary(context.elasticEnergy, context.collisionHandler, session.integrator);

    const bool dumpDeformThisFrame = runtimeConfig.dumpDeformEveryFrame || (framei % runtimeConfig.frameGap == 0);
    if (dumpDeformThisFrame)
      output.writeState(framei, session.u, session.uvel, session.uacc);

    if (runtimeConfig.outputVonMises)
      output.writeVonMisesStressJson(framei, runtimeConfig.timestep, context, session.u);

    if (framei % runtimeConfig.frameGap == 0) {
      pgo::Mesh::TriMeshGeo mesh = context.surfaceMesh;
      ES::mv(context.surfaceFromSimulationDispMap, session.u, session.usurf);
      const ES::VXd psurf = context.surfaceRestPositions + session.usurf;
      for (int vi = 0; vi < mesh.numVertices(); ++vi)
        mesh.pos(vi) = psurf.segment<3>(vi * 3) / runtimeConfig.scale;
      output.writeSurface(framei / runtimeConfig.frameGap, mesh);
    }
  }

  if (!executedStep)
    logRunIPCSimMaxStepSummary(context.elasticEnergy, context.collisionHandler, session.integrator);
}
}  // namespace pgo::RunIPCSim
