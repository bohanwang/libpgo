#include "app/loop.h"

#include "implicitBackwardEulerTimeIntegrator.h"
#include "multiVertexPullingSoftConstraints.h"

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

    context.contactBackend->beginFrame(framei, runtimeConfig, context, session);
    context.contactBackend->addForces(framei, runtimeConfig, context, session);
    if (context.surfacePressureForceEnabled) {
      const double ramp = std::min(1.0, static_cast<double>(framei + 1) / static_cast<double>(context.surfacePressureRampSteps));
      session.fext.noalias() = session.gravityForce + ramp * context.surfacePressureSimulationForce;
      session.integrator->setExternalForce(session.fext.data());
    }
    session.integrator->setqState(session.u, session.uvel, session.uacc);
    session.integrator->doTimestep(1, 3, 1);
    executedStep = true;
    session.integrator->getq(session.u);
    session.integrator->getqvel(session.uvel);
    session.integrator->getqacc(session.uacc);
    context.contactBackend->afterStep(framei, runtimeConfig, context, session);
    context.contactBackend->logSummary(context, session);

    const bool dumpDeformThisFrame = runtimeConfig.dumpDeformEveryFrame || (framei % runtimeConfig.frameGap == 0);
    const bool dumpSurfaceThisFrame = (framei % runtimeConfig.frameGap == 0);
    output.writeStateAndSurfaceFrame(framei, framei / runtimeConfig.frameGap, context,
      session.u, session.uvel, session.uacc, runtimeConfig.scale,
      dumpDeformThisFrame, dumpSurfaceThisFrame);

    if (runtimeConfig.outputVonMises)
      output.writeVonMisesStressJson(framei, runtimeConfig.timestep, context, session.u);
  }

  if (!executedStep)
    context.contactBackend->logSummary(context, session);
}
}  // namespace pgo::RunIPCSim
