#include "app/session.h"

#include "deformationModelEnergy.h"
#include "implicitBackwardEulerTimeIntegrator.h"
#include "multiVertexPullingSoftConstraints.h"

#include <iostream>

namespace pgo::RunIPCSim
{
namespace ES = pgo::EigenSupport;

RunIPCSimSession createRunIPCSimSession(const RunIPCSimRuntimeConfig &runtimeConfig,
  const IpcSimulationContext &context)
{
  const int n3 = static_cast<int>(context.simulationRestPosition.size());
  const int surfn3 = static_cast<int>(context.surfaceRestPositions.size());
  const int n = n3 / 3;

  ES::VXd g(n3);
  for (int vi = 0; vi < n; ++vi)
    g.segment<3>(vi * 3) = runtimeConfig.gravity;

  RunIPCSimSession session;
  session.gravityForce = ES::VXd(n3);
  ES::mv(context.M, g, session.gravityForce);
  session.fext = session.gravityForce;

  session.integrator = std::make_shared<pgo::Simulation::ImplicitBackwardEulerTimeIntegrator>(context.M, context.elasticEnergy,
    runtimeConfig.dampingParams[0], runtimeConfig.dampingParams[1], runtimeConfig.timestep,
    runtimeConfig.solverMaxIter, runtimeConfig.solverEps);

  for (auto &pullingEnergy : context.pullingEnergies)
    session.integrator->addImplicitForceModel(pullingEnergy, 0, 0);

  session.integrator->setExternalForce(session.fext.data());

  session.u = ES::VXd::Zero(n3);
  session.uvel = ES::VXd::Zero(n3);
  session.uacc = ES::VXd::Zero(n3);
  session.usurf = ES::VXd::Zero(surfn3);

  for (int i = 0; i < n; ++i)
    session.uvel.segment<3>(i * 3) = runtimeConfig.initialVelocity;

  return session;
}

void restoreRestartStateIfRequested(const RunIPCSimRuntimeConfig &runtimeConfig,
  const RunIPCSimOutput &output, RunIPCSimSession &session)
{
  if (runtimeConfig.restartFromU) {
    session.frameStart = output.loadLatestRestartState(
      runtimeConfig.numSimSteps,
      static_cast<int>(session.u.size()),
      session.u, session.uvel, session.uacc);
    if (session.frameStart >= 0) {
      std::cout << "Restarting from frame " << session.frameStart << std::endl;
      return;
    }

    std::cout << "No restart state found in " << output.directories().states << ". Starting from frame 0." << std::endl;
    return;
  }

  session.frameStart = -1;
  std::cout << "Starting from frame 0." << std::endl;
}
}  // namespace pgo::RunIPCSim
