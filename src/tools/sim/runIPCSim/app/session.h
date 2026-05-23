#pragma once

#include "EigenSupport.h"
#include "app/config.h"
#include "app/output.h"
#include "setup/setup.h"

#include <memory>

namespace pgo::Simulation
{
class ImplicitBackwardEulerTimeIntegrator;
}

namespace pgo::RunIPCSim
{
struct RunIPCSimSession
{
  std::shared_ptr<pgo::Simulation::ImplicitBackwardEulerTimeIntegrator> integrator;
  EigenSupport::VXd u;
  EigenSupport::VXd uvel;
  EigenSupport::VXd uacc;
  EigenSupport::VXd usurf;
  EigenSupport::VXd gravityForce;
  EigenSupport::VXd fext;
  int frameStart = -1;
};

RunIPCSimSession createRunIPCSimSession(const RunIPCSimRuntimeConfig &runtimeConfig,
  const IpcSimulationContext &context);
void restoreRestartStateIfRequested(const RunIPCSimRuntimeConfig &runtimeConfig,
  const RunIPCSimOutput &output, RunIPCSimSession &session);
}  // namespace pgo::RunIPCSim
