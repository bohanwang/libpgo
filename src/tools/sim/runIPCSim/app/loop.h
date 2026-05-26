#pragma once

#include "app/config.h"
#include "app/output.h"
#include "app/session.h"
#include "setup/setup.h"

namespace pgo::RunIPCSim
{
void runIPCSimLoop(const RunIPCSimRuntimeConfig &runtimeConfig,
  IpcSimulationContext &context,
  RunIPCSimSession &session,
  const RunIPCSimOutput &output);
}  // namespace pgo::RunIPCSim
