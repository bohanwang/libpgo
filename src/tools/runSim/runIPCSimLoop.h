#pragma once

#include "runIPCSimConfig.h"
#include "runIPCSimOutput.h"
#include "runIPCSimSession.h"
#include "runIPCSimSetup.h"

namespace pgo::RunIPCSim
{
void runIPCSimLoop(const RunIPCSimRuntimeConfig &runtimeConfig,
  IpcSimulationContext &context,
  RunIPCSimSession &session,
  const RunIPCSimOutput &output);
}  // namespace pgo::RunIPCSim
