#pragma once

namespace pgo::RunIPCSim
{
struct IpcSimulationContext;
struct RunIPCSimRuntimeConfig;
class RunIPCSimOutput;

void runIPCSimStaticSolve(
  const RunIPCSimRuntimeConfig &runtimeConfig,
  IpcSimulationContext &context,
  const RunIPCSimOutput &output);
}  // namespace pgo::RunIPCSim
