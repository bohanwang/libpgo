#pragma once

namespace pgo
{
class ConfigFileJSON;

namespace SimulationRunner
{
struct IPCSimulationRunOptions;
int runCudaIPCSimulation(const ConfigFileJSON &config, bool enableCliLog,
  const IPCSimulationRunOptions &options);
}  // namespace SimulationRunner
}  // namespace pgo
