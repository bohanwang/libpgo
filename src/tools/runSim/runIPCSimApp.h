#pragma once

#include "runIPCSimSetup.h"

#include <filesystem>

namespace pgo
{
class ConfigFileJSON;
}

namespace pgo::RunIPCSim
{
struct RunIPCSimOptions
{
  bool enableCliLog = false;
};

bool openRunIPCSimConfig(const std::filesystem::path &configPath, pgo::ConfigFileJSON &config);
IpcSimulationContext buildIpcSimulation(const pgo::ConfigFileJSON &config);
int runFromConfig(const std::filesystem::path &configPath, const RunIPCSimOptions &options);
}  // namespace pgo::RunIPCSim
