#pragma once

#include "runIPCSimContactBackend.h"
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
  ContactBackendKind contactBackendKind = ContactBackendKind::Ipc;
};

bool openRunIPCSimConfig(const std::filesystem::path &configPath, pgo::ConfigFileJSON &config);
IpcSimulationContext buildIpcSimulation(const pgo::ConfigFileJSON &config);
IpcSimulationContext buildRunIPCSimSimulation(const pgo::ConfigFileJSON &config, const RunIPCSimOptions &options);
int runFromConfig(const std::filesystem::path &configPath, const RunIPCSimOptions &options);
}  // namespace pgo::RunIPCSim
