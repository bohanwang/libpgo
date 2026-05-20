#pragma once

#include "runIPCSimApp.h"

#include <filesystem>

namespace pgo::RunIPCSim
{
struct RunIPCSimCliOptions
{
  std::filesystem::path configPath;
  RunIPCSimOptions runOptions;
};

RunIPCSimCliOptions parseRunIPCSimCli(int argc, char *argv[]);
int runCli(int argc, char *argv[]);
}  // namespace pgo::RunIPCSim
