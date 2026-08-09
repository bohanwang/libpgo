#pragma once

#include <filesystem>

namespace pgo::SimulationRunner
{
// These are the two existing simulation implementations extracted from the
// command-line tools. They intentionally remain separate.
int runSampledSimulationFromConfig(
  const std::filesystem::path &configFilename, bool enableCliLog = false);
int runIPCSimulationFromConfig(
  const std::filesystem::path &configFilename, bool enableCliLog = false);

// Convenience dispatcher for the config-driven C and Python APIs.
int runSimulationFromConfig(const std::filesystem::path &configFilename);
}  // namespace pgo::SimulationRunner
