#pragma once

#include <filesystem>

namespace pgo::SimulationRunner
{
// Invocation controls do not change the configured physical time or loading
// schedule. A bounded invocation can later resume from its ordinary checkpoint.
struct IPCSimulationRunOptions
{
  int maxFrames = -1;  // -1: all remaining frames; otherwise strictly positive.
  std::filesystem::path timingOutput;
};

// These are the two existing simulation implementations extracted from the
// command-line tools. They intentionally remain separate.
int runSampledSimulationFromConfig(
  const std::filesystem::path &configFilename, bool enableCliLog = false);
int runIPCSimulationFromConfig(
  const std::filesystem::path &configFilename, bool enableCliLog = false);
int runIPCSimulationFromConfig(
  const std::filesystem::path &configFilename, bool enableCliLog,
  const IPCSimulationRunOptions &options);

// Convenience dispatcher for the config-driven C and Python APIs.
int runSimulationFromConfig(const std::filesystem::path &configFilename);
}  // namespace pgo::SimulationRunner
