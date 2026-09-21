#pragma once

// Logging support shared by the sampled and IPC config runners.
#include "pgoLogging.h"

#include <filesystem>
#include <string>

namespace pgo
{
class ConfigFileJSON;
}

namespace pgo::RunSim
{
std::filesystem::path deriveDefaultLogPathFromConfig(const std::filesystem::path &configPath);
spdlog::level::level_enum resolveConfiguredLogLevel(const ConfigFileJSON &config);

// Interprets a NonlinearOptimization::NewtonSolver::SolveStatus for a runner.
// A numerical failure is fatal and returns true. A solve that stopped short of
// its gradient tolerance (iteration limit, tiny step, failed line search) logs
// a warning and returns false so the run continues with the returned state.
bool newtonStatusIsFatal(int status, const std::string &context);

class ScopedRunSimCliLogRedirect
{
public:
  explicit ScopedRunSimCliLogRedirect(const std::string &logFilename);
  ~ScopedRunSimCliLogRedirect();

  ScopedRunSimCliLogRedirect(const ScopedRunSimCliLogRedirect &) = delete;
  ScopedRunSimCliLogRedirect &operator=(const ScopedRunSimCliLogRedirect &) = delete;

private:
  int savedStdoutFd = -1;
  int savedStderrFd = -1;
};
}  // namespace pgo::RunSim
