#pragma once

#include <spdlog/common.h>

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
