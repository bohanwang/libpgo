#pragma once

#include <filesystem>
#include <string>

namespace pgo::RunSim
{
std::filesystem::path deriveDefaultLogPathFromConfig(const std::filesystem::path &configPath);

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
