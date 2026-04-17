#pragma once

#include <filesystem>
#include <fstream>
#include <iosfwd>
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
  std::ofstream logStream;
  std::streambuf *coutBuffer = nullptr;
  std::streambuf *cerrBuffer = nullptr;
};
}  // namespace pgo::RunSim
