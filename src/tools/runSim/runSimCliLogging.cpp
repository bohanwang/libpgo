#include "runSimCliLogging.h"

#include <iostream>
#include <stdexcept>

namespace pgo::RunSim
{
std::filesystem::path deriveDefaultLogPathFromConfig(const std::filesystem::path &configPath)
{
  return configPath.parent_path() / (configPath.stem().string() + ".log");
}

ScopedRunSimCliLogRedirect::ScopedRunSimCliLogRedirect(const std::string &logFilename)
{
  const std::filesystem::path logPath(logFilename);
  if (logPath.has_parent_path()) {
    std::filesystem::create_directories(logPath.parent_path());
  }

  logStream.open(logPath, std::ios::out | std::ios::trunc);
  if (!logStream.is_open()) {
    throw std::runtime_error("Failed to open runSim log file: " + logFilename);
  }

  coutBuffer = std::cout.rdbuf(logStream.rdbuf());
  cerrBuffer = std::cerr.rdbuf(logStream.rdbuf());
}

ScopedRunSimCliLogRedirect::~ScopedRunSimCliLogRedirect()
{
  if (coutBuffer) {
    std::cout.rdbuf(coutBuffer);
  }
  if (cerrBuffer) {
    std::cerr.rdbuf(cerrBuffer);
  }
}
}  // namespace pgo::RunSim
