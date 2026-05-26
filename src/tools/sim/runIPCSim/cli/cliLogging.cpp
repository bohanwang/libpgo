#include "cli/cliLogging.h"

#include "configFileJSON.h"

#include <cstdio>
#include <cerrno>
#include <iostream>
#include <stdexcept>
#include <system_error>

#if defined(_WIN32)
#include <fcntl.h>
#include <io.h>
#include <sys/stat.h>

namespace
{
constexpr int kLogOpenFlags = _O_WRONLY | _O_CREAT | _O_TRUNC | _O_BINARY;
constexpr int kLogFileMode = _S_IREAD | _S_IWRITE;

int dupFD(int fd)
{
  return _dup(fd);
}

int dup2FD(int oldfd, int newfd)
{
  return _dup2(oldfd, newfd);
}

int closeFD(int fd)
{
  return _close(fd);
}

int openFD(const char *path)
{
  return _open(path, kLogOpenFlags, kLogFileMode);
}
}  // namespace
#else
#include <fcntl.h>
#include <unistd.h>

namespace
{
constexpr int kLogOpenFlags = O_WRONLY | O_CREAT | O_TRUNC;
constexpr int kLogFileMode = 0644;

int dupFD(int fd)
{
  return dup(fd);
}

int dup2FD(int oldfd, int newfd)
{
  return dup2(oldfd, newfd);
}

int closeFD(int fd)
{
  return close(fd);
}

int openFD(const char *path)
{
  return open(path, kLogOpenFlags, kLogFileMode);
}
}  // namespace
#endif

namespace pgo::RunSim
{
std::filesystem::path deriveDefaultLogPathFromConfig(const std::filesystem::path &configPath)
{
  return configPath.parent_path() / (configPath.stem().string() + ".log");
}

spdlog::level::level_enum resolveConfiguredLogLevel(const ConfigFileJSON &config)
{
  const std::string configuredLevel = config.exist("loglevel")
    ? config.getString("loglevel")
    : "info";

  if (configuredLevel == "trace")
    return spdlog::level::trace;
  if (configuredLevel == "debug")
    return spdlog::level::debug;
  if (configuredLevel == "info")
    return spdlog::level::info;
  if (configuredLevel == "warn")
    return spdlog::level::warn;
  if (configuredLevel == "error")
    return spdlog::level::err;

  throw std::invalid_argument("Unsupported loglevel: " + configuredLevel);
}

ScopedRunSimCliLogRedirect::ScopedRunSimCliLogRedirect(const std::string &logFilename)
{
  const std::filesystem::path logPath(logFilename);
  if (logPath.has_parent_path()) {
    std::filesystem::create_directories(logPath.parent_path());
  }

  std::cout.flush();
  std::cerr.flush();
  std::fflush(stdout);
  std::fflush(stderr);

  savedStdoutFd = dupFD(fileno(stdout));
  if (savedStdoutFd < 0) {
    throw std::system_error(errno, std::generic_category(), "Failed to backup stdout fd");
  }

  savedStderrFd = dupFD(fileno(stderr));
  if (savedStderrFd < 0) {
    closeFD(savedStdoutFd);
    savedStdoutFd = -1;
    throw std::system_error(errno, std::generic_category(), "Failed to backup stderr fd");
  }

  int logFd = openFD(logPath.string().c_str());
  if (logFd < 0) {
    closeFD(savedStdoutFd);
    closeFD(savedStderrFd);
    savedStdoutFd = -1;
    savedStderrFd = -1;
    throw std::runtime_error("Failed to open runSim log file: " + logFilename);
  }

  if (dup2FD(logFd, fileno(stdout)) < 0 || dup2FD(logFd, fileno(stderr)) < 0) {
    const int dup2Errno = errno;
    closeFD(logFd);
    if (savedStdoutFd >= 0) {
      dup2FD(savedStdoutFd, fileno(stdout));
      closeFD(savedStdoutFd);
      savedStdoutFd = -1;
    }
    if (savedStderrFd >= 0) {
      dup2FD(savedStderrFd, fileno(stderr));
      closeFD(savedStderrFd);
      savedStderrFd = -1;
    }
    throw std::system_error(dup2Errno, std::generic_category(), "Failed to redirect stdout/stderr to runSim log file");
  }

  closeFD(logFd);
}

ScopedRunSimCliLogRedirect::~ScopedRunSimCliLogRedirect()
{
  std::cout.flush();
  std::cerr.flush();
  std::fflush(stdout);
  std::fflush(stderr);

  if (savedStdoutFd >= 0) {
    dup2FD(savedStdoutFd, fileno(stdout));
    closeFD(savedStdoutFd);
    savedStdoutFd = -1;
  }

  if (savedStderrFd >= 0) {
    dup2FD(savedStderrFd, fileno(stderr));
    closeFD(savedStderrFd);
    savedStderrFd = -1;
  }
}
}  // namespace pgo::RunSim
