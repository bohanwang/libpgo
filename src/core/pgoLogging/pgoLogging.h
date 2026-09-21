#pragma once

#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_TRACE
#include <spdlog/spdlog.h>
#include <spdlog/fmt/ostr.h>

#include <memory>

namespace pgo
{
class Logging
{
public:
  // Called once by the process entry point (a tool's main() or pgo_init()),
  // never by library code.
  static void init(const char *filename = nullptr, spdlog::level::level_enum level = spdlog::level::info);
  static void setLevel(spdlog::level::level_enum level);
  // Re-evaluates whether the console sink emits ANSI colors. Call after
  // redirecting stdout into a file, and again after restoring it.
  static void refreshConsoleColorMode();
  // The logger created by init(). Before init() this is spdlog's default
  // logger, so library code never logs through a null pointer.
  static inline std::shared_ptr<spdlog::logger> lgr() { return logger ? logger : spdlog::default_logger(); }

protected:
  static std::shared_ptr<spdlog::logger> logger;
  static std::shared_ptr<spdlog::sinks::sink> consoleSink;
};

#define PGO_ALOG(cond)                                                        \
  do {                                                                        \
    if ((cond) == false) {                                                    \
      SPDLOG_LOGGER_CRITICAL(pgo::Logging::lgr(), "{} failed. Abort", #cond); \
    }                                                                         \
  } while (0)
}  // namespace pgo
