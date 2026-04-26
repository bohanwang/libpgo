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
  static void init(const char *filename = nullptr, spdlog::level::level_enum level = spdlog::level::info);
  static void setLevel(spdlog::level::level_enum level);
  static inline std::shared_ptr<spdlog::logger> lgr() { return logger ? logger : spdlog::default_logger(); }

protected:
  static std::shared_ptr<spdlog::logger> logger;
};

#define PGO_ALOG(cond)                                                        \
  do {                                                                        \
    if ((cond) == false) {                                                    \
      SPDLOG_LOGGER_CRITICAL(pgo::Logging::lgr(), "{} failed. Abort", #cond); \
    }                                                                         \
  } while (0)
}  // namespace pgo
