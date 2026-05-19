#pragma once

#include <chrono>
#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

namespace pgo::Profiling
{

struct ProfileStat
{
  std::string name;
  std::uint64_t callCount = 0;
  double totalSeconds = 0.0;
  double maxSeconds = 0.0;
};

struct ProfileCounterStat
{
  std::string name;
  std::uint64_t sampleCount = 0;
  std::uint64_t total = 0;
  std::uint64_t max = 0;
};

void setProfilingEnabled(bool enabled);
bool isProfilingEnabled();

void resetProfileStatistics();
std::vector<ProfileStat> snapshotProfileStatistics();
void recordProfileCounter(std::string_view name, std::uint64_t value);
std::vector<ProfileCounterStat> snapshotProfileCounterStatistics();

class ScopedProfileSection
{
public:
  explicit ScopedProfileSection(std::string_view name);
  ~ScopedProfileSection();

  ScopedProfileSection(const ScopedProfileSection &) = delete;
  ScopedProfileSection &operator=(const ScopedProfileSection &) = delete;
  ScopedProfileSection(ScopedProfileSection &&) = delete;
  ScopedProfileSection &operator=(ScopedProfileSection &&) = delete;

private:
  using Clock = std::chrono::steady_clock;

  std::string_view name_;
  Clock::time_point start_;
  bool active_ = false;
};

}  // namespace pgo::Profiling
