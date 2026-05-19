#include "scopedProfileSection.h"

#include <algorithm>
#include <atomic>
#include <mutex>
#include <unordered_map>

namespace pgo::Profiling
{
namespace
{
struct AggregatedProfileStat
{
  std::uint64_t callCount = 0;
  double totalSeconds = 0.0;
  double maxSeconds = 0.0;
};

struct AggregatedProfileCounterStat
{
  std::uint64_t sampleCount = 0;
  std::uint64_t total = 0;
  std::uint64_t max = 0;
};

std::atomic<bool> gProfilingEnabled{ false };
std::mutex gProfileMutex;
std::unordered_map<std::string, AggregatedProfileStat> gProfileStats;
std::unordered_map<std::string, AggregatedProfileCounterStat> gProfileCounterStats;

void recordProfileSample(std::string_view name, double seconds)
{
  std::lock_guard<std::mutex> lock(gProfileMutex);
  AggregatedProfileStat &stat = gProfileStats[std::string(name)];
  stat.callCount += 1;
  stat.totalSeconds += seconds;
  stat.maxSeconds = std::max(stat.maxSeconds, seconds);
}

template<typename Stat>
void sortByName(std::vector<Stat> &stats)
{
  std::sort(stats.begin(), stats.end(),
    [](const Stat &lhs, const Stat &rhs) {
      return lhs.name < rhs.name;
    });
}
}  // namespace

void setProfilingEnabled(bool enabled)
{
  gProfilingEnabled.store(enabled, std::memory_order_relaxed);
}

bool isProfilingEnabled()
{
  return gProfilingEnabled.load(std::memory_order_relaxed);
}

void resetProfileStatistics()
{
  std::lock_guard<std::mutex> lock(gProfileMutex);
  gProfileStats.clear();
  gProfileCounterStats.clear();
}

std::vector<ProfileStat> snapshotProfileStatistics()
{
  std::lock_guard<std::mutex> lock(gProfileMutex);

  std::vector<ProfileStat> snapshot;
  snapshot.reserve(gProfileStats.size());
  for (const auto &[name, stat] : gProfileStats) {
    snapshot.push_back(ProfileStat{
      .name = name,
      .callCount = stat.callCount,
      .totalSeconds = stat.totalSeconds,
      .maxSeconds = stat.maxSeconds,
    });
  }

  sortByName(snapshot);

  return snapshot;
}

void recordProfileCounter(std::string_view name, std::uint64_t value)
{
  if (!isProfilingEnabled())
    return;

  std::lock_guard<std::mutex> lock(gProfileMutex);
  AggregatedProfileCounterStat &stat = gProfileCounterStats[std::string(name)];
  stat.sampleCount += 1;
  stat.total += value;
  stat.max = std::max(stat.max, value);
}

std::vector<ProfileCounterStat> snapshotProfileCounterStatistics()
{
  std::lock_guard<std::mutex> lock(gProfileMutex);

  std::vector<ProfileCounterStat> snapshot;
  snapshot.reserve(gProfileCounterStats.size());
  for (const auto &[name, stat] : gProfileCounterStats) {
    snapshot.push_back(ProfileCounterStat{
      .name = name,
      .sampleCount = stat.sampleCount,
      .total = stat.total,
      .max = stat.max,
    });
  }

  sortByName(snapshot);

  return snapshot;
}

ScopedProfileSection::ScopedProfileSection(std::string_view name):
  name_(name)
{
  if (!isProfilingEnabled())
    return;

  start_ = Clock::now();
  active_ = true;
}

ScopedProfileSection::~ScopedProfileSection()
{
  if (!active_)
    return;

  const std::chrono::duration<double> elapsed = Clock::now() - start_;
  recordProfileSample(name_, elapsed.count());
}

}  // namespace pgo::Profiling
