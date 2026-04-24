#include <gtest/gtest.h>

#include "scopedProfileSection.h"

#include <algorithm>
#include <chrono>
#include <string_view>
#include <thread>
#include <vector>

namespace
{
using pgo::Profiling::ProfileStat;
using pgo::Profiling::ScopedProfileSection;

const ProfileStat *findStat(const std::vector<ProfileStat> &stats, std::string_view name)
{
  const auto it = std::find_if(stats.begin(), stats.end(),
    [name](const ProfileStat &stat) { return stat.name == name; });
  return it == stats.end() ? nullptr : &(*it);
}

class ScopedProfileSectionGTest : public ::testing::Test
{
protected:
  void SetUp() override
  {
    pgo::Profiling::setProfilingEnabled(false);
    pgo::Profiling::resetProfileStatistics();
  }

  void TearDown() override
  {
    pgo::Profiling::setProfilingEnabled(false);
    pgo::Profiling::resetProfileStatistics();
  }
};
}  // namespace

TEST_F(ScopedProfileSectionGTest, DisabledProfilingLeavesNoRecords)
{
  {
    ScopedProfileSection scope("profiling.disabled");
    std::this_thread::sleep_for(std::chrono::milliseconds(1));
  }

  EXPECT_TRUE(pgo::Profiling::snapshotProfileStatistics().empty());
}

TEST_F(ScopedProfileSectionGTest, RepeatedSectionsAggregateByName)
{
  pgo::Profiling::setProfilingEnabled(true);

  {
    ScopedProfileSection scope("profiling.repeat");
    std::this_thread::sleep_for(std::chrono::milliseconds(1));
  }
  {
    ScopedProfileSection scope("profiling.repeat");
    std::this_thread::sleep_for(std::chrono::milliseconds(1));
  }

  const auto stats = pgo::Profiling::snapshotProfileStatistics();
  ASSERT_EQ(stats.size(), 1u);
  EXPECT_EQ(stats.front().name, "profiling.repeat");
  EXPECT_EQ(stats.front().callCount, 2u);
  EXPECT_GT(stats.front().totalSeconds, 0.0);
  EXPECT_GT(stats.front().maxSeconds, 0.0);
}

TEST_F(ScopedProfileSectionGTest, DistinctSectionsRemainIndependentAndSorted)
{
  pgo::Profiling::setProfilingEnabled(true);

  {
    ScopedProfileSection scope("profiling.zeta");
    std::this_thread::sleep_for(std::chrono::milliseconds(1));
  }
  {
    ScopedProfileSection scope("profiling.alpha");
    std::this_thread::sleep_for(std::chrono::milliseconds(1));
  }

  const auto stats = pgo::Profiling::snapshotProfileStatistics();
  ASSERT_EQ(stats.size(), 2u);
  EXPECT_EQ(stats[0].name, "profiling.alpha");
  EXPECT_EQ(stats[1].name, "profiling.zeta");
  EXPECT_EQ(stats[0].callCount, 1u);
  EXPECT_EQ(stats[1].callCount, 1u);
}

TEST_F(ScopedProfileSectionGTest, NestedSectionsAreTrackedInclusively)
{
  pgo::Profiling::setProfilingEnabled(true);

  {
    ScopedProfileSection outer("profiling.outer");
    std::this_thread::sleep_for(std::chrono::milliseconds(1));
    {
      ScopedProfileSection innerSameName("profiling.outer");
      std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
    {
      ScopedProfileSection innerDifferentName("profiling.inner");
      std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
  }

  const auto stats = pgo::Profiling::snapshotProfileStatistics();
  const ProfileStat *outer = findStat(stats, "profiling.outer");
  const ProfileStat *inner = findStat(stats, "profiling.inner");
  ASSERT_NE(outer, nullptr);
  ASSERT_NE(inner, nullptr);
  EXPECT_EQ(outer->callCount, 2u);
  EXPECT_EQ(inner->callCount, 1u);
  EXPECT_GE(outer->totalSeconds, inner->totalSeconds);
}

TEST_F(ScopedProfileSectionGTest, ResetClearsCollectedStatistics)
{
  pgo::Profiling::setProfilingEnabled(true);

  {
    ScopedProfileSection scope("profiling.reset");
    std::this_thread::sleep_for(std::chrono::milliseconds(1));
  }

  ASSERT_FALSE(pgo::Profiling::snapshotProfileStatistics().empty());
  pgo::Profiling::resetProfileStatistics();
  EXPECT_TRUE(pgo::Profiling::snapshotProfileStatistics().empty());
}
