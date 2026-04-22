#include <gtest/gtest.h>

#include "CIPC.h"
#include "scopedProfileSection.h"
#include "surfaceIPCProfiling.h"

#include <algorithm>
#include <string_view>
#include <vector>

namespace
{
using pgo::Contact::CIPC::CIPCPotentialEnergy;
using pgo::EigenSupport::MXd;
using pgo::EigenSupport::MXi;
using pgo::EigenSupport::VXd;
using pgo::Profiling::ProfileStat;

const ProfileStat *findStat(const std::vector<ProfileStat> &stats, std::string_view name)
{
  const auto it = std::find_if(stats.begin(), stats.end(),
    [name](const ProfileStat &stat) { return stat.name == name; });
  return it == stats.end() ? nullptr : &(*it);
}
}  // namespace

TEST(CIPCProfilingGTest, FuncRecordsSurfaceSectionsButNotWrapperSections)
{
  MXd V(3, 3);
  V <<
    0.0, 0.0, 0.0,
    1.0, 0.0, 0.0,
    0.0, 1.0, 0.0;

  MXi F(1, 3);
  F << 0, 1, 2;

  CIPCPotentialEnergy energy(/*dhat=*/0.1, /*kappa=*/1.0, /*isInputDisp=*/false);
  energy.setMesh(V, F);

  VXd x(9);
  x <<
    0.0, 0.0, 0.0,
    1.0, 0.0, 0.0,
    0.0, 1.0, 0.0;

  pgo::Profiling::setProfilingEnabled(true);
  pgo::Profiling::resetProfileStatistics();

  const double value = energy.func(x);
  (void)value;

  const auto stats = pgo::Profiling::snapshotProfileStatistics();
  EXPECT_NE(findStat(stats, pgo::Contact::SurfaceIPCProfileSections::kPairBuildStatic), nullptr);
  EXPECT_NE(findStat(stats, pgo::Contact::SurfaceIPCProfileSections::kEnergy), nullptr);
  EXPECT_EQ(findStat(stats, pgo::Contact::SurfaceIPCProfileSections::kWrapperSync), nullptr);
  EXPECT_EQ(findStat(stats, pgo::Contact::SurfaceIPCProfileSections::kFloorPostPass), nullptr);

  pgo::Profiling::setProfilingEnabled(false);
  pgo::Profiling::resetProfileStatistics();
}
