#include <gtest/gtest.h>

#include "CIPC.h"
#include "scopedProfileSection.h"
#include "ipc/profiling/surfaceIPCProfiling.h"

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

TEST(CIPCProfilingGTest, FuncRecordsWrapperSyncAndSurfaceSectionsWithoutFloorPostPass)
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
  EXPECT_NE(findStat(stats, pgo::Contact::SurfaceIPCProfileSections::kWrapperSync), nullptr);
  EXPECT_NE(findStat(stats, pgo::Contact::SurfaceIPCProfileSections::kBuildActiveSet), nullptr);
  EXPECT_NE(findStat(stats, pgo::Contact::SurfaceIPCProfileSections::kPairBuildStatic), nullptr);
  EXPECT_NE(findStat(stats, pgo::Contact::SurfaceIPCProfileSections::kActiveSetEnergy), nullptr);
  EXPECT_EQ(findStat(stats, pgo::Contact::SurfaceIPCProfileSections::kFloorPostPass), nullptr);

  pgo::Profiling::setProfilingEnabled(false);
  pgo::Profiling::resetProfileStatistics();
}

TEST(CIPCProfilingGTest, FuncRecordsFloorPostPassWhenFloorIsActive)
{
  MXd V(6, 3);
  V <<
    0.0, 0.0, 0.00,
    1.0, 0.0, 0.00,
    0.0, 1.0, 0.00,
    0.2, 0.2, 0.05,
    1.2, 0.2, 0.05,
    0.2, 1.2, 0.05;

  MXi F(2, 3);
  F <<
    0, 1, 2,
    3, 4, 5;

  CIPCPotentialEnergy energy(/*dhat=*/0.1, /*kappa=*/1.0, /*isInputDisp=*/false,
    /*eps_ee=*/0.0, /*useFloor=*/true, /*floorHeight=*/0.08, /*floorKappa=*/2.0);
  energy.setMesh(V, F);

  VXd x(18);
  for (int vi = 0; vi < V.rows(); ++vi)
    x.segment<3>(3 * vi) = V.row(vi).transpose();

  pgo::Profiling::setProfilingEnabled(true);
  pgo::Profiling::resetProfileStatistics();

  (void)energy.func(x);

  const auto stats = pgo::Profiling::snapshotProfileStatistics();
  EXPECT_NE(findStat(stats, pgo::Contact::SurfaceIPCProfileSections::kWrapperSync), nullptr);
  EXPECT_NE(findStat(stats, pgo::Contact::SurfaceIPCProfileSections::kBuildActiveSet), nullptr);
  EXPECT_NE(findStat(stats, pgo::Contact::SurfaceIPCProfileSections::kPairBuildStatic), nullptr);
  EXPECT_NE(findStat(stats, pgo::Contact::SurfaceIPCProfileSections::kActiveSetEnergy), nullptr);
  EXPECT_NE(findStat(stats, pgo::Contact::SurfaceIPCProfileSections::kFloorPostPass), nullptr);

  pgo::Profiling::setProfilingEnabled(false);
  pgo::Profiling::resetProfileStatistics();
}
