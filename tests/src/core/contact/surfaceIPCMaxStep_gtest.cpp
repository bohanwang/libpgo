#include <gtest/gtest.h>

#include "ipc/core/surfaceIPCMaxStep.h"
#include "ipc/topology/surfaceIPCTopology.h"
#include "ipc/core/surfaceIPCCore.h"
#include "ipc/profiling/surfaceIPCProfiling.h"
#include "scopedProfileSection.h"

#include "testCIPCHelpers.h"

#include <algorithm>
#include <string_view>
#include <vector>

namespace
{
namespace ES = pgo::EigenSupport;
using pgo::Contact::CIPC::SurfaceIPCCore;
using pgo::Contact::CIPC::SurfaceIPCTopology;
using pgo::Contact::CIPCTest::flattenPositions;
using pgo::Contact::CIPCTest::makeTwoTriangleMesh;
using pgo::Profiling::ProfileCounterStat;

const ProfileCounterStat *findCounterStat(const std::vector<ProfileCounterStat> &stats, std::string_view name)
{
  const auto it = std::find_if(stats.begin(), stats.end(),
    [name](const ProfileCounterStat &stat) { return stat.name == name; });
  return it == stats.end() ? nullptr : &(*it);
}
}  // namespace

TEST(SurfaceIPCMaxStepGTest, HelperMatchesSurfaceIPCCoreMaxStep)
{
  const auto [V, F] = makeTwoTriangleMesh();
  const ES::VXd x = flattenPositions(V);
  ES::VXd dx = ES::VXd::Zero(x.size());
  for (int vi = 3; vi < 6; ++vi)
    dx[3 * vi + 2] = -0.1;

  SurfaceIPCTopology topology;
  topology.setMesh(V, F);

  SurfaceIPCCore core;
  SurfaceIPCCore::Parameters params;
  params.dhat = 0.1;
  params.kappa = 1.0;
  params.eps_ee = 0.0;
  params.slackness = 0.9;
  core.setParameters(params);
  core.setMesh(V, F);

  const double helperAlpha = computeSelfMaxStep(topology, x, dx, params.dhat, params.slackness);
  const double coreAlpha = core.computeMaxStepLimit(x, dx).alpha;

  EXPECT_GT(helperAlpha, 0.0);
  EXPECT_LT(helperAlpha, 1.0);
  EXPECT_NEAR(helperAlpha, coreAlpha, 1e-12);
}

TEST(SurfaceIPCMaxStepGTest, ProfilingRecordsSelfSweptCandidateCounters)
{
  const auto [V, F] = makeTwoTriangleMesh();
  const ES::VXd x = flattenPositions(V);
  ES::VXd dx = ES::VXd::Zero(x.size());
  for (int vi = 3; vi < 6; ++vi)
    dx[3 * vi + 2] = -0.1;

  SurfaceIPCTopology topology;
  topology.setMesh(V, F);

  pgo::Profiling::setProfilingEnabled(true);
  pgo::Profiling::resetProfileStatistics();

  const double alpha = computeSelfMaxStep(topology, x, dx, 0.1, 0.9);
  (void)alpha;

  const auto stats = pgo::Profiling::snapshotProfileCounterStatistics();
  const auto *ptHashCandidates = findCounterStat(stats, pgo::Contact::SurfaceIPCProfileSections::kMaxStepSelfPTHashCandidates);
  const auto *ptCCDTests = findCounterStat(stats, pgo::Contact::SurfaceIPCProfileSections::kMaxStepSelfPTCCDTests);
  const auto *eeHashCandidates = findCounterStat(stats, pgo::Contact::SurfaceIPCProfileSections::kMaxStepSelfEEHashCandidates);
  const auto *eeCCDTests = findCounterStat(stats, pgo::Contact::SurfaceIPCProfileSections::kMaxStepSelfEECCDTests);

  ASSERT_NE(ptHashCandidates, nullptr);
  ASSERT_NE(ptCCDTests, nullptr);
  ASSERT_NE(eeHashCandidates, nullptr);
  ASSERT_NE(eeCCDTests, nullptr);
  EXPECT_GT(ptHashCandidates->total, 0u);
  EXPECT_GT(ptCCDTests->total, 0u);
  EXPECT_GT(eeHashCandidates->total, 0u);
  EXPECT_GT(eeCCDTests->total, 0u);

  pgo::Profiling::setProfilingEnabled(false);
  pgo::Profiling::resetProfileStatistics();
}
