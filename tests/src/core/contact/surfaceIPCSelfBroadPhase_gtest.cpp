#include <gtest/gtest.h>

#include "ipc/broadPhase/surfaceIPCBroadPhase.h"
#include "ipc/topology/surfaceIPCTopology.h"
#include "ipc/core/surfaceIPCCore.h"
#include "ipc/geometry/ipcDistancePrimitives.h"
#include "ipc/profiling/surfaceIPCProfiling.h"
#include "scopedProfileSection.h"

#include "testCIPCHelpers.h"

#include <algorithm>
#include <string_view>
#include <tuple>
#include <vector>

namespace
{
namespace ES = pgo::EigenSupport;
using pgo::Contact::CIPC::EEPair;
using pgo::Contact::CIPC::PTPair;
using pgo::Contact::CIPC::SelfPairSet;
using pgo::Contact::CIPC::SurfaceIPCCore;
using pgo::Contact::CIPC::SurfaceIPCTopology;
namespace distance = pgo::Contact::CIPC::distance;
using pgo::Contact::CIPCTest::flattenPositions;
using pgo::Contact::CIPCTest::makeTwoTriangleMesh;
using pgo::Profiling::ProfileCounterStat;

std::vector<std::tuple<int, int, int, int>> canonicalPT(const std::vector<PTPair> &pairs)
{
  std::vector<std::tuple<int, int, int, int>> keys;
  keys.reserve(pairs.size());
  for (const auto &pair : pairs)
    keys.emplace_back(pair.p, pair.t0, pair.t1, pair.t2);
  std::sort(keys.begin(), keys.end());
  return keys;
}

std::vector<std::tuple<int, int, int, int>> canonicalEE(const std::vector<EEPair> &pairs)
{
  std::vector<std::tuple<int, int, int, int>> keys;
  keys.reserve(pairs.size());
  for (const auto &pair : pairs)
    keys.emplace_back(pair.ea0, pair.ea1, pair.eb0, pair.eb1);
  std::sort(keys.begin(), keys.end());
  return keys;
}

const ProfileCounterStat *findCounterStat(const std::vector<ProfileCounterStat> &stats, std::string_view name)
{
  const auto it = std::find_if(stats.begin(), stats.end(),
    [name](const ProfileCounterStat &stat) { return stat.name == name; });
  return it == stats.end() ? nullptr : &(*it);
}

ES::V3d vertex(const ES::VXd &x, int vi)
{
  return x.segment<3>(3 * vi);
}

SelfPairSet buildBruteForcePairs(const SurfaceIPCTopology &topology, const ES::VXd &x, double dhat)
{
  const double dhat2 = dhat * dhat;
  SelfPairSet pairs;

  for (int vi = 0; vi < topology.numVerts; ++vi) {
    for (int fi = 0; fi < static_cast<int>(topology.triangles.size()); ++fi) {
      const auto &tri = topology.triangles[fi];
      if (vi == tri[0] || vi == tri[1] || vi == tri[2])
        continue;

      const double d2 = distance::computePTSqDist(
        vertex(x, vi),
        vertex(x, tri[0]),
        vertex(x, tri[1]),
        vertex(x, tri[2]));
      if (d2 < dhat2)
        pairs.ptPairs.push_back({ vi, tri[0], tri[1], tri[2], topology.vertexArea[vi] * topology.triArea[fi] });
    }
  }

  for (int ei = 0; ei < static_cast<int>(topology.edges.size()); ++ei) {
    const int a0 = topology.edges[ei][0];
    const int a1 = topology.edges[ei][1];
    for (int ej = ei + 1; ej < static_cast<int>(topology.edges.size()); ++ej) {
      const int b0 = topology.edges[ej][0];
      const int b1 = topology.edges[ej][1];
      if (a0 == b0 || a0 == b1 || a1 == b0 || a1 == b1)
        continue;

      const double d2 = distance::computeEESqDist(
        vertex(x, a0), vertex(x, a1),
        vertex(x, b0), vertex(x, b1));
      if (d2 < dhat2)
        pairs.eePairs.push_back({ a0, a1, b0, b1, topology.edgeLength[ei] * topology.edgeLength[ej] });
    }
  }

  return pairs;
}
}  // namespace

TEST(SurfaceIPCSelfBroadPhaseGTest, BuilderMatchesSurfaceIPCCorePairSet)
{
  const auto [V, F] = makeTwoTriangleMesh();
  const ES::VXd x = flattenPositions(V);

  SurfaceIPCTopology topology;
  topology.setMesh(V, F);

  SelfPairSet broadPhasePairs;
  buildSelfPairs(topology, x, 0.1, broadPhasePairs);

  SurfaceIPCCore core;
  SurfaceIPCCore::Parameters params;
  params.dhat = 0.1;
  params.kappa = 1.0;
  params.eps_ee = 0.0;
  params.slackness = 0.9;
  core.setParameters(params);
  core.setMesh(V, F);
  const auto activeSet = core.buildActiveSet(x);

  EXPECT_EQ(canonicalPT(broadPhasePairs.ptPairs), canonicalPT(activeSet.selfPairs.ptPairs));
  EXPECT_EQ(canonicalEE(broadPhasePairs.eePairs), canonicalEE(activeSet.selfPairs.eePairs));
}

TEST(SurfaceIPCSelfBroadPhaseGTest, BuilderMatchesBruteForceNearThreshold)
{
  ES::MXd V(9, 3);
  V << 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0,
    0.0, 1.0, 0.0,
    0.2, 0.2, 0.09,
    1.2, 0.2, 0.09,
    0.2, 1.2, 0.09,
    3.0, 3.0, 3.0,
    4.0, 3.0, 3.0,
    3.0, 4.0, 3.0;

  ES::MXi F(3, 3);
  F << 0, 1, 2,
    3, 4, 5,
    6, 7, 8;

  const double dhat = 0.1;
  const ES::VXd x = flattenPositions(V);

  SurfaceIPCTopology topology;
  topology.setMesh(V, F);

  SelfPairSet broadPhasePairs;
  buildSelfPairs(topology, x, dhat, broadPhasePairs);

  const SelfPairSet bruteForcePairs = buildBruteForcePairs(topology, x, dhat);

  EXPECT_EQ(canonicalPT(broadPhasePairs.ptPairs), canonicalPT(bruteForcePairs.ptPairs));
  EXPECT_EQ(canonicalEE(broadPhasePairs.eePairs), canonicalEE(bruteForcePairs.eePairs));
}

TEST(SurfaceIPCSelfBroadPhaseGTest, ProfilingRecordsSelfCandidateCounters)
{
  const auto [V, F] = makeTwoTriangleMesh();
  const ES::VXd x = flattenPositions(V);

  SurfaceIPCTopology topology;
  topology.setMesh(V, F);

  pgo::Profiling::setProfilingEnabled(true);
  pgo::Profiling::resetProfileStatistics();

  SelfPairSet broadPhasePairs;
  buildSelfPairs(topology, x, 0.1, broadPhasePairs);

  const auto stats = pgo::Profiling::snapshotProfileCounterStatistics();
  const auto *ptHashCandidates = findCounterStat(stats, pgo::Contact::SurfaceIPCProfileSections::kPairBuildSelfPTHashCandidates);
  const auto *ptDistanceTests = findCounterStat(stats, pgo::Contact::SurfaceIPCProfileSections::kPairBuildSelfPTDistanceTests);
  const auto *ptAcceptedPairs = findCounterStat(stats, pgo::Contact::SurfaceIPCProfileSections::kPairBuildSelfPTAcceptedPairs);
  const auto *eeHashCandidates = findCounterStat(stats, pgo::Contact::SurfaceIPCProfileSections::kPairBuildSelfEEHashCandidates);
  const auto *eeDistanceTests = findCounterStat(stats, pgo::Contact::SurfaceIPCProfileSections::kPairBuildSelfEEDistanceTests);
  const auto *eeAcceptedPairs = findCounterStat(stats, pgo::Contact::SurfaceIPCProfileSections::kPairBuildSelfEEAcceptedPairs);

  ASSERT_NE(ptHashCandidates, nullptr);
  ASSERT_NE(ptDistanceTests, nullptr);
  ASSERT_NE(ptAcceptedPairs, nullptr);
  ASSERT_NE(eeHashCandidates, nullptr);
  ASSERT_NE(eeDistanceTests, nullptr);
  ASSERT_NE(eeAcceptedPairs, nullptr);
  EXPECT_GT(ptHashCandidates->total, 0u);
  EXPECT_GT(ptDistanceTests->total, 0u);
  EXPECT_EQ(ptAcceptedPairs->total, broadPhasePairs.ptPairs.size());
  EXPECT_GT(eeHashCandidates->total, 0u);
  EXPECT_GT(eeDistanceTests->total, 0u);
  EXPECT_EQ(eeAcceptedPairs->total, broadPhasePairs.eePairs.size());

  pgo::Profiling::setProfilingEnabled(false);
  pgo::Profiling::resetProfileStatistics();
}
