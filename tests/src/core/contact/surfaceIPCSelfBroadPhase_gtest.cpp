#include <gtest/gtest.h>

#include "ipc/broadPhase/surfaceIPCSelfBroadPhase.h"
#include "ipc/topology/surfaceIPCTopology.h"
#include "ipc/core/surfaceIPCCore.h"

#include "testCIPCHelpers.h"

#include <algorithm>
#include <tuple>
#include <vector>

namespace
{
namespace ES = pgo::EigenSupport;
using pgo::Contact::CIPC::EEPair;
using pgo::Contact::CIPC::PTPair;
using pgo::Contact::CIPC::SurfaceIPCCore;
using pgo::Contact::CIPC::SurfaceIPCSelfBroadPhase;
using pgo::Contact::CIPC::SurfaceIPCTopology;
using pgo::Contact::CIPCTest::flattenPositions;
using pgo::Contact::CIPCTest::makeTwoTriangleMesh;

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
}  // namespace

TEST(SurfaceIPCSelfBroadPhaseGTest, BuilderMatchesSurfaceIPCCorePairSet)
{
  const auto [V, F] = makeTwoTriangleMesh();
  const ES::VXd x = flattenPositions(V);

  SurfaceIPCTopology topology;
  topology.setMesh(V, F);

  std::vector<PTPair> broadPhasePTPairs;
  std::vector<EEPair> broadPhaseEEPairs;
  SurfaceIPCSelfBroadPhase().buildPairs(topology, x, 0.1, broadPhasePTPairs, broadPhaseEEPairs);

  SurfaceIPCCore core;
  SurfaceIPCCore::Parameters params;
  params.dhat = 0.1;
  params.kappa = 1.0;
  params.eps_ee = 0.0;
  params.slackness = 0.9;
  core.setParameters(params);
  core.setMesh(V, F);
  core.computeEnergy(x);

  EXPECT_EQ(canonicalPT(broadPhasePTPairs), canonicalPT(core.getPTPairs()));
  EXPECT_EQ(canonicalEE(broadPhaseEEPairs), canonicalEE(core.getEEPairs()));
}

TEST(SurfaceIPCSelfBroadPhaseGTest, MixedOwnershipRetainsOnlyPairsWithDeformableGeometry)
{
  const auto [V, F] = makeTwoTriangleMesh();
  const ES::VXd x = flattenPositions(V);
  const std::vector<uint8_t> mixedMask = { 1, 1, 1, 0, 0, 0 };

  SurfaceIPCTopology topology;
  topology.setMesh(V, F, mixedMask);

  std::vector<PTPair> ptPairs;
  std::vector<EEPair> eePairs;
  SurfaceIPCSelfBroadPhase().buildPairs(topology, x, 0.25, ptPairs, eePairs);

  ASSERT_FALSE(ptPairs.empty());
  ASSERT_FALSE(eePairs.empty());
  bool sawDynamicPointStaticTriangle = false;
  bool sawStaticPointDynamicTriangle = false;
  for (const auto &pair : ptPairs) {
    const bool pointDynamic = topology.isVertexDeformable(pair.p);
    const bool triangleDynamic = topology.isVertexDeformable(pair.t0) ||
      topology.isVertexDeformable(pair.t1) || topology.isVertexDeformable(pair.t2);
    EXPECT_TRUE(pointDynamic || triangleDynamic);
    sawDynamicPointStaticTriangle |= pointDynamic && !triangleDynamic;
    sawStaticPointDynamicTriangle |= !pointDynamic && triangleDynamic;
  }
  EXPECT_TRUE(sawDynamicPointStaticTriangle);
  EXPECT_TRUE(sawStaticPointDynamicTriangle);
  for (const auto &pair : eePairs) {
    const bool edgeADynamic = topology.isVertexDeformable(pair.ea0) || topology.isVertexDeformable(pair.ea1);
    const bool edgeBDynamic = topology.isVertexDeformable(pair.eb0) || topology.isVertexDeformable(pair.eb1);
    EXPECT_TRUE(edgeADynamic || edgeBDynamic);
    EXPECT_NE(edgeADynamic, edgeBDynamic);
  }

  SurfaceIPCCore core;
  SurfaceIPCCore::Parameters params;
  params.dhat = 0.25;
  params.kappa = 1.0;
  params.slackness = 0.9;
  core.setParameters(params);
  core.setMesh(V, F, mixedMask);
  core.computeEnergy(x);
  EXPECT_EQ(canonicalPT(ptPairs), canonicalPT(core.getPTPairs()));
  EXPECT_EQ(canonicalEE(eePairs), canonicalEE(core.getEEPairs()));
}

TEST(SurfaceIPCSelfBroadPhaseGTest, FullyExternalGeometryProducesNoPairs)
{
  const auto [V, F] = makeTwoTriangleMesh();
  const ES::VXd x = flattenPositions(V);
  SurfaceIPCTopology topology;
  topology.setMesh(V, F, std::vector<uint8_t>(6, 0));

  std::vector<PTPair> ptPairs;
  std::vector<EEPair> eePairs;
  SurfaceIPCSelfBroadPhase().buildPairs(topology, x, 0.1, ptPairs, eePairs);
  EXPECT_TRUE(ptPairs.empty());
  EXPECT_TRUE(eePairs.empty());
}

TEST(SurfaceIPCSelfBroadPhaseGTest, ExplicitAllDeformableMaskMatchesDefault)
{
  const auto [V, F] = makeTwoTriangleMesh();
  const ES::VXd x = flattenPositions(V);
  SurfaceIPCTopology defaultTopology;
  defaultTopology.setMesh(V, F);
  SurfaceIPCTopology explicitTopology;
  explicitTopology.setMesh(V, F, std::vector<uint8_t>(6, 1));

  std::vector<PTPair> defaultPT, explicitPT;
  std::vector<EEPair> defaultEE, explicitEE;
  SurfaceIPCSelfBroadPhase().buildPairs(defaultTopology, x, 0.1, defaultPT, defaultEE);
  SurfaceIPCSelfBroadPhase().buildPairs(explicitTopology, x, 0.1, explicitPT, explicitEE);
  EXPECT_EQ(canonicalPT(defaultPT), canonicalPT(explicitPT));
  EXPECT_EQ(canonicalEE(defaultEE), canonicalEE(explicitEE));
}
