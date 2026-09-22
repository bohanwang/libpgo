#include <gtest/gtest.h>

#include "ipc/broadPhase/surfaceIPCSelfBroadPhase.h"
#include "ipc/topology/surfaceIPCTopology.h"
#include "ipc/core/surfaceIPCCore.h"

#include "testCIPCHelpers.h"

#include <algorithm>
#include <array>
#include <cstdint>
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

// Triangle 1 stands vertically in the plane y = 0.2 and passes through the
// interior of triangle 0, which lies in the plane z = 0. Two of its edges cross
// triangle 0 while every point-triangle and edge-edge distance is at least 0.2,
// so a distance-based check with dhat = 0.1 sees nothing.
std::tuple<ES::MXd, ES::MXi> makePiercingTriangleMesh()
{
  ES::MXd V(6, 3);
  V <<
    0.0, 0.0, 0.0,
    1.0, 0.0, 0.0,
    0.0, 1.0, 0.0,
    0.2, 0.2, -0.5,
    0.2, 0.2, 0.5,
    0.4, 0.2, 0.5;

  ES::MXi F(2, 3);
  F <<
    0, 1, 2,
    3, 4, 5;

  return { V, F };
}

int findEdgeIndex(const SurfaceIPCTopology &topology, int a, int b)
{
  for (int ei = 0; ei < static_cast<int>(topology.edges.size()); ++ei) {
    const auto &edge = topology.edges[ei];
    if ((edge[0] == a && edge[1] == b) || (edge[0] == b && edge[1] == a))
      return ei;
  }
  return -1;
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

TEST(SurfaceIPCSelfBroadPhaseGTest, FindsEdgesCrossingTriangleInteriorsThatDistancesCannotSee)
{
  const auto [V, F] = makePiercingTriangleMesh();
  const ES::VXd x = flattenPositions(V);

  SurfaceIPCTopology topology;
  topology.setMesh(V, F);

  std::vector<PTPair> ptPairs;
  std::vector<EEPair> eePairs;
  SurfaceIPCSelfBroadPhase().buildPairs(topology, x, 0.1, ptPairs, eePairs);
  EXPECT_TRUE(ptPairs.empty());
  EXPECT_TRUE(eePairs.empty());

  std::vector<std::array<int, 2>> crossings;
  SurfaceIPCSelfBroadPhase().findIntersectingEdgeTrianglePairs(topology, x, crossings);

  std::vector<std::array<int, 2>> expected = {
    { findEdgeIndex(topology, 3, 4), 0 },
    { findEdgeIndex(topology, 3, 5), 0 },
  };
  std::sort(expected.begin(), expected.end());
  ASSERT_GE(expected[0][0], 0);
  EXPECT_EQ(crossings, expected);
}

TEST(SurfaceIPCSelfBroadPhaseGTest, IntersectionQueryIgnoresSeparatedAdjacentAndStaticOnlyPairs)
{
  std::vector<std::array<int, 2>> crossings;

  // separated triangles
  {
    const auto [V, F] = makeTwoTriangleMesh();
    SurfaceIPCTopology topology;
    topology.setMesh(V, F);
    SurfaceIPCSelfBroadPhase().findIntersectingEdgeTrianglePairs(topology, flattenPositions(V), crossings);
    EXPECT_TRUE(crossings.empty());
  }

  // two coplanar triangles sharing an edge touch along it; shared vertices exclude the pair
  {
    ES::MXd V(4, 3);
    V <<
      0.0, 0.0, 0.0,
      1.0, 0.0, 0.0,
      0.0, 1.0, 0.0,
      1.0, 1.0, 0.0;
    ES::MXi F(2, 3);
    F <<
      0, 1, 2,
      1, 3, 2;
    SurfaceIPCTopology topology;
    topology.setMesh(V, F);
    SurfaceIPCSelfBroadPhase().findIntersectingEdgeTrianglePairs(topology, flattenPositions(V), crossings);
    EXPECT_TRUE(crossings.empty());
  }

  // a crossing between two static triangles is not the simulation's problem
  {
    const auto [V, F] = makePiercingTriangleMesh();
    SurfaceIPCTopology topology;
    topology.setMesh(V, F, std::vector<uint8_t>(6, 0));
    SurfaceIPCSelfBroadPhase().findIntersectingEdgeTrianglePairs(topology, flattenPositions(V), crossings);
    EXPECT_TRUE(crossings.empty());

    // a deformable triangle pierced by a static one is
    topology.setMesh(V, F, std::vector<uint8_t>{ 1, 1, 1, 0, 0, 0 });
    SurfaceIPCSelfBroadPhase().findIntersectingEdgeTrianglePairs(topology, flattenPositions(V), crossings);
    EXPECT_EQ(crossings.size(), 2u);
  }
}
