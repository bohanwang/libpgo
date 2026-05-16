#include <gtest/gtest.h>

#include "ipc/broadPhase/surfaceIPCBroadPhase.h"
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
using pgo::Contact::CIPC::SelfPairSet;
using pgo::Contact::CIPC::SurfaceIPCCore;
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
  core.computeEnergy(x);

  EXPECT_EQ(canonicalPT(broadPhasePairs.ptPairs), canonicalPT(core.preparedState().selfPairs.ptPairs));
  EXPECT_EQ(canonicalEE(broadPhasePairs.eePairs), canonicalEE(core.preparedState().selfPairs.eePairs));
}
