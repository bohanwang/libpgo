#include "ipc/broadPhase/surfaceIPCBroadPhase.h"
#include "ipc/core/surfaceIPCCore.h"
#include "ipc/external/obstacleSurface.h"
#include "ipc/topology/surfaceIPCTopology.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <tuple>
#include <utility>
#include <vector>

namespace ES = pgo::EigenSupport;
using pgo::Contact::CIPC::ExternalEEPair;
using pgo::Contact::CIPC::ExternalPairSet;
using pgo::Contact::CIPC::ExternalPTPair;
using pgo::Contact::CIPC::ExternalTPPair;
using pgo::Contact::CIPC::ObstacleSurface;
using pgo::Contact::CIPC::SurfaceIPCCore;
using pgo::Contact::CIPC::SurfaceIPCTopology;

static ES::VXd flattenRows(const ES::MXd &V)
{
  ES::VXd x(V.rows() * 3);
  for (int vi = 0; vi < V.rows(); ++vi)
    x.segment<3>(3 * vi) = V.row(vi).transpose();
  return x;
}

static std::pair<ES::MXd, ES::MXi> makeUnitSquareMesh()
{
  ES::MXd V(4, 3);
  V << 0.0, 0.0, 0.0,
       1.0, 0.0, 0.0,
       0.0, 1.0, 0.0,
       1.0, 1.0, 0.0;
  ES::MXi F(2, 3);
  F << 0, 1, 2,
       1, 3, 2;
  return { V, F };
}

static std::pair<ES::MXd, ES::MXi> makeSmallBoxObstacle()
{
  ES::MXd V(8, 3);
  V << -0.5, -0.5, -0.5,
        0.5, -0.5, -0.5,
        0.5,  0.5, -0.5,
       -0.5,  0.5, -0.5,
       -0.5, -0.5,  0.5,
        0.5, -0.5,  0.5,
        0.5,  0.5,  0.5,
       -0.5,  0.5,  0.5;
  ES::MXi F(12, 3);
  F << 0, 2, 1,  0, 3, 2,
       4, 5, 6,  4, 6, 7,
       0, 1, 5,  0, 5, 4,
       1, 2, 6,  1, 6, 5,
       2, 3, 7,  2, 7, 6,
       3, 0, 4,  3, 4, 7;
  return { V, F };
}

static long long weightKey(double weight)
{
  return std::llround(weight * 1e12);
}

static auto canonicalPT(const std::vector<ExternalPTPair> &pairs)
{
  std::vector<std::tuple<int32_t, int, int, int, int, long long>> keys;
  keys.reserve(pairs.size());
  for (const auto &pair : pairs)
    keys.emplace_back(
      pair.obstacleSlot, pair.dynVertex, pair.obsTri[0], pair.obsTri[1], pair.obsTri[2], weightKey(pair.weight));
  std::sort(keys.begin(), keys.end());
  return keys;
}

static auto canonicalTP(const std::vector<ExternalTPPair> &pairs)
{
  std::vector<std::tuple<int32_t, int, int, int, int, long long>> keys;
  keys.reserve(pairs.size());
  for (const auto &pair : pairs)
    keys.emplace_back(
      pair.obstacleSlot, pair.dynTri[0], pair.dynTri[1], pair.dynTri[2], pair.obsVertex, weightKey(pair.weight));
  std::sort(keys.begin(), keys.end());
  return keys;
}

static auto canonicalEE(const std::vector<ExternalEEPair> &pairs)
{
  std::vector<std::tuple<int32_t, int, int, int, int, long long>> keys;
  keys.reserve(pairs.size());
  for (const auto &pair : pairs)
    keys.emplace_back(
      pair.obstacleSlot, pair.dynEdge[0], pair.dynEdge[1], pair.obsEdge[0], pair.obsEdge[1], weightKey(pair.weight));
  std::sort(keys.begin(), keys.end());
  return keys;
}

TEST(SurfaceIPCExternalBroadPhaseGTest, BuilderMatchesSurfaceIPCCoreExternalPairs)
{
  auto [V, F] = makeUnitSquareMesh();
  auto [obsV, obsF] = makeSmallBoxObstacle();
  const ES::VXd obsRest = flattenRows(obsV);

  auto makeObs = [&]() {
    ObstacleSurface o(obsV, obsF,
      pgo::Contact::CIPC::makeLinearTrajectorySampler(obsRest, ES::V3d::Zero()));
    o.update(0.0, 0.0);
    return o;
  };

  SurfaceIPCCore::Parameters params;
  params.dhat_external = 1.0;

  std::vector<ObstacleSurface> coreObstacles;
  coreObstacles.emplace_back(makeObs());
  SurfaceIPCCore core(params, std::move(coreObstacles));
  core.setMesh(V, F);

  const ES::VXd x = flattenRows(V);
  const auto activeSet = core.buildActiveSet(x);

  SurfaceIPCTopology topology;
  topology.setMesh(V, F);
  std::vector<ObstacleSurface> obstacles;
  obstacles.emplace_back(makeObs());
  obstacles.front().setObjectId(0);

  ExternalPairSet pairs;
  buildExternalPairs(topology, x, obstacles, params.dhat_external, pairs);

  EXPECT_EQ(canonicalPT(pairs.ptPairs), canonicalPT(activeSet.externalPairs.ptPairs));
  EXPECT_EQ(canonicalTP(pairs.tpPairs), canonicalTP(activeSet.externalPairs.tpPairs));
  EXPECT_EQ(canonicalEE(pairs.eePairs), canonicalEE(activeSet.externalPairs.eePairs));
  EXPECT_GT(pairs.size(), 0u);
}

TEST(SurfaceIPCExternalBroadPhaseGTest, MovingObstacleProducesGoldenPairsAndWeights)
{
  ES::MXd V(3, 3);
  V << 0.0, 0.0, 0.0,
       1.0, 0.0, 0.0,
       0.0, 1.0, 0.0;
  ES::MXi F(1, 3);
  F << 0, 1, 2;

  ES::MXd obsV(3, 3);
  obsV << 0.0, 0.0, 0.2,
          1.0, 0.0, 0.2,
          0.0, 1.0, 0.2;
  ES::MXi obsF(1, 3);
  obsF << 0, 1, 2;

  ObstacleSurface obs(
    obsV, obsF,
    pgo::Contact::CIPC::makeLinearTrajectorySampler(flattenRows(obsV), ES::V3d(0.0, 0.0, -0.15)));
  obs.setObjectId(7);
  obs.update(0.0, 1.0);

  SurfaceIPCTopology topology;
  topology.setMesh(V, F);
  std::vector<ObstacleSurface> obstacles;
  obstacles.emplace_back(std::move(obs));

  ExternalPairSet pairs;
  buildExternalPairs(topology, flattenRows(V), obstacles, 0.2, pairs);

  const long long oneTwelfth = weightKey(1.0 / 12.0);
  const long long half = weightKey(0.5);
  const long long one = weightKey(1.0);
  const long long sqrtTwo = weightKey(std::sqrt(2.0));
  const long long two = weightKey(2.0);

  EXPECT_EQ(canonicalPT(pairs.ptPairs),
    (std::vector<std::tuple<int32_t, int, int, int, int, long long>>{
      { 7, 0, 0, 1, 2, oneTwelfth },
      { 7, 1, 0, 1, 2, oneTwelfth },
      { 7, 2, 0, 1, 2, oneTwelfth },
    }));
  EXPECT_EQ(canonicalTP(pairs.tpPairs),
    (std::vector<std::tuple<int32_t, int, int, int, int, long long>>{
      { 7, 0, 1, 2, 0, half },
      { 7, 0, 1, 2, 1, half },
      { 7, 0, 1, 2, 2, half },
    }));
  EXPECT_EQ(canonicalEE(pairs.eePairs),
    (std::vector<std::tuple<int32_t, int, int, int, int, long long>>{
      { 7, 0, 1, 0, 1, one },
      { 7, 0, 1, 0, 2, one },
      { 7, 0, 1, 1, 2, sqrtTwo },
      { 7, 0, 2, 0, 1, one },
      { 7, 0, 2, 0, 2, one },
      { 7, 0, 2, 1, 2, sqrtTwo },
      { 7, 1, 2, 0, 1, sqrtTwo },
      { 7, 1, 2, 0, 2, sqrtTwo },
      { 7, 1, 2, 1, 2, two },
    }));
}
