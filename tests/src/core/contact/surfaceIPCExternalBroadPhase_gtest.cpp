#include "ipc/broadPhase/surfaceIPCBroadPhase.h"
#include "ipc/core/surfaceIPCCore.h"
#include "ipc/external/obstacleSurface.h"
#include "ipc/topology/surfaceIPCTopology.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <memory>
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
      pair.obstacleObjectId, pair.dynVertex, pair.obsTri[0], pair.obsTri[1], pair.obsTri[2], weightKey(pair.weight));
  std::sort(keys.begin(), keys.end());
  return keys;
}

static auto canonicalTP(const std::vector<ExternalTPPair> &pairs)
{
  std::vector<std::tuple<int32_t, int, int, int, int, long long>> keys;
  keys.reserve(pairs.size());
  for (const auto &pair : pairs)
    keys.emplace_back(
      pair.obstacleObjectId, pair.dynTri[0], pair.dynTri[1], pair.dynTri[2], pair.obsVertex, weightKey(pair.weight));
  std::sort(keys.begin(), keys.end());
  return keys;
}

static auto canonicalEE(const std::vector<ExternalEEPair> &pairs)
{
  std::vector<std::tuple<int32_t, int, int, int, int, long long>> keys;
  keys.reserve(pairs.size());
  for (const auto &pair : pairs)
    keys.emplace_back(
      pair.obstacleObjectId, pair.dynEdge[0], pair.dynEdge[1], pair.obsEdge[0], pair.obsEdge[1], weightKey(pair.weight));
  std::sort(keys.begin(), keys.end());
  return keys;
}

TEST(SurfaceIPCExternalBroadPhaseGTest, BuilderMatchesSurfaceIPCCoreExternalPairs)
{
  auto [V, F] = makeUnitSquareMesh();
  auto [obsV, obsF] = makeSmallBoxObstacle();
  const ES::VXd obsRest = flattenRows(obsV);

  auto obs = std::make_shared<ObstacleSurface>(
    obsV, obsF,
    pgo::Contact::CIPC::makeLinearTrajectorySampler(obsRest, ES::V3d::Zero()));
  obs->update(0.0, 0.0);

  SurfaceIPCCore::Parameters params;
  params.dhat_external = 1.0;
  SurfaceIPCCore core(params);
  core.setMesh(V, F);
  core.addObstacleSurface(obs);

  const ES::VXd x = flattenRows(V);
  core.prepareForSurfacePositions(x);

  SurfaceIPCTopology topology;
  topology.setMesh(V, F);
  std::vector<std::shared_ptr<ObstacleSurface>> obstacles = { obs };

  ExternalPairSet pairs;
  buildExternalPairs(topology, x, obstacles, params.dhat_external, pairs);

  EXPECT_EQ(canonicalPT(pairs.ptPairs), canonicalPT(core.preparedState().externalPairs.ptPairs));
  EXPECT_EQ(canonicalTP(pairs.tpPairs), canonicalTP(core.preparedState().externalPairs.tpPairs));
  EXPECT_EQ(canonicalEE(pairs.eePairs), canonicalEE(core.preparedState().externalPairs.eePairs));
  EXPECT_GT(pairs.size(), 0u);
}
