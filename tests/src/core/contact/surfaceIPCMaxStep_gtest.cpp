#include <gtest/gtest.h>

#include "ipc/core/surfaceIPCMaxStep.h"
#include "ipc/topology/surfaceIPCTopology.h"
#include "ipc/core/surfaceIPCCore.h"
#include "ipc/geometry/ipcDistancePrimitives.h"

#include "testCIPCHelpers.h"

namespace
{
namespace ES = pgo::EigenSupport;
using pgo::Contact::CIPC::SurfaceIPCCore;
using pgo::Contact::CIPC::SurfaceIPCMaxStep;
using pgo::Contact::CIPC::SurfaceIPCTopology;
using pgo::Contact::CIPCTest::flattenPositions;
using pgo::Contact::CIPCTest::makeTwoTriangleMesh;

std::pair<double, double> minNonAdjacentPointTriangleAndEdgeEdgeDistancesSquared(const ES::VXd &x)
{
  const auto vertex = [&](int i) {
    return x.segment<3>(3 * i);
  };
  double minPointTriangle = std::numeric_limits<double>::infinity();
  for (int p = 0; p < 3; ++p)
    minPointTriangle = std::min(minPointTriangle, pgo::Contact::CIPC::distance::computePTSqDist(vertex(p), vertex(3), vertex(4), vertex(5)));
  for (int p = 3; p < 6; ++p)
    minPointTriangle = std::min(minPointTriangle, pgo::Contact::CIPC::distance::computePTSqDist(vertex(p), vertex(0), vertex(1), vertex(2)));

  constexpr int edges[3][2] = { { 0, 1 }, { 1, 2 }, { 2, 0 } };
  double minEdgeEdge = std::numeric_limits<double>::infinity();
  for (const auto &edgeA : edges) {
    for (const auto &edgeB : edges) {
      minEdgeEdge = std::min(minEdgeEdge, pgo::Contact::CIPC::distance::computeEESqDist(vertex(edgeA[0]), vertex(edgeA[1]), vertex(3 + edgeB[0]), vertex(3 + edgeB[1])));
    }
  }
  return { minPointTriangle, minEdgeEdge };
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

  const double helperAlpha = SurfaceIPCMaxStep().compute(topology, x, dx, params.dhat, params.slackness);
  const double coreAlpha = core.computeMaxStepSize(x, dx);

  EXPECT_GT(helperAlpha, 0.0);
  EXPECT_LT(helperAlpha, 1.0);
  EXPECT_NEAR(helperAlpha, coreAlpha, 1e-12);
}

TEST(SurfaceIPCMaxStepGTest, CCDStepKeepsPointTriangleAndEdgeEdgeDistancesPositive)
{
  const auto [V, _] = makeTwoTriangleMesh();
  const ES::VXd x = flattenPositions(V);
  ES::VXd dx = ES::VXd::Zero(x.size());
  for (int vi = 3; vi < 6; ++vi)
    dx[3 * vi + 2] = -0.1;

  SurfaceIPCTopology topology;
  const auto [meshVertices, meshTriangles] = makeTwoTriangleMesh();
  topology.setMesh(meshVertices, meshTriangles);
  constexpr double kDhat = 0.1;
  constexpr double kSlackness = 0.9;
  const double alpha = SurfaceIPCMaxStep().compute(topology, x, dx, kDhat, kSlackness);
  ASSERT_GT(alpha, 0.0);
  ASSERT_LT(alpha, 1.0);

  const ES::VXd safeState = x + alpha * dx;
  const ES::VXd collisionState = x + (alpha / kSlackness) * dx;
  const auto [safePointTriangle, safeEdgeEdge] = minNonAdjacentPointTriangleAndEdgeEdgeDistancesSquared(safeState);
  const auto [unscaledPointTriangle, unscaledEdgeEdge] = minNonAdjacentPointTriangleAndEdgeEdgeDistancesSquared(collisionState);
  EXPECT_GT(safePointTriangle, 1e-12);
  EXPECT_GT(safeEdgeEdge, 1e-12);
  EXPECT_LT(std::min(unscaledPointTriangle, unscaledEdgeEdge), std::min(safePointTriangle, safeEdgeEdge));
}

TEST(SurfaceIPCMaxStepGTest, DynamicTriangleIsLimitedByStaticTriangle)
{
  const auto [V, F] = makeTwoTriangleMesh();
  const ES::VXd x = flattenPositions(V);
  ES::VXd dx = ES::VXd::Zero(x.size());
  for (int vi = 3; vi < 6; ++vi)
    dx[3 * vi + 2] = -0.1;

  constexpr double kDhat = 0.1;
  constexpr double kSlackness = 0.9;
  const std::vector<uint8_t> mixedMask = { 0, 0, 0, 1, 1, 1 };
  SurfaceIPCTopology topology;
  topology.setMesh(V, F, mixedMask);
  const double helperAlpha = SurfaceIPCMaxStep().compute(topology, x, dx, kDhat, kSlackness);
  EXPECT_GT(helperAlpha, 0.0);
  EXPECT_LT(helperAlpha, 1.0);

  SurfaceIPCCore::Parameters params;
  params.dhat = kDhat;
  params.kappa = 1.0;
  params.slackness = kSlackness;
  SurfaceIPCCore core(params);
  core.setMesh(V, F, mixedMask);
  EXPECT_NEAR(core.computeMaxStepSize(x, dx), helperAlpha, 1e-12);
}

TEST(SurfaceIPCMaxStepGTest, FullyExternalGeometryDoesNotLimitStep)
{
  const auto [V, F] = makeTwoTriangleMesh();
  const ES::VXd x = flattenPositions(V);
  ES::VXd dx = ES::VXd::Zero(x.size());
  for (int vi = 3; vi < 6; ++vi)
    dx[3 * vi + 2] = -0.1;

  SurfaceIPCTopology topology;
  topology.setMesh(V, F, std::vector<uint8_t>(6, 0));
  EXPECT_DOUBLE_EQ(SurfaceIPCMaxStep().compute(topology, x, dx, 0.1, 0.9), 1.0);
}
