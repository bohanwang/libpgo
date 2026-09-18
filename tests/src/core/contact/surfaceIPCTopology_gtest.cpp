#include <gtest/gtest.h>

#include "ipc/topology/surfaceIPCTopology.h"

#include "testCIPCHelpers.h"

#include <algorithm>
#include <cstdint>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace
{
namespace ES = pgo::EigenSupport;
using pgo::Contact::CIPC::SurfaceIPCTopology;
using pgo::Contact::CIPCTest::makeTwoTriangleMesh;
}  // namespace

TEST(SurfaceIPCTopologyGTest, SetMeshBuildsEdgesWeightsAndDofs)
{
  const auto [V, F] = makeTwoTriangleMesh();

  SurfaceIPCTopology topology;
  topology.setMesh(V, F);

  ASSERT_EQ(topology.numVerts, 6);
  EXPECT_EQ(topology.numSurfaceDOFs(), 18);
  ASSERT_EQ(topology.triangles.size(), 2u);
  ASSERT_EQ(topology.edges.size(), 6u);
  ASSERT_EQ(topology.vertexArea.size(), 6u);
  ASSERT_EQ(topology.triArea.size(), 2u);
  ASSERT_EQ(topology.edgeLength.size(), 6u);
  ASSERT_EQ(topology.vertexIsDeformable.size(), 6u);
  ASSERT_EQ(topology.triangleHasDeformableVertex.size(), 2u);
  ASSERT_EQ(topology.edgeHasDeformableVertex.size(), 6u);
  EXPECT_TRUE(std::all_of(topology.vertexIsDeformable.begin(), topology.vertexIsDeformable.end(),
    [](uint8_t value) { return value == 1; }));

  EXPECT_EQ(topology.triangles[0], (std::array<int, 3>{ 0, 1, 2 }));
  EXPECT_EQ(topology.edges[0], (std::array<int, 2>{ 0, 1 }));
  EXPECT_EQ(topology.edges[1], (std::array<int, 2>{ 0, 2 }));
  EXPECT_EQ(topology.edges[2], (std::array<int, 2>{ 1, 2 }));
  EXPECT_NEAR(topology.triArea[0], 0.5, 1e-14);

  const double totalTriArea = std::accumulate(topology.triArea.begin(), topology.triArea.end(), 0.0);
  const double totalVertexArea = std::accumulate(topology.vertexArea.begin(), topology.vertexArea.end(), 0.0);
  EXPECT_NEAR(totalVertexArea, totalTriArea, 1e-14);
}

TEST(SurfaceIPCTopologyGTest, MixedMaskBuildsPrimitiveActivity)
{
  const auto [V, F] = makeTwoTriangleMesh();
  const std::vector<uint8_t> mask = { 1, 1, 1, 0, 0, 0 };

  SurfaceIPCTopology topology;
  topology.setMesh(V, F, mask);

  EXPECT_EQ(topology.vertexIsDeformable, mask);
  ASSERT_EQ(topology.triangleHasDeformableVertex.size(), 2u);
  EXPECT_TRUE(topology.triangleContainsDeformableVertex(0));
  EXPECT_FALSE(topology.triangleContainsDeformableVertex(1));
  ASSERT_EQ(topology.edgeHasDeformableVertex.size(), 6u);
  for (int ei = 0; ei < 3; ++ei)
    EXPECT_TRUE(topology.edgeContainsDeformableVertex(ei));
  for (int ei = 3; ei < 6; ++ei)
    EXPECT_FALSE(topology.edgeContainsDeformableVertex(ei));
}

TEST(SurfaceIPCTopologyGTest, RejectsInvalidMaskAndRebuildsActivity)
{
  const auto [V, F] = makeTwoTriangleMesh();
  SurfaceIPCTopology topology;

  EXPECT_THROW(topology.setMesh(V, F, std::vector<uint8_t>{ 1, 1 }), std::invalid_argument);
  EXPECT_THROW(topology.setMesh(V, F, std::vector<uint8_t>{ 1, 1, 1, 0, 0, 2 }), std::invalid_argument);

  topology.setMesh(V, F, std::vector<uint8_t>{ 1, 1, 1, 0, 0, 0 });
  EXPECT_FALSE(topology.triangleContainsDeformableVertex(1));
  topology.setMesh(V, F);
  EXPECT_TRUE(topology.triangleContainsDeformableVertex(0));
  EXPECT_TRUE(topology.triangleContainsDeformableVertex(1));
  EXPECT_TRUE(std::all_of(topology.edgeHasDeformableVertex.begin(), topology.edgeHasDeformableVertex.end(),
    [](uint8_t value) { return value == 1; }));
}
