#include <gtest/gtest.h>

#include "ipc/topology/surfaceIPCTopology.h"

#include "testCIPCHelpers.h"

#include <numeric>

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

  EXPECT_EQ(topology.triangles[0], (std::array<int, 3>{ 0, 1, 2 }));
  EXPECT_EQ(topology.edges[0], (std::array<int, 2>{ 0, 1 }));
  EXPECT_EQ(topology.edges[1], (std::array<int, 2>{ 0, 2 }));
  EXPECT_EQ(topology.edges[2], (std::array<int, 2>{ 1, 2 }));
  EXPECT_NEAR(topology.triArea[0], 0.5, 1e-14);

  const double totalTriArea = std::accumulate(topology.triArea.begin(), topology.triArea.end(), 0.0);
  const double totalVertexArea = std::accumulate(topology.vertexArea.begin(), topology.vertexArea.end(), 0.0);
  EXPECT_NEAR(totalVertexArea, totalTriArea, 1e-14);
}
