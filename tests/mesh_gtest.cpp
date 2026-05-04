#include <gtest/gtest.h>

#include "boundingVolumeTree.h"
#include "geometryQuery.h"
#include "triMeshGeo.h"
#include "triMeshNeighbor.h"

#include <cmath>
#include <algorithm>
#include <vector>

namespace
{
using pgo::Mesh::TriMeshBVTree;
using pgo::Mesh::TriMeshGeo;
using pgo::Mesh::TriMeshNeighbor;
using pgo::Vec3d;
using pgo::Vec3i;

TriMeshGeo makeUnitSquareMesh()
{
  std::vector<Vec3d> vertices{
    Vec3d(0.0, 0.0, 0.0),
    Vec3d(1.0, 0.0, 0.0),
    Vec3d(1.0, 1.0, 0.0),
    Vec3d(0.0, 1.0, 0.0),
  };
  std::vector<Vec3i> triangles{
    Vec3i(0, 1, 2),
    Vec3i(0, 2, 3),
  };
  return TriMeshGeo(std::move(vertices), std::move(triangles));
}
}

TEST(MeshGTest, ComputesBasicTriangleGeometry)
{
  const Vec3d v0(0.0, 0.0, 0.0);
  const Vec3d v1(1.0, 0.0, 0.0);
  const Vec3d v2(0.0, 1.0, 0.0);

  const auto scaledNormal = pgo::Mesh::getTriangleScaledNormal(v0, v1, v2);
  const auto normal = pgo::Mesh::getTriangleNormal(v0, v1, v2);
  const auto center = pgo::Mesh::getTriangleCenterOfMass(v0, v1, v2);
  const auto bary = pgo::Mesh::getBarycentricWeightProjectedOnTrianglePlane(Vec3d(0.25, 0.25, 0.0), v0, v1, v2);

  EXPECT_NEAR(pgo::Mesh::getTriangleArea(v0, v1, v2), 0.5, 1e-12);
  EXPECT_NEAR(scaledNormal[0], 0.0, 1e-12);
  EXPECT_NEAR(scaledNormal[1], 0.0, 1e-12);
  EXPECT_NEAR(scaledNormal[2], 1.0, 1e-12);
  EXPECT_NEAR(normal[0], 0.0, 1e-12);
  EXPECT_NEAR(normal[1], 0.0, 1e-12);
  EXPECT_NEAR(normal[2], 1.0, 1e-12);
  EXPECT_NEAR(center[0], 1.0 / 3.0, 1e-12);
  EXPECT_NEAR(center[1], 1.0 / 3.0, 1e-12);
  EXPECT_NEAR(center[2], 0.0, 1e-12);
  EXPECT_NEAR(bary[0], 0.5, 1e-12);
  EXPECT_NEAR(bary[1], 0.25, 1e-12);
  EXPECT_NEAR(bary[2], 0.25, 1e-12);
}

TEST(MeshGTest, BuildsNeighborInformationForUnitSquare)
{
  const auto mesh = makeUnitSquareMesh();
  TriMeshNeighbor neighbor(mesh);

  EXPECT_EQ(mesh.numVertices(), 4);
  EXPECT_EQ(mesh.numTriangles(), 2);
  EXPECT_NEAR(mesh.ref().computeSurfaceArea(), 1.0, 1e-12);

  const auto tri0Neighbors = neighbor.getTriangleNeighbors(0);
  const auto tri1Neighbors = neighbor.getTriangleNeighbors(1);

  EXPECT_EQ(tri0Neighbors[0], -1);
  EXPECT_EQ(tri0Neighbors[1], -1);
  EXPECT_EQ(tri0Neighbors[2], 1);

  EXPECT_EQ(tri1Neighbors[0], 0);
  EXPECT_EQ(tri1Neighbors[1], -1);
  EXPECT_EQ(tri1Neighbors[2], -1);
}

TEST(MeshGTest, ComputesFastEdgeConnectivityStats)
{
  std::vector<Vec3i> triangles{
    Vec3i(0, 1, 2),
    Vec3i(0, 2, 3),
    Vec3i(4, 5, 6),
    Vec3i(4, 6, 7),
    Vec3i(8, 9, 10),
  };

  const auto components = pgo::Mesh::getConnectedComponentsByEdge(pgo::BasicAlgorithms::makeArrayRef(triangles));
  std::vector<int> expectedCounts;
  for (const auto &component : components)
    expectedCounts.push_back((int)component.size());
  std::sort(expectedCounts.begin(), expectedCounts.end(), std::greater<int>());

  const auto stats = pgo::Mesh::computeTriangleEdgeConnectivityStats(pgo::BasicAlgorithms::makeArrayRef(triangles));
  EXPECT_TRUE(stats.isManifold);
  EXPECT_GT(stats.boundaryOrNonManifoldEdges, 0);
  EXPECT_EQ(stats.componentsByEdge, 3);
  EXPECT_EQ(stats.componentTriangleCountsByEdge, expectedCounts);
}

TEST(MeshGTest, FiltersSmallEdgeConnectedTriangleComponents)
{
  std::vector<Vec3d> vertices{
    Vec3d(0.0, 0.0, 0.0),
    Vec3d(1.0, 0.0, 0.0),
    Vec3d(1.0, 1.0, 0.0),
    Vec3d(0.0, 1.0, 0.0),
    Vec3d(10.0, 0.0, 0.0),
    Vec3d(11.0, 0.0, 0.0),
    Vec3d(10.0, 1.0, 0.0),
  };
  std::vector<Vec3i> triangles{
    Vec3i(0, 1, 2),
    Vec3i(0, 2, 3),
    Vec3i(4, 5, 6),
  };

  const TriMeshGeo mesh(std::move(vertices), std::move(triangles));
  const TriMeshGeo filtered =
    pgo::Mesh::filterSmallTriangleComponentsByEdge(mesh, /*minTriangleCount=*/2);

  EXPECT_EQ(filtered.numVertices(), 4);
  EXPECT_EQ(filtered.numTriangles(), 2);

  const auto stats = pgo::Mesh::computeTriangleEdgeConnectivityStats(filtered.triangles());
  EXPECT_EQ(stats.componentsByEdge, 1);
  ASSERT_EQ(stats.componentTriangleCountsByEdge.size(), 1u);
  EXPECT_EQ(stats.componentTriangleCountsByEdge[0], 2);
}

TEST(MeshGTest, QueriesClosestTriangleWithBVTree)
{
  const auto mesh = makeUnitSquareMesh();
  TriMeshBVTree bvTree;
  bvTree.buildByInertiaPartition(mesh);

  const auto query = Vec3d(0.75, 0.25, 1.0);
  const auto result = bvTree.closestTriangleQuery(mesh, query);

  ASSERT_EQ(result.triID, 0);
  EXPECT_NEAR(result.dist2, 1.0, 1e-12);
  EXPECT_NEAR(result.closestPosition[0], 0.75, 1e-12);
  EXPECT_NEAR(result.closestPosition[1], 0.25, 1e-12);
  EXPECT_NEAR(result.closestPosition[2], 0.0, 1e-12);
  EXPECT_NEAR(result.triBaryWeight.sum(), 1.0, 1e-12);
  EXPECT_NEAR(result.triBaryWeight[0], 0.25, 1e-12);
  EXPECT_NEAR(result.triBaryWeight[1], 0.5, 1e-12);
  EXPECT_NEAR(result.triBaryWeight[2], 0.25, 1e-12);
}
