#include <gtest/gtest.h>

#include "tetgenInterface.h"
#include "triMeshGeo.h"

#include <cmath>

namespace
{
namespace ES = pgo::EigenSupport;

pgo::Mesh::TriMeshGeo makeCubeSurface()
{
  const double vertices[] = { 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1 };
  const int triangles[] = {
    0, 2, 1, 0, 3, 2,  // bottom (z = 0), outward normal -z
    4, 5, 6, 4, 6, 7,  // top
    0, 1, 5, 0, 5, 4,  // front (y = 0)
    1, 2, 6, 1, 6, 5,  // right (x = 1)
    2, 3, 7, 2, 7, 6,  // back (y = 1)
    3, 0, 4, 3, 4, 7   // left (x = 0)
  };
  return pgo::Mesh::TriMeshGeo(8, vertices, 12, triangles);
}
}  // namespace

TEST(TetgenInterfaceGTest, LinearTetrahedralizationYieldsValidConnectivity)
{
  ES::MXd vertices;
  ES::MXi tets;
  ASSERT_EQ(pgo::TetgenInterface::computeTetMesh(makeCubeSurface(), "pq1.414a0.05Q", vertices, tets), 0);
  ASSERT_GT(tets.rows(), 0);
  ASSERT_EQ(tets.cols(), 4);
  EXPECT_GE(tets.minCoeff(), 0);
  EXPECT_LT(tets.maxCoeff(), vertices.rows());

  double volume = 0.0;
  for (int t = 0; t < tets.rows(); ++t) {
    const ES::V3d a = vertices.row(tets(t, 0)), b = vertices.row(tets(t, 1)), c = vertices.row(tets(t, 2)), d = vertices.row(tets(t, 3));
    volume += std::abs((b - a).dot((c - a).cross(d - a))) / 6.0;
  }
  EXPECT_NEAR(volume, 1.0, 1e-9);
}

TEST(TetgenInterfaceGTest, SecondOrderOutputIsRejectedInsteadOfMisread)
{
  ES::MXd vertices;
  ES::MXi tets;
  // With o2 TetGen emits ten node indices per tetrahedron; reading them with a
  // stride of four would silently produce garbage connectivity.
  EXPECT_NE(pgo::TetgenInterface::computeTetMesh(makeCubeSurface(), "pq1.414a0.05o2Q", vertices, tets), 0);
}
