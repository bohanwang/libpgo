#include <gtest/gtest.h>

#include "field/gridSpec.h"
#include "field/denseGrid.h"
#include "operations/booleanOps.h"
#include "extraction/marchingCubesExtractor.h"
#include "geometry/meshDistance.h"
#include "geometry/shellThickening.h"
#include "geometry/sphereField.h"

#ifdef PGO_HAS_OPENVDB
#include "extraction/openVDBExtractor.h"
#include <openvdb/openvdb.h>
#endif

#include <cmath>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <vector>
#include <cstdio>

namespace IS = pgo::ImplicitSurface;
namespace ES = pgo::EigenSupport;

// =========================================================================
// GridSpec
// =========================================================================

TEST(GridSpecTest, ValidGridSpec)
{
  IS::GridSpec spec;
  spec.bmin = ES::V3d(0.0, 0.0, 0.0);
  spec.bmax = ES::V3d(1.0, 1.0, 1.0);
  spec.resolution = 4;
  IS::validateGridSpec(spec);
  EXPECT_EQ(spec, spec);
}

TEST(GridSpecTest, ResolutionTooSmall)
{
  IS::GridSpec spec;
  spec.bmin = ES::V3d(0.0, 0.0, 0.0);
  spec.bmax = ES::V3d(1.0, 1.0, 1.0);
  spec.resolution = 1;
  EXPECT_THROW(IS::validateGridSpec(spec), std::runtime_error);
}

TEST(GridSpecTest, BmaxNotGreaterThanBmin)
{
  IS::GridSpec spec;
  spec.bmin = ES::V3d(0.0, 0.0, 0.0);
  spec.bmax = ES::V3d(0.0, 1.0, 1.0);
  spec.resolution = 4;
  EXPECT_THROW(IS::validateGridSpec(spec), std::runtime_error);
}

TEST(GridSpecTest, NonFiniteBounds)
{
  IS::GridSpec spec;
  spec.bmin = ES::V3d(0.0, 0.0, 0.0);
  spec.bmax = ES::V3d(std::numeric_limits<double>::infinity(), 1.0, 1.0);
  spec.resolution = 4;
  EXPECT_THROW(IS::validateGridSpec(spec), std::runtime_error);
}

TEST(GridSpecTest, DifferentSpecsNotEqual)
{
  IS::GridSpec a;
  a.bmin = ES::V3d(0.0, 0.0, 0.0);
  a.bmax = ES::V3d(1.0, 1.0, 1.0);
  a.resolution = 4;

  IS::GridSpec b = a;
  b.resolution = 8;
  EXPECT_FALSE(a == b);
}

// =========================================================================
// DenseGrid
// =========================================================================

TEST(DenseGridTest, ConstructAndAccess)
{
  IS::GridSpec spec;
  spec.bmin = ES::V3d(0.0, 0.0, 0.0);
  spec.bmax = ES::V3d(1.0, 1.0, 1.0);
  spec.resolution = 4;

  IS::DenseGrid grid(spec);
  EXPECT_EQ(grid.resolution(), 4);
  EXPECT_EQ(grid.size(), 64);

  grid.fill(0.0);
  EXPECT_DOUBLE_EQ(grid[0], 0.0);

  grid.at(0, 0, 0) = 1.0;
  EXPECT_DOUBLE_EQ(grid.at(0, 0, 0), 1.0);
  EXPECT_DOUBLE_EQ(grid[0], 1.0);
  EXPECT_DOUBLE_EQ(grid.at(1, 0, 0), 0.0);
  EXPECT_DOUBLE_EQ(grid.at(0, 1, 0), 0.0);
  EXPECT_DOUBLE_EQ(grid.at(0, 0, 1), 0.0);
}

TEST(DenseGridTest, InvalidGridSpecThrows)
{
  IS::GridSpec spec;
  spec.bmin = ES::V3d(0.0, 0.0, 0.0);
  spec.bmax = ES::V3d(1.0, 1.0, 1.0);
  spec.resolution = 1;
  EXPECT_THROW(IS::DenseGrid grid(spec), std::runtime_error);
}

TEST(DenseGridTest, LinearIndexOrder)
{
  IS::GridSpec spec;
  spec.bmin = ES::V3d(0.0, 0.0, 0.0);
  spec.bmax = ES::V3d(1.0, 1.0, 1.0);
  spec.resolution = 4;

  IS::DenseGrid grid(spec);
  for (int z = 0; z < 4; ++z)
    for (int y = 0; y < 4; ++y)
      for (int x = 0; x < 4; ++x)
        grid.at(x, y, z) = static_cast<double>(IS::linearIndex(x, y, z, 4));

  EXPECT_DOUBLE_EQ(grid[0], 0.0);
  EXPECT_DOUBLE_EQ(grid[1], 1.0);
  EXPECT_DOUBLE_EQ(grid[4], 4.0);
  EXPECT_DOUBLE_EQ(grid[16], 16.0);
  EXPECT_DOUBLE_EQ(grid[63], 63.0);
}

TEST(DenseGridTest, SetZeroClearsAllValues)
{
  IS::GridSpec spec;
  spec.bmin = ES::V3d(0.0, 0.0, 0.0);
  spec.bmax = ES::V3d(1.0, 1.0, 1.0);
  spec.resolution = 4;

  IS::DenseGrid grid(spec);
  grid.fill(3.14);
  grid.setZero();
  for (int i = 0; i < grid.size(); ++i)
    EXPECT_DOUBLE_EQ(grid[i], 0.0);
}

TEST(DenseGridTest, DataPointerAccess)
{
  IS::GridSpec spec;
  spec.bmin = ES::V3d(0.0, 0.0, 0.0);
  spec.bmax = ES::V3d(1.0, 1.0, 1.0);
  spec.resolution = 3;

  IS::DenseGrid grid(spec);
  grid.fill(7.0);
  const double *data = grid.data();
  EXPECT_DOUBLE_EQ(data[0], 7.0);
  EXPECT_DOUBLE_EQ(data[grid.size() - 1], 7.0);

  double *mutableData = grid.data();
  mutableData[0] = 42.0;
  EXPECT_DOUBLE_EQ(grid[0], 42.0);
}

// =========================================================================
// BooleanOps
// =========================================================================

TEST(BooleanOpsTest, Union)
{
  IS::GridSpec spec;
  spec.bmin = ES::V3d(0.0, 0.0, 0.0);
  spec.bmax = ES::V3d(1.0, 1.0, 1.0);
  spec.resolution = 3;

  IS::DenseGrid a(spec), b(spec), out(spec);
  a.fill(-1.0);
  b.fill(1.0);
  IS::applyBoolean(a, b, IS::BooleanOp::Union, out);
  for (int i = 0; i < out.size(); ++i)
    EXPECT_DOUBLE_EQ(out[i], -1.0);
}

TEST(BooleanOpsTest, Intersection)
{
  IS::GridSpec spec;
  spec.bmin = ES::V3d(0.0, 0.0, 0.0);
  spec.bmax = ES::V3d(1.0, 1.0, 1.0);
  spec.resolution = 3;

  IS::DenseGrid a(spec), b(spec), out(spec);
  a.fill(-1.0);
  b.fill(1.0);
  IS::applyBoolean(a, b, IS::BooleanOp::Intersection, out);
  for (int i = 0; i < out.size(); ++i)
    EXPECT_DOUBLE_EQ(out[i], 1.0);
}

TEST(BooleanOpsTest, Difference)
{
  IS::GridSpec spec;
  spec.bmin = ES::V3d(0.0, 0.0, 0.0);
  spec.bmax = ES::V3d(1.0, 1.0, 1.0);
  spec.resolution = 3;

  IS::DenseGrid a(spec), b(spec), out(spec);
  a.fill(-1.0);
  b.fill(1.0);
  IS::applyBoolean(a, b, IS::BooleanOp::Difference, out);
  for (int i = 0; i < out.size(); ++i)
    EXPECT_DOUBLE_EQ(out[i], -1.0);
}

TEST(BooleanOpsTest, MismatchedSpecThrows)
{
  IS::GridSpec spec1;
  spec1.bmin = ES::V3d(0.0, 0.0, 0.0);
  spec1.bmax = ES::V3d(1.0, 1.0, 1.0);
  spec1.resolution = 3;

  IS::GridSpec spec2;
  spec2.bmin = ES::V3d(0.0, 0.0, 0.0);
  spec2.bmax = ES::V3d(2.0, 1.0, 1.0);
  spec2.resolution = 3;

  IS::DenseGrid a(spec1), b(spec2), out(spec1);
  EXPECT_THROW(IS::applyBoolean(a, b, IS::BooleanOp::Union, out), std::runtime_error);
}

TEST(BooleanOpsTest, UnionPointwiseCheck)
{
  IS::GridSpec spec;
  spec.bmin = ES::V3d(0.0, 0.0, 0.0);
  spec.bmax = ES::V3d(1.0, 1.0, 1.0);
  spec.resolution = 3;

  IS::DenseGrid a(spec), b(spec), out(spec);
  a.fill(0.0);
  b.fill(0.0);

  // a: alternating -1 and 1, b: all 0 → Union = min
  for (int i = 0; i < a.size(); ++i)
    a[i] = (i % 2 == 0) ? -1.0 : 1.0;

  IS::applyBoolean(a, b, IS::BooleanOp::Union, out);
  EXPECT_DOUBLE_EQ(out[0], -1.0);   // min(-1, 0)
  EXPECT_DOUBLE_EQ(out[1], 0.0);    // min(1, 0)
}

// =========================================================================
// MarchingCubes
// =========================================================================

TEST(MarchingCubesTest, SphereFieldExtraction)
{
  const int res = 32;
  IS::GridSpec spec;
  spec.bmin = ES::V3d(-1.5, -1.5, -1.5);
  spec.bmax = ES::V3d(1.5, 1.5, 1.5);
  spec.resolution = res;

  IS::DenseGrid field(spec);
  const ES::V3d center(0.0, 0.0, 0.0);
  const double radius = 1.0;

  const ES::V3d delta = (spec.bmax - spec.bmin) / static_cast<double>(res - 1);
  for (int z = 0; z < res; ++z)
    for (int y = 0; y < res; ++y)
      for (int x = 0; x < res; ++x) {
        const ES::V3d p = spec.bmin + delta.cwiseProduct(ES::V3d(x, y, z).cast<double>());
        field.at(x, y, z) = (p - center).norm() - radius;
      }

  IS::MarchingCubesOptions mcOpts;
  mcOpts.isoOffset = 0.0;
  pgo::Mesh::TriMeshGeo outMesh;
  IS::extractMarchingCubes(field, mcOpts, outMesh);

  EXPECT_GT(outMesh.numVertices(), 0);
  EXPECT_GT(outMesh.numTriangles(), 0);
}

TEST(MarchingCubesTest, NonZeroIsoOffsetShiftsSurface)
{
  const int res = 32;
  IS::GridSpec spec;
  spec.bmin = ES::V3d(-1.5, -1.5, -1.5);
  spec.bmax = ES::V3d(1.5, 1.5, 1.5);
  spec.resolution = res;

  IS::DenseGrid field(spec);
  const ES::V3d center(0.0, 0.0, 0.0);
  const double radius = 1.0;

  const ES::V3d delta = (spec.bmax - spec.bmin) / static_cast<double>(res - 1);
  for (int z = 0; z < res; ++z)
    for (int y = 0; y < res; ++y)
      for (int x = 0; x < res; ++x) {
        const ES::V3d p = spec.bmin + delta.cwiseProduct(ES::V3d(x, y, z).cast<double>());
        field.at(x, y, z) = (p - center).norm() - radius;
      }

  pgo::Mesh::TriMeshGeo mesh0, meshOffset;
  IS::MarchingCubesOptions mcOpts0;
  mcOpts0.isoOffset = 0.0;
  IS::extractMarchingCubes(field, mcOpts0, mesh0);

  IS::MarchingCubesOptions mcOptsOffset;
  mcOptsOffset.isoOffset = 0.1;
  IS::extractMarchingCubes(field, mcOptsOffset, meshOffset);

  // Both should produce non-empty meshes
  EXPECT_GT(mesh0.numTriangles(), 0);
  EXPECT_GT(meshOffset.numTriangles(), 0);
  // Different iso offsets should produce different meshes
  EXPECT_NE(mesh0.numVertices(), meshOffset.numVertices());
}

// =========================================================================
// MeshVolume
// =========================================================================

TEST(MeshVolumeTest, UnitCubeVolume)
{
  std::vector<ES::V3d> positions = {
    ES::V3d(0.0, 0.0, 0.0), ES::V3d(1.0, 0.0, 0.0),
    ES::V3d(0.0, 1.0, 0.0), ES::V3d(0.0, 0.0, 1.0),
    ES::V3d(1.0, 1.0, 0.0), ES::V3d(1.0, 0.0, 1.0),
    ES::V3d(0.0, 1.0, 1.0), ES::V3d(1.0, 1.0, 1.0),
  };
  std::vector<ES::V3i> triangles = {
    ES::V3i(0, 2, 1), ES::V3i(1, 2, 4),  // front
    ES::V3i(0, 1, 3), ES::V3i(1, 5, 3),  // bottom
    ES::V3i(0, 3, 2), ES::V3i(2, 3, 6),  // left
    ES::V3i(1, 4, 5), ES::V3i(4, 7, 5),  // right
    ES::V3i(2, 6, 4), ES::V3i(4, 6, 7),  // top
    ES::V3i(3, 5, 6), ES::V3i(5, 7, 6),  // back
  };
  pgo::Mesh::TriMeshGeo cube(std::move(positions), std::move(triangles));
  EXPECT_NEAR(pgo::Mesh::computeMeshVolume(cube), 1.0, 1e-10);
}

TEST(MeshVolumeTest, TranslatedCubeSameVolume)
{
  std::vector<ES::V3d> positions = {
    ES::V3d(5.0, 5.0, 5.0), ES::V3d(6.0, 5.0, 5.0),
    ES::V3d(5.0, 6.0, 5.0), ES::V3d(5.0, 5.0, 6.0),
    ES::V3d(6.0, 6.0, 5.0), ES::V3d(6.0, 5.0, 6.0),
    ES::V3d(5.0, 6.0, 6.0), ES::V3d(6.0, 6.0, 6.0),
  };
  std::vector<ES::V3i> triangles = {
    ES::V3i(0, 2, 1), ES::V3i(1, 2, 4),
    ES::V3i(0, 1, 3), ES::V3i(1, 5, 3),
    ES::V3i(0, 3, 2), ES::V3i(2, 3, 6),
    ES::V3i(1, 4, 5), ES::V3i(4, 7, 5),
    ES::V3i(2, 6, 4), ES::V3i(4, 6, 7),
    ES::V3i(3, 5, 6), ES::V3i(5, 7, 6),
  };
  pgo::Mesh::TriMeshGeo cube(std::move(positions), std::move(triangles));
  EXPECT_NEAR(pgo::Mesh::computeMeshVolume(cube), 1.0, 1e-10);
}

TEST(MeshVolumeTest, ScaledCubeVolume)
{
  std::vector<ES::V3d> positions = {
    ES::V3d(0.0, 0.0, 0.0), ES::V3d(2.0, 0.0, 0.0),
    ES::V3d(0.0, 3.0, 0.0), ES::V3d(0.0, 0.0, 4.0),
    ES::V3d(2.0, 3.0, 0.0), ES::V3d(2.0, 0.0, 4.0),
    ES::V3d(0.0, 3.0, 4.0), ES::V3d(2.0, 3.0, 4.0),
  };
  std::vector<ES::V3i> triangles = {
    ES::V3i(0, 2, 1), ES::V3i(1, 2, 4),
    ES::V3i(0, 1, 3), ES::V3i(1, 5, 3),
    ES::V3i(0, 3, 2), ES::V3i(2, 3, 6),
    ES::V3i(1, 4, 5), ES::V3i(4, 7, 5),
    ES::V3i(2, 6, 4), ES::V3i(4, 6, 7),
    ES::V3i(3, 5, 6), ES::V3i(5, 7, 6),
  };
  pgo::Mesh::TriMeshGeo box(std::move(positions), std::move(triangles));
  EXPECT_NEAR(pgo::Mesh::computeMeshVolume(box), 24.0, 1e-10);
}

// =========================================================================
// ShellThickening (thickenMeshShell) + sphereField evaluation
// =========================================================================

TEST(ShellThickeningTest, ThickenMeshShellProducesNegativeInside)
{
  const int res = 16;
  IS::GridSpec spec;
  spec.bmin = ES::V3d(-2.0, -2.0, -2.0);
  spec.bmax = ES::V3d(2.0, 2.0, 2.0);
  spec.resolution = res;

  IS::DenseGrid distance(spec);
  IS::DenseGrid outField(spec);

  const ES::V3d delta = (spec.bmax - spec.bmin) / static_cast<double>(res - 1);
  for (int z = 0; z < res; ++z)
    for (int y = 0; y < res; ++y)
      for (int x = 0; x < res; ++x) {
        const ES::V3d p = spec.bmin + delta.cwiseProduct(ES::V3d(x, y, z).cast<double>());
        distance.at(x, y, z) = p.norm();
      }

  IS::thickenMeshShell(distance, 2.0, outField);

  // With thickness=2.0, half=1.0. Near origin distance≈0 → shell≈-1 (inside)
  int negativeCount = 0;
  for (int i = 0; i < outField.size(); ++i)
    if (outField[i] < 0.0) ++negativeCount;
  EXPECT_GT(negativeCount, 0);
}

TEST(SphereFieldTest, ThickenSphereShellProducesShellAroundSphere)
{
  const int res = 16;
  IS::GridSpec spec;
  spec.bmin = ES::V3d(-2.0, -2.0, -2.0);
  spec.bmax = ES::V3d(2.0, 2.0, 2.0);
  spec.resolution = res;

  IS::DenseGrid outField(spec);
  IS::SphereField sphere;
  sphere.center = ES::V3d(0.0, 0.0, 0.0);
  sphere.radius = 1.0;

  IS::thickenSphereShell(sphere, 0.2, spec, outField);

  // At grid center (8,8,8): p≈(0.133,0.133,0.133), |p|≈0.231
  // sphereShell = |0.231-1.0|-0.1 = 0.769-0.1 = 0.669
  EXPECT_NEAR(outField.at(8, 8, 8), 0.669, 0.01);
}

TEST(SphereFieldTest, EvaluateBallSDFIsNegativeInside)
{
  const int res = 16;
  IS::GridSpec spec;
  spec.bmin = ES::V3d(-2.0, -2.0, -2.0);
  spec.bmax = ES::V3d(2.0, 2.0, 2.0);
  spec.resolution = res;

  IS::DenseGrid outField(spec);
  IS::SphereField sphere;
  sphere.center = ES::V3d(0.0, 0.0, 0.0);
  sphere.radius = 1.0;

  IS::evaluateBallSDF(sphere, spec, outField);

  // Center (r=0): 0 - 1 = -1 (inside ball)
  EXPECT_LT(outField.at(8, 8, 8), 0.0);
  // Far outside: should be positive
  EXPECT_GT(outField.at(0, 0, 0), 0.0);
}

TEST(SphereFieldTest, TruncationViaBooleanIntersection)
{
  const int res = 32;
  IS::GridSpec spec;
  spec.bmin = ES::V3d(-2.0, -2.0, -2.0);
  spec.bmax = ES::V3d(2.0, 2.0, 2.0);
  spec.resolution = res;

  // Create a synthetic mesh distance field (zero at center, grows outward)
  IS::DenseGrid distance(spec);
  const ES::V3d delta = (spec.bmax - spec.bmin) / static_cast<double>(res - 1);
  for (int z = 0; z < res; ++z)
    for (int y = 0; y < res; ++y)
      for (int x = 0; x < res; ++x) {
        const ES::V3d p = spec.bmin + delta.cwiseProduct(ES::V3d(x, y, z).cast<double>());
        distance.at(x, y, z) = p.norm();
      }

  IS::SphereField sphere;
  sphere.center = ES::V3d(0.0, 0.0, 0.0);
  sphere.radius = 1.0;

  // Build mesh shell
  IS::DenseGrid shellGrid(spec);
  IS::thickenMeshShell(distance, 0.2, shellGrid);

  // Build ball SDF
  IS::DenseGrid ballGrid(spec);
  IS::evaluateBallSDF(sphere, spec, ballGrid);

  // Truncate: shell ∩ ball
  IS::DenseGrid truncated(spec);
  IS::applyBoolean(shellGrid, ballGrid, IS::BooleanOp::Intersection, truncated);

  // Far outside the ball (>1.5), truncated field should be positive
  int farPositiveCount = 0;
  for (int z = 0; z < res; ++z)
    for (int y = 0; y < res; ++y)
      for (int x = 0; x < res; ++x) {
        const ES::V3d p = spec.bmin + delta.cwiseProduct(ES::V3d(x, y, z).cast<double>());
        if (p.norm() > 1.5 && truncated.at(x, y, z) > 0.0)
          ++farPositiveCount;
      }
  EXPECT_GT(farPositiveCount, 0);
}

// =========================================================================
// SphereField
// =========================================================================

TEST(SphereFieldTest, ComputeFromBBox)
{
  std::vector<ES::V3d> positions = {
    ES::V3d(1.0, 0.0, 0.0), ES::V3d(0.0, 1.0, 0.0),
    ES::V3d(0.0, 0.0, 1.0), ES::V3d(-1.0, 0.0, 0.0),
    ES::V3d(0.0, -1.0, 0.0), ES::V3d(0.0, 0.0, -1.0),
  };
  std::vector<ES::V3i> triangles = {
    ES::V3i(0, 1, 2), ES::V3i(0, 2, 3),
    ES::V3i(0, 3, 4), ES::V3i(0, 4, 1),
  };
  pgo::Mesh::TriMeshGeo mesh(std::move(positions), std::move(triangles));

  IS::SphereField sphere = IS::computeSphereFieldFromBBox(mesh);
  EXPECT_NEAR(sphere.center[0], 0.0, 0.5);
  EXPECT_NEAR(sphere.center[1], 0.0, 0.5);
  EXPECT_NEAR(sphere.center[2], 0.0, 0.5);
  EXPECT_NEAR(sphere.radius, 1.0, 0.2);
}

TEST(SphereFieldTest, EmptyMeshThrows)
{
  pgo::Mesh::TriMeshGeo mesh;
  EXPECT_THROW(IS::computeSphereFieldFromBBox(mesh), std::runtime_error);
}

TEST(SphereFieldTest, ProjectOpenBoundaryToSphere)
{
  std::vector<ES::V3d> positions = {
    ES::V3d(0.5, 0.0, 0.0),   // 0: interior, edge(0,1) and (0,2) are boundary
    ES::V3d(2.0, 0.0, 0.0),   // 1: boundary
    ES::V3d(0.0, 2.0, 0.0),   // 2: boundary
    ES::V3d(0.5, 0.5, 0.0),   // 3: interior shared by two tris
  };
  std::vector<ES::V3i> triangles = {
    ES::V3i(0, 1, 3),
    ES::V3i(0, 3, 2),
  };
  pgo::Mesh::TriMeshGeo mesh(std::move(positions), std::move(triangles));

  IS::SphereField sphere;
  sphere.center = ES::V3d(0.0, 0.0, 0.0);
  sphere.radius = 1.0;

  const int n = IS::projectOpenBoundaryToSphere(mesh, sphere);
  EXPECT_GT(n, 0);

  // Boundary vertices (1 and 2) should now be on the sphere
  EXPECT_NEAR(mesh.pos(1).norm(), 1.0, 1e-10);
  EXPECT_NEAR(mesh.pos(2).norm(), 1.0, 1e-10);
}

// =========================================================================
// ComponentFilter
// =========================================================================

TEST(ComponentFilterTest, SingleLargeComponentIsKept)
{
  std::vector<ES::V3d> positions = {
    ES::V3d(0.0, 0.0, 0.0), ES::V3d(1.0, 0.0, 0.0), ES::V3d(0.0, 1.0, 0.0),
    ES::V3d(0.5, 0.5, 0.0), ES::V3d(1.0, 1.0, 0.0), ES::V3d(1.5, 0.5, 0.0),
  };
  std::vector<ES::V3i> triangles = {
    ES::V3i(0, 1, 2),
    ES::V3i(3, 4, 5),
  };
  pgo::Mesh::TriMeshGeo mesh(std::move(positions), std::move(triangles));

  pgo::Mesh::filterSmallComponents(mesh, /*minComponentTriangles=*/1, /*keepLargestComponents=*/-1);
  EXPECT_EQ(mesh.numTriangles(), 2);
}

TEST(ComponentFilterTest, TinyComponentsAreRemoved)
{
  std::vector<ES::V3d> positions = {
    ES::V3d(0.0, 0.0, 0.0), ES::V3d(1.0, 0.0, 0.0),
    ES::V3d(0.0, 1.0, 0.0), ES::V3d(1.0, 1.0, 0.0),
    ES::V3d(100.0, 100.0, 0.0), ES::V3d(101.0, 100.0, 0.0),
    ES::V3d(100.0, 101.0, 0.0),
  };
  std::vector<ES::V3i> triangles = {
    ES::V3i(0, 1, 2), ES::V3i(1, 3, 2),
    ES::V3i(4, 5, 6),
  };
  pgo::Mesh::TriMeshGeo mesh(std::move(positions), std::move(triangles));

  pgo::Mesh::filterSmallComponents(mesh, /*minComponentTriangles=*/2, /*keepLargestComponents=*/-1);
  EXPECT_EQ(mesh.numTriangles(), 2);
}

// =========================================================================
// MeshDistance (basic smoke tests — full coverage needs real meshes)
// =========================================================================

TEST(MeshDistanceTest, UnionBBoxProducesValidBounds)
{
  std::vector<ES::V3d> posA = {
    ES::V3d(0.0, 0.0, 0.0), ES::V3d(1.0, 0.0, 0.0), ES::V3d(0.0, 1.0, 0.0),
  };
  std::vector<ES::V3i> triA = { ES::V3i(0, 1, 2) };
  pgo::Mesh::TriMeshGeo meshA(std::move(posA), std::move(triA));

  std::vector<ES::V3d> posB = {
    ES::V3d(2.0, 2.0, 2.0), ES::V3d(3.0, 2.0, 2.0), ES::V3d(2.0, 3.0, 2.0),
  };
  std::vector<ES::V3i> triB = { ES::V3i(0, 1, 2) };
  pgo::Mesh::TriMeshGeo meshB(std::move(posB), std::move(triB));

  ES::V3d bmin, bmax;
  pgo::Mesh::computeUnionBBox(meshA, meshB, 0.1, 0.1, 0.05, bmin, bmax);

  // bbox should cover both meshes + expansion
  EXPECT_LT(bmin[0], 0.0);  // expanded below meshA min x (0.0)
  EXPECT_GT(bmax[0], 3.0);  // expanded above meshB max x (3.0)
  EXPECT_GT(bmax[0] - bmin[0], 0.0);  // positive side lengths

  // Also validate as a GridSpec
  IS::GridSpec spec;
  spec.bmin = bmin;
  spec.bmax = bmax;
  spec.resolution = 4;
  IS::validateGridSpec(spec);  // should not throw
}

// =========================================================================
// OpenVDB Extraction (conditional on PGO_HAS_OPENVDB)
// =========================================================================

#ifdef PGO_HAS_OPENVDB

namespace {

pgo::Mesh::TriMeshGeo makeUnitSphereMesh(int uSteps = 8, int vSteps = 8)
{
  std::vector<ES::V3d> positions;
  std::vector<ES::V3i> triangles;
  const double pi = 3.14159265358979323846;
  // latitude rings from pole to pole
  for (int j = 0; j <= vSteps; ++j) {
    const double theta = pi * static_cast<double>(j) / static_cast<double>(vSteps);
    const double sinTheta = std::sin(theta);
    const double cosTheta = std::cos(theta);
    const int nLong = (j == 0 || j == vSteps) ? 1 : uSteps;  // single vertex at poles
    for (int i = 0; i < nLong; ++i) {
      const double phi = 2.0 * pi * static_cast<double>(i) / static_cast<double>(uSteps);
      positions.push_back(ES::V3d(std::cos(phi) * sinTheta, cosTheta, std::sin(phi) * sinTheta));
    }
  }
  // Build triangles
  int prevRingStart = 0;
  for (int j = 1; j <= vSteps; ++j) {
    const bool prevIsPole = (j - 1 == 0);
    const bool currIsPole = (j == vSteps);
    const int prevRingSize = prevIsPole ? 1 : uSteps;

    for (int i = 0; i < (prevIsPole ? uSteps : prevRingSize); ++i) {
      const int pi0 = prevRingStart + (prevIsPole ? 0 : i);
      const int pi1 = prevRingStart + (prevIsPole ? 0 : (i + 1) % uSteps);
      const int ci0 = prevRingStart + prevRingSize + (currIsPole ? 0 : i);
      const int ci1 = prevRingStart + prevRingSize + (currIsPole ? 0 : (i + 1) % uSteps);

      if (!prevIsPole)
        triangles.push_back(ES::V3i(pi0, pi1, ci0));
      if (!currIsPole)
        triangles.push_back(ES::V3i(pi1, ci1, ci0));
    }
    prevRingStart += prevRingSize;
  }
  return pgo::Mesh::TriMeshGeo(std::move(positions), std::move(triangles));
}

}  // namespace

TEST(OpenVDBExtractorTest, SphereShellMeshIsNonEmpty)
{
  openvdb::initialize();

  IS::SphereField sphere;
  sphere.center = ES::V3d(0.0, 0.0, 0.0);
  sphere.radius = 1.0;

  IS::OpenVDBOptions opts;
  opts.voxelSize = 0.05;
  opts.halfWidth = 3.0;
  opts.adaptivity = 0.0;
  opts.smoothSteps = 0;
  IS::validateOpenVDBOptions(opts);

  auto shell = IS::buildOpenVDBSphereShell(sphere, 0.2, opts);
  ASSERT_NE(shell, nullptr);
  EXPECT_NE(shell->grid, nullptr);

  pgo::Mesh::TriMeshGeo outMesh;
  IS::extractOpenVDBLevelSet(*shell, opts, outMesh);
  EXPECT_GT(outMesh.numVertices(), 0);
  EXPECT_GT(outMesh.numTriangles(), 0);
}

TEST(OpenVDBExtractorTest, BallLevelSetMeshIsNonEmpty)
{
  openvdb::initialize();

  IS::SphereField sphere;
  sphere.center = ES::V3d(0.0, 0.0, 0.0);
  sphere.radius = 0.5;

  IS::OpenVDBOptions opts;
  opts.voxelSize = 0.05;
  opts.halfWidth = 3.0;
  IS::validateOpenVDBOptions(opts);

  auto ball = IS::buildOpenVDBBallLevelSet(sphere, opts);
  ASSERT_NE(ball, nullptr);

  pgo::Mesh::TriMeshGeo outMesh;
  IS::extractOpenVDBLevelSet(*ball, opts, outMesh);
  EXPECT_GT(outMesh.numVertices(), 0);
  EXPECT_GT(outMesh.numTriangles(), 0);
}

TEST(OpenVDBExtractorTest, MeshShellExtractionIsNonEmpty)
{
  openvdb::initialize();

  pgo::Mesh::TriMeshGeo sphereMesh = makeUnitSphereMesh(16, 16);

  IS::OpenVDBOptions opts;
  opts.voxelSize = 0.05;
  opts.halfWidth = 3.0;
  IS::validateOpenVDBOptions(opts);

  auto shell = IS::buildOpenVDBShellFromMesh(sphereMesh, 0.1, opts);
  ASSERT_NE(shell, nullptr);

  pgo::Mesh::TriMeshGeo outMesh;
  IS::extractOpenVDBLevelSet(*shell, opts, outMesh);
  EXPECT_GT(outMesh.numVertices(), 0);
  EXPECT_GT(outMesh.numTriangles(), 0);
}

TEST(OpenVDBExtractorTest, CSGUnionOfTwoShells)
{
  openvdb::initialize();

  IS::SphereField sphere1;
  sphere1.center = ES::V3d(-0.3, 0.0, 0.0);
  sphere1.radius = 0.5;

  IS::SphereField sphere2;
  sphere2.center = ES::V3d(0.3, 0.0, 0.0);
  sphere2.radius = 0.5;

  IS::OpenVDBOptions opts;
  opts.voxelSize = 0.04;
  opts.halfWidth = 3.0;
  IS::validateOpenVDBOptions(opts);

  auto shell1 = IS::buildOpenVDBSphereShell(sphere1, 0.15, opts);
  auto shell2 = IS::buildOpenVDBSphereShell(sphere2, 0.15, opts);
  ASSERT_NE(shell1, nullptr);
  ASSERT_NE(shell2, nullptr);

  auto combined = IS::combineOpenVDBLevelSets(*shell1, *shell2, IS::BooleanOp::Union);
  ASSERT_NE(combined, nullptr);

  pgo::Mesh::TriMeshGeo outMesh;
  IS::extractOpenVDBLevelSet(*combined, opts, outMesh);
  EXPECT_GT(outMesh.numVertices(), 0);
  EXPECT_GT(outMesh.numTriangles(), 0);

  // Union of two separate shells should have more volume than either alone
  pgo::Mesh::TriMeshGeo mesh1, mesh2;
  IS::extractOpenVDBLevelSet(*shell1, opts, mesh1);
  IS::extractOpenVDBLevelSet(*shell2, opts, mesh2);
  const double volUnion = pgo::Mesh::computeMeshVolume(outMesh);
  const double vol1 = pgo::Mesh::computeMeshVolume(mesh1);
  const double vol2 = pgo::Mesh::computeMeshVolume(mesh2);
  // Union volume should be roughly sum of individual volumes (slightly less if they overlap)
  EXPECT_GT(volUnion, std::max(vol1, vol2));
}

TEST(OpenVDBExtractorTest, CSGIntersectionReducesVolume)
{
  openvdb::initialize();

  // Two overlapping spheres
  IS::SphereField sphere;
  sphere.center = ES::V3d(0.0, 0.0, 0.0);
  sphere.radius = 0.5;

  IS::OpenVDBOptions opts;
  opts.voxelSize = 0.04;
  opts.halfWidth = 3.0;
  IS::validateOpenVDBOptions(opts);

  auto ball = IS::buildOpenVDBBallLevelSet(sphere, opts);
  ASSERT_NE(ball, nullptr);

  pgo::Mesh::TriMeshGeo meshBall;
  IS::extractOpenVDBLevelSet(*ball, opts, meshBall);
  const double ballVol = pgo::Mesh::computeMeshVolume(meshBall);

  // Intersection of ball with itself = same ball
  auto intersection = IS::combineOpenVDBLevelSets(*ball, *ball, IS::BooleanOp::Intersection);
  ASSERT_NE(intersection, nullptr);

  pgo::Mesh::TriMeshGeo meshIntersection;
  IS::extractOpenVDBLevelSet(*intersection, opts, meshIntersection);
  const double intersectVol = pgo::Mesh::computeMeshVolume(meshIntersection);

  EXPECT_NEAR(intersectVol, ballVol, ballVol * 0.1);  // ~same volume
}

TEST(OpenVDBExtractorTest, SmoothingChangesMesh)
{
  openvdb::initialize();

  IS::SphereField sphere;
  sphere.center = ES::V3d(0.0, 0.0, 0.0);
  sphere.radius = 0.5;

  IS::OpenVDBOptions optsNoSmooth;
  optsNoSmooth.voxelSize = 0.05;
  optsNoSmooth.halfWidth = 3.0;
  optsNoSmooth.smoothSteps = 0;
  IS::validateOpenVDBOptions(optsNoSmooth);

  IS::OpenVDBOptions optsSmooth = optsNoSmooth;
  optsSmooth.smoothSteps = 3;

  auto shell1 = IS::buildOpenVDBSphereShell(sphere, 0.15, optsNoSmooth);
  auto shell2 = IS::buildOpenVDBSphereShell(sphere, 0.15, optsSmooth);

  pgo::Mesh::TriMeshGeo meshNoSmooth, meshSmooth;
  IS::extractOpenVDBLevelSet(*shell1, optsNoSmooth, meshNoSmooth);
  IS::extractOpenVDBLevelSet(*shell2, optsSmooth, meshSmooth);

  EXPECT_GT(meshNoSmooth.numTriangles(), 0);
  EXPECT_GT(meshSmooth.numTriangles(), 0);
  // Smoothing should change the mesh
  EXPECT_NE(meshNoSmooth.numVertices(), meshSmooth.numVertices());
}

#endif  // PGO_HAS_OPENVDB

// =========================================================================
// OpenVDB stub tests (when OpenVDB is not available)
// =========================================================================

#ifndef PGO_HAS_OPENVDB

TEST(OpenVDBStubTest, BuildFunctionsThrowWhenDisabled)
{
  IS::OpenVDBOptions opts;
  opts.voxelSize = 0.1;
  opts.halfWidth = 3.0;
  IS::validateOpenVDBOptions(opts);

  pgo::Mesh::TriMeshGeo dummyMesh;
  IS::SphereField dummySphere;

  EXPECT_THROW(IS::buildOpenVDBShellFromMesh(dummyMesh, 0.1, opts), std::runtime_error);
  EXPECT_THROW(IS::buildOpenVDBSphereShell(dummySphere, 0.1, opts), std::runtime_error);
  EXPECT_THROW(IS::buildOpenVDBBallLevelSet(dummySphere, opts), std::runtime_error);

  IS::OpenVDBLevelSet dummyLS;
  EXPECT_THROW(IS::combineOpenVDBLevelSets(dummyLS, dummyLS, IS::BooleanOp::Union), std::runtime_error);

  pgo::Mesh::TriMeshGeo outMesh;
  EXPECT_THROW(IS::extractOpenVDBLevelSet(dummyLS, opts, outMesh), std::runtime_error);
}

#endif  // !PGO_HAS_OPENVDB
