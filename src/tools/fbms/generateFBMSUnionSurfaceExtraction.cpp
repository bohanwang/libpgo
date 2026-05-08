#include "generateFBMSUnionSurfaceExtraction.h"

#include "generateFBMSUnionSurfaceVolume.h"
#include "libiglInterface.h"

#ifdef PGO_HAS_OPENVDB
#include <openvdb/openvdb.h>
#include <openvdb/tools/Composite.h>
#include <openvdb/tools/LevelSetFilter.h>
#include <openvdb/tools/LevelSetSphere.h>
#include <openvdb/tools/MeshToVolume.h>
#include <openvdb/tools/VolumeToMesh.h>
#endif

#include <iostream>
#include <stdexcept>
#include <utility>
#include <vector>

namespace pgo::Tools::FBMSUnionSurface
{

void extractSurfaceMarchingCubes(const EigenSupport::V3d &bmin, const EigenSupport::V3d &bmax,
  int resolution, const EigenSupport::VXd &field, double isoOffset, Mesh::TriMeshGeo &outMesh)
{
  libiglInterface::computeMarchingCubes(bmin, bmax, resolution, field, outMesh, isoOffset);
}

#ifdef PGO_HAS_OPENVDB
namespace
{

std::vector<openvdb::Vec3s> meshPointsForOpenVDB(const Mesh::TriMeshGeo &mesh)
{
  std::vector<openvdb::Vec3s> points;
  points.reserve(mesh.numVertices());
  for (int i = 0; i < mesh.numVertices(); ++i) {
    const Vec3d &p = mesh.pos(i);
    points.emplace_back(static_cast<float>(p[0]), static_cast<float>(p[1]), static_cast<float>(p[2]));
  }
  return points;
}

std::vector<openvdb::Vec3I> meshTrianglesForOpenVDB(const Mesh::TriMeshGeo &mesh)
{
  std::vector<openvdb::Vec3I> triangles;
  triangles.reserve(mesh.numTriangles());
  for (int i = 0; i < mesh.numTriangles(); ++i) {
    const Vec3i &tri = mesh.tri(i);
    triangles.emplace_back(tri[0], tri[1], tri[2]);
  }
  return triangles;
}

void offsetOpenVDBGrid(openvdb::FloatGrid &grid, float offset)
{
  for (auto iter = grid.beginValueOn(); iter; ++iter)
    iter.setValue(iter.getValue() + offset);
  grid.tree().root().setBackground(grid.background() + offset, /*updateChildNodes=*/false);
}

openvdb::FloatGrid::Ptr createFBMSShellGrid(const Mesh::TriMeshGeo &fbmsMesh,
  double fbmsThickness, double voxelSize, double halfWidth)
{
  const std::vector<openvdb::Vec3s> points = meshPointsForOpenVDB(fbmsMesh);
  const std::vector<openvdb::Vec3I> triangles = meshTrianglesForOpenVDB(fbmsMesh);
  const std::vector<openvdb::Vec4I> quads;
  const openvdb::math::Transform::Ptr transform =
    openvdb::math::Transform::createLinearTransform(voxelSize);

  const double shellHalfThickness = 0.5 * fbmsThickness;
  const float bandWidth = static_cast<float>(halfWidth + shellHalfThickness / voxelSize);
  openvdb::FloatGrid::Ptr grid = openvdb::tools::meshToUnsignedDistanceField<openvdb::FloatGrid>(
    *transform, points, triangles, quads, bandWidth);

  offsetOpenVDBGrid(*grid, static_cast<float>(-shellHalfThickness));
  grid->setGridClass(openvdb::GRID_LEVEL_SET);
  grid->setName("fbms_shell");
  return grid;
}

openvdb::FloatGrid::Ptr createLevelSetSphereGrid(const SphereParameters &sphere,
  double radius, double voxelSize, double halfWidth)
{
  return openvdb::tools::createLevelSetSphere<openvdb::FloatGrid>(
    static_cast<float>(radius),
    openvdb::Vec3f(
      static_cast<float>(sphere.center[0]),
      static_cast<float>(sphere.center[1]),
      static_cast<float>(sphere.center[2])),
    static_cast<float>(voxelSize),
    static_cast<float>(halfWidth));
}

openvdb::FloatGrid::Ptr createSphereShellGrid(const SphereParameters &sphere,
  double sphereThickness, double voxelSize, double halfWidth)
{
  const double shellHalfThickness = 0.5 * sphereThickness;
  const double innerRadius = sphere.radius - shellHalfThickness;
  const double outerRadius = sphere.radius + shellHalfThickness;
  if (innerRadius <= 0.0)
    throw std::runtime_error("Sphere shell inner radius is non-positive");

  openvdb::FloatGrid::Ptr outer = createLevelSetSphereGrid(sphere, outerRadius, voxelSize, halfWidth);
  openvdb::FloatGrid::Ptr inner = createLevelSetSphereGrid(sphere, innerRadius, voxelSize, halfWidth);
  openvdb::tools::csgDifference(*outer, *inner);
  outer->setName("sphere_shell");
  return outer;
}

Mesh::TriMeshGeo openVDBMeshToTriMeshGeo(const std::vector<openvdb::Vec3s> &points,
  const std::vector<openvdb::Vec3I> &triangles, const std::vector<openvdb::Vec4I> &quads)
{
  std::vector<Vec3d> positions;
  positions.reserve(points.size());
  for (const openvdb::Vec3s &p : points)
    positions.emplace_back(p[0], p[1], p[2]);

  std::vector<Vec3i> outTriangles;
  outTriangles.reserve(triangles.size() + quads.size() * 2);
  for (const openvdb::Vec3I &tri : triangles)
    outTriangles.emplace_back(tri[0], tri[1], tri[2]);
  for (const openvdb::Vec4I &quad : quads) {
    outTriangles.emplace_back(quad[0], quad[1], quad[2]);
    outTriangles.emplace_back(quad[0], quad[2], quad[3]);
  }

  return Mesh::TriMeshGeo(std::move(positions), std::move(outTriangles));
}

void meshOpenVDBGrid(openvdb::FloatGrid &grid, double adaptivity, int smoothSteps,
  Mesh::TriMeshGeo &outMesh)
{
  for (int i = 0; i < smoothSteps; ++i) {
    openvdb::tools::LevelSetFilter<openvdb::FloatGrid> filter(grid);
    filter.meanCurvature();
  }

  std::vector<openvdb::Vec3s> points;
  std::vector<openvdb::Vec3I> triangles;
  std::vector<openvdb::Vec4I> quads;
  openvdb::tools::volumeToMesh(grid, points, triangles, quads,
    /*isovalue=*/0.0, adaptivity, /*relaxDisorientedTriangles=*/true);

  outMesh = openVDBMeshToTriMeshGeo(points, triangles, quads);
}

}  // namespace

void extractSurfaceOpenVDB(const Mesh::TriMeshGeo &fbmsMesh, const SphereParameters &sphere,
  double fbmsThickness, double sphereThickness, bool enableTruncating, const std::string &debugFieldMode,
  double voxelSize, double halfWidth, double adaptivity, int smoothSteps,
  Mesh::TriMeshGeo &outMesh)
{
  openvdb::initialize();

  openvdb::FloatGrid::Ptr grid;
  if (debugFieldMode != "sphere") {
    grid = createFBMSShellGrid(fbmsMesh, fbmsThickness, voxelSize, halfWidth);
    if (enableTruncating) {
      openvdb::FloatGrid::Ptr truncationBall = createLevelSetSphereGrid(sphere, sphere.radius, voxelSize, halfWidth);
      openvdb::tools::csgIntersection(*grid, *truncationBall);
    }
  }

  if (debugFieldMode == "sphere") {
    grid = createSphereShellGrid(sphere, sphereThickness, voxelSize, halfWidth);
  }
  else if (debugFieldMode == "union") {
    openvdb::FloatGrid::Ptr sphereShell = createSphereShellGrid(sphere, sphereThickness, voxelSize, halfWidth);
    openvdb::tools::csgUnion(*grid, *sphereShell);
  }
  else if (debugFieldMode == "union-minus-sphere") {
    openvdb::FloatGrid::Ptr sphereShellForUnion = createSphereShellGrid(sphere, sphereThickness, voxelSize, halfWidth);
    openvdb::tools::csgUnion(*grid, *sphereShellForUnion);

    openvdb::FloatGrid::Ptr sphereShellForDifference = createSphereShellGrid(sphere, sphereThickness, voxelSize, halfWidth);
    openvdb::tools::csgDifference(*grid, *sphereShellForDifference);
  }

  if (!grid)
    throw std::runtime_error("OpenVDB backend failed to create a level-set grid");

  meshOpenVDBGrid(*grid, adaptivity, smoothSteps, outMesh);
}

ThicknessSearchResult findThicknessForVolumeBudgetOpenVDB(
  const Mesh::TriMeshGeo &fbmsMesh, const SphereParameters &sphere,
  double sphereThickness, bool enableTruncating,
  double voxelSize, double halfWidth, double adaptivity, int smoothSteps,
  double targetVolume, double tLo, double tHi, int maxIterations, double relTol)
{
  openvdb::initialize();

  auto evalAt = [&](double t, Mesh::TriMeshGeo &outMesh) -> double {
    extractSurfaceOpenVDB(fbmsMesh, sphere,
      t, sphereThickness, enableTruncating, "union",
      voxelSize, halfWidth, adaptivity, smoothSteps,
      outMesh);
    if (outMesh.numTriangles() == 0)
      return 0.0;
    return computeMeshVolume(outMesh);
  };

  ThicknessSearchResult result;
  Mesh::TriMeshGeo meshLo, meshHi, meshMid;

  const double vLo = evalAt(tLo, meshLo);
  std::cout << "[volume-search] tLo=" << tLo << " vol=" << vLo << " (target=" << targetVolume << ")" << std::endl;

  if (vLo > targetVolume) {
    result.budgetExceededAtMin = true;
    result.thickness = tLo;
    result.volume = vLo;
    result.mesh = std::move(meshLo);
    return result;
  }

  const double vHi = evalAt(tHi, meshHi);
  std::cout << "[volume-search] tHi=" << tHi << " vol=" << vHi << std::endl;

  if (vHi <= targetVolume) {
    result.budgetUnreachedAtMax = true;
    result.thickness = tHi;
    result.volume = vHi;
    result.mesh = std::move(meshHi);
    return result;
  }

  result.thickness = tLo;
  result.volume = vLo;
  result.mesh = std::move(meshLo);

  for (int it = 0; it < maxIterations; ++it) {
    const double mid = 0.5 * (tLo + tHi);
    const double vMid = evalAt(mid, meshMid);
    std::cout << "[volume-search] iter=" << it << " t=" << mid << " vol=" << vMid << std::endl;

    if (vMid <= targetVolume) {
      tLo = mid;
      result.thickness = mid;
      result.volume = vMid;
      result.mesh = std::move(meshMid);
    }
    else {
      tHi = mid;
    }
    result.iterations = it + 1;
    if (tHi <= 0.0 || (tHi - tLo) / tHi < relTol)
      break;
  }
  return result;
}
#else
void extractSurfaceOpenVDB(const Mesh::TriMeshGeo &, const SphereParameters &,
  double, double, bool, const std::string &, double, double, double, int,
  Mesh::TriMeshGeo &)
{
  throw std::runtime_error("OpenVDB backend is unavailable; configure with -DPGO_ENABLE_OPENVDB=ON");
}

ThicknessSearchResult findThicknessForVolumeBudgetOpenVDB(
  const Mesh::TriMeshGeo &, const SphereParameters &,
  double, bool, double, double, double, int,
  double, double, double, int, double)
{
  throw std::runtime_error("OpenVDB backend is unavailable; configure with -DPGO_ENABLE_OPENVDB=ON");
}
#endif

}  // namespace pgo::Tools::FBMSUnionSurface
