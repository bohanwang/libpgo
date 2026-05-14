#include "extraction/openVDBExtractor.h"

#ifdef PGO_HAS_OPENVDB
#include <openvdb/openvdb.h>
#include <openvdb/tools/Composite.h>
#include <openvdb/tools/LevelSetFilter.h>
#include <openvdb/tools/LevelSetSphere.h>
#include <openvdb/tools/MeshToVolume.h>
#include <openvdb/tools/VolumeToMesh.h>
#endif

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace pgo::ImplicitSurface {

void validateOpenVDBOptions(const OpenVDBOptions &options)
{
  if (options.voxelSize <= 0.0 || !std::isfinite(options.voxelSize))
    throw std::runtime_error("OpenVDBOptions: voxelSize must be positive");
  if (options.halfWidth <= 0.0 || !std::isfinite(options.halfWidth))
    throw std::runtime_error("OpenVDBOptions: halfWidth must be positive");
  if (options.adaptivity < 0.0 || !std::isfinite(options.adaptivity))
    throw std::runtime_error("OpenVDBOptions: adaptivity must be non-negative");
  if (options.smoothSteps < 0)
    throw std::runtime_error("OpenVDBOptions: smoothSteps must be non-negative");
}

#ifdef PGO_HAS_OPENVDB
namespace {

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

}  // namespace

std::unique_ptr<OpenVDBLevelSet> buildOpenVDBShellFromMesh(
  const Mesh::TriMeshGeo &surfaceMesh, double shellThickness,
  const OpenVDBOptions &options)
{
  validateOpenVDBOptions(options);

  const double shellHalfThickness = 0.5 * shellThickness;

  const std::vector<openvdb::Vec3s> points = meshPointsForOpenVDB(surfaceMesh);
  const std::vector<openvdb::Vec3I> triangles = meshTrianglesForOpenVDB(surfaceMesh);
  const std::vector<openvdb::Vec4I> quads;
  const openvdb::math::Transform::Ptr transform =
    openvdb::math::Transform::createLinearTransform(options.voxelSize);

  const float bandWidth = static_cast<float>(options.halfWidth + shellHalfThickness / options.voxelSize);
  openvdb::FloatGrid::Ptr grid = openvdb::tools::meshToUnsignedDistanceField<openvdb::FloatGrid>(
    *transform, points, triangles, quads, bandWidth);

  offsetOpenVDBGrid(*grid, static_cast<float>(-shellHalfThickness));
  grid->setGridClass(openvdb::GRID_LEVEL_SET);
  grid->setName("shell");

  return std::make_unique<OpenVDBLevelSet>(std::move(grid));
}

std::unique_ptr<OpenVDBLevelSet> buildOpenVDBSphereShell(
  const SphereField &sphere, double sphereShellThickness,
  const OpenVDBOptions &options)
{
  validateOpenVDBOptions(options);

  const double shellHalfThickness = 0.5 * sphereShellThickness;
  const double innerRadius = sphere.radius - shellHalfThickness;
  const double outerRadius = sphere.radius + shellHalfThickness;

  if (innerRadius <= 0.0)
    throw std::runtime_error("Sphere shell inner radius is non-positive");

  auto outer = openvdb::tools::createLevelSetSphere<openvdb::FloatGrid>(
    static_cast<float>(outerRadius),
    openvdb::Vec3f(
      static_cast<float>(sphere.center[0]),
      static_cast<float>(sphere.center[1]),
      static_cast<float>(sphere.center[2])),
    static_cast<float>(options.voxelSize),
    static_cast<float>(options.halfWidth));

  auto inner = openvdb::tools::createLevelSetSphere<openvdb::FloatGrid>(
    static_cast<float>(innerRadius),
    openvdb::Vec3f(
      static_cast<float>(sphere.center[0]),
      static_cast<float>(sphere.center[1]),
      static_cast<float>(sphere.center[2])),
    static_cast<float>(options.voxelSize),
    static_cast<float>(options.halfWidth));

  openvdb::tools::csgDifference(*outer, *inner);
  outer->setName("sphere_shell");

  return std::make_unique<OpenVDBLevelSet>(std::move(outer));
}

std::unique_ptr<OpenVDBLevelSet> buildOpenVDBBallLevelSet(
  const SphereField &sphere, const OpenVDBOptions &options)
{
  validateOpenVDBOptions(options);

  auto grid = openvdb::tools::createLevelSetSphere<openvdb::FloatGrid>(
    static_cast<float>(sphere.radius),
    openvdb::Vec3f(
      static_cast<float>(sphere.center[0]),
      static_cast<float>(sphere.center[1]),
      static_cast<float>(sphere.center[2])),
    static_cast<float>(options.voxelSize),
    static_cast<float>(options.halfWidth));

  return std::make_unique<OpenVDBLevelSet>(std::move(grid));
}

std::unique_ptr<OpenVDBLevelSet> combineOpenVDBLevelSets(
  const OpenVDBLevelSet &a, const OpenVDBLevelSet &b, BooleanOp op)
{
  if (!a.grid || !b.grid)
    throw std::runtime_error("combineOpenVDBLevelSets: one or both level sets are null");

  openvdb::FloatGrid::Ptr result = a.grid->deepCopy();

  switch (op) {
    case BooleanOp::Union:
      openvdb::tools::csgUnion(*result, *b.grid);
      break;
    case BooleanOp::Intersection:
      openvdb::tools::csgIntersection(*result, *b.grid);
      break;
    case BooleanOp::Difference:
      openvdb::tools::csgDifference(*result, *b.grid);
      break;
  }

  return std::make_unique<OpenVDBLevelSet>(std::move(result));
}

void extractOpenVDBLevelSet(const OpenVDBLevelSet &levelSet,
  const OpenVDBOptions &options, Mesh::TriMeshGeo &outMesh)
{
  if (!levelSet.grid)
    throw std::runtime_error("extractOpenVDBLevelSet: level set grid is null");

  openvdb::FloatGrid &grid = *levelSet.grid;

  for (int i = 0; i < options.smoothSteps; ++i) {
    openvdb::tools::LevelSetFilter<openvdb::FloatGrid> filter(grid);
    filter.meanCurvature();
  }

  std::vector<openvdb::Vec3s> points;
  std::vector<openvdb::Vec3I> triangles;
  std::vector<openvdb::Vec4I> quads;
  openvdb::tools::volumeToMesh(grid, points, triangles, quads,
    /*isovalue=*/0.0, options.adaptivity, /*relaxDisorientedTriangles=*/true);

  outMesh = openVDBMeshToTriMeshGeo(points, triangles, quads);
}

#else  // !PGO_HAS_OPENVDB

void validateOpenVDBOptions(const OpenVDBOptions &options)
{
  // Validating options without OpenVDB is a no-op; actual usage will throw.
}

std::unique_ptr<OpenVDBLevelSet> buildOpenVDBShellFromMesh(
  const Mesh::TriMeshGeo &, double, const OpenVDBOptions &)
{
  throw std::runtime_error("OpenVDB backend is unavailable; configure with -DPGO_ENABLE_OPENVDB=ON");
}

std::unique_ptr<OpenVDBLevelSet> buildOpenVDBSphereShell(
  const SphereField &, double, const OpenVDBOptions &)
{
  throw std::runtime_error("OpenVDB backend is unavailable; configure with -DPGO_ENABLE_OPENVDB=ON");
}

std::unique_ptr<OpenVDBLevelSet> buildOpenVDBBallLevelSet(
  const SphereField &, const OpenVDBOptions &)
{
  throw std::runtime_error("OpenVDB backend is unavailable; configure with -DPGO_ENABLE_OPENVDB=ON");
}

std::unique_ptr<OpenVDBLevelSet> combineOpenVDBLevelSets(
  const OpenVDBLevelSet &, const OpenVDBLevelSet &, BooleanOp)
{
  throw std::runtime_error("OpenVDB backend is unavailable; configure with -DPGO_ENABLE_OPENVDB=ON");
}

void extractOpenVDBLevelSet(const OpenVDBLevelSet &,
  const OpenVDBOptions &, Mesh::TriMeshGeo &)
{
  throw std::runtime_error("OpenVDB backend is unavailable; configure with -DPGO_ENABLE_OPENVDB=ON");
}

#endif

}  // namespace pgo::ImplicitSurface
