#include "triangleMeshVoxelizer.h"

#include "arrayRef.h"
#include "boundingVolumeTree.h"
#include "geometryQuery.h"
#include "pointInsideOutsideQuery.h"
#include "triMeshGeo.h"
#include "triMeshNeighbor.h"
#include "triple.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cctype>
#include <filesystem>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace cubic_mesher
{
namespace
{

using GridIndex = pgo::triple<int, int, int>;

struct VoxelGrid
{
  pgo::Mesh::BoundingBox inputBoundingBox;
  pgo::Vec3d gridMin = pgo::Vec3d::Zero();
  double cellSize = 0.0;
  int nx = 0;
  int ny = 0;
  int nz = 0;

  pgo::Vec3d voxelMin(int i, int j, int k) const
  {
    return gridMin + cellSize * pgo::Vec3d((double)i, (double)j, (double)k);
  }

  pgo::Vec3d voxelCenter(int i, int j, int k) const
  {
    return voxelMin(i, j, k) + pgo::Vec3d::Constant(0.5 * cellSize);
  }

  pgo::Mesh::BoundingBox voxelBoundingBox(int i, int j, int k) const
  {
    const pgo::Vec3d bmin = voxelMin(i, j, k);
    return pgo::Mesh::BoundingBox(bmin, bmin + pgo::Vec3d::Constant(cellSize));
  }
};

std::string toLower(std::string value)
{
  std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
    return (char)std::tolower(c);
  });
  return value;
}

void validateOptions(const TriangleMeshVoxelizerOptions &options)
{
  if (options.inputMesh.empty())
    throw std::runtime_error("Triangle mesh voxelizer requires a non-empty --input-mesh");

  if (options.resolution <= 0)
    throw std::runtime_error("Resolution must be positive");

  if (options.E <= 0.0)
    throw std::runtime_error("Young's modulus must be positive");

  if (options.nu <= -1.0 || options.nu >= 0.5)
    throw std::runtime_error("Poisson ratio must be in (-1, 0.5)");

  if (options.density <= 0.0)
    throw std::runtime_error("Density must be positive");
}

pgo::Mesh::TriMeshGeo loadTriangleMesh(const std::string &inputMesh)
{
  const std::filesystem::path inputPath(inputMesh);
  if (toLower(inputPath.extension().string()) != ".obj")
    throw std::runtime_error("Phase 1 cubicMesher only supports .obj input");

  pgo::Mesh::TriMeshGeo mesh;
  if (mesh.load(inputMesh) != true)
    throw std::runtime_error("Failed to load triangle mesh: " + inputMesh);

  if (mesh.numTriangles() == 0)
    throw std::runtime_error("Input triangle mesh contains no triangles");

  return mesh;
}

void validateTriangleMesh(const pgo::Mesh::TriMeshGeo &mesh)
{
  pgo::Mesh::TriMeshBVTree bvTree;
  bvTree.buildByInertiaPartition(mesh.ref());

  std::vector<std::pair<int, int>> selfIntersections;
  bvTree.selfIntersectionExact(mesh.ref(), selfIntersections);
  if (selfIntersections.empty() == false)
    throw std::runtime_error("Input triangle mesh has self intersections");

  if (pgo::Mesh::areTrianglesManifold(pgo::BasicAlgorithms::makeArrayRef(mesh.triangles())) == false)
    throw std::runtime_error("Input triangle mesh is not manifold");

  if (pgo::Mesh::getExteriorEdges(pgo::BasicAlgorithms::makeArrayRef(mesh.triangles())).empty() == false)
    throw std::runtime_error("Input triangle mesh is not closed");
}

VoxelGrid buildVoxelGrid(const pgo::Mesh::TriMeshGeo &mesh, int resolution)
{
  VoxelGrid grid;
  grid.inputBoundingBox = mesh.ref().computeTriangleBoundingBox();

  const pgo::Vec3d sides = grid.inputBoundingBox.sides();
  const double smin = sides.minCoeff();
  if (smin <= 0.0)
    throw std::runtime_error("Input triangle mesh has a degenerate bounding box");

  grid.cellSize = smin / (double)resolution;
  grid.nx = (int)std::ceil(sides[0] / grid.cellSize);
  grid.ny = (int)std::ceil(sides[1] / grid.cellSize);
  grid.nz = (int)std::ceil(sides[2] / grid.cellSize);

  const pgo::Vec3d gridSideLengths(grid.nx * grid.cellSize, grid.ny * grid.cellSize, grid.nz * grid.cellSize);
  grid.gridMin = grid.inputBoundingBox.center() - 0.5 * gridSideLengths;
  return grid;
}

bool triangleOverlapsVoxel(const pgo::Mesh::TriMeshGeo &mesh, const pgo::Mesh::BoundingBox &voxelBox,
  const std::vector<pgo::Mesh::BoundingBox> &triangleBoundingBoxes, const pgo::Mesh::BoundingBoxBVTree &triangleBoxTree,
  std::vector<int> &candidateTriangleIDs)
{
  candidateTriangleIDs.clear();
  triangleBoxTree.intersectAABB(pgo::BasicAlgorithms::makeArrayRef(triangleBoundingBoxes), voxelBox, candidateTriangleIDs);

  for (int triID : candidateTriangleIDs) {
    if (pgo::Mesh::whetherTriangleIntersectBoundingBox(mesh.pos(triID, 0), mesh.pos(triID, 1), mesh.pos(triID, 2),
        voxelBox.bmin(), voxelBox.bmax())) {
      return true;
    }
  }

  return false;
}

std::vector<GridIndex> collectOccupiedVoxels(const pgo::Mesh::TriMeshGeo &mesh, const VoxelGrid &grid)
{
  pgo::Mesh::PointInsideOutsideQuery insideOutside(mesh);
  const std::vector<pgo::Mesh::BoundingBox> triangleBoundingBoxes = mesh.ref().getTriangleBoundingBoxes();

  pgo::Mesh::BoundingBoxBVTree triangleBoxTree;
  triangleBoxTree.buildByInertiaPartition(pgo::BasicAlgorithms::makeArrayRef(triangleBoundingBoxes));

  std::vector<GridIndex> occupied;
  std::vector<int> candidateTriangleIDs;

  for (int i = 0; i < grid.nx; ++i) {
    for (int j = 0; j < grid.ny; ++j) {
      for (int k = 0; k < grid.nz; ++k) {
        const pgo::Vec3d center = grid.voxelCenter(i, j, k);
        bool occupiedVoxel = insideOutside.isPointInside(center);

        if (!occupiedVoxel) {
          const pgo::Mesh::BoundingBox voxelBox = grid.voxelBoundingBox(i, j, k);
          occupiedVoxel = triangleOverlapsVoxel(mesh, voxelBox, triangleBoundingBoxes, triangleBoxTree, candidateTriangleIDs);
        }

        if (occupiedVoxel)
          occupied.emplace_back(i, j, k);
      }
    }
  }

  if (occupied.empty())
    throw std::runtime_error("Triangle mesh voxelization produced no occupied voxels");

  return occupied;
}

std::unique_ptr<pgo::VolumetricMeshes::CubicMesh> buildCubicMesh(
  const std::vector<GridIndex> &occupiedVoxels, const VoxelGrid &grid, const TriangleMeshVoxelizerOptions &options)
{
  constexpr std::array<int, 8> vtxI{ 0, 1, 1, 0, 0, 1, 1, 0 };
  constexpr std::array<int, 8> vtxJ{ 0, 0, 1, 1, 0, 0, 1, 1 };
  constexpr std::array<int, 8> vtxK{ 0, 0, 0, 0, 1, 1, 1, 1 };

  std::set<GridIndex> vertexSet;
  for (const GridIndex &voxel : occupiedVoxels) {
    for (int corner = 0; corner < 8; ++corner) {
      vertexSet.emplace(voxel.first + vtxI[corner], voxel.second + vtxJ[corner], voxel.third + vtxK[corner]);
    }
  }

  std::vector<double> vertices;
  vertices.resize(vertexSet.size() * 3);
  std::map<GridIndex, int> vertexMap;

  int vertexID = 0;
  for (const GridIndex &gridVertex : vertexSet) {
    const pgo::Vec3d worldPos = grid.gridMin +
      grid.cellSize * pgo::Vec3d((double)gridVertex.first, (double)gridVertex.second, (double)gridVertex.third);
    vertices[3 * vertexID + 0] = worldPos[0];
    vertices[3 * vertexID + 1] = worldPos[1];
    vertices[3 * vertexID + 2] = worldPos[2];
    vertexMap.emplace(gridVertex, vertexID);
    ++vertexID;
  }

  std::vector<int> elements;
  elements.resize(occupiedVoxels.size() * 8);

  int elementID = 0;
  for (const GridIndex &voxel : occupiedVoxels) {
    for (int corner = 0; corner < 8; ++corner) {
      const GridIndex cornerIndex(voxel.first + vtxI[corner], voxel.second + vtxJ[corner], voxel.third + vtxK[corner]);
      elements[8 * elementID + corner] = vertexMap.at(cornerIndex);
    }
    ++elementID;
  }

  return std::make_unique<pgo::VolumetricMeshes::CubicMesh>((int)vertexSet.size(), vertices.data(),
    (int)occupiedVoxels.size(), elements.data(), options.E, options.nu, options.density);
}

}  // namespace

std::unique_ptr<pgo::VolumetricMeshes::CubicMesh> createTriangleMeshCubicMesh(
  const TriangleMeshVoxelizerOptions &options)
{
  validateOptions(options);

  const pgo::Mesh::TriMeshGeo mesh = loadTriangleMesh(options.inputMesh);
  validateTriangleMesh(mesh);

  const VoxelGrid grid = buildVoxelGrid(mesh, options.resolution);
  const std::vector<GridIndex> occupiedVoxels = collectOccupiedVoxels(mesh, grid);

  return buildCubicMesh(occupiedVoxels, grid, options);
}

}  // namespace cubic_mesher
