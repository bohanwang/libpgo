#include "tetMesherBackend.h"

#include "generateSurfaceMesh.h"
#include "triMeshGeo.h"

#include <iostream>
#include <stdexcept>
#include <utility>
#include <vector>

namespace tet_mesher
{
namespace
{

pgo::Mesh::BoundingBox computeBoundingBox(const pgo::VolumetricMeshes::TetMesh &tetMesh)
{
  return pgo::Mesh::BoundingBox(tetMesh.getNumVertices(), tetMesh.getVertices());
}

void writeTetSurfaceMesh(const pgo::VolumetricMeshes::TetMesh &tetMesh, const std::string &outputSurface)
{
  std::vector<pgo::EigenSupport::V3d> vertices;
  std::vector<std::vector<int>> faces;
  pgo::VolumetricMeshes::GenerateSurfaceMesh::computeMesh(&tetMesh, vertices, faces, true, false);

  std::vector<pgo::Vec3i> triangles;
  triangles.reserve(faces.size());
  for (const auto &face : faces) {
    if (face.size() != 3)
      throw std::runtime_error("Tet surface extraction returned a non-triangle face");

    triangles.emplace_back(face[0], face[1], face[2]);
  }

  pgo::Mesh::TriMeshGeo rawSurface(std::move(vertices), std::move(triangles));
  pgo::Mesh::TriMeshGeo surface = pgo::Mesh::removeIsolatedVertices(rawSurface.ref());
  if (surface.save(outputSurface) != true)
    throw std::runtime_error("Failed to save tet mesh surface to " + outputSurface);
}

void printStats(const pgo::VolumetricMeshes::TetMesh &tetMesh)
{
  const pgo::Mesh::BoundingBox bbox = computeBoundingBox(tetMesh);
  std::cout << "Generated tet mesh with " << tetMesh.getNumVertices() << " vertices and "
            << tetMesh.getNumElements() << " elements." << std::endl;
  std::cout << "BBox min: " << bbox.bmin()[0] << " " << bbox.bmin()[1] << " " << bbox.bmin()[2]
            << std::endl;
  std::cout << "BBox max: " << bbox.bmax()[0] << " " << bbox.bmax()[1] << " " << bbox.bmax()[2]
            << std::endl;
}

}  // namespace

void saveTetMeshOutputs(const pgo::VolumetricMeshes::TetMesh &tetMesh, const CommonOptions &options)
{
  if (tetMesh.getNumVertices() <= 0)
    throw std::runtime_error("Tet mesher generated an empty vertex set");

  if (tetMesh.getNumElements() <= 0)
    throw std::runtime_error("Tet mesher generated zero tetrahedra");

  if (tetMesh.save(options.outputMesh.c_str()) != 0)
    throw std::runtime_error("Failed to save tet mesh to " + options.outputMesh);

  pgo::VolumetricMeshes::TetMesh reloadedMesh(options.outputMesh.c_str());
  if (reloadedMesh.getNumVertices() != tetMesh.getNumVertices() ||
    reloadedMesh.getNumElements() != tetMesh.getNumElements()) {
    throw std::runtime_error("Saved tet mesh reload check failed for " + options.outputMesh);
  }

  if (options.outputSurface.empty() == false)
    writeTetSurfaceMesh(tetMesh, options.outputSurface);

  if (options.printStats)
    printStats(tetMesh);
}

}  // namespace tet_mesher
