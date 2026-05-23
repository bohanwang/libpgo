#include "cubicMesherIO.h"

#include "generateSurfaceMesh.h"
#include "triMeshGeo.h"

#include <stdexcept>
#include <utility>
#include <vector>

namespace cubic_mesher
{

void saveCubicMesh(const pgo::VolumetricMeshes::CubicMesh &cubicMesh, const std::string &outputMesh)
{
  if (cubicMesh.save(outputMesh.c_str()) != 0)
    throw std::runtime_error("Failed to save cubic mesh to " + outputMesh);

  pgo::VolumetricMeshes::CubicMesh reloadedMesh(outputMesh.c_str());
  if (reloadedMesh.getNumVertices() != cubicMesh.getNumVertices() ||
    reloadedMesh.getNumElements() != cubicMesh.getNumElements()) {
    throw std::runtime_error("Saved cubic mesh reload check failed for " + outputMesh);
  }
}

void writeSurfaceMesh(const pgo::VolumetricMeshes::CubicMesh &cubicMesh, const std::string &outputSurface)
{
  std::vector<pgo::EigenSupport::V3d> vertices;
  std::vector<std::vector<int>> faces;
  pgo::VolumetricMeshes::GenerateSurfaceMesh::computeMesh(&cubicMesh, vertices, faces, true, false);

  std::vector<pgo::Vec3i> triangles;
  triangles.reserve(faces.size());
  for (const auto &face : faces) {
    if (face.size() != 3)
      throw std::runtime_error("Surface extraction returned a non-triangle face");

    triangles.emplace_back(face[0], face[1], face[2]);
  }

  pgo::Mesh::TriMeshGeo rawSurface(std::move(vertices), std::move(triangles));
  pgo::Mesh::TriMeshGeo surface = pgo::Mesh::removeIsolatedVertices(rawSurface.ref());
  if (surface.save(outputSurface) != true)
    throw std::runtime_error("Failed to save cubic mesh surface to " + outputSurface);
}

}  // namespace cubic_mesher
