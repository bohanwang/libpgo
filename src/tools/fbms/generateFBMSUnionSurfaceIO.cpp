#include "generateFBMSUnionSurfaceIO.h"

#include "triMeshNeighbor.h"

#include <filesystem>
#include <iostream>
#include <stdexcept>

namespace pgo::Tools::FBMSUnionSurface
{

void loadMeshOrThrow(const std::string &path, Mesh::TriMeshGeo &mesh, const std::string &label)
{
  if (!mesh.load(path))
    throw std::runtime_error("Failed to load " + label + " mesh: " + path);
  if (mesh.numVertices() == 0)
    throw std::runtime_error(label + " mesh has no vertices: " + path);
}

void printVector(const char *label, const EigenSupport::V3d &v)
{
  std::cout << label << " = (" << v[0] << ", " << v[1] << ", " << v[2] << ")" << std::endl;
}

void saveSurfaceOrThrow(const Mesh::TriMeshGeo &mesh, const std::string &path)
{
  const std::filesystem::path outputPath(path);
  if (outputPath.has_parent_path()) {
    std::error_code ec;
    std::filesystem::create_directories(outputPath.parent_path(), ec);
    if (ec)
      throw std::runtime_error("Failed to create output directory: " + outputPath.parent_path().string());
  }

  if (!mesh.save(path))
    throw std::runtime_error("Failed to save output surface: " + path);
}

void applySmallComponentFilterIfRequested(Mesh::TriMeshGeo &rawSurface, const Options &options)
{
  if (!options.filterSmallComponents)
    return;

  const Mesh::TriangleEdgeConnectivityStats beforeStats =
    Mesh::computeTriangleEdgeConnectivityStats(rawSurface.triangles());
  Mesh::TriMeshGeo filteredSurface =
    Mesh::filterSmallTriangleComponentsByEdge(
      rawSurface, options.minComponentTriangles, options.keepLargestComponents);
  const Mesh::TriangleEdgeConnectivityStats afterStats =
    Mesh::computeTriangleEdgeConnectivityStats(filteredSurface.triangles());

  std::cout << "Small component filter = enabled" << std::endl;
  std::cout << "Min component triangles = " << options.minComponentTriangles << std::endl;
  if (options.keepLargestComponents > 0)
    std::cout << "Keep largest components = " << options.keepLargestComponents << std::endl;
  else
    std::cout << "Keep largest components = all above threshold" << std::endl;
  std::cout << "Components before filter = " << beforeStats.componentsByEdge << std::endl;
  std::cout << "Components after filter = " << afterStats.componentsByEdge << std::endl;
  std::cout << "Vertices removed by component filter = "
            << rawSurface.numVertices() - filteredSurface.numVertices() << std::endl;
  std::cout << "Triangles removed by component filter = "
            << rawSurface.numTriangles() - filteredSurface.numTriangles() << std::endl;

  if (filteredSurface.numVertices() == 0 || filteredSurface.numTriangles() == 0)
    throw std::runtime_error("Small component filter removed the entire surface");

  rawSurface = std::move(filteredSurface);
  std::cout << "Filtered vertex count = " << rawSurface.numVertices() << std::endl;
  std::cout << "Filtered face count = " << rawSurface.numTriangles() << std::endl;
}

}  // namespace pgo::Tools::FBMSUnionSurface
