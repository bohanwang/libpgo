#include "geometry/sphereField.h"

#include "boundingBox.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>

namespace pgo::ImplicitSurface {

SphereField computeSphereFieldFromBBox(const Mesh::TriMeshGeo &sphereMesh)
{
  if (sphereMesh.numVertices() == 0)
    throw std::runtime_error("Cannot derive sphere parameters from an empty sphere mesh");

  const Mesh::BoundingBox bbox(sphereMesh.positions());
  SphereField field;
  field.center = bbox.center();
  field.radius = 0.0;
  for (int i = 0; i < sphereMesh.numVertices(); ++i) {
    field.radius = std::max(field.radius, (sphereMesh.pos(i) - field.center).norm());
  }

  if (field.radius <= 0.0 || !std::isfinite(field.radius))
    throw std::runtime_error("Derived sphere radius is not positive");

  return field;
}

int projectOpenBoundaryToSphere(Mesh::TriMeshGeo &mesh, const SphereField &sphere)
{
  std::map<std::pair<int, int>, int> edgeIncidentCounts;
  for (int triID = 0; triID < mesh.numTriangles(); ++triID) {
    const Vec3i &tri = mesh.tri(triID);
    for (int i = 0; i < 3; ++i) {
      int a = tri[i];
      int b = tri[(i + 1) % 3];
      if (a > b)
        std::swap(a, b);
      ++edgeIncidentCounts[std::make_pair(a, b)];
    }
  }

  std::vector<char> isBoundaryVertex(mesh.numVertices(), 0);
  for (const auto &entry : edgeIncidentCounts) {
    if (entry.second != 1)
      continue;
    isBoundaryVertex[entry.first.first] = 1;
    isBoundaryVertex[entry.first.second] = 1;
  }

  int projectedVertices = 0;
  for (int vertexID = 0; vertexID < mesh.numVertices(); ++vertexID) {
    if (!isBoundaryVertex[vertexID])
      continue;

    const EigenSupport::V3d radial = mesh.pos(vertexID) - sphere.center;
    const double radialNorm = radial.norm();
    if (radialNorm <= 0.0)
      continue;

    mesh.pos(vertexID) = sphere.center + radial * (sphere.radius / radialNorm);
    ++projectedVertices;
  }

  return projectedVertices;
}

void thickenSphereShell(const SphereField &sphere, double thickness,
  const GridSpec &grid, DenseGrid &out)
{
  if (out.gridSpec() != grid)
    throw std::runtime_error("thickenSphereShell: grid spec mismatch");

  const int resolution = grid.resolution;
  const EigenSupport::V3d delta = (grid.bmax - grid.bmin) / static_cast<double>(resolution - 1);
  const double halfThickness = 0.5 * thickness;

  for (int z = 0; z < resolution; ++z) {
    for (int y = 0; y < resolution; ++y) {
      for (int x = 0; x < resolution; ++x) {
        const EigenSupport::V3d p =
          grid.bmin + delta.cwiseProduct(EigenSupport::V3d(x, y, z).cast<double>());
        const double radial = (p - sphere.center).norm();
        out.at(x, y, z) = std::abs(radial - sphere.radius) - halfThickness;
      }
    }
  }
}

void evaluateBallSDF(const SphereField &sphere, const GridSpec &grid, DenseGrid &out)
{
  if (out.gridSpec() != grid)
    throw std::runtime_error("evaluateBallSDF: grid spec mismatch");

  const int resolution = grid.resolution;
  const EigenSupport::V3d delta = (grid.bmax - grid.bmin) / static_cast<double>(resolution - 1);

  for (int z = 0; z < resolution; ++z) {
    for (int y = 0; y < resolution; ++y) {
      for (int x = 0; x < resolution; ++x) {
        const EigenSupport::V3d p =
          grid.bmin + delta.cwiseProduct(EigenSupport::V3d(x, y, z).cast<double>());
        out.at(x, y, z) = (p - sphere.center).norm() - sphere.radius;
      }
    }
  }
}

}  // namespace pgo::ImplicitSurface
