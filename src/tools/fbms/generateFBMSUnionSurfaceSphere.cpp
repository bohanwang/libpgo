#include "generateFBMSUnionSurfaceSphere.h"

#include "boundingBox.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <optional>
#include <regex>
#include <stdexcept>

namespace pgo::Tools::FBMSUnionSurface
{
namespace
{

std::optional<EigenSupport::V3d> parseCenterLine(const std::string &line)
{
  static const std::regex centerRegex(
    R"(^\s*#\s*center\s*=\s*\(\s*([^,\s]+)\s*,\s*([^,\s]+)\s*,\s*([^,\s\)]+)\s*\)\s*$)",
    std::regex::icase);

  std::smatch match;
  if (!std::regex_match(line, match, centerRegex))
    return std::nullopt;

  return EigenSupport::V3d(std::stod(match[1].str()), std::stod(match[2].str()), std::stod(match[3].str()));
}

std::optional<double> parseRadiusLine(const std::string &line)
{
  static const std::regex radiusRegex(
    R"(^\s*#\s*radius\s*=\s*([^\s#]+).*$)",
    std::regex::icase);

  std::smatch match;
  if (!std::regex_match(line, match, radiusRegex))
    return std::nullopt;

  return std::stod(match[1].str());
}

std::optional<SphereParameters> parseSphereHeader(const std::string &spherePath)
{
  std::ifstream in(spherePath);
  if (!in.is_open())
    throw std::runtime_error("Failed to open sphere OBJ: " + spherePath);

  std::optional<EigenSupport::V3d> center;
  std::optional<double> radius;
  std::string line;
  while (std::getline(in, line)) {
    if (!line.empty() && line[0] != '#')
      break;

    if (!center)
      center = parseCenterLine(line);
    if (!radius)
      radius = parseRadiusLine(line);
  }

  if (!center || !radius || *radius <= 0.0 || !std::isfinite(*radius))
    return std::nullopt;

  return SphereParameters{ *center, *radius, true };
}

SphereParameters computeSphereParametersFromBBox(const Mesh::TriMeshGeo &sphereMesh)
{
  if (sphereMesh.numVertices() == 0)
    throw std::runtime_error("Cannot derive sphere parameters from an empty sphere mesh");

  const Mesh::BoundingBox bbox(sphereMesh.positions());
  SphereParameters params;
  params.center = bbox.center();
  params.radius = 0.0;
  for (const Vec3d &pos : sphereMesh.positions()) {
    params.radius = std::max(params.radius, (pos - params.center).norm());
  }

  if (params.radius <= 0.0 || !std::isfinite(params.radius))
    throw std::runtime_error("Derived sphere radius is not positive");

  return params;
}

}  // namespace

SphereParameters loadSphereParameters(const std::string &spherePath, const Mesh::TriMeshGeo &sphereMesh)
{
  std::optional<SphereParameters> parsed = parseSphereHeader(spherePath);
  if (parsed)
    return *parsed;

  return computeSphereParametersFromBBox(sphereMesh);
}

int projectBoundaryVerticesToSphere(Mesh::TriMeshGeo &mesh, const SphereParameters &sphere)
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

}  // namespace pgo::Tools::FBMSUnionSurface
