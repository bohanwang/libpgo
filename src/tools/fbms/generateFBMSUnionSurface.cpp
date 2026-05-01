#include "boundingBox.h"
#include "libiglInterface.h"
#include "pgoLogging.h"
#include "triMeshGeo.h"

#include <argparse/argparse.hpp>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <optional>
#include <regex>
#include <stdexcept>
#include <string>

namespace
{

using pgo::EigenSupport::V3d;

struct Options
{
  std::string fbmsPath;
  std::string spherePath;
  std::string outputSurfacePath;
  double fbmsThickness = 0.0;
  double sphereThickness = 0.0;
  double paddingRatio = 0.0;
  int resolution = 0;
  bool enableTruncating = false;
};

struct SphereParameters
{
  V3d center = V3d::Zero();
  double radius = 0.0;
  bool fromHeader = false;
};

Options parseOptions(int argc, char *argv[])
{
  argparse::ArgumentParser program("generateFBMSUnionSurface");

  program.add_argument("--fbms")
    .help("Input FBMS triangle surface OBJ")
    .required()
    .metavar("PATH");
  program.add_argument("--sphere")
    .help("Input bounding sphere OBJ")
    .required()
    .metavar("PATH");
  program.add_argument("--fbms-thickness")
    .help("FBMS unsigned-distance thickening width. Positive: use as-is. "
          "Negative V: treat |V| as a target volume budget and binary-search "
          "the largest thickness whose union mesh has volume <= |V|.")
    .required()
    .metavar("FLOAT")
    .scan<'g', double>();
  program.add_argument("--sphere-thickness")
    .help("Bounding sphere shell thickening width")
    .required()
    .metavar("FLOAT")
    .scan<'g', double>();
  program.add_argument("--resolution")
    .help("Uniform SDF grid resolution per axis")
    .required()
    .metavar("INT")
    .scan<'i', int>();
  program.add_argument("--padding-ratio")
    .help("BBox diagonal ratio added to each side after thickness expansion")
    .required()
    .metavar("FLOAT")
    .scan<'g', double>();
  program.add_argument("--output-surface")
    .help("Output raw union surface OBJ")
    .required()
    .metavar("PATH");
  program.add_argument("--enable-truncating")
    .help("Truncate the FBMS shell by the outer surface of the thickened sphere; "
          "only FBMS structure inside that surface is kept.")
    .default_value(false)
    .implicit_value(true);

  program.parse_args(argc, argv);

  Options options;
  options.fbmsPath = program.get<std::string>("--fbms");
  options.spherePath = program.get<std::string>("--sphere");
  options.fbmsThickness = program.get<double>("--fbms-thickness");
  options.sphereThickness = program.get<double>("--sphere-thickness");
  options.resolution = program.get<int>("--resolution");
  options.paddingRatio = program.get<double>("--padding-ratio");
  options.outputSurfacePath = program.get<std::string>("--output-surface");
  options.enableTruncating = program.get<bool>("--enable-truncating");
  return options;
}

void validateOptions(const Options &options)
{
  if (options.fbmsThickness == 0.0 || !std::isfinite(options.fbmsThickness))
    throw std::runtime_error("--fbms-thickness must be a finite non-zero value "
                             "(positive = explicit thickness; negative = volume budget)");
  if (options.sphereThickness <= 0.0 || !std::isfinite(options.sphereThickness))
    throw std::runtime_error("--sphere-thickness must be positive");
  if (options.resolution < 2)
    throw std::runtime_error("--resolution must be at least 2");
  if (options.paddingRatio < 0.0 || !std::isfinite(options.paddingRatio))
    throw std::runtime_error("--padding-ratio must be finite and non-negative");

  const long long gridCount = 1LL * options.resolution * options.resolution * options.resolution;
  if (gridCount > std::numeric_limits<int>::max())
    throw std::runtime_error("--resolution is too large for the current marching cubes interface");
}

std::optional<V3d> parseCenterLine(const std::string &line)
{
  static const std::regex centerRegex(
    R"(^\s*#\s*center\s*=\s*\(\s*([^,\s]+)\s*,\s*([^,\s]+)\s*,\s*([^,\s\)]+)\s*\)\s*$)",
    std::regex::icase);

  std::smatch match;
  if (!std::regex_match(line, match, centerRegex))
    return std::nullopt;

  return V3d(std::stod(match[1].str()), std::stod(match[2].str()), std::stod(match[3].str()));
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

  std::optional<V3d> center;
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

SphereParameters computeSphereParametersFromBBox(const pgo::Mesh::TriMeshGeo &sphereMesh)
{
  if (sphereMesh.numVertices() == 0)
    throw std::runtime_error("Cannot derive sphere parameters from an empty sphere mesh");

  const pgo::Mesh::BoundingBox bbox(sphereMesh.positions());
  SphereParameters params;
  params.center = bbox.center();
  params.radius = 0.0;
  for (const pgo::Vec3d &pos : sphereMesh.positions()) {
    params.radius = std::max(params.radius, (pos - params.center).norm());
  }

  if (params.radius <= 0.0 || !std::isfinite(params.radius))
    throw std::runtime_error("Derived sphere radius is not positive");

  return params;
}

SphereParameters loadSphereParameters(const std::string &spherePath, const pgo::Mesh::TriMeshGeo &sphereMesh)
{
  std::optional<SphereParameters> parsed = parseSphereHeader(spherePath);
  if (parsed)
    return *parsed;

  return computeSphereParametersFromBBox(sphereMesh);
}

void loadMeshOrThrow(const std::string &path, pgo::Mesh::TriMeshGeo &mesh, const std::string &label)
{
  if (!mesh.load(path))
    throw std::runtime_error("Failed to load " + label + " mesh: " + path);
  if (mesh.numVertices() == 0)
    throw std::runtime_error(label + " mesh has no vertices: " + path);
}

void printVector(const char *label, const V3d &v)
{
  std::cout << label << " = (" << v[0] << ", " << v[1] << ", " << v[2] << ")" << std::endl;
}

void computeUnionBBox(const pgo::Mesh::TriMeshGeo &fbmsMesh, const pgo::Mesh::TriMeshGeo &sphereMesh,
  double fbmsThickness, double sphereThickness, double paddingRatio,
  V3d &bmin, V3d &bmax)
{
  const pgo::Mesh::BoundingBox fbmsBBox(fbmsMesh.positions());
  const pgo::Mesh::BoundingBox sphereBBox(sphereMesh.positions());

  bmin = fbmsBBox.bmin().cwiseMin(sphereBBox.bmin());
  bmax = fbmsBBox.bmax().cwiseMax(sphereBBox.bmax());

  const V3d baseSides = bmax - bmin;
  const double baseDiag = baseSides.norm();
  const double expansion = std::max(fbmsThickness, sphereThickness) + paddingRatio * baseDiag;

  bmin.array() -= expansion;
  bmax.array() += expansion;

  const V3d sides = bmax - bmin;
  for (int axis = 0; axis < 3; ++axis) {
    if (sides[axis] <= 0.0 || !std::isfinite(sides[axis]))
      throw std::runtime_error("Expanded grid bbox has a non-positive or invalid side length");
  }
}

void computeFBMSDistance(const pgo::Mesh::TriMeshGeo &fbmsMesh,
  const V3d &bmin, const V3d &bmax, int resolution,
  pgo::EigenSupport::VXd &fbmsDistance)
{
  pgo::libiglInterface::computeDistanceField(
    fbmsMesh, bmin, bmax, resolution,
    /*robust=*/1, /*sign=*/0, fbmsDistance);
}

void assembleUnionField(const pgo::EigenSupport::VXd &fbmsDistance, const SphereParameters &sphere,
  const V3d &bmin, const V3d &bmax, int resolution,
  double fbmsThickness, double sphereThickness, bool enableTruncating,
  pgo::EigenSupport::VXd &unionField, double &fieldMin, double &fieldMax)
{
  unionField.resize(fbmsDistance.size());
  const V3d delta = (bmax - bmin) / static_cast<double>(resolution - 1);

  // Outer surface of the thickened sphere shell. When --enable-truncating is on,
  // the FBMS shell is CSG-intersected with the solid ball of this radius, so any
  // FBMS structure outside the sphere's outer surface is dropped.
  const double truncationRadius = sphere.radius + sphereThickness * 0.5;

  fieldMin = std::numeric_limits<double>::infinity();
  fieldMax = -std::numeric_limits<double>::infinity();

  for (int z = 0; z < resolution; ++z) {
    for (int y = 0; y < resolution; ++y) {
      for (int x = 0; x < resolution; ++x) {
        const int index = z * resolution * resolution + y * resolution + x;
        const V3d p = bmin + delta.cwiseProduct(V3d(x, y, z).cast<double>());
        const double radial = (p - sphere.center).norm();
        double fbmsField = fbmsDistance[index] - fbmsThickness * 0.5;
        if (enableTruncating) {
          const double ballSDF = radial - truncationRadius;
          fbmsField = std::max(fbmsField, ballSDF);
        }
        const double sphereShellField = std::abs(radial - sphere.radius) - sphereThickness * 0.5;
        const double value = std::min(fbmsField, sphereShellField);
        unionField[index] = value;
        fieldMin = std::min(fieldMin, value);
        fieldMax = std::max(fieldMax, value);
      }
    }
  }
}

double computeMeshVolume(const pgo::Mesh::TriMeshGeo &mesh)
{
  // Signed volume of a closed polyhedron via the divergence theorem on triangles:
  // V = (1/6) * sum_t (a . (b x c)). Sign depends on triangle winding; absolute
  // value gives the enclosed volume for a watertight mesh.
  double v6 = 0.0;
  const int numTri = mesh.numTriangles();
  for (int t = 0; t < numTri; ++t) {
    const pgo::Vec3d &a = mesh.pos(t, 0);
    const pgo::Vec3d &b = mesh.pos(t, 1);
    const pgo::Vec3d &c = mesh.pos(t, 2);
    v6 += a.dot(b.cross(c));
  }
  return std::abs(v6) / 6.0;
}

struct ThicknessSearchResult
{
  double thickness = 0.0;
  double volume = 0.0;
  pgo::Mesh::TriMeshGeo mesh;
  bool budgetExceededAtMin = false;
  bool budgetUnreachedAtMax = false;
  int iterations = 0;
};

// Finds the largest fbmsThickness in [tLo, tHi] whose union mesh has volume <= targetVolume.
// fbmsDistance is the precomputed unsigned distance field on the fixed grid; the grid bbox
// must be sized to accommodate tHi.
ThicknessSearchResult findThicknessForVolumeBudget(
  const pgo::EigenSupport::VXd &fbmsDistance, const SphereParameters &sphere,
  const V3d &bmin, const V3d &bmax, int resolution,
  double sphereThickness, bool enableTruncating, double targetVolume,
  double tLo, double tHi, int maxIterations, double relTol)
{
  pgo::EigenSupport::VXd unionField;
  double fmin = 0.0;
  double fmax = 0.0;

  auto evalAt = [&](double t, pgo::Mesh::TriMeshGeo &outMesh) -> double {
    assembleUnionField(fbmsDistance, sphere, bmin, bmax, resolution,
      t, sphereThickness, enableTruncating, unionField, fmin, fmax);
    if (!(fmin <= 0.0 && fmax >= 0.0)) {
      // No zero crossing: nothing to mesh; treat as volume 0.
      outMesh = pgo::Mesh::TriMeshGeo();
      return 0.0;
    }
    pgo::libiglInterface::computeMarchingCubes(bmin, bmax, resolution, unionField, outMesh);
    if (outMesh.numTriangles() == 0)
      return 0.0;
    return computeMeshVolume(outMesh);
  };

  ThicknessSearchResult result;
  pgo::Mesh::TriMeshGeo meshLo, meshHi, meshMid;

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
  result.mesh = meshLo;

  for (int it = 0; it < maxIterations; ++it) {
    const double mid = 0.5 * (tLo + tHi);
    const double vMid = evalAt(mid, meshMid);
    std::cout << "[volume-search] iter=" << it << " t=" << mid << " vol=" << vMid << std::endl;

    if (vMid <= targetVolume) {
      tLo = mid;
      result.thickness = mid;
      result.volume = vMid;
      result.mesh = meshMid;
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

void saveSurfaceOrThrow(const pgo::Mesh::TriMeshGeo &mesh, const std::string &path)
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

}  // namespace

int main(int argc, char *argv[])
{
  try {
    const Options options = parseOptions(argc, argv);
    validateOptions(options);

    pgo::Logging::init();

    pgo::Mesh::TriMeshGeo fbmsMesh;
    pgo::Mesh::TriMeshGeo sphereMesh;
    loadMeshOrThrow(options.fbmsPath, fbmsMesh, "FBMS");
    loadMeshOrThrow(options.spherePath, sphereMesh, "sphere");

    if (fbmsMesh.numTriangles() == 0)
      throw std::runtime_error("FBMS mesh has no triangles: " + options.fbmsPath);

    const SphereParameters sphere = loadSphereParameters(options.spherePath, sphereMesh);
    std::cout << "Sphere parameter source: " << (sphere.fromHeader ? "header" : "bbox fallback") << std::endl;
    printVector("Sphere center", sphere.center);
    std::cout << "Sphere radius = " << sphere.radius << std::endl;

    if (options.enableTruncating)
      std::cout << "Truncation enabled: FBMS clipped to ball of radius "
                << (sphere.radius + options.sphereThickness * 0.5) << std::endl;

    pgo::Mesh::TriMeshGeo rawSurface;
    V3d bmin, bmax;

    if (options.fbmsThickness > 0.0) {
      computeUnionBBox(fbmsMesh, sphereMesh, options.fbmsThickness,
        options.sphereThickness, options.paddingRatio, bmin, bmax);

      pgo::EigenSupport::VXd fbmsDistance;
      computeFBMSDistance(fbmsMesh, bmin, bmax, options.resolution, fbmsDistance);

      pgo::EigenSupport::VXd unionField;
      double fieldMin = 0.0;
      double fieldMax = 0.0;
      assembleUnionField(fbmsDistance, sphere, bmin, bmax, options.resolution,
        options.fbmsThickness, options.sphereThickness, options.enableTruncating,
        unionField, fieldMin, fieldMax);

      if (!(fieldMin <= 0.0 && fieldMax >= 0.0))
        throw std::runtime_error("Union field does not cross the zero isosurface");

      pgo::libiglInterface::computeMarchingCubes(bmin, bmax, options.resolution, unionField, rawSurface);

      std::cout << "Grid resolution = " << options.resolution << std::endl;
      printVector("BBox min", bmin);
      printVector("BBox max", bmax);
      std::cout << "Field min = " << fieldMin << std::endl;
      std::cout << "Field max = " << fieldMax << std::endl;
    }
    else {
      // Volume-budget mode: |options.fbmsThickness| = target volume.
      const double targetVolume = -options.fbmsThickness;

      const pgo::Mesh::BoundingBox fbmsBBox(fbmsMesh.positions());
      const pgo::Mesh::BoundingBox sphereBBox(sphereMesh.positions());
      const V3d baseSides =
        fbmsBBox.bmax().cwiseMax(sphereBBox.bmax())
        - fbmsBBox.bmin().cwiseMin(sphereBBox.bmin());
      const double baseDiag = baseSides.norm();

      // Search ceiling. With truncating, volume saturates once t exceeds the
      // ball diameter; without truncating, t can grow until the union fills
      // a sizable fraction of the union bbox.
      const double tHiCap = options.enableTruncating
        ? std::max(2.0 * sphere.radius, 4.0 * options.sphereThickness)
        : std::max(0.5 * baseDiag, 4.0 * options.sphereThickness);

      // Size the grid so the precomputed FBMS UDF is valid across the whole
      // search range. With truncating, the meshable region is clipped to the
      // ball, so we only need padding for the sphere shell. Otherwise we must
      // accommodate the search ceiling.
      const double bboxThickness = options.enableTruncating ? options.sphereThickness : tHiCap;
      computeUnionBBox(fbmsMesh, sphereMesh, bboxThickness,
        options.sphereThickness, options.paddingRatio, bmin, bmax);

      pgo::EigenSupport::VXd fbmsDistance;
      computeFBMSDistance(fbmsMesh, bmin, bmax, options.resolution, fbmsDistance);

      const V3d voxel = (bmax - bmin) / static_cast<double>(options.resolution - 1);
      const double tLo = 2.0 * voxel.maxCoeff();

      std::cout << "[volume-search] target volume = " << targetVolume
                << ", search range t in [" << tLo << ", " << tHiCap << "]" << std::endl;

      const ThicknessSearchResult result = findThicknessForVolumeBudget(
        fbmsDistance, sphere, bmin, bmax, options.resolution,
        options.sphereThickness, options.enableTruncating, targetVolume,
        tLo, tHiCap, /*maxIterations=*/20, /*relTol=*/1e-3);

      if (result.budgetExceededAtMin)
        std::cout << "[volume-search] WARNING: even t=" << result.thickness
                  << " produces volume " << result.volume
                  << " > target " << targetVolume
                  << "; saving the minimum-thickness mesh anyway." << std::endl;
      else if (result.budgetUnreachedAtMax)
        std::cout << "[volume-search] WARNING: top of search range t=" << result.thickness
                  << " still has volume " << result.volume
                  << " <= target " << targetVolume
                  << "; budget is not exhausted." << std::endl;

      std::cout << "[volume-search] selected thickness = " << result.thickness
                << " volume = " << result.volume
                << " iterations = " << result.iterations << std::endl;

      rawSurface = result.mesh;

      std::cout << "Grid resolution = " << options.resolution << std::endl;
      printVector("BBox min", bmin);
      printVector("BBox max", bmax);
    }

    std::cout << "Raw vertex count = " << rawSurface.numVertices() << std::endl;
    std::cout << "Raw face count = " << rawSurface.numTriangles() << std::endl;

    if (rawSurface.numVertices() == 0 || rawSurface.numTriangles() == 0)
      throw std::runtime_error("Marching cubes produced an empty raw surface");

    saveSurfaceOrThrow(rawSurface, options.outputSurfacePath);
    return 0;
  }
  catch (const std::exception &err) {
    std::cerr << "Error: " << err.what() << std::endl;
    return 1;
  }
}
