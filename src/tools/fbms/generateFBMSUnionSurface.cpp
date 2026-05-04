#include "boundingBox.h"
#include "libiglInterface.h"
#include "pgoLogging.h"
#include "triMeshGeo.h"
#include "triMeshNeighbor.h"

#include <argparse/argparse.hpp>

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
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
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
  std::string debugFieldMode = "union";
  std::string extractionBackend = "marching-cubes";
  std::string isoOffsetMode = "zero";
  double fbmsThickness = 0.0;
  double sphereThickness = 0.0;
  double isoOffset = 0.0;
  double paddingRatio = 0.0;
  double vdbVoxelSize = 0.0;
  double vdbHalfWidth = 3.0;
  double vdbAdaptivity = 0.0;
  int resolution = 0;
  int vdbSmoothSteps = 0;
  int minComponentTriangles = 100;
  int keepLargestComponents = -1;
  bool enableTruncating = false;
  bool filterSmallComponents = false;
  bool projectFBMSBoundaryToSphere = false;
};

struct SphereParameters
{
  V3d center = V3d::Zero();
  double radius = 0.0;
  bool fromHeader = false;
};

void addCommonOptions(argparse::ArgumentParser &program)
{
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
    .help("Truncate the FBMS shell by the bounding sphere mid-surface; "
          "only FBMS structure inside the sphere radius is kept.")
    .default_value(false)
    .implicit_value(true);
  program.add_argument("--project-fbms-boundary-to-sphere")
    .help("Before thickening, project open FBMS boundary vertices to the "
          "bounding sphere radius so the sphere shell consumes the thickened "
          "boundary side wall.")
    .default_value(false)
    .implicit_value(true);
  program.add_argument("--debug-field-mode")
    .help("Diagnostic field to polygonize: union, fbms, or sphere")
    .default_value(std::string("union"))
    .metavar("MODE");
  program.add_argument("--filter-small-components")
    .help("Remove tiny edge-connected triangle components after surface extraction")
    .default_value(false)
    .implicit_value(true);
  program.add_argument("--min-component-triangles")
    .help("Minimum triangle count kept by --filter-small-components")
    .default_value(100)
    .metavar("INT")
    .scan<'i', int>();
  program.add_argument("--keep-largest-components")
    .help("With --filter-small-components, optionally keep only the largest N components after thresholding; -1 keeps all components above threshold")
    .default_value(-1)
    .metavar("INT")
    .scan<'i', int>();
}

void addMarchingCubesOptions(argparse::ArgumentParser &program)
{
  program.add_argument("--iso-offset-mode")
    .help("Iso-surface extraction mode: zero, fixed, or auto. zero keeps the historical 0-level extraction. fixed uses --iso-offset. auto chooses a small nearby iso level away from sampled grid values.")
    .default_value(std::string("zero"))
    .metavar("MODE");
  program.add_argument("--iso-offset")
    .help("Iso level used when --iso-offset-mode fixed is selected")
    .default_value(0.0)
    .metavar("FLOAT")
    .scan<'g', double>();
}

void addOpenVDBOptions(argparse::ArgumentParser &program)
{
  program.add_argument("--vdb-voxel-size")
    .help("OpenVDB voxel size. If 0, derive it from --resolution and the expanded bbox.")
    .default_value(0.0)
    .metavar("FLOAT")
    .scan<'g', double>();
  program.add_argument("--vdb-half-width")
    .help("OpenVDB level-set half width in voxels")
    .default_value(3.0)
    .metavar("FLOAT")
    .scan<'g', double>();
  program.add_argument("--vdb-adaptivity")
    .help("OpenVDB volumeToMesh adaptivity in [0, 1]")
    .default_value(0.0)
    .metavar("FLOAT")
    .scan<'g', double>();
  program.add_argument("--vdb-smooth-steps")
    .help("Number of OpenVDB mean-curvature smoothing iterations before meshing")
    .default_value(0)
    .metavar("INT")
    .scan<'i', int>();
}

void readCommonOptions(const argparse::ArgumentParser &program, Options &options)
{
  options.fbmsPath = program.get<std::string>("--fbms");
  options.spherePath = program.get<std::string>("--sphere");
  options.fbmsThickness = program.get<double>("--fbms-thickness");
  options.sphereThickness = program.get<double>("--sphere-thickness");
  options.resolution = program.get<int>("--resolution");
  options.paddingRatio = program.get<double>("--padding-ratio");
  options.outputSurfacePath = program.get<std::string>("--output-surface");
  options.enableTruncating = program.get<bool>("--enable-truncating");
  options.projectFBMSBoundaryToSphere = program.get<bool>("--project-fbms-boundary-to-sphere");
  options.debugFieldMode = program.get<std::string>("--debug-field-mode");
  options.filterSmallComponents = program.get<bool>("--filter-small-components");
  options.minComponentTriangles = program.get<int>("--min-component-triangles");
  options.keepLargestComponents = program.get<int>("--keep-largest-components");
}

void readMarchingCubesOptions(const argparse::ArgumentParser &program, Options &options)
{
  options.isoOffsetMode = program.get<std::string>("--iso-offset-mode");
  options.isoOffset = program.get<double>("--iso-offset");
}

void readOpenVDBOptions(const argparse::ArgumentParser &program, Options &options)
{
  options.vdbVoxelSize = program.get<double>("--vdb-voxel-size");
  options.vdbHalfWidth = program.get<double>("--vdb-half-width");
  options.vdbAdaptivity = program.get<double>("--vdb-adaptivity");
  options.vdbSmoothSteps = program.get<int>("--vdb-smooth-steps");
}

bool isBackendSubcommand(const char *arg)
{
  const std::string value(arg == nullptr ? "" : arg);
  return value == "marching-cubes" || value == "openvdb";
}

Options parseSubcommandOptions(int argc, char *argv[])
{
  argparse::ArgumentParser program("generateFBMSUnionSurface");

  argparse::ArgumentParser marchingCubesCommand("marching-cubes");
  marchingCubesCommand.add_description("Extract the union surface with the historical dense grid marching-cubes backend");
  addCommonOptions(marchingCubesCommand);
  addMarchingCubesOptions(marchingCubesCommand);

  argparse::ArgumentParser openVDBCommand("openvdb");
  openVDBCommand.add_description("Extract the union surface with OpenVDB sparse level-set meshing");
  addCommonOptions(openVDBCommand);
  addOpenVDBOptions(openVDBCommand);

  program.add_subparser(marchingCubesCommand);
  program.add_subparser(openVDBCommand);
  program.parse_args(argc, argv);

  Options options;
  if (program.is_subcommand_used(marchingCubesCommand)) {
    readCommonOptions(marchingCubesCommand, options);
    readMarchingCubesOptions(marchingCubesCommand, options);
    options.extractionBackend = "marching-cubes";
  }
  else if (program.is_subcommand_used(openVDBCommand)) {
    readCommonOptions(openVDBCommand, options);
    readOpenVDBOptions(openVDBCommand, options);
    options.extractionBackend = "openvdb";
  }
  else {
    throw std::runtime_error("Expected subcommand: marching-cubes or openvdb");
  }

  return options;
}

Options parseLegacyOptions(int argc, char *argv[])
{
  argparse::ArgumentParser program("generateFBMSUnionSurface");
  addCommonOptions(program);
  program.add_argument("--extraction-backend")
    .help("Surface extraction backend: marching-cubes or openvdb")
    .default_value(std::string("marching-cubes"))
    .metavar("BACKEND");
  addMarchingCubesOptions(program);
  addOpenVDBOptions(program);

  program.parse_args(argc, argv);

  Options options;
  readCommonOptions(program, options);
  options.extractionBackend = program.get<std::string>("--extraction-backend");
  readMarchingCubesOptions(program, options);
  readOpenVDBOptions(program, options);
  return options;
}

Options parseOptions(int argc, char *argv[])
{
  if (argc > 1 && isBackendSubcommand(argv[1]))
    return parseSubcommandOptions(argc, argv);

  return parseLegacyOptions(argc, argv);
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
  if (options.debugFieldMode != "union" && options.debugFieldMode != "fbms" && options.debugFieldMode != "sphere")
    throw std::runtime_error("--debug-field-mode must be one of: union, fbms, sphere");
  if (options.minComponentTriangles < 1)
    throw std::runtime_error("--min-component-triangles must be at least 1");
  if (options.keepLargestComponents != -1 && options.keepLargestComponents < 1)
    throw std::runtime_error("--keep-largest-components must be -1 or a positive integer");
  if (options.extractionBackend != "marching-cubes" && options.extractionBackend != "openvdb")
    throw std::runtime_error("--extraction-backend must be one of: marching-cubes, openvdb");
  if (options.isoOffsetMode != "zero" && options.isoOffsetMode != "fixed" && options.isoOffsetMode != "auto")
    throw std::runtime_error("--iso-offset-mode must be one of: zero, fixed, auto");
  if (!std::isfinite(options.isoOffset))
    throw std::runtime_error("--iso-offset must be finite");
  if (options.isoOffsetMode != "fixed" && options.isoOffset != 0.0)
    throw std::runtime_error("--iso-offset is only used with --iso-offset-mode fixed");
  if (options.extractionBackend == "openvdb" && options.isoOffsetMode != "zero")
    throw std::runtime_error("--iso-offset-mode is only supported by the marching-cubes backend");
  if (options.extractionBackend == "openvdb" && options.fbmsThickness < 0.0)
    throw std::runtime_error("Volume-budget mode is only supported by the marching-cubes backend");
  if (options.vdbVoxelSize < 0.0 || !std::isfinite(options.vdbVoxelSize))
    throw std::runtime_error("--vdb-voxel-size must be finite and non-negative");
  if (options.vdbHalfWidth <= 0.0 || !std::isfinite(options.vdbHalfWidth))
    throw std::runtime_error("--vdb-half-width must be positive");
  if (options.vdbAdaptivity < 0.0 || options.vdbAdaptivity > 1.0 || !std::isfinite(options.vdbAdaptivity))
    throw std::runtime_error("--vdb-adaptivity must be finite and in [0, 1]");
  if (options.vdbSmoothSteps < 0)
    throw std::runtime_error("--vdb-smooth-steps must be non-negative");

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

int projectBoundaryVerticesToSphere(pgo::Mesh::TriMeshGeo &mesh, const SphereParameters &sphere)
{
  std::map<std::pair<int, int>, int> edgeIncidentCounts;
  for (int triID = 0; triID < mesh.numTriangles(); ++triID) {
    const pgo::Vec3i &tri = mesh.tri(triID);
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

    const V3d radial = mesh.pos(vertexID) - sphere.center;
    const double radialNorm = radial.norm();
    if (radialNorm <= 0.0)
      continue;

    mesh.pos(vertexID) = sphere.center + radial * (sphere.radius / radialNorm);
    ++projectedVertices;
  }

  return projectedVertices;
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
  double fbmsThickness, double sphereThickness, bool enableTruncating, const std::string &debugFieldMode,
  pgo::EigenSupport::VXd &unionField, double &fieldMin, double &fieldMax)
{
  unionField.resize(fbmsDistance.size());
  const V3d delta = (bmax - bmin) / static_cast<double>(resolution - 1);

  // When --enable-truncating is on, the FBMS shell is CSG-intersected with the
  // solid ball at the sphere mid-surface radius. Keeping the clipping cap inside
  // the thickened sphere shell avoids placing it exactly on the shell's outer
  // zero surface, which tends to produce degenerate/non-manifold marching-cubes
  // output.
  const double truncationRadius = sphere.radius;

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
        double value = std::min(fbmsField, sphereShellField);
        if (debugFieldMode == "fbms")
          value = fbmsField;
        else if (debugFieldMode == "sphere")
          value = sphereShellField;
        unionField[index] = value;
        fieldMin = std::min(fieldMin, value);
        fieldMax = std::max(fieldMax, value);
      }
    }
  }
}

double computeAutoIsoOffsetLimit(double fbmsThickness, double sphereThickness,
  const V3d &bmin, const V3d &bmax, int resolution)
{
  const V3d voxel = (bmax - bmin) / static_cast<double>(resolution - 1);
  const double voxelLimit = 0.25 * voxel.minCoeff();
  const double thicknessLimit = 0.1 * std::min(fbmsThickness, sphereThickness);
  return std::max(0.0, std::min(voxelLimit, thicknessLimit));
}

double selectAutoIsoOffset(const pgo::EigenSupport::VXd &field, double limit)
{
  if (limit <= 0.0 || !std::isfinite(limit))
    return 0.0;

  std::vector<double> nearValues;
  for (int i = 0; i < field.size(); ++i) {
    const double value = field[i];
    if (std::isfinite(value) && value >= -limit && value <= limit)
      nearValues.push_back(value);
  }

  if (nearValues.empty())
    return 0.0;

  std::sort(nearValues.begin(), nearValues.end());
  nearValues.erase(std::unique(nearValues.begin(), nearValues.end()), nearValues.end());

  double bestIso = 0.0;
  double bestScore = -1.0;
  double bestClearance = -1.0;
  auto considerGap = [&](double lo, double hi) {
    if (!(lo < hi))
      return;

    const double candidate = 0.5 * (lo + hi);
    const double clearance = 0.5 * (hi - lo);
    const double normalizedDistance = std::abs(candidate) / limit;
    const double score = clearance / (1.0 + normalizedDistance);
    if (score > bestScore ||
      (score == bestScore && std::abs(candidate) < std::abs(bestIso))) {
      bestScore = score;
      bestClearance = clearance;
      bestIso = candidate;
    }
  };

  considerGap(-limit, nearValues.front());
  for (size_t i = 1; i < nearValues.size(); ++i)
    considerGap(nearValues[i - 1], nearValues[i]);
  considerGap(nearValues.back(), limit);

  if (bestClearance <= 0.0)
    return 0.0;

  return bestIso;
}

double selectIsoOffset(const std::string &isoOffsetMode, double fixedIsoOffset,
  const pgo::EigenSupport::VXd &field, double fbmsThickness, double sphereThickness,
  const V3d &bmin, const V3d &bmax, int resolution)
{
  if (isoOffsetMode == "zero")
    return 0.0;
  if (isoOffsetMode == "fixed")
    return fixedIsoOffset;

  const double limit = computeAutoIsoOffsetLimit(fbmsThickness, sphereThickness, bmin, bmax, resolution);
  return selectAutoIsoOffset(field, limit);
}

#ifdef PGO_HAS_OPENVDB
std::vector<openvdb::Vec3s> meshPointsForOpenVDB(const pgo::Mesh::TriMeshGeo &mesh)
{
  std::vector<openvdb::Vec3s> points;
  points.reserve(mesh.numVertices());
  for (int i = 0; i < mesh.numVertices(); ++i) {
    const pgo::Vec3d &p = mesh.pos(i);
    points.emplace_back(static_cast<float>(p[0]), static_cast<float>(p[1]), static_cast<float>(p[2]));
  }
  return points;
}

std::vector<openvdb::Vec3I> meshTrianglesForOpenVDB(const pgo::Mesh::TriMeshGeo &mesh)
{
  std::vector<openvdb::Vec3I> triangles;
  triangles.reserve(mesh.numTriangles());
  for (int i = 0; i < mesh.numTriangles(); ++i) {
    const pgo::Vec3i &tri = mesh.tri(i);
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

openvdb::FloatGrid::Ptr createFBMSShellGrid(const pgo::Mesh::TriMeshGeo &fbmsMesh,
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

pgo::Mesh::TriMeshGeo openVDBMeshToTriMeshGeo(const std::vector<openvdb::Vec3s> &points,
  const std::vector<openvdb::Vec3I> &triangles, const std::vector<openvdb::Vec4I> &quads)
{
  std::vector<pgo::Vec3d> positions;
  positions.reserve(points.size());
  for (const openvdb::Vec3s &p : points)
    positions.emplace_back(p[0], p[1], p[2]);

  std::vector<pgo::Vec3i> outTriangles;
  outTriangles.reserve(triangles.size() + quads.size() * 2);
  for (const openvdb::Vec3I &tri : triangles)
    outTriangles.emplace_back(tri[0], tri[1], tri[2]);
  for (const openvdb::Vec4I &quad : quads) {
    outTriangles.emplace_back(quad[0], quad[1], quad[2]);
    outTriangles.emplace_back(quad[0], quad[2], quad[3]);
  }

  return pgo::Mesh::TriMeshGeo(std::move(positions), std::move(outTriangles));
}

void extractSurfaceOpenVDB(const pgo::Mesh::TriMeshGeo &fbmsMesh, const SphereParameters &sphere,
  double fbmsThickness, double sphereThickness, bool enableTruncating, const std::string &debugFieldMode,
  double voxelSize, double halfWidth, double adaptivity, int smoothSteps,
  pgo::Mesh::TriMeshGeo &outMesh)
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

  if (!grid)
    throw std::runtime_error("OpenVDB backend failed to create a level-set grid");

  for (int i = 0; i < smoothSteps; ++i) {
    openvdb::tools::LevelSetFilter<openvdb::FloatGrid> filter(*grid);
    filter.meanCurvature();
  }

  std::vector<openvdb::Vec3s> points;
  std::vector<openvdb::Vec3I> triangles;
  std::vector<openvdb::Vec4I> quads;
  openvdb::tools::volumeToMesh(*grid, points, triangles, quads,
    /*isovalue=*/0.0, adaptivity, /*relaxDisorientedTriangles=*/true);

  outMesh = openVDBMeshToTriMeshGeo(points, triangles, quads);
}
#else
void extractSurfaceOpenVDB(const pgo::Mesh::TriMeshGeo &, const SphereParameters &,
  double, double, bool, const std::string &, double, double, double, int,
  pgo::Mesh::TriMeshGeo &)
{
  throw std::runtime_error("OpenVDB backend is unavailable; configure with -DPGO_ENABLE_OPENVDB=ON");
}
#endif

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
  double sphereThickness, bool enableTruncating, const std::string &isoOffsetMode, double fixedIsoOffset,
  double targetVolume,
  double tLo, double tHi, int maxIterations, double relTol)
{
  pgo::EigenSupport::VXd unionField;
  double fmin = 0.0;
  double fmax = 0.0;

  auto evalAt = [&](double t, pgo::Mesh::TriMeshGeo &outMesh) -> double {
    assembleUnionField(fbmsDistance, sphere, bmin, bmax, resolution,
      t, sphereThickness, enableTruncating, "union", unionField, fmin, fmax);
    const double isoOffset = selectIsoOffset(isoOffsetMode, fixedIsoOffset,
      unionField, t, sphereThickness, bmin, bmax, resolution);
    if (!(fmin <= isoOffset && fmax >= isoOffset)) {
      // No selected isosurface crossing: nothing to mesh; treat as volume 0.
      outMesh = pgo::Mesh::TriMeshGeo();
      return 0.0;
    }
    pgo::libiglInterface::computeMarchingCubes(bmin, bmax, resolution, unionField, outMesh, isoOffset);
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

    if (options.projectFBMSBoundaryToSphere) {
      const int projectedVertices = projectBoundaryVerticesToSphere(fbmsMesh, sphere);
      std::cout << "Projected FBMS boundary vertices to sphere = " << projectedVertices << std::endl;
    }

    if (options.enableTruncating)
      std::cout << "Truncation enabled: FBMS clipped to ball of radius "
                << sphere.radius << std::endl;

    pgo::Mesh::TriMeshGeo rawSurface;
    V3d bmin, bmax;

    if (options.fbmsThickness > 0.0) {
      computeUnionBBox(fbmsMesh, sphereMesh, options.fbmsThickness,
        options.sphereThickness, options.paddingRatio, bmin, bmax);

      if (options.extractionBackend == "openvdb") {
        const double derivedVoxelSize = (bmax - bmin).maxCoeff() / static_cast<double>(options.resolution - 1);
        const double voxelSize = options.vdbVoxelSize > 0.0 ? options.vdbVoxelSize : derivedVoxelSize;
        extractSurfaceOpenVDB(fbmsMesh, sphere,
          options.fbmsThickness, options.sphereThickness, options.enableTruncating, options.debugFieldMode,
          voxelSize, options.vdbHalfWidth, options.vdbAdaptivity, options.vdbSmoothSteps,
          rawSurface);

        std::cout << "Extraction backend = openvdb" << std::endl;
        std::cout << "Reference grid resolution = " << options.resolution << std::endl;
        printVector("Reference bbox min", bmin);
        printVector("Reference bbox max", bmax);
        std::cout << "OpenVDB voxel size = " << voxelSize << std::endl;
        std::cout << "OpenVDB half width = " << options.vdbHalfWidth << std::endl;
        std::cout << "OpenVDB adaptivity = " << options.vdbAdaptivity << std::endl;
        std::cout << "OpenVDB smooth steps = " << options.vdbSmoothSteps << std::endl;
      }
      else {
        pgo::EigenSupport::VXd fbmsDistance;
        computeFBMSDistance(fbmsMesh, bmin, bmax, options.resolution, fbmsDistance);

        pgo::EigenSupport::VXd unionField;
        double fieldMin = 0.0;
        double fieldMax = 0.0;
        assembleUnionField(fbmsDistance, sphere, bmin, bmax, options.resolution,
          options.fbmsThickness, options.sphereThickness, options.enableTruncating, options.debugFieldMode,
          unionField, fieldMin, fieldMax);

        const double isoOffset = selectIsoOffset(options.isoOffsetMode, options.isoOffset,
          unionField, options.fbmsThickness, options.sphereThickness, bmin, bmax, options.resolution);

        if (!(fieldMin <= isoOffset && fieldMax >= isoOffset))
          throw std::runtime_error("Union field does not cross the selected isosurface");

        pgo::libiglInterface::computeMarchingCubes(bmin, bmax, options.resolution, unionField, rawSurface, isoOffset);

        std::cout << "Extraction backend = marching-cubes" << std::endl;
        std::cout << "Grid resolution = " << options.resolution << std::endl;
        printVector("BBox min", bmin);
        printVector("BBox max", bmax);
        std::cout << "Field min = " << fieldMin << std::endl;
        std::cout << "Field max = " << fieldMax << std::endl;
        std::cout << "Iso offset mode = " << options.isoOffsetMode << std::endl;
        std::cout << "Selected iso offset = " << isoOffset << std::endl;
      }
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
        options.sphereThickness, options.enableTruncating, options.isoOffsetMode, options.isoOffset, targetVolume,
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
      throw std::runtime_error("Surface extraction produced an empty raw surface");

    if (options.filterSmallComponents) {
      const pgo::Mesh::TriangleEdgeConnectivityStats beforeStats =
        pgo::Mesh::computeTriangleEdgeConnectivityStats(rawSurface.triangles());
      pgo::Mesh::TriMeshGeo filteredSurface =
        pgo::Mesh::filterSmallTriangleComponentsByEdge(
          rawSurface, options.minComponentTriangles, options.keepLargestComponents);
      const pgo::Mesh::TriangleEdgeConnectivityStats afterStats =
        pgo::Mesh::computeTriangleEdgeConnectivityStats(filteredSurface.triangles());

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

    saveSurfaceOrThrow(rawSurface, options.outputSurfacePath);
    return 0;
  }
  catch (const std::exception &err) {
    std::cerr << "Error: " << err.what() << std::endl;
    return 1;
  }
}
