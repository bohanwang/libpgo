#include "generateFBMSUnionSurfaceAdapter.h"

#include "field/gridSpec.h"
#include "field/denseGrid.h"
#include "geometry/meshDistance.h"
#include "geometry/shellThickening.h"
#include "geometry/sphereField.h"
#include "operations/booleanOps.h"
#include "extraction/marchingCubesExtractor.h"

#ifdef PGO_HAS_OPENVDB
#include "extraction/openVDBExtractor.h"
#include <openvdb/openvdb.h>
#endif

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <optional>
#include <regex>
#include <stdexcept>
#include <vector>

namespace pgo::Tools::FBMSUnionSurface
{

namespace {

// FBMS-specific OBJ header parsing. The FBMS pipeline embeds sphere metadata as
// OBJ comments:  # center = (x, y, z)  and  # radius = value.
// This format is an FBMS convention and does not belong in the neutral core.

std::optional<EigenSupport::V3d> parseCenterLine(const std::string &line)
{
  static const std::regex centerRegex(
    R"(^\s*#\s*center\s*=\s*\(\s*([^,\s]+)\s*,\s*([^,\s]+)\s*,\s*([^,\s\)]+)\s*\)\s*$)",
    std::regex::icase);

  std::smatch match;
  if (!std::regex_match(line, match, centerRegex))
    return std::nullopt;

  return EigenSupport::V3d(
    std::stod(match[1].str()),
    std::stod(match[2].str()),
    std::stod(match[3].str()));
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

struct FBMSSphereHeader {
  EigenSupport::V3d center = EigenSupport::V3d::Zero();
  double radius = 0.0;
  bool valid = false;
};

FBMSSphereHeader parseFBMSHeader(const std::string &spherePath)
{
  std::ifstream in(spherePath);
  if (!in.is_open())
    return {};

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

  if (center && radius && *radius > 0.0 && std::isfinite(*radius)) {
    FBMSSphereHeader h;
    h.center = *center;
    h.radius = *radius;
    h.valid = true;
    return h;
  }
  return {};
}

}  // anonymous namespace

// ---- Sphere operations ----

SphereParameters loadSphereParameters(const std::string &spherePath, const Mesh::TriMeshGeo &sphereMesh)
{
  FBMSSphereHeader header = parseFBMSHeader(spherePath);

  SphereParameters params;
  if (header.valid) {
    params.center = header.center;
    params.radius = header.radius;
    params.fromHeader = true;
  }
  else {
    pgo::ImplicitSurface::SphereField sf =
      pgo::ImplicitSurface::computeSphereFieldFromBBox(sphereMesh);
    params.center = sf.center;
    params.radius = sf.radius;
    params.fromHeader = false;
  }
  return params;
}

int projectBoundaryVerticesToSphere(Mesh::TriMeshGeo &mesh, const SphereParameters &sphere)
{
  pgo::ImplicitSurface::SphereField sf;
  sf.center = sphere.center;
  sf.radius = sphere.radius;
  return pgo::ImplicitSurface::projectOpenBoundaryToSphere(mesh, sf);
}

// ---- Field operations ----

void computeUnionBBox(const Mesh::TriMeshGeo &fbmsMesh, const Mesh::TriMeshGeo &sphereMesh,
  double fbmsThickness, double sphereThickness, double paddingRatio,
  EigenSupport::V3d &bmin, EigenSupport::V3d &bmax)
{
  pgo::Mesh::computeUnionBBox(
    fbmsMesh, sphereMesh, fbmsThickness, sphereThickness, paddingRatio, bmin, bmax);
}

void computeFBMSDistance(const Mesh::TriMeshGeo &fbmsMesh,
  const EigenSupport::V3d &bmin, const EigenSupport::V3d &bmax, int resolution,
  EigenSupport::VXd &fbmsDistance)
{

  pgo::ImplicitSurface::GridSpec spec;
  spec.bmin = bmin;
  spec.bmax = bmax;
  spec.resolution = resolution;

  pgo::ImplicitSurface::DenseGrid grid(spec);

  pgo::ImplicitSurface::computeMeshUnsignedDistance(fbmsMesh, spec, grid);

  fbmsDistance.resize(grid.size());
  for (int i = 0; i < grid.size(); ++i)
    fbmsDistance[i] = grid[i];
}

void assembleUnionField(const EigenSupport::VXd &fbmsDistance, const SphereParameters &sphere,
  const EigenSupport::V3d &bmin, const EigenSupport::V3d &bmax, int resolution,
  double fbmsThickness, double sphereThickness, bool enableTruncating, const std::string &debugFieldMode,
  EigenSupport::VXd &unionField, double &fieldMin, double &fieldMax)
{

  pgo::ImplicitSurface::GridSpec spec;
  spec.bmin = bmin;
  spec.bmax = bmax;
  spec.resolution = resolution;

  // Copy FBMS distance data into a DenseGrid
  pgo::ImplicitSurface::DenseGrid meshDistGrid(spec);
  for (int i = 0; i < fbmsDistance.size(); ++i)
    meshDistGrid[i] = fbmsDistance[i];

  // Convert SphereParameters to core SphereField
  pgo::ImplicitSurface::SphereField sphereField;
  sphereField.center = sphere.center;
  sphereField.radius = sphere.radius;

  // Output grid
  pgo::ImplicitSurface::DenseGrid outGrid(spec);

  // Build mesh shell field, with optional ball-SDF truncation
  pgo::ImplicitSurface::DenseGrid shellGrid(spec);
  pgo::ImplicitSurface::thickenMeshShell(meshDistGrid, fbmsThickness, shellGrid);

  if (enableTruncating) {
    pgo::ImplicitSurface::DenseGrid ballGrid(spec);
    pgo::ImplicitSurface::evaluateBallSDF(sphereField, spec, ballGrid);
    pgo::ImplicitSurface::DenseGrid truncatedGrid(spec);
    pgo::ImplicitSurface::applyBoolean(shellGrid, ballGrid, pgo::ImplicitSurface::BooleanOp::Intersection, truncatedGrid);
    shellGrid = std::move(truncatedGrid);
  }

  // Build sphere shell field
  pgo::ImplicitSurface::DenseGrid sphereShellGrid(spec);
  pgo::ImplicitSurface::thickenSphereShell(sphereField, sphereThickness, spec, sphereShellGrid);

  // Compose according to the FBMS mode
  if (debugFieldMode == "fbms") {
    outGrid = std::move(shellGrid);
  }
  else if (debugFieldMode == "sphere") {
    outGrid = std::move(sphereShellGrid);
  }
  else if (debugFieldMode == "union-minus-sphere") {
    pgo::ImplicitSurface::DenseGrid tmp(spec);
    pgo::ImplicitSurface::applyBoolean(shellGrid, sphereShellGrid, pgo::ImplicitSurface::BooleanOp::Union, tmp);
    pgo::ImplicitSurface::applyBoolean(tmp, sphereShellGrid, pgo::ImplicitSurface::BooleanOp::Difference, outGrid);
  }
  else {
    // "union" (and default)
    pgo::ImplicitSurface::applyBoolean(shellGrid, sphereShellGrid, pgo::ImplicitSurface::BooleanOp::Union, outGrid);
  }

  // Copy output back to VXd and compute min/max
  unionField.resize(outGrid.size());
  fieldMin = std::numeric_limits<double>::infinity();
  fieldMax = -std::numeric_limits<double>::infinity();
  for (int i = 0; i < outGrid.size(); ++i) {
    const double val = outGrid[i];
    unionField[i] = val;
    if (val < fieldMin) fieldMin = val;
    if (val > fieldMax) fieldMax = val;
  }
}

namespace
{

double computeAutoIsoOffsetLimit(double fbmsThickness, double sphereThickness,
  const EigenSupport::V3d &bmin, const EigenSupport::V3d &bmax, int resolution)
{
  const EigenSupport::V3d voxel = (bmax - bmin) / static_cast<double>(resolution - 1);
  const double voxelLimit = 0.25 * voxel.minCoeff();
  const double thicknessLimit = 0.1 * std::min(fbmsThickness, sphereThickness);
  return std::max(0.0, std::min(voxelLimit, thicknessLimit));
}

double selectAutoIsoOffset(const EigenSupport::VXd &field, double limit)
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

}  // anonymous namespace

double selectIsoOffset(const std::string &isoOffsetMode, double fixedIsoOffset,
  const EigenSupport::VXd &field, double fbmsThickness, double sphereThickness,
  const EigenSupport::V3d &bmin, const EigenSupport::V3d &bmax, int resolution)
{
  if (isoOffsetMode == "zero")
    return 0.0;
  if (isoOffsetMode == "fixed")
    return fixedIsoOffset;

  const double limit = computeAutoIsoOffsetLimit(fbmsThickness, sphereThickness, bmin, bmax, resolution);
  return selectAutoIsoOffset(field, limit);
}

// ---- Extraction operations ----

void extractSurfaceMarchingCubes(const EigenSupport::V3d &bmin, const EigenSupport::V3d &bmax,
  int resolution, const EigenSupport::VXd &field, double isoOffset, Mesh::TriMeshGeo &outMesh)
{

  pgo::ImplicitSurface::GridSpec spec;
  spec.bmin = bmin;
  spec.bmax = bmax;
  spec.resolution = resolution;

  pgo::ImplicitSurface::DenseGrid grid(spec);
  for (int i = 0; i < field.size(); ++i)
    grid[i] = field[i];

  pgo::ImplicitSurface::MarchingCubesOptions mcOptions;
  mcOptions.isoOffset = isoOffset;

  pgo::ImplicitSurface::extractMarchingCubes(grid, mcOptions, outMesh);
}

#ifdef PGO_HAS_OPENVDB
void extractSurfaceOpenVDB(const Mesh::TriMeshGeo &fbmsMesh, const SphereParameters &sphere,
  double fbmsThickness, double sphereThickness, bool enableTruncating, const std::string &debugFieldMode,
  double voxelSize, double halfWidth, double adaptivity, int smoothSteps,
  Mesh::TriMeshGeo &outMesh)
{
  openvdb::initialize();

  pgo::ImplicitSurface::OpenVDBOptions options;
  options.voxelSize = voxelSize;
  options.halfWidth = halfWidth;
  options.adaptivity = adaptivity;
  options.smoothSteps = smoothSteps;

  pgo::ImplicitSurface::SphereField sphereField;
  sphereField.center = sphere.center;
  sphereField.radius = sphere.radius;

  std::unique_ptr<pgo::ImplicitSurface::OpenVDBLevelSet> grid;

  if (debugFieldMode != "sphere") {
    grid = pgo::ImplicitSurface::buildOpenVDBShellFromMesh(fbmsMesh, fbmsThickness, options);
    if (enableTruncating) {
      auto truncationBall = pgo::ImplicitSurface::buildOpenVDBBallLevelSet(sphereField, options);
      grid = pgo::ImplicitSurface::combineOpenVDBLevelSets(
        *grid, *truncationBall, pgo::ImplicitSurface::BooleanOp::Intersection);
    }
  }

  if (debugFieldMode == "sphere") {
    grid = pgo::ImplicitSurface::buildOpenVDBSphereShell(sphereField, sphereThickness, options);
  }
  else if (debugFieldMode == "union") {
    auto sphereShell = pgo::ImplicitSurface::buildOpenVDBSphereShell(sphereField, sphereThickness, options);
    grid = pgo::ImplicitSurface::combineOpenVDBLevelSets(
      *grid, *sphereShell, pgo::ImplicitSurface::BooleanOp::Union);
  }
  else if (debugFieldMode == "union-minus-sphere") {
    auto sphereShellForUnion = pgo::ImplicitSurface::buildOpenVDBSphereShell(sphereField, sphereThickness, options);
    grid = pgo::ImplicitSurface::combineOpenVDBLevelSets(
      *grid, *sphereShellForUnion, pgo::ImplicitSurface::BooleanOp::Union);

    auto sphereShellForDiff = pgo::ImplicitSurface::buildOpenVDBSphereShell(sphereField, sphereThickness, options);
    grid = pgo::ImplicitSurface::combineOpenVDBLevelSets(
      *grid, *sphereShellForDiff, pgo::ImplicitSurface::BooleanOp::Difference);
  }

  if (!grid)
    throw std::runtime_error("OpenVDB backend failed to create a level-set grid");

  pgo::ImplicitSurface::extractOpenVDBLevelSet(*grid, options, outMesh);
}
#else
void extractSurfaceOpenVDB(const Mesh::TriMeshGeo &, const SphereParameters &,
  double, double, bool, const std::string &, double, double, double, int,
  Mesh::TriMeshGeo &)
{
  throw std::runtime_error("OpenVDB backend is unavailable; configure with -DPGO_ENABLE_OPENVDB=ON");
}
#endif

// ---- Volume operations ----

double computeMeshVolume(const Mesh::TriMeshGeo &mesh)
{
  return pgo::Mesh::computeMeshVolume(mesh);
}

ThicknessSearchResult findThicknessForVolumeBudget(
  const EigenSupport::VXd &fbmsDistance, const SphereParameters &sphere,
  const EigenSupport::V3d &bmin, const EigenSupport::V3d &bmax, int resolution,
  double sphereThickness, bool enableTruncating, const std::string &isoOffsetMode, double fixedIsoOffset,
  double targetVolume,
  double tLo, double tHi, int maxIterations, double relTol)
{
  EigenSupport::VXd unionField;
  double fmin = 0.0;
  double fmax = 0.0;

  auto evalAt = [&](double t, Mesh::TriMeshGeo &outMesh) -> double {
    assembleUnionField(fbmsDistance, sphere, bmin, bmax, resolution,
      t, sphereThickness, enableTruncating, "union", unionField, fmin, fmax);
    const double isoOffset = selectIsoOffset(isoOffsetMode, fixedIsoOffset,
      unionField, t, sphereThickness, bmin, bmax, resolution);
    if (!(fmin <= isoOffset && fmax >= isoOffset)) {
      outMesh = Mesh::TriMeshGeo();
      return 0.0;
    }
    extractSurfaceMarchingCubes(bmin, bmax, resolution, unionField, isoOffset, outMesh);
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

#ifdef PGO_HAS_OPENVDB
ThicknessSearchResult findThicknessForVolumeBudgetOpenVDB(
  const Mesh::TriMeshGeo &fbmsMesh, const SphereParameters &sphere,
  double sphereThickness, bool enableTruncating,
  double voxelSize, double halfWidth, double adaptivity, int smoothSteps,
  double targetVolume, double tLo, double tHi, int maxIterations, double relTol)
{
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
ThicknessSearchResult findThicknessForVolumeBudgetOpenVDB(
  const Mesh::TriMeshGeo &, const SphereParameters &,
  double, bool, double, double, double, int,
  double, double, double, int, double)
{
  throw std::runtime_error("OpenVDB backend is unavailable; configure with -DPGO_ENABLE_OPENVDB=ON");
}
#endif

}  // namespace pgo::Tools::FBMSUnionSurface
