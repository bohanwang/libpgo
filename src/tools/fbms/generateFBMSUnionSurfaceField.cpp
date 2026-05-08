#include "generateFBMSUnionSurfaceField.h"

#include "boundingBox.h"
#include "libiglInterface.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace pgo::Tools::FBMSUnionSurface
{

void computeUnionBBox(const Mesh::TriMeshGeo &fbmsMesh, const Mesh::TriMeshGeo &sphereMesh,
  double fbmsThickness, double sphereThickness, double paddingRatio,
  EigenSupport::V3d &bmin, EigenSupport::V3d &bmax)
{
  const Mesh::BoundingBox fbmsBBox(fbmsMesh.positions());
  const Mesh::BoundingBox sphereBBox(sphereMesh.positions());

  bmin = fbmsBBox.bmin().cwiseMin(sphereBBox.bmin());
  bmax = fbmsBBox.bmax().cwiseMax(sphereBBox.bmax());

  const EigenSupport::V3d baseSides = bmax - bmin;
  const double baseDiag = baseSides.norm();
  const double expansion = std::max(fbmsThickness, sphereThickness) + paddingRatio * baseDiag;

  bmin.array() -= expansion;
  bmax.array() += expansion;

  const EigenSupport::V3d sides = bmax - bmin;
  for (int axis = 0; axis < 3; ++axis) {
    if (sides[axis] <= 0.0 || !std::isfinite(sides[axis]))
      throw std::runtime_error("Expanded grid bbox has a non-positive or invalid side length");
  }
}

void computeFBMSDistance(const Mesh::TriMeshGeo &fbmsMesh,
  const EigenSupport::V3d &bmin, const EigenSupport::V3d &bmax, int resolution,
  EigenSupport::VXd &fbmsDistance)
{
  libiglInterface::computeDistanceField(
    fbmsMesh, bmin, bmax, resolution,
    /*robust=*/1, /*sign=*/0, fbmsDistance);
}

void assembleUnionField(const EigenSupport::VXd &fbmsDistance, const SphereParameters &sphere,
  const EigenSupport::V3d &bmin, const EigenSupport::V3d &bmax, int resolution,
  double fbmsThickness, double sphereThickness, bool enableTruncating, const std::string &debugFieldMode,
  EigenSupport::VXd &unionField, double &fieldMin, double &fieldMax)
{
  unionField.resize(fbmsDistance.size());
  const EigenSupport::V3d delta = (bmax - bmin) / static_cast<double>(resolution - 1);

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
        const EigenSupport::V3d p = bmin + delta.cwiseProduct(EigenSupport::V3d(x, y, z).cast<double>());
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
        else if (debugFieldMode == "union-minus-sphere")
          value = std::max(value, -sphereShellField);
        unionField[index] = value;
        fieldMin = std::min(fieldMin, value);
        fieldMax = std::max(fieldMax, value);
      }
    }
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

}  // namespace

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

}  // namespace pgo::Tools::FBMSUnionSurface
