/*
copyright to Bohan Wang
*/

#include "ipcCCD.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace pgo {
namespace Contact {
namespace CIPC {
// =========================================================================
//  CCD — Additive Continuous Collision Detection (ACCD)
//
//  Based on Codim-IPC [Li et al. 2021]. Uses conservative advancement:
//  iteratively steps forward by a safe amount based on current distance
//  and maximum displacement magnitude. Much more robust than cubic-solver
//  CCD in floating point.
// =========================================================================
namespace ccd {

using EigenSupport::V3d;
using EigenSupport::V12d;

static constexpr long ACCD_MAX_ITER = 1000000;
static constexpr double ACCD_CONSERVATIVE_RESCALING = 0.1;

// Closest-feature residuals only propose separating axes. The lower bound
// comes from convex support projections, so a misclassified or rounded
// closest feature can never replace the distance by an upper bound.
using WideVector = Eigen::Matrix<long double, 3, 1>;
using WideState = Eigen::Matrix<long double, 12, 1>;
static constexpr long double wideEpsilon = std::numeric_limits<long double>::epsilon();
static constexpr long double poseEpsilon = std::numeric_limits<double>::epsilon();

static WideVector pointSegmentResidual(const WideVector &p, const WideVector &a, const WideVector &b)
{
  const WideVector edge = b - a, relative = p - a;
  const long double length2 = edge.squaredNorm();
  const long double parameter = length2 > 0 ? std::clamp(relative.dot(edge) / length2, 0.0L, 1.0L) : 0;
  return relative - parameter * edge;
}

static long double axisLowerBound(const WideState &positions, const WideState &errors,
  const WideVector &candidate, bool pointTriangle,
  const WideState *endPositions = nullptr, const WideState *endErrors = nullptr)
{
  const long double scale = candidate.cwiseAbs().maxCoeff();
  if (!(scale > 0) || !std::isfinite(scale))
    return 0;
  // Axis construction need not be exact: every represented nonzero direction
  // is a valid support axis. Scaling keeps the subsequent norm well conditioned.
  const WideVector axis = candidate / scale;
  const long double normUpper = std::nextafter(axis.norm() * (1 + 8 * wideEpsilon),
    std::numeric_limits<long double>::infinity());
  long double positive = std::numeric_limits<long double>::infinity();
  long double negative = positive;
  const int firstEnd = pointTriangle ? 1 : 2;
  for (int a = 0; a < firstEnd; ++a) {
    for (int b = firstEnd; b < 4; ++b) {
      for (int ta = 0; ta < (endPositions ? 2 : 1); ++ta) {
        for (int tb = 0; tb < (endPositions ? 2 : 1); ++tb) {
          const WideState &aPositions = ta ? *endPositions : positions, &bPositions = tb ? *endPositions : positions;
          const WideState &aErrors = ta ? *endErrors : errors, &bErrors = tb ? *endErrors : errors;
          const WideVector pa = aPositions.segment<3>(3 * a), pb = bPositions.segment<3>(3 * b);
          const long double projection = (pa - pb).dot(axis);
          // Eight long-double operations cover coordinate subtraction, three
          // products and their sum. Pose radii also cover double interpolation
          // during original-coordinate reconstruction.
          const long double arithmetic = 16 * wideEpsilon *
            ((pa.cwiseAbs() + pb.cwiseAbs()).cwiseProduct(axis.cwiseAbs())).sum();
          const long double uncertainty = arithmetic +
            ((aErrors.segment<3>(3 * a) + bErrors.segment<3>(3 * b)).cwiseProduct(axis.cwiseAbs())).sum() +
            std::numeric_limits<long double>::denorm_min();
          if (!std::isfinite(projection) || !std::isfinite(uncertainty))
            return 0;
          positive = std::min(positive, std::nextafter(projection - uncertainty,
            -std::numeric_limits<long double>::infinity()));
          negative = std::min(negative, std::nextafter(-projection - uncertainty,
            -std::numeric_limits<long double>::infinity()));
        }
      }
    }
  }
  const long double gap = std::max({ positive, negative, 0.0L });
  return gap > 0 ? std::max(0.0L, std::nextafter(gap / normUpper, 0.0L)) : 0;
}

static long double distanceLowerBound(const WideState &positions, const WideState &errors, bool pointTriangle,
  const WideState *endPositions = nullptr, const WideState *endErrors = nullptr)
{
  const WideVector a = positions.segment<3>(0), b = positions.segment<3>(3);
  const WideVector c = positions.segment<3>(6), d = positions.segment<3>(9);
  long double bound;
  if (pointTriangle) {
    bound = axisLowerBound(positions, errors, (c - b).cross(d - b), true, endPositions, endErrors);
    bound = std::max(bound, axisLowerBound(positions, errors, pointSegmentResidual(a, b, c), true, endPositions, endErrors));
    bound = std::max(bound, axisLowerBound(positions, errors, pointSegmentResidual(a, c, d), true, endPositions, endErrors));
    bound = std::max(bound, axisLowerBound(positions, errors, pointSegmentResidual(a, d, b), true, endPositions, endErrors));
  }
  else {
    bound = axisLowerBound(positions, errors, (b - a).cross(d - c), false, endPositions, endErrors);
    bound = std::max(bound, axisLowerBound(positions, errors, pointSegmentResidual(a, c, d), false, endPositions, endErrors));
    bound = std::max(bound, axisLowerBound(positions, errors, pointSegmentResidual(b, c, d), false, endPositions, endErrors));
    bound = std::max(bound, axisLowerBound(positions, errors, pointSegmentResidual(c, a, b), false, endPositions, endErrors));
    bound = std::max(bound, axisLowerBound(positions, errors, pointSegmentResidual(d, a, b), false, endPositions, endErrors));
  }
  return bound;
}

static double additiveCCD(const V12d &original, const V12d &motion,
  bool pointTriangle, double minDistance, double tmax)
{
  if (!original.allFinite() || !motion.allFinite() ||
    !std::isfinite(minDistance) || minDistance < 0 || !std::isfinite(tmax) || tmax < 0)
    return 0;
  const WideState initial = original.cast<long double>();
  const WideState displacement = motion.cast<long double>();
  // The maximum relative endpoint speed bounds every convex-combination
  // relative speed, and removes common translation without modifying poses.
  long double speed = 0;
  const int firstEnd = pointTriangle ? 1 : 2;
  for (int a = 0; a < firstEnd; ++a)
    for (int b = firstEnd; b < 4; ++b) {
      const WideVector relative = displacement.segment<3>(3 * a) - displacement.segment<3>(3 * b);
      const long double magnitude2 = relative.squaredNorm();
      if (!std::isfinite(static_cast<double>(magnitude2)))
        return 0;
      speed = std::max(speed, std::sqrt(magnitude2));
    }
  speed = std::nextafter(speed * (1 + 16 * wideEpsilon), std::numeric_limits<long double>::infinity());
  auto separation = [&](double time) {
    // Reconstruct from the original pose. Incremental x += dt*dx can drift
    // across contact during repeated conservative advancement.
    const WideState positions = initial + static_cast<long double>(time) * displacement;
    const WideState radius = 8 * poseEpsilon * (initial.cwiseAbs() + static_cast<long double>(time) * displacement.cwiseAbs());
    return distanceLowerBound(positions, radius, pointTriangle);
  };
  long double distance = separation(0);
  if (!(distance > minDistance))
    return 0;
  if (speed <= std::numeric_limits<long double>::denorm_min())
    return tmax;
  const long double eta = ACCD_CONSERVATIVE_RESCALING;
  const long double initialGapUpper = std::nextafter(distance - minDistance,
    std::numeric_limits<long double>::infinity());
  const long double reserve = std::nextafter(eta * initialGapUpper,
    std::numeric_limits<long double>::infinity());
  double time = 0, safePrefix = 0;
  for (long iteration = 0; iteration < ACCD_MAX_ITER; ++iteration) {
    const WideState current = initial + static_cast<long double>(time) * displacement;
    const WideState final = initial + static_cast<long double>(tmax) * displacement;
    const WideState currentError = 8 * poseEpsilon * (initial.cwiseAbs() + static_cast<long double>(time) * displacement.cwiseAbs());
    const WideState finalError = 8 * poseEpsilon * (initial.cwiseAbs() + static_cast<long double>(tmax) * displacement.cwiseAbs());
    // A fixed axis separating both primitives throughout the remaining
    // endpoint hull proves the entire remaining linear motion clear.
    if (distanceLowerBound(current, currentError, pointTriangle, &final, &finalError) > minDistance)
      return tmax;
    // Keep advancing internally to prove the entire original chord. If that
    // proof stalls near contact, return a previously certified step whose
    // endpoint retains a finite gap, rather than the near-contact proof pose.
    // The reserve is relative to this query's initial certified clear gap.
    const long double available = std::nextafter(
      std::nextafter(distance - minDistance, 0.0L) - reserve, 0.0L);
    if (available > 0) {
      const long double scaled = std::nextafter((1 - eta) * available, 0.0L);
      const long double advance = std::nextafter(scaled / speed, 0.0L);
      const long double sum = std::nextafter(static_cast<long double>(time) + advance, 0.0L);
      const double candidate = std::nextafter(static_cast<double>(sum), time);
      if (candidate > safePrefix && candidate < tmax)
        safePrefix = candidate;
    }
    const long double lower = std::nextafter((1 - eta) * (distance - minDistance) / speed, 0.0L);
    if (!(lower > 0) || !std::isfinite(lower))
      return safePrefix;
    if (lower > static_cast<long double>(tmax) - time)
      return tmax;
    const double next = std::nextafter(static_cast<double>(static_cast<long double>(time) + lower), time);
    if (!(next > time))
      return safePrefix;
    distance = separation(next);
    if (!(distance > minDistance))
      return safePrefix;
    time = next;
  }
  return safePrefix;
}

double pointTriangleCCD(const V3d &p, const V3d &t0,
  const V3d &t1, const V3d &t2,
  const V3d &dp, const V3d &dt0,
  const V3d &dt1, const V3d &dt2,
  double thickness, double tMax)
{
  V12d positions, displacements;
  positions << p, t0, t1, t2;
  displacements << dp, dt0, dt1, dt2;
  return additiveCCD(positions, displacements, true, thickness, tMax);
}

double edgeEdgeCCD(const V3d &ea0, const V3d &ea1,
  const V3d &eb0, const V3d &eb1,
  const V3d &dea0, const V3d &dea1,
  const V3d &deb0, const V3d &deb1,
  double thickness, double tMax)
{
  V12d positions, displacements;
  positions << ea0, ea1, eb0, eb1;
  displacements << dea0, dea1, deb0, deb1;
  return additiveCCD(positions, displacements, false, thickness, tMax);
}

}  // namespace ccd

}  // namespace CIPC
}  // namespace Contact
}  // namespace pgo
