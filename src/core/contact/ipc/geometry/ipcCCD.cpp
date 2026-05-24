/*
copyright to Bohan Wang
*/

#include "ipcCCD.h"

#include "ipcDistancePrimitives.h"

#include <algorithm>
#include <cmath>
#include <functional>

namespace pgo {
namespace Contact {
namespace IPC {
// =========================================================================
//  CCD — Additive Continuous Collision Detection (ACCD)
//
//  Based on Codim-IPC [Li et al. 2021]. Uses conservative advancement:
//  iteratively steps forward by a safe amount based on current distance
//  and maximum displacement magnitude. Much more robust than cubic-solver
//  CCD in floating point.
// =========================================================================
namespace ccd {

static constexpr long ACCD_MAX_ITER = 1000000;
static constexpr double ACCD_CONSERVATIVE_RESCALING = 0.1;

// Core ACCD algorithm: given stacked positions x, displacement dx,
// a distance function, and max displacement magnitude, find TOI.
// Returns true if collision found, false otherwise.
static bool additiveCCD(
  V12d x,  // mutable copy of stacked positions
  const V12d &dx,
  const std::function<double(const V12d &)> &distanceSquared,
  double maxDispMag,
  double &toi,
  double minDistance,
  double tmax)
{
  const double eta = ACCD_CONSERVATIVE_RESCALING;
  const double minDistSq = minDistance * minDistance;

  double d_sq = distanceSquared(x);
  double d = std::sqrt(d_sq);
  if (d <= minDistance) {
    toi = 0.0;
    return true;
  }

  double d_func = d_sq - minDistSq;
  // gap = η * (d² - ξ²) / (d + ξ)
  const double gap = eta * d_func / (d + minDistance);

  toi = 0.0;
  for (long i = 0; i < ACCD_MAX_ITER; ++i) {
    // Conservative lower bound on step:
    //   t_l = (1 - η) * (d² - ξ²) / ((d + ξ) * l_p)
    double toi_lower = (1.0 - eta) * d_func / ((d + minDistance) * maxDispMag);

    x += toi_lower * dx;

    d = std::sqrt(d_sq = distanceSquared(x));
    d_func = d_sq - minDistSq;

    if (d_func <= 0.0) {
      // Shouldn't happen with conservative stepping, but handle it
      break;
    }

    if (toi > 0.0 && d_func / (d + minDistance) < gap) {
      break;  // distance has decreased below gap threshold
    }

    toi += toi_lower;
    if (toi > tmax) {
      return false;  // no collision within tmax
    }
  }

  return true;
}

// Helper: subtract mean displacement for translation invariance
template<typename... Args>
static void subtractMean(Args &...args)
{
  constexpr int n = sizeof...(args);
  V3d mean = V3d::Zero();
  for (const V3d &v : { args... })
    mean += v;
  mean /= n;
  for (V3d *v : { &args... })
    *v -= mean;
}

double pointTriangleCCD(const V3d &p, const V3d &t0,
  const V3d &t1, const V3d &t2,
  const V3d &dp, const V3d &dt0,
  const V3d &dt1, const V3d &dt2,
  double thickness, double tMax)
{
  // Check initial distance
  double initDist = distance::computePTSqDist(p, t0, t1, t2);
  if (initDist <= thickness * thickness) {
    return 0.0;
  }

  // Compute displacements and subtract mean for translation invariance
  V3d ddp = dp;
  V3d ddt0 = dt0;
  V3d ddt1 = dt1;
  V3d ddt2 = dt2;
  subtractMean(ddp, ddt0, ddt1, ddt2);

  // Max displacement magnitude
  double maxDisp = ddp.norm() + std::sqrt(std::max({ ddt0.squaredNorm(), ddt1.squaredNorm(), ddt2.squaredNorm() }));
  if (maxDisp < 1e-30)
    return tMax;

  // Stack positions and displacements
  V12d x, dx;
  x.segment<3>(0) = p;
  x.segment<3>(3) = t0;
  x.segment<3>(6) = t1;
  x.segment<3>(9) = t2;

  dx.segment<3>(0) = ddp;
  dx.segment<3>(3) = ddt0;
  dx.segment<3>(6) = ddt1;
  dx.segment<3>(9) = ddt2;

  auto distSq = [](const V12d &xx) -> double {
    return distance::computePTSqDist(
      xx.segment<3>(0), xx.segment<3>(3),
      xx.segment<3>(6), xx.segment<3>(9));
  };

  double toi;
  if (additiveCCD(x, dx, distSq, maxDisp, toi, thickness, tMax)) {
    return toi;
  }
  return tMax;
}

double edgeEdgeCCD(const V3d &ea0, const V3d &ea1,
  const V3d &eb0, const V3d &eb1,
  const V3d &dea0, const V3d &dea1,
  const V3d &deb0, const V3d &deb1,
  double thickness, double tMax)
{
  // Check initial distance
  double initDist = distance::computeEESqDist(ea0, ea1, eb0, eb1);
  if (initDist <= thickness * thickness) {
    return 0.0;
  }

  // Compute displacements and subtract mean
  V3d ddea0 = dea0, ddea1 = dea1;
  V3d ddeb0 = deb0, ddeb1 = deb1;
  subtractMean(ddea0, ddea1, ddeb0, ddeb1);

  double maxDisp =
    std::sqrt(std::max(ddea0.squaredNorm(), ddea1.squaredNorm())) + std::sqrt(std::max(ddeb0.squaredNorm(), ddeb1.squaredNorm()));
  if (maxDisp < 1e-30)
    return tMax;

  double minDistSq = thickness * thickness;
  V12d x, dx;
  x.segment<3>(0) = ea0;
  x.segment<3>(3) = ea1;
  x.segment<3>(6) = eb0;
  x.segment<3>(9) = eb1;
  dx.segment<3>(0) = ddea0;
  dx.segment<3>(3) = ddea1;
  dx.segment<3>(6) = ddeb0;
  dx.segment<3>(9) = ddeb1;

  auto distSq = [minDistSq](const V12d &xx) -> double {
    double d_sq = distance::computeEESqDist(
      xx.segment<3>(0), xx.segment<3>(3),
      xx.segment<3>(6), xx.segment<3>(9));

    if (d_sq - minDistSq <= 0) {
      // Near-parallel edges: fall back to minimum PP distance
      const auto &a0 = xx.segment<3>(0);
      const auto &a1 = xx.segment<3>(3);
      const auto &b0 = xx.segment<3>(6);
      const auto &b1 = xx.segment<3>(9);
      d_sq = std::min({ (a0 - b0).squaredNorm(), (a0 - b1).squaredNorm(),
        (a1 - b0).squaredNorm(), (a1 - b1).squaredNorm() });
    }
    return d_sq;
  };

  double toi;
  if (additiveCCD(x, dx, distSq, maxDisp, toi, thickness, tMax)) {
    return toi;
  }
  return tMax;
}

}  // namespace ccd

}  // namespace IPC
}  // namespace Contact
}  // namespace pgo
