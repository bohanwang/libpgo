/*
copyright to Bohan Wang
*/
#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <cuda_runtime.h>
#include <math_constants.h>

namespace pgo::Contact::CIPC::CudaIPCContinuousCollision
{
using Vec3 = Eigen::Matrix<double, 3, 1>;

// These intervals contain the exact translated path x_i+t*(v_i-v_0), enlarged
// by the caller's surface-mapping error. Reconstruct from the original path at
// every time: accumulated floating-point pose updates can step across a root.
struct Pose
{
  Vec3 point[4], lower[4], upper[4];
};

__device__ inline Pose evaluate(const Vec3 (&x)[4], const Vec3 (&v)[4],
  const Vec3 (&error)[4], double time)
{
  Pose pose;
  for (int i = 0; i < 4; ++i) {
    for (int k = 0; k < 3; ++k) {
      const double lo = __dsub_rd(v[i][k], v[0][k]);
      const double hi = __dsub_ru(v[i][k], v[0][k]);
      pose.point[i][k] = fma(time, v[i][k] - v[0][k], x[i][k]);
      pose.lower[i][k] = __dsub_rd(__dadd_rd(x[i][k], __dmul_rd(time, lo)), error[i][k]);
      pose.upper[i][k] = __dadd_ru(__dadd_ru(x[i][k], __dmul_ru(time, hi)), error[i][k]);
    }
  }
  return pose;
}

// The axis is only a proposal. Its support intervals, rather than the
// floating-point closest-point parameters used to propose it, certify the gap.
__device__ inline double projectedSeparation(const Pose &pose, Vec3 axis, int split)
{
  const double scale = axis.cwiseAbs().maxCoeff();
  if (!(scale > 0.0) || !isfinite(scale))
    return 0.0;
  axis /= scale;
  double normSquared = 0.0;
  for (int k = 0; k < 3; ++k)
    normSquared = __dadd_ru(normSquared, __dmul_ru(axis[k], axis[k]));
  const double norm = __dsqrt_ru(normSquared);
  double minimum[2] = { CUDART_INF, CUDART_INF };
  double maximum[2] = { -CUDART_INF, -CUDART_INF };
  for (int i = 0; i < 4; ++i) {
    double lower = 0.0, upper = 0.0;
    for (int k = 0; k < 3; ++k) {
      // A common represented origin reduces cancellation without changing
      // either primitive's support separation.
      const double lo = __dsub_rd(pose.lower[i][k], pose.point[0][k]);
      const double hi = __dsub_ru(pose.upper[i][k], pose.point[0][k]);
      lower = __dadd_rd(lower, __dmul_rd(axis[k], axis[k] >= 0 ? lo : hi));
      upper = __dadd_ru(upper, __dmul_ru(axis[k], axis[k] >= 0 ? hi : lo));
    }
    if (!isfinite(lower) || !isfinite(upper))
      return 0.0;
    const int group = i < split ? 0 : 1;
    minimum[group] = fmin(minimum[group], lower);
    maximum[group] = fmax(maximum[group], upper);
  }
  const double gap = fmax(0.0, fmax(__dsub_rd(minimum[0], maximum[1]), __dsub_rd(minimum[1], maximum[0])));
  return __ddiv_rd(gap, norm);
}

__device__ inline Vec3 scaledDirection(Vec3 v)
{
  const double scale = v.cwiseAbs().maxCoeff();
  if (scale > 0.0 && isfinite(scale))
    v /= scale;
  return v;
}

__device__ inline Vec3 pointSegmentAxis(const Vec3 &p, const Vec3 &a, const Vec3 &b)
{
  const Vec3 edge = b - a;
  const double scale = edge.cwiseAbs().maxCoeff();
  if (!(scale > 0.0))
    return p - a;
  const Vec3 direction = edge / scale;
  double t = ((p - a) / scale).dot(direction) / direction.squaredNorm();
  if (!isfinite(t))
    return Vec3::Zero();
  t = fmin(1.0, fmax(0.0, t));
  return p - (a + t * edge);
}

__device__ inline double separationLowerBound(const Pose &pose, int kind)
{
  const Vec3 *p = pose.point;
  const int split = kind == 0 ? 1 : 2;
  Vec3 normal = kind == 0 ? scaledDirection(p[2] - p[1]).cross(scaledDirection(p[3] - p[1])).eval() :
                            scaledDirection(p[1] - p[0]).cross(scaledDirection(p[3] - p[2])).eval();
  double lower = projectedSeparation(pose, normal, split);
  if (kind == 0) {
    for (int j = 1; j < 4; ++j)
      lower = fmax(lower, projectedSeparation(pose, pointSegmentAxis(p[0], p[j], p[j == 3 ? 1 : j + 1]), split));
  }
  else {
    for (int j = 0; j < 4; ++j) {
      const int first = j < 2 ? 2 : 0;
      lower = fmax(lower, projectedSeparation(pose, pointSegmentAxis(p[j], p[first], p[first + 1]), split));
    }
  }
  return lower;
}

__device__ inline double speedUpperBound(const Vec3 (&motion)[4], int kind)
{
  const int split = kind == 0 ? 1 : 2;
  double speed[2] = {};
  for (int i = 0; i < 4; ++i) {
    double normSquared = 0.0;
    for (int k = 0; k < 3; ++k) {
      const double lo = __dsub_rd(motion[i][k], motion[0][k]);
      const double hi = __dsub_ru(motion[i][k], motion[0][k]);
      const double magnitude = fmax(fabs(lo), fabs(hi));
      normSquared = __dadd_ru(normSquared, __dmul_ru(magnitude, magnitude));
    }
    const int group = i < split ? 0 : 1;
    speed[group] = fmax(speed[group], __dsqrt_ru(normSquared));
  }
  return __dadd_ru(speed[0], speed[1]);
}

// Every vertex follows an affine trajectory in the original time parameter.
// A fixed-axis projection therefore attains its extrema at an interval's
// endpoints. Hull the signed supports before choosing a separation direction;
// hulling endpoint distances would incorrectly certify a through-plane motion.
__device__ inline double projectedIntervalSeparation(const Pose &first,
  const Pose &last, Vec3 axis, int split)
{
  const double scale = axis.cwiseAbs().maxCoeff();
  if (!(scale > 0.0) || !isfinite(scale))
    return 0.0;
  axis /= scale;
  double normSquared = 0.0;
  for (int k = 0; k < 3; ++k)
    normSquared = __dadd_ru(normSquared, __dmul_ru(axis[k], axis[k]));
  const double norm = __dsqrt_ru(normSquared);
  double minimum[2] = { CUDART_INF, CUDART_INF };
  double maximum[2] = { -CUDART_INF, -CUDART_INF };
  for (int endpoint = 0; endpoint < 2; ++endpoint) {
    const Pose &pose = endpoint == 0 ? first : last;
    for (int i = 0; i < 4; ++i) {
      double lower = 0.0, upper = 0.0;
      for (int k = 0; k < 3; ++k) {
        // Use one common origin for both endpoint support intervals.
        const double lo = __dsub_rd(pose.lower[i][k], first.point[0][k]);
        const double hi = __dsub_ru(pose.upper[i][k], first.point[0][k]);
        lower = __dadd_rd(lower, __dmul_rd(axis[k], axis[k] >= 0 ? lo : hi));
        upper = __dadd_ru(upper, __dmul_ru(axis[k], axis[k] >= 0 ? hi : lo));
      }
      if (!isfinite(lower) || !isfinite(upper))
        return 0.0;
      const int group = i < split ? 0 : 1;
      minimum[group] = fmin(minimum[group], lower);
      maximum[group] = fmax(maximum[group], upper);
    }
  }
  const double gap = fmax(0.0, fmax(__dsub_rd(minimum[0], maximum[1]), __dsub_rd(minimum[1], maximum[0])));
  return __ddiv_rd(gap, norm);
}

__device__ inline double intervalSeparationLowerBound(const Pose &first,
  const Pose &last, const Pose &midpoint, int kind)
{
  const Vec3 *p = midpoint.point;
  const int split = kind == 0 ? 1 : 2;
  const Vec3 normal = kind == 0 ? scaledDirection(p[2] - p[1]).cross(scaledDirection(p[3] - p[1])).eval() :
                                  scaledDirection(p[1] - p[0]).cross(scaledDirection(p[3] - p[2])).eval();
  double lower = projectedIntervalSeparation(first, last, normal, split);
  for (int j = kind == 0 ? 1 : 0; j < 4; ++j) {
    const int a = kind == 0 ? j : j < 2 ? 2 :
                                          0;
    const int b = kind == 0 ? j == 3 ? 1 : j + 1 : a + 1;
    lower = fmax(lower, projectedIntervalSeparation(first, last, pointSegmentAxis(p[kind == 0 ? 0 : j], p[a], p[b]), split));
  }
  return lower;
}

// Slow fallback for a clear, nearly grazing rotating primitive that exhausts
// additive advancement. Each closed interval is evaluated from the original
// inputs; rounded intermediate poses never become new exact path endpoints.
// Bounded work/depth or uncertain separation returns false, never a certificate.
__device__ inline bool certifyRemainingInterval(const Vec3 (&x)[4],
  const Vec3 (&motion)[4], const Vec3 (&error)[4], int kind, double firstTime)
{
  struct Interval
  {
    double first, last;
    int depth;
  };
  constexpr int maxDepth = 32, maxWork = 4096;
  Interval stack[maxDepth + 1];
  int size = 1;
  stack[0] = { firstTime, 1.0, 0 };
  for (int work = 0; work < maxWork && size > 0; ++work) {
    const Interval interval = stack[--size];
    const double middle = interval.first + 0.5 * (interval.last - interval.first);
    const Pose first = evaluate(x, motion, error, interval.first);
    const Pose last = evaluate(x, motion, error, interval.last);
    const Pose midpoint = evaluate(x, motion, error, middle);
    if (intervalSeparationLowerBound(first, last, midpoint, kind) > 0.0)
      continue;
    if (interval.depth == maxDepth || !(middle > interval.first && middle < interval.last) || size + 2 > maxDepth + 1)
      return false;
    // Sharing the same represented midpoint covers the complete parent.
    stack[size++] = { middle, interval.last, interval.depth + 1 };
    stack[size++] = { interval.first, middle, interval.depth + 1 };
  }
  return size == 0;
}

// Additive conservative advancement following the selected libuipc CCD design.
// A separating-axis gap is a lower distance bound regardless of feature
// classification; directed arithmetic includes projection and pose roundoff.
// A return value of one certifies the whole interval. Every other return is a
// conservative prefix, including numerical uncertainty or iteration exhaustion.
__device__ inline double maximumStep(const Vec3 (&x)[4], const Vec3 (&motion)[4],
  const Vec3 (&error)[4], int kind)
{
  const double speed = speedUpperBound(motion, kind);
  if (!isfinite(speed))
    return 0.0;
  const Pose endpoint = evaluate(x, motion, error, 1.0);
  double time = 0.0, retainedPrefix = 0.0, retainedGap = 0.0;
  for (int iteration = 0; iteration < 1000; ++iteration) {
    const Pose pose = evaluate(x, motion, error, time);
    const double gap = separationLowerBound(pose, kind);
    if (!(gap > 0.0) || !isfinite(gap))
      return retainedPrefix;
    if (iteration == 0)
      retainedGap = __dmul_ru(0.1, gap);
    // A fixed axis separating both endpoint hulls separates every linear
    // intermediate pose. This certifies long near-miss tangential motion
    // without taking thousands of steps limited by its small normal gap.
    Pose swept = pose;
    for (int i = 0; i < 4; ++i)
      for (int k = 0; k < 3; ++k) {
        swept.lower[i][k] = fmin(pose.lower[i][k], endpoint.lower[i][k]);
        swept.upper[i][k] = fmax(pose.upper[i][k], endpoint.upper[i][k]);
      }
    if (separationLowerBound(swept, kind) > 0.0)
      return 1.0;
    if (speed == 0.0)
      return 1.0;
    // Keep Newton away from the barrier singularity. Continue proving the
    // complete original chord, but on a partial result return a prefix whose
    // Lipschitz gap still retains 10% of the initial conservative separation.
    // Returning the last advancement pose instead can approach collision to
    // a few ULPs and create enormous barrier gradients when slackness is one.
    if (gap > retainedGap) {
      const double reserveAdvance = __dmul_rd(0.9, __ddiv_rd(__dsub_rd(gap, retainedGap), speed));
      retainedPrefix = fmax(retainedPrefix, __dadd_rd(time, reserveAdvance));
    }
    const double advance = __dmul_rd(0.9, __ddiv_rd(gap, speed));
    const double next = __dadd_rd(time, advance);
    if (next >= 1.0)
      return 1.0;
    if (!(next > time))
      return retainedPrefix;
    time = next;
  }
  return certifyRemainingInterval(x, motion, error, kind, time) ? 1.0 : retainedPrefix;
}
}  // namespace pgo::Contact::CIPC::CudaIPCContinuousCollision
