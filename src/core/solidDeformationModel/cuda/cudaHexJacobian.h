#pragma once

#include <cuda_runtime.h>
#include <math_constants.h>
#include <cfloat>

// Whole-cell feasibility for trilinear hexes. This is deliberately independent
// of the CPU trajectory checker's Bernstein implementation. A Jacobian column
// is multiaffine in the other two material coordinates and linear time, so its
// range on a parameter box lies in the hull of its eight corner values.
// Directed arithmetic encloses those values and the determinant. Subdivision
// retains the original displacement chord; no rounded midpoint is a new pose.
namespace pgo::SolidDeformationModel::HexJacobian
{
struct Interval
{
  double lo, hi;
  __device__ Interval(double value = 0):
    lo(value), hi(value) {}
  __device__ Interval(double lower, double upper):
    lo(lower), hi(upper) {}
};
__device__ inline Interval operator+(Interval a, Interval b)
{
  return { __dadd_rd(a.lo, b.lo), __dadd_ru(a.hi, b.hi) };
}
__device__ inline Interval operator-(Interval a, Interval b)
{
  return { __dsub_rd(a.lo, b.hi), __dsub_ru(a.hi, b.lo) };
}
__device__ inline Interval operator*(Interval a, Interval b)
{
  if (!isfinite(a.lo) || !isfinite(a.hi) || !isfinite(b.lo) || !isfinite(b.hi))
    return { -CUDART_INF, CUDART_INF };
  const double lower[4]{ __dmul_rd(a.lo, b.lo), __dmul_rd(a.lo, b.hi), __dmul_rd(a.hi, b.lo), __dmul_rd(a.hi, b.hi) };
  const double upper[4]{ __dmul_ru(a.lo, b.lo), __dmul_ru(a.lo, b.hi), __dmul_ru(a.hi, b.lo), __dmul_ru(a.hi, b.hi) };
  for (int i = 0; i < 4; ++i)
    if (!isfinite(lower[i]) || !isfinite(upper[i]))
      return { -CUDART_INF, CUDART_INF };
  return { fmin(fmin(lower[0], lower[1]), fmin(lower[2], lower[3])),
    fmax(fmax(upper[0], upper[1]), fmax(upper[2], upper[3])) };
}
__device__ inline bool finite(Interval a)
{
  return isfinite(a.lo) && isfinite(a.hi) && a.lo <= a.hi;
}
__device__ inline Interval determinant(const Interval m[3][3])
{
  return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
    m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
    m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
}

struct Element
{
  int vertices[8];
  double rest[8][3];
  double orientation;
};
struct Motion
{
  Interval from[8][3], delta[8][3];
};
struct Box
{
  double lo[4]{ 0, 0, 0, 0 }, hi[4]{ 1, 1, 1, 1 };
  int depth = 0;
};
__device__ inline int vertex(int bits)
{
  // SimulationMesh's perimeter ordering differs from binary xyz ordering.
  const int ids[8]{ 0, 1, 3, 2, 4, 5, 7, 6 };
  return ids[bits];
}

__device__ inline bool loadMotion(const Element &element, const double *x, const double *dx, Motion &motion)
{
  for (int v = 0; v < 8; ++v) {
    for (int c = 0; c < 3; ++c) {
      const double rest = element.rest[v][c];
      const double u = x ? x[3 * element.vertices[v] + c] : 0;
      const double d = dx ? dx[3 * element.vertices[v] + c] : 0;
      if (!isfinite(rest) || !isfinite(u) || !isfinite(d))
        return false;
      motion.from[v][c] = Interval(rest) + Interval(u);
      // certifyChord supplies fl(trial-base). Adjacent doubles enclose the
      // exact difference, including subnormal values. Zero differences between
      // represented inputs are exact. Never rebase this interval at a split.
      motion.delta[v][c] = d == 0 ? Interval(0) : Interval(nextafter(d, -CUDART_INF), nextafter(d, CUDART_INF));
      if (dx) {
        // The original rest/displacement scale also bounds reconstruction of
        // rounded x+t*dx and rest+(x+t*dx), and the exported endpoint chord.
        // Four DBL_EPSILON exceeds the gamma_3 bound for these operations;
        // retain an outward subnormal allowance as well. No world-unit floor.
        const double magnitude = __dadd_ru(fabs(rest), __dadd_ru(fabs(u), fabs(d)));
        const double radius = __dadd_ru(__dmul_ru(4 * DBL_EPSILON, magnitude), 4 * 0x0.0000000000001p-1022);
        motion.from[v][c] = motion.from[v][c] + Interval(-radius, radius);
      }
      if (!finite(motion.from[v][c]) || !finite(motion.delta[v][c]))
        return false;
    }
  }
  return true;
}

__device__ inline Interval edge(const Element &element, const Motion &motion,
  int bits, int axis, int row, double time, bool reference)
{
  const int a = vertex(bits), b = vertex(bits | (1 << axis));
  if (reference)
    return Interval(element.rest[b][row]) - Interval(element.rest[a][row]);
  return (motion.from[b][row] - motion.from[a][row]) +
    Interval(time) * (motion.delta[b][row] - motion.delta[a][row]);
}

__device__ inline bool jacobianRange(const Element &element, const Motion &motion,
  const Box &box, bool reference, Interval matrix[3][3], double variation[4])
{
  for (int axis = 0; axis < 3; ++axis) {
    const int other[2]{ (axis + 1) % 3, (axis + 2) % 3 };
    for (int row = 0; row < 3; ++row) {
      Interval values[8];
      for (int corner = 0; corner < 8; ++corner) {
        const Interval u((corner & 1) ? box.hi[other[0]] : box.lo[other[0]]);
        const Interval v((corner & 2) ? box.hi[other[1]] : box.lo[other[1]]);
        const double time = (corner & 4) ? box.hi[3] : box.lo[3];
        Interval value;
        for (int bits = 0; bits < 4; ++bits) {
          const auto weight = ((bits & 1) ? u : Interval(1) - u) * ((bits & 2) ? v : Interval(1) - v);
          const int location = ((bits & 1) ? 1 << other[0] : 0) | ((bits & 2) ? 1 << other[1] : 0);
          value = value + weight * edge(element, motion, location, axis, row, time, reference);
        }
        if (!finite(value))
          return false;
        values[corner] = value;
        if (corner == 0)
          matrix[row][axis] = value;
        else {
          matrix[row][axis].lo = fmin(matrix[row][axis].lo, value.lo);
          matrix[row][axis].hi = fmax(matrix[row][axis].hi, value.hi);
        }
      }
      if (!reference) {
        for (int direction = 0; direction < 3; ++direction) {
          const int parameter = direction == 2 ? 3 : other[direction];
          for (int corner = 0; corner < 8; ++corner) {
            if (corner & (1 << direction))
              continue;
            const auto difference = values[corner | (1 << direction)] - values[corner];
            variation[parameter] = fmax(variation[parameter], fmax(fabs(difference.lo), fabs(difference.hi)));
          }
        }
      }
    }
  }
  return true;
}

__device__ inline bool certify(const Element &element, const Motion &motion, double limit, double determinantFloor)
{
  // Samples only reject. They cannot certify a whole cell, but cheaply avoid
  // subdividing many positive boxes before finding a known inverted corner.
  for (int time = 0; time < 3; ++time) {
    const double t = time == 0 ? 0 : time == 1 ? 0.5 * limit :
                                                 limit;
    for (int corner = 0; corner < 8; ++corner) {
      Interval current[3][3], reference[3][3];
      for (int axis = 0; axis < 3; ++axis)
        for (int row = 0; row < 3; ++row) {
          // At a material corner each derivative is its incident parallel edge.
          const int bits = corner & ~(1 << axis);
          current[row][axis] = edge(element, motion, bits, axis, row, t, false);
          reference[row][axis] = edge(element, motion, bits, axis, row, t, true);
        }
      const auto bound = determinant(current) * Interval(element.orientation) -
        Interval(determinantFloor) * determinant(reference) * Interval(element.orientation);
      if (!finite(bound) || bound.hi <= 0)
        return false;
    }
  }
  // A warped but positive corner can require refinement in all four
  // parameters. Real dragon chords with det(F) >= 2e-8 require depth94;
  // a shallower limit imposed an artificial boundary well above the 1e-8
  // material floor. The independent work bound still fails closed.
  constexpr int maxDepth = 96, maxWork = 4096;
  Box stack[maxDepth + 1];
  stack[0].hi[3] = limit;
  int pending = 1;
  for (int work = 0; work < maxWork && pending; ++work) {
    const Box box = stack[--pending];
    Interval current[3][3], reference[3][3];
    double variation[4]{};
    if (!jacobianRange(element, motion, box, false, current, variation) ||
      !jacobianRange(element, motion, box, true, reference, variation))
      return false;
    const auto restDet = determinant(reference) * Interval(element.orientation);
    const auto bound = determinant(current) * Interval(element.orientation) - Interval(determinantFloor) * restDet;
    if (!finite(bound) || !finite(restDet))
      return false;
    if (bound.lo > 0)
      continue;
    if (bound.hi <= 0 || box.depth == maxDepth)
      return false;
    int axis = -1;
    double maximum = -1;
    for (int candidate = 0; candidate < 4; ++candidate) {
      if (box.lo[candidate] == box.hi[candidate])
        continue;
      if (variation[candidate] > maximum) {
        maximum = variation[candidate];
        axis = candidate;
      }
    }
    if (axis < 0 || pending + 2 > maxDepth + 1)
      return false;
    const double middle = box.lo[axis] + 0.5 * (box.hi[axis] - box.lo[axis]);
    if (!(middle > box.lo[axis] && middle < box.hi[axis]))
      return false;
    Box left = box, right = box;
    left.hi[axis] = middle;
    right.lo[axis] = middle;
    left.depth = right.depth = box.depth + 1;
    stack[pending++] = right;
    stack[pending++] = left;
  }
  return pending == 0;
}

// Most cells are slightly strained. Every column is then certified directly
// from its four parallel edges at the original chord's two time endpoints.
// This does not need local subdivision storage in the common kernel path.
__device__ inline bool certifyEdgeHull(const Element &element, const Motion &motion, double limit, double determinantFloor)
{
  Interval current[3][3], reference[3][3];
  for (int axis = 0; axis < 3; ++axis) {
    for (int row = 0; row < 3; ++row) {
      current[row][axis] = reference[row][axis] = Interval(CUDART_INF, -CUDART_INF);
      for (int bits = 0; bits < 8; ++bits) {
        if (bits & (1 << axis))
          continue;
        const auto rest = edge(element, motion, bits, axis, row, 0, true);
        if (!finite(rest))
          return false;
        reference[row][axis].lo = fmin(reference[row][axis].lo, rest.lo);
        reference[row][axis].hi = fmax(reference[row][axis].hi, rest.hi);
        for (int endpoint = 0; endpoint < 2; ++endpoint) {
          const auto value = edge(element, motion, bits, axis, row, endpoint ? limit : 0, false);
          if (!finite(value))
            return false;
          current[row][axis].lo = fmin(current[row][axis].lo, value.lo);
          current[row][axis].hi = fmax(current[row][axis].hi, value.hi);
        }
      }
    }
  }
  const auto restDet = determinant(reference) * Interval(element.orientation);
  const auto bound = determinant(current) * Interval(element.orientation) - Interval(determinantFloor) * restDet;
  return finite(bound) && bound.lo > 0;
}
}  // namespace pgo::SolidDeformationModel::HexJacobian
