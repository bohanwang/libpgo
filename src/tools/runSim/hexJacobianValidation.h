#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

// Independent CPU certificate for a trilinear hexahedron. Its Jacobian
// determinant has tensor degree (2,2,2), and degree 3 along linear time motion.
// All bounds are rounded outward. An inconclusive certificate is a failure,
// never evidence that a cell is valid. No CUDA/FEM implementation is reused.
namespace pgo::RunSim::HexValidation
{
struct Interval
{
  double lo = 0, hi = 0;
  Interval() = default;
  Interval(double value):
    lo(value), hi(value) {}
  Interval(double lower, double upper):
    lo(lower), hi(upper) {}
};
inline double down(double value)
{
  return std::nextafter(value, -std::numeric_limits<double>::infinity());
}
inline double up(double value)
{
  return std::nextafter(value, std::numeric_limits<double>::infinity());
}
inline Interval operator+(Interval a, Interval b)
{
  return { down(a.lo + b.lo), up(a.hi + b.hi) };
}
inline Interval operator-(Interval a, Interval b)
{
  return { down(a.lo - b.hi), up(a.hi - b.lo) };
}
inline Interval operator*(Interval a, Interval b)
{
  const std::array<double, 4> products{ a.lo * b.lo, a.lo * b.hi, a.hi * b.lo, a.hi * b.hi };
  return { down(*std::min_element(products.begin(), products.end())), up(*std::max_element(products.begin(), products.end())) };
}
inline Interval fraction(int numerator, int denominator)
{
  const double value = double(numerator) / denominator;
  return { down(value), up(value) };
}
using Point = std::array<Interval, 3>;
using Positions = std::array<Point, 8>;
using Matrix = std::array<std::array<Interval, 3>, 3>;
constexpr std::array<int, 8> bitToVertex{ 0, 1, 3, 2, 4, 5, 7, 6 };
inline Interval determinant(const Matrix &m)
{
  return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
    m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
    m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
}
inline Interval jacobian(const Positions &positions, std::array<Interval, 3> parameter)
{
  Matrix matrix{};
  for (int axis = 0; axis < 3; ++axis) {
    for (int bits = 0; bits < 8; ++bits) {
      if (bits & (1 << axis))
        continue;
      Interval weight(1);
      for (int other = 0; other < 3; ++other)
        if (other != axis)
          weight = weight * ((bits & (1 << other)) ? parameter[other] : Interval(1) - parameter[other]);
      for (int row = 0; row < 3; ++row)
        matrix[row][axis] = matrix[row][axis] + weight * (positions[bitToVertex[bits | (1 << axis)]][row] - positions[bitToVertex[bits]][row]);
    }
  }
  return determinant(matrix);
}

struct Certificate
{
  bool certified = false;
  bool nonpositiveSample = false;
  double lower = -std::numeric_limits<double>::infinity();
  double upper = std::numeric_limits<double>::infinity();
  int subdivisions = 0;
};
using Coefficients = std::array<Interval, 108>;
constexpr std::array<int, 4> strides{ 1, 3, 9, 27 };

inline Certificate certifyCoefficients(const Coefficients &coefficients, int times, int depth = 0)
{
  Certificate result;
  result.lower = std::numeric_limits<double>::infinity();
  result.upper = -result.lower;
  const int size = 27 * times;
  for (int i = 0; i < size; ++i) {
    result.lower = std::min(result.lower, coefficients[i].lo);
    result.upper = std::max(result.upper, coefficients[i].hi);
  }
  if (std::isfinite(result.lower) && std::isfinite(result.upper) && result.lower > 0) {
    result.certified = true;
    return result;
  }
  if (!std::isfinite(result.lower) || !std::isfinite(result.upper) || result.upper <= 0 || depth == 12)
    return result;
  // Split the axis with the greatest variation in control coefficients.
  int axis = 0;
  double variation = -1;
  for (int candidate = 0; candidate < (times == 1 ? 3 : 4); ++candidate) {
    const int count = candidate == 3 ? times : 3, stride = strides[candidate];
    double maximum = 0;
    for (int i = 0; i < size; ++i)
      if ((i / stride) % count < count - 1)
        maximum = std::max(maximum, std::abs(coefficients[i + stride].lo - coefficients[i].lo));
    if (maximum > variation) {
      variation = maximum;
      axis = candidate;
    }
  }
  const int count = axis == 3 ? times : 3, stride = strides[axis];
  Coefficients left = coefficients, right = coefficients;
  for (int i = 0; i < size; ++i) {
    if ((i / stride) % count != 0)
      continue;
    std::array<Interval, 4> work{};
    for (int j = 0; j < count; ++j)
      work[j] = coefficients[i + j * stride];
    left[i] = work[0];
    right[i + (count - 1) * stride] = work[count - 1];
    for (int level = 1; level < count; ++level) {
      for (int j = 0; j < count - level; ++j)
        work[j] = Interval(0.5) * (work[j] + work[j + 1]);
      left[i + level * stride] = work[0];
      right[i + (count - 1 - level) * stride] = work[count - 1 - level];
    }
  }
  const auto a = certifyCoefficients(left, times, depth + 1);
  if (!a.certified)
    return a;
  const auto b = certifyCoefficients(right, times, depth + 1);
  return { a.certified && b.certified, false, std::min(a.lower, b.lower), std::max(a.upper, b.upper),
    1 + a.subdivisions + b.subdivisions };
}

inline Certificate certify(const Positions &from, const Positions *to = nullptr, double orientation = 1)
{
  if (orientation != 1 && orientation != -1)
    return {};
  for (int endpoint = 0; endpoint < (to ? 2 : 1); ++endpoint)
    for (const auto &vertex : endpoint == 0 ? from : *to)
      for (const auto value : vertex)
        if (!std::isfinite(value.lo) || !std::isfinite(value.hi) || value.lo > value.hi)
          return {};
  // Each Jacobian entry is a convex combination of the four parallel
  // hexahedron edges, also between the two time endpoints. This cheap interval
  // bound certifies most slightly strained cells without polynomial conversion.
  Matrix bounds;
  for (int axis = 0; axis < 3; ++axis) {
    for (int row = 0; row < 3; ++row) {
      bounds[row][axis] = { std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity() };
      for (int endpoint = 0; endpoint < (to ? 2 : 1); ++endpoint) {
        const auto &positions = endpoint == 0 ? from : *to;
        for (int bits = 0; bits < 8; ++bits) {
          if (bits & (1 << axis))
            continue;
          const auto edge = positions[bitToVertex[bits | (1 << axis)]][row] - positions[bitToVertex[bits]][row];
          if (!std::isfinite(edge.lo) || !std::isfinite(edge.hi))
            return {};
          bounds[row][axis].lo = std::min(bounds[row][axis].lo, edge.lo);
          bounds[row][axis].hi = std::max(bounds[row][axis].hi, edge.hi);
        }
      }
    }
  }
  const auto simple = determinant(bounds) * Interval(orientation);
  if (std::isfinite(simple.lo) && std::isfinite(simple.hi) && simple.lo > 0)
    return { true, false, simple.lo, simple.hi, 0 };

  const int times = to ? 4 : 1, size = 27 * times;
  Coefficients coefficients{};
  for (int time = 0; time < times; ++time) {
    Positions positions = from;
    if (to) {
      const auto t = time == 0 ? Interval(0) : time == 3 ? Interval(1) :
                                                           fraction(time, 3);
      for (int v = 0; v < 8; ++v)
        for (int row = 0; row < 3; ++row)
          positions[v][row] = from[v][row] * (Interval(1) - t) + (*to)[v][row] * t;
    }
    for (int z = 0; z < 3; ++z)
      for (int y = 0; y < 3; ++y)
        for (int x = 0; x < 3; ++x) {
          const auto value = jacobian(positions, { Interval(0.5 * x), Interval(0.5 * y), Interval(0.5 * z) }) * Interval(orientation);
          if (value.hi <= 0 || !std::isfinite(value.lo) || !std::isfinite(value.hi))
            return { false, value.hi <= 0, value.lo, value.hi, 0 };
          coefficients[x + 3 * y + 9 * z + 27 * time] = value;
        }
  }
  // Convert tensor nodal values at equally spaced points to Bernstein form.
  for (int axis = 0; axis < (to ? 4 : 3); ++axis) {
    const int count = axis == 3 ? times : 3, stride = strides[axis];
    for (int i = 0; i < size; ++i) {
      if ((i / stride) % count != 0)
        continue;
      if (count == 3)
        coefficients[i + stride] = Interval(2) * coefficients[i + stride] -
          Interval(0.5) * (coefficients[i] + coefficients[i + 2 * stride]);
      else {
        const auto a = Interval(27) * coefficients[i + stride] - Interval(8) * coefficients[i] - coefficients[i + 3 * stride];
        const auto b = Interval(27) * coefficients[i + 2 * stride] - coefficients[i] - Interval(8) * coefficients[i + 3 * stride];
        coefficients[i + stride] = (Interval(2) * a - b) * fraction(1, 18);
        coefficients[i + 2 * stride] = (Interval(2) * b - a) * fraction(1, 18);
      }
    }
  }
  return certifyCoefficients(coefficients, times);
}
}  // namespace pgo::RunSim::HexValidation
