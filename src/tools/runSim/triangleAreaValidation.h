#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

// Independent CPU certificate for a triangle's nonzero area. Endpoints are
// represented world coordinates; the path is their exact linear interpolation.
// The cross product has Bernstein degree2 in time and its squared norm degree4.
// Every arithmetic bound is rounded outwards. An inconclusive bound fails.
// No normal orientation, constitutive model, or GPU routine is involved.
namespace pgo::RunSim::TriangleAreaValidation
{
using Point = std::array<double, 3>;
using Positions = std::array<Point, 3>;

struct Interval
{
  double lo = 0, hi = 0;
  Interval() = default;
  Interval(double value): lo(value), hi(value) {}
  Interval(double lower, double upper): lo(lower), hi(upper) {}
};

inline double down(double value)
{
  return std::nextafter(value, -std::numeric_limits<double>::infinity());
}

inline double up(double value)
{
  return std::nextafter(value, std::numeric_limits<double>::infinity());
}

inline Interval uncertain()
{
  return { -std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity() };
}

inline bool finite(Interval value)
{
  return std::isfinite(value.lo) && std::isfinite(value.hi) && value.lo <= value.hi;
}

inline Interval operator+(Interval a, Interval b)
{
  if (!finite(a) || !finite(b))
    return uncertain();
  return { down(a.lo + b.lo), up(a.hi + b.hi) };
}

inline Interval operator-(Interval a, Interval b)
{
  if (!finite(a) || !finite(b))
    return uncertain();
  return { down(a.lo - b.hi), up(a.hi - b.lo) };
}

inline Interval operator*(Interval a, Interval b)
{
  if (!finite(a) || !finite(b))
    return uncertain();
  const std::array<double, 4> products{ a.lo * b.lo, a.lo * b.hi, a.hi * b.lo, a.hi * b.hi };
  for (double value : products)
    if (!std::isfinite(value))
      return uncertain();
  return { down(*std::min_element(products.begin(), products.end())), up(*std::max_element(products.begin(), products.end())) };
}

inline Interval square(Interval a)
{
  if (!finite(a))
    return uncertain();
  const double low = a.lo <= 0 && a.hi >= 0 ? 0 : std::min(a.lo * a.lo, a.hi * a.hi);
  const double high = std::max(a.lo * a.lo, a.hi * a.hi);
  if (!std::isfinite(low) || !std::isfinite(high))
    return uncertain();
  return { std::max(0.0, down(low)), up(high) };
}

inline Interval fraction(int numerator, int denominator)
{
  const double value = double(numerator) / denominator;
  return { down(value), up(value) };
}

// A shared power-of-two edge scale prevents tiny rest triangles from
// underflowing when their area is squared. It cancels in the rest-area ratio.
inline int edgeScalePower(const Positions &rest)
{
  long double maximum = 0;
  for (int v = 1; v < 3; ++v)
    for (int c = 0; c < 3; ++c)
      maximum = std::max(maximum, std::abs(static_cast<long double>(rest[v][c]) - rest[0][c]));
  if (!(maximum > 0) || !std::isfinite(maximum))
    return 0;
  int exponent;
  std::frexp(maximum, &exponent);
  return -exponent;
}

using Vector = std::array<Interval, 3>;
inline Vector edge(const Positions &positions, int vertex, int power)
{
  Vector result;
  for (int c = 0; c < 3; ++c) {
    const auto value = Interval(positions[vertex][c]) - Interval(positions[0][c]);
    result[c] = finite(value) ? Interval(down(std::scalbn(value.lo, power)), up(std::scalbn(value.hi, power))) : uncertain();
  }
  return result;
}

inline Vector cross(const Vector &a, const Vector &b)
{
  Vector result;
  for (int c = 0; c < 3; ++c)
    result[c] = a[(c + 1) % 3] * b[(c + 2) % 3] - a[(c + 2) % 3] * b[(c + 1) % 3];
  return result;
}

inline Interval dot(const Vector &a, const Vector &b)
{
  Interval result;
  for (int c = 0; c < 3; ++c)
    result = result + a[c] * b[c];
  return result;
}

struct Certificate
{
  bool certified = false;
  double lower = 0, upper = std::numeric_limits<double>::infinity();
  int subdivisions = 0;
};
using Coefficients = std::array<Interval, 5>;

inline Certificate certifyCoefficients(const Coefficients &coefficients, int depth = 0)
{
  Certificate result;
  double minimum = std::numeric_limits<double>::infinity(), maximum = -minimum;
  for (auto value : coefficients) {
    if (!finite(value))
      return result;
    minimum = std::min(minimum, value.lo);
    maximum = std::max(maximum, value.hi);
  }
  result.lower = std::max(0.0, minimum);
  result.upper = maximum;
  if (minimum > 0) {
    result.certified = true;
    return result;
  }
  if (maximum <= 0 || depth == 32)
    return result;
  // de Casteljau subdivision carries the original coefficient intervals;
  // no rounded midpoint is ever promoted to a new exact trajectory endpoint.
  Coefficients left, right, work = coefficients;
  left[0] = work[0];
  right[4] = work[4];
  for (int level = 1; level <= 4; ++level) {
    for (int i = 0; i <= 4 - level; ++i)
      work[i] = Interval(0.5) * (work[i] + work[i + 1]);
    left[level] = work[0];
    right[4 - level] = work[4 - level];
  }
  const auto a = certifyCoefficients(left, depth + 1);
  if (!a.certified)
    return { false, 0, maximum, 1 + a.subdivisions };
  const auto b = certifyCoefficients(right, depth + 1);
  return { a.certified && b.certified, std::min(a.lower, b.lower), std::max(a.upper, b.upper), 1 + a.subdivisions + b.subdivisions };
}

inline Certificate certify(const Positions &from, const Positions *to = nullptr, int edgePower = 0)
{
  for (int endpoint = 0; endpoint < (to ? 2 : 1); ++endpoint)
    for (const auto &vertex : endpoint == 0 ? from : *to)
      for (double value : vertex)
        if (!std::isfinite(value))
          return {};
  const auto a0 = edge(from, 1, edgePower), b0 = edge(from, 2, edgePower);
  const auto n0 = cross(a0, b0);
  if (!to) {
    Interval value;
    for (auto component : n0)
      value = value + square(component);
    return finite(value) ? Certificate{ value.lo > 0, std::max(0.0, value.lo), value.hi, 0 } : Certificate{};
  }
  const auto a1 = edge(*to, 1, edgePower), b1 = edge(*to, 2, edgePower);
  const auto mixed0 = cross(a0, b1), mixed1 = cross(a1, b0);
  std::array<Vector, 3> normals{ n0, {}, cross(a1, b1) };
  for (int c = 0; c < 3; ++c)
    normals[1][c] = Interval(0.5) * (mixed0[c] + mixed1[c]);
  constexpr int choose2[3] = { 1, 2, 1 }, choose4[5] = { 1, 4, 6, 4, 1 };
  Coefficients coefficients{};
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      coefficients[i + j] = coefficients[i + j] + fraction(choose2[i] * choose2[j], choose4[i + j]) * dot(normals[i], normals[j]);
  return certifyCoefficients(coefficients);
}
}  // namespace pgo::RunSim::TriangleAreaValidation
