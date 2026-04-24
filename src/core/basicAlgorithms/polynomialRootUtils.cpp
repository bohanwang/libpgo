

#include "polynomialRootUtils.h"

#include <algorithm>
#include <cmath>

namespace pgo::BasicAlgorithms
{
namespace
{
constexpr double kBoundaryEvalTol = 1e-12;
constexpr int kBoundaryBisectionIterations = 64;

void insertCriticalPointIfInUnitInterval(std::vector<double> &roots, double value)
{
  if (value > 0.0 && value < 1.0) {
    roots.push_back(value);
  }
}
}

double CubicPolynomial::eval(double alpha) const
{
  return ((c3 * alpha + c2) * alpha + c1) * alpha + c0;
}

std::vector<double> findCriticalPointsInUnitInterval(const CubicPolynomial &poly)
{
  const double qa = 3.0 * poly.c3;
  const double qb = 2.0 * poly.c2;
  const double qc = poly.c1;

  std::vector<double> roots;
  if (qa == 0.0) {
    if (qb == 0.0) {
      return roots;
    }

    insertCriticalPointIfInUnitInterval(roots, -qc / qb);
    return roots;
  }

  const double disc = qb * qb - 4.0 * qa * qc;
  if (disc < 0.0) {
    return roots;
  }

  const double sqrtDisc = std::sqrt(std::max(0.0, disc));
  const double denom = 2.0 * qa;
  insertCriticalPointIfInUnitInterval(roots, (-qb - sqrtDisc) / denom);
  insertCriticalPointIfInUnitInterval(roots, (-qb + sqrtDisc) / denom);

  std::sort(roots.begin(), roots.end());
  roots.erase(std::unique(roots.begin(), roots.end(), [](double lhs, double rhs) {
    return std::abs(lhs - rhs) <= kBoundaryEvalTol;
  }), roots.end());
  return roots;
}

double findFirstBoundaryRootByBisection(const CubicPolynomial &poly, double left, double right)
{
  double feasible = left;
  double infeasible = right;

  if (std::abs(poly.eval(infeasible)) <= kBoundaryEvalTol) {
    return infeasible;
  }

  for (int iter = 0; iter < kBoundaryBisectionIterations; iter++) {
    const double mid = 0.5 * (feasible + infeasible);
    const double gmid = poly.eval(mid);
    if (gmid > 0.0) {
      feasible = mid;
    }
    else {
      infeasible = mid;
    }
  }

  return feasible;
}
}
