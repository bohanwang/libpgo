

#pragma once

#include <vector>

namespace pgo::BasicAlgorithms
{
struct CubicPolynomial
{
  double c0 = 0.0;
  double c1 = 0.0;
  double c2 = 0.0;
  double c3 = 0.0;

  double eval(double alpha) const;
};

std::vector<double> findCriticalPointsInUnitInterval(const CubicPolynomial &poly);
double findFirstBoundaryRootByBisection(const CubicPolynomial &poly, double left, double right);
}
