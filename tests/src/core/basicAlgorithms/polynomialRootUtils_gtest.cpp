#include <gtest/gtest.h>

#include "polynomialRootUtils.h"

#include <cmath>

namespace
{
using pgo::BasicAlgorithms::CubicPolynomial;
using pgo::BasicAlgorithms::findCriticalPointsInUnitInterval;
using pgo::BasicAlgorithms::findFirstBoundaryRootByBisection;
}

TEST(PolynomialRootUtilsGTest, DegreeZeroPolynomialHasNoCriticalPoints)
{
  const CubicPolynomial poly{ 3.0, 0.0, 0.0, 0.0 };
  EXPECT_TRUE(findCriticalPointsInUnitInterval(poly).empty());
  EXPECT_DOUBLE_EQ(poly.eval(0.25), 3.0);
}

TEST(PolynomialRootUtilsGTest, DegreeOnePolynomialHasNoCriticalPointsAndBisectsRoot)
{
  const CubicPolynomial poly{ 1.0, -2.0, 0.0, 0.0 };
  EXPECT_TRUE(findCriticalPointsInUnitInterval(poly).empty());
  EXPECT_NEAR(findFirstBoundaryRootByBisection(poly, 0.0, 1.0), 0.5, 1e-12);
}

TEST(PolynomialRootUtilsGTest, DegreeTwoPolynomialReturnsSingleCriticalPointInUnitInterval)
{
  const CubicPolynomial poly{ 0.0, -1.0, 1.0, 0.0 };
  const std::vector<double> roots = findCriticalPointsInUnitInterval(poly);
  ASSERT_EQ(roots.size(), 1u);
  EXPECT_NEAR(roots[0], 0.5, 1e-12);
}

TEST(PolynomialRootUtilsGTest, DegreeThreePolynomialReturnsTwoCriticalPoints)
{
  const CubicPolynomial poly{ 1.0, 0.5, -1.5, 1.0 };
  const std::vector<double> roots = findCriticalPointsInUnitInterval(poly);
  ASSERT_EQ(roots.size(), 2u);
  EXPECT_NEAR(roots[0], 0.2113248654051871, 1e-12);
  EXPECT_NEAR(roots[1], 0.7886751345948129, 1e-12);
}

TEST(PolynomialRootUtilsGTest, EndpointRootReturnsRightEndpointWithoutIterationDrift)
{
  const CubicPolynomial poly{ 1.0, -1.0, 0.0, 0.0 };
  EXPECT_DOUBLE_EQ(findFirstBoundaryRootByBisection(poly, 0.0, 1.0), 1.0);
}
