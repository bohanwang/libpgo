#include "triangleAreaValidation.h"

#include <gtest/gtest.h>

namespace
{
namespace Area = pgo::RunSim::TriangleAreaValidation;
Area::Positions triangle()
{
  return { Area::Point{ 0, 0, 0 }, Area::Point{ 1, 0, 0 }, Area::Point{ 0, 1, 0 } };
}

TEST(TriangleAreaValidation, CertifiesChangingNormalWithoutOrientationConstraint)
{
  const auto rest = triangle();
  auto quarter = rest, reversed = rest;
  quarter[2] = { 0, 0, 1 };
  reversed[2] = { 0, -1, 0 };
  EXPECT_TRUE(Area::certify(reversed).certified);
  const auto first = Area::certify(rest, &quarter);
  const auto second = Area::certify(quarter, &reversed);
  ASSERT_TRUE(first.certified);
  ASSERT_TRUE(second.certified);
  EXPECT_GT(first.lower, 0);
  EXPECT_LE(first.lower, 0.5);
  EXPECT_LE(second.lower, 0.5);
  EXPECT_GE(first.upper, 1);
  EXPECT_FALSE(Area::certify(rest, &reversed).certified);
}

TEST(TriangleAreaValidation, RejectsQuarticInteriorZerosDespitePositiveEndpointNormals)
{
  const auto rest = triangle();
  auto to = rest;
  to[1][0] = -1.2;
  to[2][1] = -2.4;
  ASSERT_TRUE(Area::certify(rest).certified);
  ASSERT_TRUE(Area::certify(to).certified);
  // Both endpoint normals point +Z. Distinct interior zeros occur at
  // t=1/2.2 and t=1/3.4, neither an endpoint or an equally spaced sample.
  const auto swept = Area::certify(rest, &to);
  EXPECT_FALSE(swept.certified);
  EXPECT_EQ(swept.lower, 0);
}

TEST(TriangleAreaValidation, SubdivisionCertifiesPositiveNearCollapse)
{
  auto from = triangle(), to = from;
  from[2][2] = to[2][2] = 1e-5;
  to[2][1] = -1;
  const auto result = Area::certify(from, &to);
  ASSERT_TRUE(result.certified);
  EXPECT_GT(result.subdivisions, 0);
  EXPECT_GT(result.lower, 0);
  EXPECT_LE(result.lower, 1e-10);
}

TEST(TriangleAreaValidation, SharedPowerOfTwoScaleRetainsTinyAndLargeAreaRatios)
{
  for (double scale : { 1e-80, 1.0, 1e80 }) {
    auto from = triangle();
    for (auto &vertex : from)
      for (double &value : vertex)
        value *= scale;
    auto to = from;
    to[2] = { 0, 0, scale };
    const int power = Area::edgeScalePower(from);
    const auto rest = Area::certify(from, nullptr, power);
    const auto swept = Area::certify(from, &to, power);
    ASSERT_TRUE(rest.certified);
    ASSERT_TRUE(swept.certified);
    const double ratio = Area::down(swept.lower / rest.upper);
    EXPECT_GT(ratio, 0);
    EXPECT_LE(ratio, 0.5);
  }
}

TEST(TriangleAreaValidation, RepresentedLargeOriginRetainsEdgeDifferences)
{
  auto from = triangle();
  for (auto &vertex : from)
    for (double &value : vertex)
      value += std::ldexp(1.0, 40);
  auto to = from;
  to[2][1] -= 1;
  to[2][2] += 1;
  EXPECT_TRUE(Area::certify(from, &to, Area::edgeScalePower(from)).certified);
}

TEST(TriangleAreaValidation, RejectsCollapsedNonfiniteAndOverflowingGeometry)
{
  auto positions = triangle();
  positions[2] = positions[0];
  EXPECT_FALSE(Area::certify(positions).certified);
  for (double bad : { std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity() }) {
    positions = triangle();
    positions[2][0] = bad;
    EXPECT_FALSE(Area::certify(positions).certified);
  }
  positions = triangle();
  positions[0][0] = -std::numeric_limits<double>::max();
  positions[1][0] = std::numeric_limits<double>::max();
  EXPECT_FALSE(Area::certify(positions).certified);
}
}  // namespace
