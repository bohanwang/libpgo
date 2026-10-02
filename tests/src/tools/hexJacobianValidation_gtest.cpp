#include "hexJacobianValidation.h"

#include <gtest/gtest.h>

namespace
{
namespace Hex = pgo::RunSim::HexValidation;
Hex::Positions cube()
{
  Hex::Positions positions{};
  for (int bits = 0; bits < 8; ++bits)
    for (int axis = 0; axis < 3; ++axis)
      positions[Hex::bitToVertex[bits]][axis] = Hex::Interval((bits >> axis) & 1);
  return positions;
}

TEST(HexJacobianValidation, CertifiesRigidRotationAndDistortion)
{
  const auto rest = cube();
  auto positions = rest;
  const double c = std::cos(0.85), s = std::sin(0.85);
  for (auto &vertex : positions) {
    const auto x = vertex[0], y = vertex[1];
    vertex[0] = x * c - y * s;
    vertex[1] = x * s + y * c;
  }
  const auto rotated = Hex::certify(positions);
  EXPECT_TRUE(rotated.certified);
  EXPECT_LE(rotated.lower, 1);
  EXPECT_GE(rotated.upper, 1);
  EXPECT_TRUE(Hex::certify(rest, &positions).certified);
  positions[6][0] = positions[6][0] + 0.2;
  positions[5][2] = positions[5][2] * 0.85;
  EXPECT_TRUE(Hex::certify(positions).certified);
  EXPECT_TRUE(Hex::certify(rest, &positions).certified);
}

TEST(HexJacobianValidation, RejectsCornerInversionMissedByGaussPoints)
{
  auto positions = cube();
  for (int axis = 0; axis < 3; ++axis)
    positions[6][axis] = 0.5;
  const double g = (1 - 1 / std::sqrt(3.0)) / 2;
  for (double x : { g, 1 - g })
    for (double y : { g, 1 - g })
      for (double z : { g, 1 - g })
        EXPECT_GT(Hex::jacobian(positions, { x, y, z }).lo, 0);
  const auto result = Hex::certify(positions);
  EXPECT_FALSE(result.certified);
  EXPECT_TRUE(result.nonpositiveSample);
}

TEST(HexJacobianValidation, RejectsInversionBetweenPositiveEndpoints)
{
  const auto from = cube();
  auto to = from;
  for (auto &vertex : to) {
    vertex[0] = vertex[0] * -2;
    vertex[1] = vertex[1] * -0.5;
  }
  EXPECT_TRUE(Hex::certify(from).certified);
  EXPECT_TRUE(Hex::certify(to).certified);
  EXPECT_FALSE(Hex::certify(from, &to).certified);
}

TEST(HexJacobianValidation, PreservesRestOrientationAndRejectsSingularity)
{
  auto positions = cube();
  for (auto &vertex : positions)
    vertex[0] = vertex[0] * -1;
  EXPECT_TRUE(Hex::certify(positions, nullptr, -1).certified);
  EXPECT_FALSE(Hex::certify(positions).certified);
  for (auto &vertex : positions)
    vertex[2] = 0;
  EXPECT_FALSE(Hex::certify(positions).certified);
}

TEST(HexJacobianValidation, RejectsNonfiniteAndOverflowingGeometry)
{
  for (double bad : { std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity() }) {
    auto positions = cube();
    positions[6][0] = bad;
    EXPECT_FALSE(Hex::certify(positions).certified);
  }
  auto positions = cube();
  positions[0][0] = -std::numeric_limits<double>::max();
  positions[1][0] = std::numeric_limits<double>::max();
  EXPECT_FALSE(Hex::certify(positions).certified);
}
}  // namespace
