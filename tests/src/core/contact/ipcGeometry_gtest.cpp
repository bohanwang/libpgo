#include <gtest/gtest.h>

#include "ipc/geometry/ipcBarrier.h"
#include "ipc/geometry/ipcCCD.h"
#include "ipc/geometry/ipcDistancePrimitives.h"
#include "ipc/geometry/ipcHessianProjection.h"

#include "testCIPCHelpers.h"

#include <array>

namespace
{
namespace ES = pgo::EigenSupport;
namespace barrier = pgo::Contact::CIPC::barrier;
namespace ccd = pgo::Contact::CIPC::ccd;
namespace distance = pgo::Contact::CIPC::distance;
using pgo::Contact::CIPC::projectToPSD;
using pgo::Contact::CIPCTest::finiteDifferenceGradient;
using pgo::Contact::CIPCTest::finiteDifferenceHessian;
using pgo::Contact::CIPCTest::relativeError;

ES::V12d makePointTriangleState(const ES::V3d &p)
{
  ES::V12d x;
  x << p,
    0.0, 0.0, 0.0,
    1.0, 0.0, 0.0,
    0.0, 1.0, 0.0;
  return x;
}

ES::V12d makeEdgeEdgeState(double x, double y)
{
  ES::V12d state;
  state << 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0,
    x, y, 0.1,
    x, y + 1.0, 0.1;
  return state;
}

void expectPointTriangleDerivatives(const ES::V12d &x)
{
  const auto energy = [](const ES::VXd &state) {
    return distance::computePTSqDist(state.segment<3>(0), state.segment<3>(3),
      state.segment<3>(6), state.segment<3>(9));
  };
  const auto gradient = [](const ES::VXd &state) {
    return ES::VXd(distance::computePTSqDistGrad(state.segment<3>(0), state.segment<3>(3),
      state.segment<3>(6), state.segment<3>(9)));
  };
  const ES::VXd analyticGradient = gradient(x);
  const ES::VXd fdGradient = finiteDifferenceGradient(energy, x, 1e-6);
  const ES::MXd analyticHessian = distance::computePTSqDistHess(x.segment<3>(0), x.segment<3>(3),
    x.segment<3>(6), x.segment<3>(9));
  const ES::MXd fdHessian = finiteDifferenceHessian(gradient, x, 1e-6);
  EXPECT_LT(relativeError(analyticGradient, fdGradient), 2e-5);
  EXPECT_LT(relativeError(analyticHessian, fdHessian), 2e-5);
}

void expectEdgeEdgeDerivatives(const ES::V12d &x)
{
  const auto energy = [](const ES::VXd &state) {
    return distance::computeEESqDist(state.segment<3>(0), state.segment<3>(3),
      state.segment<3>(6), state.segment<3>(9));
  };
  const auto gradient = [](const ES::VXd &state) {
    return ES::VXd(distance::computeEESqDistGrad(state.segment<3>(0), state.segment<3>(3),
      state.segment<3>(6), state.segment<3>(9)));
  };
  const ES::VXd analyticGradient = gradient(x);
  const ES::VXd fdGradient = finiteDifferenceGradient(energy, x, 1e-6);
  const ES::MXd analyticHessian = distance::computeEESqDistHess(x.segment<3>(0), x.segment<3>(3),
    x.segment<3>(6), x.segment<3>(9));
  const ES::MXd fdHessian = finiteDifferenceHessian(gradient, x, 1e-6);
  EXPECT_LT(relativeError(analyticGradient, fdGradient), 2e-5);
  EXPECT_LT(relativeError(analyticHessian, fdHessian), 2e-5);
}
}  // namespace

TEST(IPCGeometryGTest, BarrierMatchesFiniteDifferenceDerivatives)
{
  constexpr double shat = 0.04;
  constexpr double s = 0.01;
  constexpr double h = 1e-6;

  const double fd1 = (barrier::b(s + h, shat) - barrier::b(s - h, shat)) / (2.0 * h);
  const double fd2 = (barrier::dbds(s + h, shat) - barrier::dbds(s - h, shat)) / (2.0 * h);

  EXPECT_GT(barrier::b(s, shat), 0.0);
  EXPECT_NEAR(barrier::dbds(s, shat), fd1, 1e-6);
  EXPECT_NEAR(barrier::d2bds2(s, shat), fd2, 1e-3);
  EXPECT_DOUBLE_EQ(barrier::b(shat, shat), 0.0);
}

TEST(IPCGeometryGTest, PointTriangleDistanceUsesFaceDistanceForInteriorProjection)
{
  const ES::V3d p(0.25, 0.25, 0.1);
  const ES::V3d t0(0.0, 0.0, 0.0);
  const ES::V3d t1(1.0, 0.0, 0.0);
  const ES::V3d t2(0.0, 1.0, 0.0);

  EXPECT_NEAR(distance::computePTSqDist(p, t0, t1, t2), 0.01, 1e-14);
  EXPECT_EQ(distance::classifyPT(p, t0, t1, t2), pgo::Contact::CIPC::PTDistType::PT);
}

TEST(IPCGeometryGTest, PointTriangleDispatchModesMatchFiniteDifferences)
{
  const std::array<std::pair<ES::V3d, pgo::Contact::CIPC::PTDistType>, 7> cases = {
    std::pair { ES::V3d(-0.2, -0.2, 0.1), pgo::Contact::CIPC::PTDistType::PP_PT0 },
    std::pair { ES::V3d(1.2, -0.2, 0.1), pgo::Contact::CIPC::PTDistType::PP_PT1 },
    std::pair { ES::V3d(-0.2, 1.2, 0.1), pgo::Contact::CIPC::PTDistType::PP_PT2 },
    std::pair { ES::V3d(0.5, -0.2, 0.1), pgo::Contact::CIPC::PTDistType::PE_PT0T1 },
    std::pair { ES::V3d(0.6, 0.6, 0.1), pgo::Contact::CIPC::PTDistType::PE_PT1T2 },
    std::pair { ES::V3d(-0.2, 0.5, 0.1), pgo::Contact::CIPC::PTDistType::PE_PT2T0 },
    std::pair { ES::V3d(0.25, 0.25, 0.1), pgo::Contact::CIPC::PTDistType::PT },
  };

  for (const auto &[point, type] : cases) {
    const ES::V12d x = makePointTriangleState(point);
    EXPECT_EQ(distance::classifyPT(x.segment<3>(0), x.segment<3>(3), x.segment<3>(6), x.segment<3>(9)), type);
    expectPointTriangleDerivatives(x);
  }
}

TEST(IPCGeometryGTest, EdgeEdgeDispatchModesMatchFiniteDifferences)
{
  using Type = pgo::Contact::CIPC::EEDistType;
  const std::array<std::tuple<double, double, Type>, 9> cases = {
    std::tuple { -0.2, 0.2, Type::PP_Ea0Eb0 },
    std::tuple { -0.2, -1.2, Type::PP_Ea0Eb1 },
    std::tuple { 1.2, 0.2, Type::PP_Ea1Eb0 },
    std::tuple { 1.2, -1.2, Type::PP_Ea1Eb1 },
    std::tuple { -0.2, -0.5, Type::PE_Ea0_Eb },
    std::tuple { 1.2, -0.5, Type::PE_Ea1_Eb },
    std::tuple { 0.5, 0.2, Type::PE_Eb0_Ea },
    std::tuple { 0.5, -1.2, Type::PE_Eb1_Ea },
    std::tuple { 0.5, -0.5, Type::EE },
  };

  for (const auto &[edgeX, edgeY, type] : cases) {
    const ES::V12d x = makeEdgeEdgeState(edgeX, edgeY);
    EXPECT_EQ(distance::classifyEE(x.segment<3>(0), x.segment<3>(3), x.segment<3>(6), x.segment<3>(9)), type);
    expectEdgeEdgeDerivatives(x);
  }
}

TEST(IPCGeometryGTest, PointTriangleCCDDetectsCrossing)
{
  const ES::V3d p(0.25, 0.25, 0.1);
  const ES::V3d t0(0.0, 0.0, 0.0);
  const ES::V3d t1(1.0, 0.0, 0.0);
  const ES::V3d t2(0.0, 1.0, 0.0);
  const ES::V3d dp(0.0, 0.0, -0.2);
  const ES::V3d dz = ES::V3d::Zero();

  const double toi = ccd::pointTriangleCCD(p, t0, t1, t2, dp, dz, dz, dz);

  EXPECT_GT(toi, 0.0);
  EXPECT_LT(toi, 1.0);
}

TEST(IPCGeometryGTest, HessianProjectionClampsNegativeEigenvalues)
{
  ES::M12d H = ES::M12d::Identity();
  H(0, 0) = -2.0;

  const ES::M12d projected = projectToPSD(H);
  const Eigen::SelfAdjointEigenSolver<ES::M12d> eig(projected);

  ASSERT_EQ(eig.info(), Eigen::Success);
  EXPECT_GE(eig.eigenvalues().minCoeff(), -1e-12);
}
