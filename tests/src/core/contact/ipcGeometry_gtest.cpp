#include <gtest/gtest.h>

#include "ipc/geometry/ipcBarrier.h"
#include "ipc/geometry/ipcCCD.h"
#include "ipc/geometry/ipcDistancePrimitives.h"
#include "ipc/geometry/ipcHessianProjection.h"

#include "testCIPCHelpers.h"

#include <array>
#include <limits>

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

TEST(IPCGeometryGTest, EdgeCCDRejectsCoplanarInteriorCrossingAndParallelOverlap)
{
  const ES::V3d a0(-1, 0, 0), a1(1, 0, 0), zero = ES::V3d::Zero();
  for (const auto &edge : { std::array<ES::V3d, 2>{ ES::V3d(0, -1, 0), ES::V3d(0, 1, 0) },
         std::array<ES::V3d, 2>{ ES::V3d(-0.5, 0, 0), ES::V3d(0.5, 0, 0) } }) {
    for (const auto &a : { a0, a1 })
      for (const auto &b : edge)
        ASSERT_GT((a - b).norm(), 0.0);
    EXPECT_DOUBLE_EQ(ccd::edgeEdgeCCD(a0, a1, edge[0], edge[1], zero, zero, zero, zero), 0.0);
    // Even a separating direction cannot certify an initially invalid state.
    const ES::V3d away(0, 0, 1);
    EXPECT_DOUBLE_EQ(ccd::edgeEdgeCCD(a0, a1, edge[0], edge[1], zero, zero, away, away), 0.0);
  }
}

TEST(IPCGeometryGTest, CCDDetectsInteriorCrossingsAtSmallScales)
{
  const ES::V3d zero = ES::V3d::Zero();
  for (double scale : { 1.0, 1e-12, 1e-32 }) {
    const ES::V3d a0(-scale, 0, 0), a1(scale, 0, 0);
    const ES::V3d b0(0, -scale, 0.25 * scale), b1(0, scale, 0.25 * scale);
    const ES::V3d motion(0, 0, -0.5 * scale);
    const double edgeStep = ccd::edgeEdgeCCD(a0, a1, b0, b1, zero, zero, motion, motion);
    EXPECT_GT(edgeStep, 0.0);
    EXPECT_LT(edgeStep, 0.5);

    const ES::V3d p(0.25 * scale, 0.25 * scale, 0.25 * scale);
    const double pointStep = ccd::pointTriangleCCD(p, zero, ES::V3d(scale, 0, 0), ES::V3d(0, scale, 0),
      motion, zero, zero, zero);
    EXPECT_GT(pointStep, 0.0);
    EXPECT_LT(pointStep, 0.5);
    EXPECT_DOUBLE_EQ(ccd::pointTriangleCCD(p, zero, ES::V3d(scale, 0, 0), ES::V3d(0, scale, 0),
      -motion, zero, zero, zero), 1.0);
  }
}

TEST(IPCGeometryGTest, CCDHandlesCollapsedAndCollapsingPrimitives)
{
  const ES::V3d zero = ES::V3d::Zero(), motion(0, 0, -0.2);
  // A collapsed edge is a point; its distance to an edge interior still matters.
  const double edgeStep = ccd::edgeEdgeCCD(zero, zero, ES::V3d(-1, 0, 0.1), ES::V3d(1, 0, 0.1),
    zero, zero, motion, motion);
  EXPECT_GT(edgeStep, 0.0);
  EXPECT_LT(edgeStep, 0.5);
  const ES::V3d p(0.5, 0, 0.1), a(0, 0, 0), c(2, 0, 0);
  for (const ES::V3d &b : { ES::V3d(1, 0, 0), ES::V3d(1, 1, 0) }) {
    const ES::V3d db(0, -2 * b[1], 0);
    // The second triangle collapses exactly when the point reaches its line.
    const double pointStep = ccd::pointTriangleCCD(p, a, b, c, motion, zero, db, zero);
    EXPECT_GT(pointStep, 0.0);
    EXPECT_LT(pointStep, 0.5);
  }
}

TEST(IPCGeometryGTest, CCDRejectsNonfiniteAndOverflowingInputs)
{
  const ES::V3d zero = ES::V3d::Zero(), x(1, 0, 0), y(0, 1, 0), p(0.25, 0.25, 0.1);
  for (double invalid : { std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity() }) {
    ES::V3d bad = zero;
    bad[0] = invalid;
    EXPECT_DOUBLE_EQ(ccd::pointTriangleCCD(bad, zero, x, y, zero, zero, zero, zero), 0.0);
    EXPECT_DOUBLE_EQ(ccd::pointTriangleCCD(p, zero, x, y, bad, zero, zero, zero), 0.0);
    EXPECT_DOUBLE_EQ(ccd::edgeEdgeCCD(zero, x, y, p, zero, zero, bad, zero), 0.0);
  }
  const ES::V3d huge(1e308, 0, 0);
  EXPECT_DOUBLE_EQ(ccd::pointTriangleCCD(p, zero, x, y, huge, zero, zero, zero), 0.0);
  EXPECT_DOUBLE_EQ(ccd::edgeEdgeCCD(zero, x, y, p, huge, zero, zero, zero), 0.0);
}

TEST(IPCGeometryGTest, CCDFailsClosedAtWorldPoseRoundingScale)
{
  const ES::V3d zero = ES::V3d::Zero();
  for (double origin : { 1.0, 1e8 }) {
    const double adjacent = std::nextafter(origin, std::numeric_limits<double>::infinity());
    const ES::V3d a(origin, origin, origin), b = a + ES::V3d(1, 0, 0), c = a + ES::V3d(0, 1, 0);
    const ES::V3d p(origin + 0.25, origin + 0.25, adjacent);
    const ES::V3d inward(0, 0, -2 * (adjacent - origin));
    // Distinct doubles alone do not prove a positive gap after pose
    // reconstruction. Repeated suffix advancement must stop at uncertainty.
    EXPECT_DOUBLE_EQ(ccd::pointTriangleCCD(p, a, b, c, inward, zero, zero, zero), 0.0);
    const ES::V3d e0(origin, origin + 0.5, adjacent), e1(origin + 1, origin + 0.5, adjacent);
    const ES::V3d f0(origin + 0.5, origin, origin), f1(origin + 0.5, origin + 1, origin);
    EXPECT_DOUBLE_EQ(ccd::edgeEdgeCCD(e0, e1, f0, f1, inward, inward, zero, zero), 0.0);

    const ES::V3d clear = p + ES::V3d(0, 0, 0.125), motion(0, 0, -0.25);
    const double bound = ccd::pointTriangleCCD(clear, a, b, c, motion, zero, zero, zero);
    EXPECT_GT(bound, 0.0);
    EXPECT_LT(bound, 0.501);
    EXPECT_DOUBLE_EQ(ccd::pointTriangleCCD(clear, a, b, c, -motion, zero, zero, zero), 1.0);
  }
}

TEST(IPCGeometryGTest, CCDRepeatedSuffixAdvancementCannotStepAcrossGrazingContact)
{
  const ES::V3d zero = ES::V3d::Zero();
  for (double scale : { 1.0, 1e-12, 1e-32 }) {
    const ES::V3d p(0.25 * scale, 0.25 * scale, 1e-8 * scale);
    const ES::V3d a = zero, b(scale, 0, 0), c(0, scale, 0), motion(1e-7 * scale, 0, -2e-8 * scale);
    double progress = 0;
    bool falselyCertified = false;
    for (int i = 0; i < 128; ++i) {
      const ES::V3d at = p + progress * motion;
      const ES::V3d remaining = p + motion - at;
      const double alpha = ccd::pointTriangleCCD(at, a, b, c, remaining, zero, zero, zero);
      if (alpha == 1) {
        falselyCertified = true;
        break;
      }
      const double next = std::nextafter(std::fma(alpha, 1 - progress, progress), progress);
      if (!(next > progress))
        break;
      progress = next;
    }
    EXPECT_FALSE(falselyCertified);
    EXPECT_LT(progress, 0.5);
  }
}

TEST(IPCGeometryGTest, CCDRepeatedSuffixPreservesDragonInteriorCrossing)
{
  // A real frame-431 dragon candidate crosses the moving triangle twice;
  // checking only the two endpoints or rounded closest-feature distances
  // incorrectly certified it. The first root is interior at t=.1162743857.
  ES::V12d initial, motion;
  initial << -1.267687175222354, -0.27383300068135347, 0.15184045565443113,
    -1.272781006433854, -0.27206091024380585, 0.1576676778371424,
    -1.2550081601772165, -0.2813765309687974, 0.15556259714647316,
    -1.2675543945013559, -0.2719720561437912, 0.14050790887572215;
  motion << 0.0002529060242464354, -0.002683482965247508, -0.0013030811841920065,
    -0.0009267287516814271, -0.0018603615414542474, -0.0015532419446279144,
    -0.0014158529183845303, -0.0020562010885734616, -0.0017674687776073061,
    -0.0010244743031080361, -0.0017888450058192173, -0.0016116669923125215;
  double progress = 0;
  for (int i = 0; i < 128; ++i) {
    const ES::V12d at = initial + progress * motion;
    const ES::V12d remaining = initial + motion - at;
    const double alpha = ccd::pointTriangleCCD(at.segment<3>(0), at.segment<3>(3), at.segment<3>(6), at.segment<3>(9),
      remaining.segment<3>(0), remaining.segment<3>(3), remaining.segment<3>(6), remaining.segment<3>(9));
    ASSERT_LT(alpha, 1.0);
    const double next = std::nextafter(std::fma(alpha, 1 - progress, progress), progress);
    if (!(next > progress))
      break;
    progress = next;
  }
  EXPECT_GT(progress, 0.0);
  EXPECT_LT(progress, 0.11627438567556836);
}

TEST(IPCGeometryGTest, CCDFullChordCertifiesLongTangentialNearMiss)
{
  const ES::V3d zero = ES::V3d::Zero();
  const ES::V3d a(0, 0, 0), b(1, 0, 0), c(0, 1, 0);
  const ES::V3d p(-2, 0.25, 1e-8), motion(4, 0, 0);
  // Speed-only advancement would require hundreds of millions of tiny
  // steps; a fixed separating plane certifies the complete linear path.
  EXPECT_DOUBLE_EQ(ccd::pointTriangleCCD(p, a, b, c, motion, zero, zero, zero), 1.0);
  EXPECT_DOUBLE_EQ(ccd::edgeEdgeCCD(ES::V3d(-2, -1, 1e-8), ES::V3d(-2, 1, 1e-8),
    ES::V3d(0, 0, 0), ES::V3d(1, 0, 0), motion, motion, zero, zero), 1.0);
}

TEST(IPCGeometryGTest, CCDPartialBoundRetainsFiniteInitialGapReserve)
{
  const ES::V3d zero = ES::V3d::Zero();
  for (double scale : { 1.0, 1e-12, 1e-32 }) {
    const ES::V3d p(0.25 * scale, 0.25 * scale, 0.1 * scale), motion(0, 0, -0.2 * scale);
    const double thickness = 0.01 * scale;
    const double pointStep = ccd::pointTriangleCCD(p, zero, ES::V3d(scale, 0, 0), ES::V3d(0, scale, 0),
      motion, zero, zero, zero, thickness);
    EXPECT_GT(pointStep, 0.3);
    EXPECT_LT(pointStep, 0.45);
    const double remainingPointGap = p[2] + pointStep * motion[2] - thickness;
    EXPECT_GE(remainingPointGap, 0.1 * (p[2] - thickness));

    const ES::V3d a(-scale, 0, 0), b(scale, 0, 0);
    const ES::V3d c(0, -scale, p[2]), d(0, scale, p[2]);
    const double edgeStep = ccd::edgeEdgeCCD(a, b, c, d, zero, zero, motion, motion, thickness);
    EXPECT_GT(edgeStep, 0.3);
    EXPECT_LT(edgeStep, 0.45);
    EXPECT_GE(c[2] + edgeStep * motion[2] - thickness, 0.1 * (c[2] - thickness));

    // The reserve is only a fallback for an incomplete proof. A complete
    // original chord remains certifiable, including a shorter requested time.
    EXPECT_DOUBLE_EQ(ccd::pointTriangleCCD(p, zero, ES::V3d(scale, 0, 0), ES::V3d(0, scale, 0),
      motion, zero, zero, zero, thickness, 0.25), 0.25);
  }
}
