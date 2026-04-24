#include <gtest/gtest.h>

#include "ipc/geometry/ipcBarrier.h"
#include "ipc/geometry/ipcCCD.h"
#include "ipc/geometry/ipcDistancePrimitives.h"
#include "ipc/geometry/ipcHessianProjection.h"

namespace
{
namespace ES = pgo::EigenSupport;
namespace barrier = pgo::Contact::CIPC::barrier;
namespace ccd = pgo::Contact::CIPC::ccd;
namespace distance = pgo::Contact::CIPC::distance;
using pgo::Contact::CIPC::projectToPSD;
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
