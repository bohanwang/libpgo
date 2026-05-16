#include <gtest/gtest.h>

#include "pgoLogging.h"
#include "ipc/core/surfaceIPCCore.h"
#include "ipc/external/obstacleSurface.h"

#include "testCIPCHelpers.h"

#include <stdexcept>
#include <type_traits>

namespace
{
namespace ES = pgo::EigenSupport;
using pgo::Contact::CIPC::SurfaceIPCCore;
using pgo::Contact::CIPCTest::computeFloorEnergy;
using pgo::Contact::CIPCTest::computeFloorGradient;
using pgo::Contact::CIPCTest::computeFloorHessian;
using pgo::Contact::CIPCTest::finiteDifferenceGradient;
using pgo::Contact::CIPCTest::finiteDifferenceHessian;
using pgo::Contact::CIPCTest::flattenPositions;
using pgo::Contact::CIPCTest::makeTwoTriangleMesh;
using pgo::Contact::CIPCTest::relativeError;
using pgo::Contact::CIPCTest::sparseToDense;
using pgo::NonlinearOptimization::MaxStepResult;
using pgo::NonlinearOptimization::SolveDiagnostics;

constexpr double kFDStep = 1e-5;
constexpr double kGradTol = 1e-4;
constexpr double kProjectedHessFdTol = 1e-1;
constexpr double kProjectedHessMinEigenTol = 1e-8;

SurfaceIPCCore makeConfiguredCore()
{
  SurfaceIPCCore core;
  SurfaceIPCCore::Parameters params;
  params.dhat = 0.1;
  params.kappa = 1.0;
  params.eps_ee = 0.0;
  params.slackness = 0.9;
  core.setParameters(params);

  auto [V, F] = makeTwoTriangleMesh();
  core.setMesh(V, F);
  return core;
}

void initializeLogging()
{
  static const bool initialized = []() {
    pgo::Logging::init();
    return true;
  }();
  (void)initialized;
}
}  // namespace

TEST(SurfaceIPCCoreGTest, EnergyGradientMatchesFiniteDifference)
{
  SurfaceIPCCore core = makeConfiguredCore();
  const auto [V, _] = makeTwoTriangleMesh();
  const ES::VXd x = flattenPositions(V);

  ES::VXd analyticGrad(x.size());
  core.computeGradient(x, analyticGrad);

  const ES::VXd fdGrad = finiteDifferenceGradient(
    [&core](const ES::VXd &state) { return core.computeEnergy(state); },
    x,
    kFDStep);

  EXPECT_LT(relativeError(analyticGrad, fdGrad), kGradTol);
}

TEST(SurfaceIPCCoreGTest, ProjectedHessianRemainsSymmetricPSDAndTracksFiniteDifference)
{
  SurfaceIPCCore core = makeConfiguredCore();
  const auto [V, _] = makeTwoTriangleMesh();
  const ES::VXd x = flattenPositions(V);

  ES::SpMatD H;
  core.computeHessian(x, H);
  const ES::MXd analyticH = sparseToDense(H);

  const ES::MXd fdH = finiteDifferenceHessian(
    [&core](const ES::VXd &state) {
      ES::VXd grad(state.size());
      core.computeGradient(state, grad);
      return grad;
    },
    x,
    kFDStep);

  const ES::MXd symAnalyticH = 0.5 * (analyticH + analyticH.transpose());
  const ES::MXd symFdH = 0.5 * (fdH + fdH.transpose());
  const Eigen::SelfAdjointEigenSolver<ES::MXd> eig(symAnalyticH);

  // The returned Hessian is blockwise PSD-projected by design, so exact
  // equality to the raw finite-difference Jacobian of gradient is not
  // expected. We still require it to remain symmetric, PSD, and close in
  // overall scale to the unprojected finite-difference reference.
  EXPECT_LT(relativeError(analyticH, analyticH.transpose()), 1e-12);
  ASSERT_EQ(eig.info(), Eigen::Success);
  EXPECT_GE(eig.eigenvalues().minCoeff(), -kProjectedHessMinEigenTol);
  EXPECT_LT(relativeError(symAnalyticH, symFdH), kProjectedHessFdTol);
}

TEST(SurfaceIPCCoreGTest, ComputeMaxStepLimitDetectsImpendingCollision)
{
  initializeLogging();
  SurfaceIPCCore core = makeConfiguredCore();
  const auto [V, _] = makeTwoTriangleMesh();
  const ES::VXd x = flattenPositions(V);

  ES::VXd dx = ES::VXd::Zero(x.size());
  for (int vi = 3; vi < 6; ++vi)
    dx[3 * vi + 2] = -0.1;

  const MaxStepResult result = core.computeMaxStepLimit(x, dx);
  const double alpha = result.alpha;
  EXPECT_GT(alpha, 0.0);
  EXPECT_LT(alpha, 1.0);
  EXPECT_DOUBLE_EQ(result.contactAlpha, alpha);
  EXPECT_TRUE(result.contactClamped);
}

TEST(SurfaceIPCCoreGTest, ComputeMaxStepLimitCanBeRecordedInDiagnostics)
{
  initializeLogging();
  SurfaceIPCCore core = makeConfiguredCore();
  const auto [V, _] = makeTwoTriangleMesh();
  const ES::VXd x = flattenPositions(V);

  ES::VXd dx = ES::VXd::Zero(x.size());
  for (int vi = 3; vi < 6; ++vi)
    dx[3 * vi + 2] = -0.1;

  const MaxStepResult result = core.computeMaxStepLimit(x, dx);
  const double alpha = result.alpha;
  EXPECT_GT(alpha, 0.0);
  EXPECT_LT(alpha, 1.0);
  EXPECT_TRUE(result.contactClamped);

  SolveDiagnostics diagnostics;
  diagnostics.recordMaxStep(result);
  EXPECT_EQ(diagnostics.contactClampCount, 1);
  EXPECT_DOUBLE_EQ(diagnostics.minContactFeasibleAlpha, alpha);
}

TEST(SurfaceIPCCoreGTest, SmallContactAlphaWarnsAndDiagnosticsResetClearsStats)
{
  initializeLogging();
  SurfaceIPCCore core = makeConfiguredCore();
  const auto [V, _] = makeTwoTriangleMesh();
  const ES::VXd x = flattenPositions(V);

  ES::VXd dx = ES::VXd::Zero(x.size());
  for (int vi = 3; vi < 6; ++vi)
    dx[3 * vi + 2] = -10.0;

  testing::internal::CaptureStdout();
  const MaxStepResult result = core.computeMaxStepLimit(x, dx);
  const std::string logOutput = testing::internal::GetCapturedStdout();

  const double alpha = result.alpha;
  EXPECT_GT(alpha, 0.0);
  EXPECT_LT(alpha, 0.01);
  EXPECT_TRUE(result.contactClamped);
  EXPECT_NE(logOutput.find("contactFeasibleAlpha"), std::string::npos);

  SolveDiagnostics diagnostics;
  diagnostics.recordMaxStep(result);
  ASSERT_EQ(diagnostics.contactClampCount, 1);
  ASSERT_DOUBLE_EQ(diagnostics.minContactFeasibleAlpha, alpha);

  diagnostics.reset();
  EXPECT_EQ(diagnostics.contactClampCount, 0);
  EXPECT_DOUBLE_EQ(diagnostics.minContactFeasibleAlpha, 1.0);
}

TEST(SurfaceIPCCoreGTest, PairAccessorsRemainReadableAcrossComputes)
{
  SurfaceIPCCore core = makeConfiguredCore();
  const auto [V, _] = makeTwoTriangleMesh();
  const ES::VXd x = flattenPositions(V);

  const double energy = core.computeEnergy(x);
  EXPECT_GT(energy, 0.0);
  ASSERT_FALSE(core.preparedState().selfPairs.ptPairs.empty());

  const std::size_t ptCount = core.preparedState().selfPairs.ptPairs.size();
  const std::size_t eeCount = core.preparedState().selfPairs.eePairs.size();

  ES::VXd grad(x.size());
  core.computeGradient(x, grad);
  EXPECT_EQ(core.preparedState().selfPairs.ptPairs.size(), ptCount);
  EXPECT_EQ(core.preparedState().selfPairs.eePairs.size(), eeCount);

  ES::SpMatD H;
  core.computeHessian(x, H);
  EXPECT_EQ(core.preparedState().selfPairs.ptPairs.size(), ptCount);
  EXPECT_EQ(core.preparedState().selfPairs.eePairs.size(), eeCount);
}

TEST(SurfaceIPCCoreGTest, PreparedPairsMatchDirectEnergyGradientHessian)
{
  SurfaceIPCCore core = makeConfiguredCore();
  const auto [V, _] = makeTwoTriangleMesh();
  const ES::VXd x = flattenPositions(V);

  const double directEnergy = core.computeEnergy(x);
  ES::VXd directGradient = ES::VXd::Zero(x.size());
  core.computeGradient(x, directGradient);
  ES::SpMatD directHessian;
  core.computeHessian(x, directHessian);

  core.prepareForSurfacePositions(x);
  EXPECT_TRUE(core.preparedState().isPreparedFor(x));

  const double preparedEnergy = core.computeEnergyWithPreparedPairs();
  ES::VXd preparedGradient = ES::VXd::Zero(x.size());
  core.computeGradientWithPreparedPairs(preparedGradient);
  ES::SpMatD preparedHessian;
  core.computeHessianWithPreparedPairs(preparedHessian);

  EXPECT_NEAR(preparedEnergy, directEnergy, 1e-12);
  EXPECT_LT(relativeError(preparedGradient, directGradient), 1e-12);
  EXPECT_LT(relativeError(sparseToDense(preparedHessian), sparseToDense(directHessian)), 1e-12);
}

TEST(SurfaceIPCCoreGTest, PreparedStateAccessorIsReadOnlyAndExplicitlyInvalidated)
{
  static_assert(
    std::is_same_v<decltype(std::declval<const SurfaceIPCCore &>().preparedState()),
      const pgo::Contact::CIPC::SurfaceIPCPreparedState &>);

  SurfaceIPCCore core = makeConfiguredCore();
  const auto [V, _] = makeTwoTriangleMesh();
  const ES::VXd x = flattenPositions(V);

  core.prepareForSurfacePositions(x);
  ASSERT_TRUE(core.preparedState().hasState);

  const SurfaceIPCCore &constCore = core;
  constCore.invalidatePreparedState();
  EXPECT_FALSE(core.preparedState().hasState);
}

TEST(SurfaceIPCCoreGTest, ConstructorInjectedObstaclesAssignSequentialSlots)
{
  using pgo::Contact::CIPC::ObstacleSurface;

  ES::MXd dynV(4, 3);
  dynV << 0.0, 0.0, 0.0,
          1.0, 0.0, 0.0,
          0.0, 1.0, 0.0,
          1.0, 1.0, 0.0;
  ES::MXi dynF(2, 3);
  dynF << 0, 1, 2,
          1, 3, 2;

  auto buildPlaneObstacle = [](double zOffset) {
    ES::MXd obsV(3, 3);
    obsV << 0.0, 0.0, zOffset,
            1.0, 0.0, zOffset,
            0.0, 1.0, zOffset;
    ES::MXi obsF(1, 3);
    obsF << 0, 1, 2;
    ES::VXd rest(obsV.rows() * 3);
    for (int vi = 0; vi < obsV.rows(); ++vi)
      rest.segment<3>(3 * vi) = obsV.row(vi).transpose();
    return ObstacleSurface(obsV, obsF,
      pgo::Contact::CIPC::makeLinearTrajectorySampler(rest, ES::V3d::Zero()));
  };

  std::vector<ObstacleSurface> obstacles;
  obstacles.emplace_back(buildPlaneObstacle( 0.3));
  obstacles.emplace_back(buildPlaneObstacle(-0.3));
  for (auto &obs : obstacles)
    obs.update(0.0, 0.0);

  SurfaceIPCCore::Parameters params;
  params.dhat = 0.1;
  params.dhat_external = 1.0;
  params.kappa = 1.0;
  params.eps_ee = 0.0;
  params.slackness = 0.9;
  SurfaceIPCCore core(params, std::move(obstacles));
  core.setMesh(dynV, dynF);

  const ES::VXd x = flattenPositions(dynV);
  core.prepareForSurfacePositions(x);

  bool sawSlot0 = false;
  bool sawSlot1 = false;
  auto scanSlots = [&](auto &&vec) {
    for (const auto &pair : vec) {
      if (pair.obstacleSlot == 0) sawSlot0 = true;
      if (pair.obstacleSlot == 1) sawSlot1 = true;
    }
  };
  scanSlots(core.preparedState().externalPairs.ptPairs);
  scanSlots(core.preparedState().externalPairs.tpPairs);
  scanSlots(core.preparedState().externalPairs.eePairs);
  EXPECT_TRUE(sawSlot0);
  EXPECT_TRUE(sawSlot1);
}

TEST(SurfaceIPCCoreGTest, ObstacleSurfaceEmptySamplerThrows)
{
  ES::MXd V(3, 3);
  V.setZero();
  ES::MXi F(1, 3);
  F << 0, 1, 2;

  EXPECT_THROW(
    pgo::Contact::CIPC::ObstacleSurface(V, F, pgo::Contact::CIPC::ObstacleSurface::TrajectorySampler{}),
    std::invalid_argument);
}

TEST(SurfaceIPCCoreGTest, ObstacleSurfaceStoresRestPositionsRowWise)
{
  using pgo::Contact::CIPC::ObstacleSurface;

  ES::MXd V(2, 3);
  V << 1.0, 2.0, 3.0,
       4.0, 5.0, 6.0;
  ES::MXi F(0, 3);

  auto sampler = [](double, ES::RefVecXd out) { out.setZero(); };
  ObstacleSurface obs(V, F, sampler);

  ES::VXd expected(6);
  expected << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0;
  EXPECT_TRUE(obs.restPositions().isApprox(expected));
}

TEST(SurfaceIPCCoreGTest, ObstacleSurfaceInvalidVertexColumnCountThrows)
{
  using pgo::Contact::CIPC::ObstacleSurface;

  ES::MXd V(2, 4);
  V.setZero();
  ES::MXi F(0, 3);

  auto sampler = [](double, ES::RefVecXd out) { out.setZero(); };
  EXPECT_THROW(ObstacleSurface(V, F, sampler), std::invalid_argument);
}

TEST(SurfaceIPCCoreGTest, PreparedPairConsumersRequirePreparedState)
{
  SurfaceIPCCore core = makeConfiguredCore();
  EXPECT_THROW(core.computeEnergyWithPreparedPairs(), std::logic_error);
}
