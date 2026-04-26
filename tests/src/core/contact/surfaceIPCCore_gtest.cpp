#include <gtest/gtest.h>

#include "pgoLogging.h"
#include "ipc/core/surfaceIPCCore.h"

#include "testCIPCHelpers.h"

#include <stdexcept>

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
  ASSERT_FALSE(core.getPTPairs().empty());

  const std::size_t ptCount = core.getPTPairs().size();
  const std::size_t eeCount = core.getEEPairs().size();

  ES::VXd grad(x.size());
  core.computeGradient(x, grad);
  EXPECT_EQ(core.getPTPairs().size(), ptCount);
  EXPECT_EQ(core.getEEPairs().size(), eeCount);

  ES::SpMatD H;
  core.computeHessian(x, H);
  EXPECT_EQ(core.getPTPairs().size(), ptCount);
  EXPECT_EQ(core.getEEPairs().size(), eeCount);
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
  EXPECT_TRUE(core.isPreparedFor(x));

  const double preparedEnergy = core.computeEnergyWithPreparedPairs();
  ES::VXd preparedGradient = ES::VXd::Zero(x.size());
  core.computeGradientWithPreparedPairs(preparedGradient);
  ES::SpMatD preparedHessian;
  core.computeHessianWithPreparedPairs(preparedHessian);

  EXPECT_NEAR(preparedEnergy, directEnergy, 1e-12);
  EXPECT_LT(relativeError(preparedGradient, directGradient), 1e-12);
  EXPECT_LT(relativeError(sparseToDense(preparedHessian), sparseToDense(directHessian)), 1e-12);
}

TEST(SurfaceIPCCoreGTest, PreparedPairConsumersRequirePreparedState)
{
  SurfaceIPCCore core = makeConfiguredCore();
  EXPECT_THROW(core.computeEnergyWithPreparedPairs(), std::logic_error);
}
