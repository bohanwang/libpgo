#include <gtest/gtest.h>

#include "pgoLogging.h"
#include "ipc/core/surfaceIPCCore.h"

#include "testCIPCHelpers.h"

#include <limits>
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

constexpr double kFDStep = 1e-5;
constexpr double kGradTol = 1e-4;
constexpr double kProjectedHessFdTol = 1e-1;
constexpr double kProjectedHessMinEigenTol = 1e-8;

std::tuple<ES::MXd, ES::MXi> makeParallelTriangleMesh(double separation)
{
  ES::MXd V(6, 3);
  V <<
    0.0, 0.0, 0.0,
    1.0, 0.0, 0.0,
    0.0, 1.0, 0.0,
    0.2, 0.2, separation,
    1.2, 0.2, separation,
    0.2, 1.2, separation;
  ES::MXi F(2, 3);
  F << 0, 1, 2,
    3, 4, 5;
  return { V, F };
}

SurfaceIPCCore makeConfiguredCore(bool projectHessianToPSD = true)
{
  SurfaceIPCCore core;
  SurfaceIPCCore::Parameters params;
  params.dhat = 0.1;
  params.kappa = 1.0;
  params.eps_ee = 0.0;
  params.slackness = 0.9;
  params.projectHessianToPSD = projectHessianToPSD;
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

TEST(SurfaceIPCCoreGTest, RejectsInvalidParameters)
{
  SurfaceIPCCore::Parameters params;
  SurfaceIPCCore core;

  params.dhat = 0.0;
  EXPECT_THROW(core.setParameters(params), std::invalid_argument);
  params.dhat = std::numeric_limits<double>::infinity();
  EXPECT_THROW(core.setParameters(params), std::invalid_argument);
  params.dhat = 0.1;

  params.kappa = -1.0;
  EXPECT_THROW(core.setParameters(params), std::invalid_argument);
  params.kappa = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(core.setParameters(params), std::invalid_argument);
  params.kappa = 1.0;

  params.eps_ee = -1.0;
  EXPECT_THROW(core.setParameters(params), std::invalid_argument);
  params.eps_ee = std::numeric_limits<double>::infinity();
  EXPECT_THROW(core.setParameters(params), std::invalid_argument);
  params.eps_ee = 0.0;

  params.slackness = 0.0;
  EXPECT_THROW(core.setParameters(params), std::invalid_argument);
  params.slackness = 1.1;
  EXPECT_THROW(core.setParameters(params), std::invalid_argument);
}

TEST(SurfaceIPCCoreGTest, RejectsMalformedMeshAndStateInputs)
{
  const auto [V, F] = makeTwoTriangleMesh();
  SurfaceIPCCore core;
  EXPECT_THROW(core.computeEnergy(ES::VXd()), std::logic_error);

  ES::MXd wrongVertexShape(V.rows(), 2);
  wrongVertexShape.setZero();
  EXPECT_THROW(core.setMesh(wrongVertexShape, F), std::invalid_argument);

  ES::MXd nonFiniteVertices = V;
  nonFiniteVertices(0, 0) = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(core.setMesh(nonFiniteVertices, F), std::invalid_argument);

  ES::MXi wrongTriangleShape(F.rows(), 2);
  wrongTriangleShape.setZero();
  EXPECT_THROW(core.setMesh(V, wrongTriangleShape), std::invalid_argument);

  ES::MXi outOfRangeTriangles = F;
  outOfRangeTriangles(0, 0) = V.rows();
  EXPECT_THROW(core.setMesh(V, outOfRangeTriangles), std::invalid_argument);

  ES::MXi repeatedVertexTriangle = F;
  repeatedVertexTriangle(0, 2) = repeatedVertexTriangle(0, 0);
  EXPECT_THROW(core.setMesh(V, repeatedVertexTriangle), std::invalid_argument);

  core.setMesh(V, F);
  const ES::VXd x = flattenPositions(V);
  EXPECT_THROW(core.computeEnergy(x.head(x.size() - 1)), std::invalid_argument);
  ES::VXd nonFiniteState = x;
  nonFiniteState[0] = std::numeric_limits<double>::infinity();
  EXPECT_THROW(core.computeEnergy(nonFiniteState), std::invalid_argument);

  ES::VXd wrongGradient = ES::VXd::Zero(x.size() - 1);
  EXPECT_THROW(core.computeGradient(x, wrongGradient), std::invalid_argument);
  ES::VXd nonFiniteDisplacement = ES::VXd::Zero(x.size());
  nonFiniteDisplacement[0] = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(core.computeMaxStepSize(x, nonFiniteDisplacement), std::invalid_argument);
}

TEST(SurfaceIPCCoreGTest, ProjectedHessianRemainsSymmetricPSDAndTracksFiniteDifference)
{
  SurfaceIPCCore core = makeConfiguredCore();
  EXPECT_TRUE(core.getParameters().projectHessianToPSD);
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

TEST(SurfaceIPCCoreGTest, UnprojectedHessianMatchesGradientFiniteDifference)
{
  SurfaceIPCCore core = makeConfiguredCore(false);
  EXPECT_FALSE(core.getParameters().projectHessianToPSD);
  const auto [V, _] = makeTwoTriangleMesh();
  const ES::VXd x = flattenPositions(V);

  ES::SpMatD H;
  core.computeHessian(x, H);
  const ES::MXd analyticH = sparseToDense(H);
  const ES::MXd fdH = finiteDifferenceHessian(
    [&core](const ES::VXd &state) {
      ES::VXd gradient = ES::VXd::Zero(state.size());
      core.computeGradient(state, gradient);
      return gradient;
    },
    x,
    kFDStep);

  EXPECT_LT(relativeError(analyticH, analyticH.transpose()), 1e-12);
  EXPECT_LT(relativeError(analyticH, fdH), 2e-5);
}

TEST(SurfaceIPCCoreGTest, ProjectionOptionChangesHessianAndComputeAllMatchesSeparateCalls)
{
  const auto [V, _] = makeTwoTriangleMesh();
  const ES::VXd x = flattenPositions(V);
  SurfaceIPCCore projectedCore = makeConfiguredCore(true);
  SurfaceIPCCore rawCore = makeConfiguredCore(false);

  ES::SpMatD projectedHessian;
  ES::SpMatD rawHessian;
  projectedCore.computeHessian(x, projectedHessian);
  rawCore.computeHessian(x, rawHessian);
  const ES::MXd projectedDense = sparseToDense(projectedHessian);
  const ES::MXd rawDense = sparseToDense(rawHessian);
  EXPECT_GT(relativeError(projectedDense, rawDense), 1e-6);

  for (SurfaceIPCCore *core : { &projectedCore, &rawCore }) {
    const double energy = core->computeEnergy(x);
    ES::VXd gradient = ES::VXd::Zero(x.size());
    core->computeGradient(x, gradient);
    ES::SpMatD hessian;
    core->computeHessian(x, hessian);

    double allEnergy = 0.0;
    ES::VXd allGradient;
    ES::SpMatD allHessian;
    core->computeAll(x, allEnergy, allGradient, allHessian);
    EXPECT_NEAR(allEnergy, energy, 1e-12);
    EXPECT_LT(relativeError(allGradient, gradient), 1e-12);
    EXPECT_LT(relativeError(sparseToDense(allHessian), sparseToDense(hessian)), 1e-12);
  }
}

TEST(SurfaceIPCCoreGTest, RebuildsPairsAndHessianPatternAcrossBarrierCutoff)
{
  const auto [V, F] = makeParallelTriangleMesh(0.11);
  SurfaceIPCCore::Parameters params;
  params.dhat = 0.1;
  params.kappa = 1.0;
  params.slackness = 0.9;
  SurfaceIPCCore core(params);
  core.setMesh(V, F);

  const ES::VXd outside = flattenPositions(V);
  ES::VXd inside = outside;
  for (int vi = 3; vi < 6; ++vi)
    inside[3 * vi + 2] -= 0.02;

  ES::VXd cutoffOutside = outside;
  ES::VXd cutoffInside = outside;
  for (int vi = 3; vi < 6; ++vi) {
    cutoffOutside[3 * vi + 2] -= 0.009999;
    cutoffInside[3 * vi + 2] -= 0.010001;
  }

  EXPECT_DOUBLE_EQ(core.computeEnergy(cutoffOutside), 0.0);
  ES::VXd cutoffOutsideGradient = ES::VXd::Zero(cutoffOutside.size());
  core.computeGradient(cutoffOutside, cutoffOutsideGradient);
  EXPECT_EQ(cutoffOutsideGradient.norm(), 0.0);
  const double cutoffInsideEnergy = core.computeEnergy(cutoffInside);
  ES::VXd cutoffInsideGradient = ES::VXd::Zero(cutoffInside.size());
  core.computeGradient(cutoffInside, cutoffInsideGradient);
  EXPECT_GT(cutoffInsideEnergy, 0.0);
  EXPECT_LT(cutoffInsideEnergy, 1e-9);
  EXPECT_LT(cutoffInsideGradient.norm(), 1e-5);

  const double outsideEnergy = core.computeEnergy(outside);
  EXPECT_DOUBLE_EQ(outsideEnergy, 0.0);
  EXPECT_EQ(core.getPTPairs().size() + core.getEEPairs().size(), 0);
  ES::SpMatD outsideHessian;
  core.computeHessian(outside, outsideHessian);
  EXPECT_EQ(outsideHessian.nonZeros(), 0);

  const double insideEnergy = core.computeEnergy(inside);
  EXPECT_GT(insideEnergy, 0.0);
  const std::size_t insidePairCount = core.getPTPairs().size() + core.getEEPairs().size();
  EXPECT_GT(insidePairCount, 0);
  ES::SpMatD insideHessian;
  core.computeHessian(inside, insideHessian);
  EXPECT_GT(insideHessian.nonZeros(), 0);

  ES::VXd outsideGradient = ES::VXd::Zero(outside.size());
  core.computeGradient(outside, outsideGradient);
  ES::SpMatD rebuiltOutsideHessian;
  core.computeHessian(outside, rebuiltOutsideHessian);
  EXPECT_EQ(core.getPTPairs().size() + core.getEEPairs().size(), 0);
  EXPECT_EQ(rebuiltOutsideHessian.nonZeros(), 0);
  EXPECT_EQ(outsideGradient.norm(), 0.0);

  ES::SpMatD rebuiltInsideHessian;
  core.computeHessian(inside, rebuiltInsideHessian);
  EXPECT_EQ(core.getPTPairs().size() + core.getEEPairs().size(), insidePairCount);
  EXPECT_LT(relativeError(sparseToDense(rebuiltInsideHessian), sparseToDense(insideHessian)), 1e-12);
}

TEST(SurfaceIPCCoreGTest, ComputeMaxStepSizeDetectsImpendingCollision)
{
  initializeLogging();
  SurfaceIPCCore core = makeConfiguredCore();
  const auto [V, _] = makeTwoTriangleMesh();
  const ES::VXd x = flattenPositions(V);

  ES::VXd dx = ES::VXd::Zero(x.size());
  for (int vi = 3; vi < 6; ++vi)
    dx[3 * vi + 2] = -0.1;

  const double alpha = core.computeMaxStepSize(x, dx);
  EXPECT_GT(alpha, 0.0);
  EXPECT_LT(alpha, 1.0);
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
