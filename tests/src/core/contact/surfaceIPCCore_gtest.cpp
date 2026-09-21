#include <gtest/gtest.h>

#include "pgoLogging.h"
#include "embeddedSurfaceIPCPotentialEnergy.h"
#include "ipc/core/surfaceIPCCore.h"

#include "testCIPCHelpers.h"

#include <algorithm>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <tuple>

namespace
{
namespace ES = pgo::EigenSupport;
using pgo::Contact::CIPC::SurfaceIPCCore;
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

// Two n x n triangulated sheets a small distance apart: many point-triangle
// and edge-edge pairs, enough for parallel scheduling to matter.
std::tuple<ES::MXd, ES::MXi> makeSheetPairMesh(int n, double separation)
{
  const int verticesPerSheet = (n + 1) * (n + 1);
  ES::MXd V(2 * verticesPerSheet, 3);
  ES::MXi F(2 * 2 * n * n, 3);
  int f = 0;
  for (int sheet = 0; sheet < 2; ++sheet) {
    const int base = sheet * verticesPerSheet;
    const double shift = sheet * 0.37 / n;  // avoid exactly aligned vertices
    for (int j = 0; j <= n; ++j)
      for (int i = 0; i <= n; ++i)
        V.row(base + j * (n + 1) + i) << static_cast<double>(i) / n + shift,
          static_cast<double>(j) / n + shift, sheet * separation;
    for (int j = 0; j < n; ++j)
      for (int i = 0; i < n; ++i) {
        const int v00 = base + j * (n + 1) + i, v10 = v00 + 1, v01 = v00 + n + 1, v11 = v01 + 1;
        F.row(f++) << v00, v10, v11;
        F.row(f++) << v00, v11, v01;
      }
  }
  return { V, F };
}

bool bitwiseEqual(const ES::VXd &a, const ES::VXd &b)
{
  return a.size() == b.size() && std::memcmp(a.data(), b.data(), sizeof(double) * a.size()) == 0;
}

bool bitwiseEqual(const ES::SpMatD &a, const ES::SpMatD &b)
{
  return a.nonZeros() == b.nonZeros() &&
    std::memcmp(a.innerIndexPtr(), b.innerIndexPtr(), sizeof(ES::SpMatD::StorageIndex) * a.nonZeros()) == 0 &&
    std::memcmp(a.valuePtr(), b.valuePtr(), sizeof(double) * a.nonZeros()) == 0;
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

TEST(SurfaceIPCCoreGTest, DeterministicAssemblyMatchesDefaultSortsPairsAndReproducesBitwise)
{
  initializeLogging();
  const auto [V, F] = makeSheetPairMesh(12, 0.05);
  const ES::VXd x = flattenPositions(V);

  SurfaceIPCCore::Parameters params;
  params.dhat = 0.1;
  params.kappa = 1.0;
  params.frictionCoeff = 0.3;
  SurfaceIPCCore core(params);
  core.setMesh(V, F);
  EXPECT_FALSE(core.isDeterministic());

  struct Result
  {
    double energy = 0.0;
    ES::VXd gradient;
    ES::SpMatD hessian;
    double allEnergy = 0.0;
    ES::VXd allGradient;
    ES::SpMatD allHessian;
  };
  auto evaluate = [&]() {
    core.updateFriction(x, x, 0.01);
    Result r;
    r.energy = core.computeEnergy(x);
    r.gradient = ES::VXd::Zero(x.size());
    core.computeGradient(x, r.gradient);
    core.computeHessian(x, r.hessian);
    r.allGradient = ES::VXd::Zero(x.size());
    core.computeAll(x, r.allEnergy, r.allGradient, r.allHessian);
    return r;
  };

  const Result base = evaluate();
  ASSERT_GT(core.getPTPairs().size() + core.getEEPairs().size(), 100u);
  EXPECT_GT(core.getNumFrictionPairs(), 0u);
  EXPECT_GT(base.gradient.norm(), 0.0);

  core.setDeterministic(true);
  EXPECT_TRUE(core.isDeterministic());
  const Result det = evaluate();
  const Result again = evaluate();

  // canonical pair order after every broad phase
  const auto &pt = core.getPTPairs();
  const auto &ee = core.getEEPairs();
  EXPECT_TRUE(std::is_sorted(pt.begin(), pt.end(), [](const pgo::Contact::CIPC::PTPair &a, const pgo::Contact::CIPC::PTPair &b) {
    return std::tie(a.p, a.t0, a.t1, a.t2) < std::tie(b.p, b.t0, b.t1, b.t2);
  }));
  EXPECT_TRUE(std::is_sorted(ee.begin(), ee.end(), [](const pgo::Contact::CIPC::EEPair &a, const pgo::Contact::CIPC::EEPair &b) {
    return std::tie(a.ea0, a.ea1, a.eb0, a.eb1) < std::tie(b.ea0, b.ea1, b.eb0, b.eb1);
  }));

  // same values as the default path up to summation order
  EXPECT_NEAR(det.energy, base.energy, 1e-10 * std::max(1.0, std::abs(base.energy)));
  EXPECT_LT(relativeError(det.gradient, base.gradient), 1e-10);
  EXPECT_LT(relativeError(sparseToDense(det.hessian), sparseToDense(base.hessian)), 1e-10);
  EXPECT_NEAR(det.allEnergy, det.energy, 1e-12 * std::max(1.0, std::abs(det.energy)));
  EXPECT_LT(relativeError(det.allGradient, det.gradient), 1e-12);
  EXPECT_LT(relativeError(sparseToDense(det.allHessian), sparseToDense(det.hessian)), 1e-12);

  // bit-reproducible across repeated broad phases and assemblies
  EXPECT_EQ(det.energy, again.energy);
  EXPECT_TRUE(bitwiseEqual(det.gradient, again.gradient));
  EXPECT_TRUE(bitwiseEqual(det.hessian, again.hessian));
  EXPECT_EQ(det.allEnergy, again.allEnergy);
  EXPECT_TRUE(bitwiseEqual(det.allGradient, again.allGradient));
  EXPECT_TRUE(bitwiseEqual(det.allHessian, again.allHessian));

  // the embedded energy forwards the switch
  const ES::SpMatD W = ES::MXd::Identity(x.size(), x.size()).sparseView();
  pgo::Contact::CIPC::EmbeddedSurfaceIPCPotentialEnergy embedded(V, F, W, params);
  EXPECT_FALSE(embedded.isDeterministic());
  embedded.setDeterministic(true);
  EXPECT_TRUE(embedded.isDeterministic());
}

TEST(SurfaceIPCCoreGTest, PreparedPairConsumersRequirePreparedState)
{
  SurfaceIPCCore core = makeConfiguredCore();
  EXPECT_THROW(core.computeEnergyWithPreparedPairs(), std::logic_error);
}
