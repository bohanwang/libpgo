#include <gtest/gtest.h>

#include "pgoLogging.h"
#include "ipc/core/surfaceIPCCore.h"
#include "ipc/external/obstacleSurface.h"
#include "ipc/profiling/surfaceIPCProfiling.h"
#include "scopedProfileSection.h"

#include "testCIPCHelpers.h"

#include <algorithm>
#include <array>
#include <stdexcept>
#include <string_view>
#include <vector>

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
using pgo::Profiling::ProfileStat;

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

const ProfileStat *findStat(const std::vector<ProfileStat> &stats, std::string_view name)
{
  const auto it = std::find_if(stats.begin(), stats.end(),
    [name](const ProfileStat &stat) { return stat.name == name; });
  return it == stats.end() ? nullptr : &(*it);
}

std::array<int, 2> edgeKey(int a, int b)
{
  if (b < a)
    std::swap(a, b);
  return { a, b };
}

std::array<int, 3> triKey(int a, int b, int c)
{
  std::array<int, 3> key = { a, b, c };
  std::sort(key.begin(), key.end());
  return key;
}

std::array<int, 4> selfPTKey(const pgo::Contact::CIPC::PTPair &pair)
{
  const auto tri = triKey(pair.t0, pair.t1, pair.t2);
  return { pair.p, tri[0], tri[1], tri[2] };
}

std::array<int, 4> selfEEKey(const pgo::Contact::CIPC::EEPair &pair)
{
  auto edgeA = edgeKey(pair.ea0, pair.ea1);
  auto edgeB = edgeKey(pair.eb0, pair.eb1);
  if (edgeB < edgeA)
    std::swap(edgeA, edgeB);
  return { edgeA[0], edgeA[1], edgeB[0], edgeB[1] };
}

std::array<int, 5> externalPTKey(const pgo::Contact::CIPC::ExternalPTPair &pair)
{
  const auto tri = triKey(pair.obsTri[0], pair.obsTri[1], pair.obsTri[2]);
  return { static_cast<int>(pair.obstacleSlot), pair.dynVertex, tri[0], tri[1], tri[2] };
}

std::array<int, 5> externalTPKey(const pgo::Contact::CIPC::ExternalTPPair &pair)
{
  const auto tri = triKey(pair.dynTri[0], pair.dynTri[1], pair.dynTri[2]);
  return { static_cast<int>(pair.obstacleSlot), tri[0], tri[1], tri[2], pair.obsVertex };
}

std::array<int, 5> externalEEKey(const pgo::Contact::CIPC::ExternalEEPair &pair)
{
  const auto dynEdge = edgeKey(pair.dynEdge[0], pair.dynEdge[1]);
  const auto obsEdge = edgeKey(pair.obsEdge[0], pair.obsEdge[1]);
  return { static_cast<int>(pair.obstacleSlot), dynEdge[0], dynEdge[1], obsEdge[0], obsEdge[1] };
}

template<class Pair, class KeyFn>
bool containsPair(const std::vector<Pair> &pairs, const Pair &target, KeyFn keyFn)
{
  const auto targetKey = keyFn(target);
  return std::any_of(pairs.begin(), pairs.end(),
    [&](const Pair &pair) { return keyFn(pair) == targetKey; });
}

void expectActiveSetSubset(
  const pgo::Contact::CIPC::SurfaceIPCActiveSet &exact,
  const pgo::Contact::CIPC::SurfaceIPCActiveSet &superset)
{
  for (const auto &pair : exact.selfPairs.ptPairs)
    EXPECT_TRUE(containsPair(superset.selfPairs.ptPairs, pair, selfPTKey));
  for (const auto &pair : exact.selfPairs.eePairs)
    EXPECT_TRUE(containsPair(superset.selfPairs.eePairs, pair, selfEEKey));
  for (const auto &pair : exact.externalPairs.ptPairs)
    EXPECT_TRUE(containsPair(superset.externalPairs.ptPairs, pair, externalPTKey));
  for (const auto &pair : exact.externalPairs.tpPairs)
    EXPECT_TRUE(containsPair(superset.externalPairs.tpPairs, pair, externalTPKey));
  for (const auto &pair : exact.externalPairs.eePairs)
    EXPECT_TRUE(containsPair(superset.externalPairs.eePairs, pair, externalEEKey));
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

TEST(SurfaceIPCCoreGTest, BuildActiveSetCapturesPositionsAndPairs)
{
  SurfaceIPCCore core = makeConfiguredCore();
  const auto [V, F] = makeTwoTriangleMesh();
  (void)F;
  const ES::VXd x = flattenPositions(V);

  const auto activeSet = core.buildActiveSet(x);

  EXPECT_TRUE(activeSet.positions.isApprox(x));
  ASSERT_FALSE(activeSet.selfPairs.ptPairs.empty());
  EXPECT_GT(activeSet.selfPairs.size(), 0u);
  EXPECT_EQ(activeSet.externalPairs.size(), 0u);
}

TEST(SurfaceIPCCoreGTest, ActiveSetConsumersMatchStatelessEnergyGradientHessian)
{
  SurfaceIPCCore core = makeConfiguredCore();
  const auto [V, F] = makeTwoTriangleMesh();
  (void)F;
  const ES::VXd x = flattenPositions(V);

  const double directEnergy = core.computeEnergy(x);
  ES::VXd directGradient = ES::VXd::Zero(x.size());
  core.computeGradient(x, directGradient);
  ES::SpMatD directHessian;
  core.computeHessian(x, directHessian);

  const auto activeSet = core.buildActiveSet(x);

  const double activeSetEnergy = core.computeEnergy(activeSet);
  ES::VXd activeSetGradient = ES::VXd::Zero(x.size());
  core.computeGradient(activeSet, activeSetGradient);
  ES::SpMatD activeSetHessian;
  core.computeHessian(activeSet, activeSetHessian);

  EXPECT_NEAR(activeSetEnergy, directEnergy, 1e-12);
  EXPECT_LT(relativeError(activeSetGradient, directGradient), 1e-12);
  EXPECT_LT(relativeError(sparseToDense(activeSetHessian), sparseToDense(directHessian)), 1e-12);
}

TEST(SurfaceIPCCoreGTest, ComputeAllBuildsActiveSetOnce)
{
  SurfaceIPCCore core = makeConfiguredCore();
  const auto [V, F] = makeTwoTriangleMesh();
  (void)F;
  const ES::VXd x = flattenPositions(V);

  pgo::Profiling::setProfilingEnabled(true);
  pgo::Profiling::resetProfileStatistics();

  double energy = 0.0;
  ES::VXd gradient = ES::VXd::Zero(x.size());
  ES::SpMatD hessian;
  core.computeAll(x, energy, gradient, hessian);

  const auto stats = pgo::Profiling::snapshotProfileStatistics();
  const ProfileStat *activeSetBuild = findStat(stats, pgo::Contact::SurfaceIPCProfileSections::kBuildActiveSet);

  pgo::Profiling::setProfilingEnabled(false);
  pgo::Profiling::resetProfileStatistics();

  ASSERT_NE(activeSetBuild, nullptr);
  EXPECT_EQ(activeSetBuild->callCount, 1u);
  EXPECT_GT(energy, 0.0);
  EXPECT_EQ(gradient.size(), x.size());
  EXPECT_EQ(hessian.rows(), x.size());
}

TEST(SurfaceIPCCoreGTest, LineSearchActiveSetSupersetContainsExactSampledAlphas)
{
  using pgo::Contact::CIPC::ObstacleSurface;

  ES::MXd dynV(4, 3);
  dynV << 0.0, 0.0, 0.32,
          1.0, 0.0, 0.32,
          0.0, 1.0, 0.32,
          1.0, 1.0, 0.32;
  ES::MXi dynF(2, 3);
  dynF << 0, 1, 2,
          1, 3, 2;

  ES::MXd obsV(4, 3);
  obsV << 0.0, 0.0, 0.0,
          1.0, 0.0, 0.0,
          0.0, 1.0, 0.0,
          1.0, 1.0, 0.0;
  ES::MXi obsF(2, 3);
  obsF << 0, 1, 2,
          1, 3, 2;

  const ES::VXd obsRest = flattenPositions(obsV);
  std::vector<ObstacleSurface> obstacles;
  obstacles.emplace_back(obsV, obsF,
    pgo::Contact::CIPC::makeLinearTrajectorySampler(obsRest, ES::V3d::Zero()));
  obstacles.front().update(0.0);

  SurfaceIPCCore::Parameters params;
  params.dhat = 0.1;
  params.dhat_external = 0.35;
  params.kappa = 1.0;
  params.eps_ee = 0.0;
  params.slackness = 0.9;

  SurfaceIPCCore core(params, std::move(obstacles));
  core.setMesh(dynV, dynF);

  const ES::VXd x = flattenPositions(dynV);
  ES::VXd dx = ES::VXd::Zero(x.size());
  for (int vi = 0; vi < dynV.rows(); ++vi)
    dx[3 * vi + 2] = -0.28;

  const auto lineSearchSuperset = core.buildLineSearchActiveSetSuperset(x, dx);
  bool sawExternalPairs = false;

  for (double alpha : { 0.0, 0.25, 0.5, 1.0 }) {
    const ES::VXd trial = x + alpha * dx;
    const auto exact = core.buildActiveSet(trial);
    auto supersetAtTrial = lineSearchSuperset;
    supersetAtTrial.positions = trial;

    expectActiveSetSubset(exact, supersetAtTrial);
    EXPECT_NEAR(core.computeEnergy(supersetAtTrial), core.computeEnergy(exact), 1e-12);
    sawExternalPairs = sawExternalPairs || exact.externalPairs.size() > 0;
  }

  EXPECT_TRUE(sawExternalPairs);
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
    obs.update(0.0);

  SurfaceIPCCore::Parameters params;
  params.dhat = 0.1;
  params.dhat_external = 1.0;
  params.kappa = 1.0;
  params.eps_ee = 0.0;
  params.slackness = 0.9;
  SurfaceIPCCore core(params, std::move(obstacles));
  core.setMesh(dynV, dynF);

  const ES::VXd x = flattenPositions(dynV);
  const auto activeSet = core.buildActiveSet(x);

  bool sawSlot0 = false;
  bool sawSlot1 = false;
  auto scanSlots = [&](auto &&vec) {
    for (const auto &pair : vec) {
      if (pair.obstacleSlot == 0) sawSlot0 = true;
      if (pair.obstacleSlot == 1) sawSlot1 = true;
    }
  };
  scanSlots(activeSet.externalPairs.ptPairs);
  scanSlots(activeSet.externalPairs.tpPairs);
  scanSlots(activeSet.externalPairs.eePairs);
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

// update(t) must refresh every pose-derived cache so that downstream broad
// phase / max-step can borrow them without rebuilding. Verify shapes, tight
// AABB containment, positive cellSize, and hash hit on a query box that
// covers a known triangle.
TEST(SurfaceIPCCoreGTest, ObstacleSurfaceUpdateRefreshesBroadPhaseCache)
{
  using pgo::Contact::CIPC::ObstacleSurface;
  using pgo::Contact::CIPC::SpatialHashGrid;

  ES::MXd V(4, 3);
  V << 0.0, 0.0, 0.0,
       1.0, 0.0, 0.0,
       0.0, 1.0, 0.0,
       1.0, 1.0, 0.0;
  ES::MXi F(2, 3);
  F << 0, 1, 2,
       1, 3, 2;

  ES::VXd rest(V.rows() * 3);
  for (int vi = 0; vi < V.rows(); ++vi)
    rest.segment<3>(3 * vi) = V.row(vi).transpose();

  ObstacleSurface obs(V, F,
    pgo::Contact::CIPC::makeLinearTrajectorySampler(rest, ES::V3d(0.0, 0.0, 1.0)));
  obs.update(0.5);  // moves obstacle by +z 0.5
  const auto &cache = obs.cache();

  // Vector sizes match primitive counts.
  ASSERT_EQ(static_cast<int>(cache.vertBoxes.size()), V.rows());
  ASSERT_EQ(static_cast<int>(cache.triBoxes.size()), F.rows());
  ASSERT_EQ(static_cast<int>(cache.edgeBoxes.size()), obs.contactEdges().rows());
  EXPECT_GT(cache.cellSize, 0.0);

  // Each tri box exactly contains its 3 sampled vertices (un-inflated).
  for (int fi = 0; fi < F.rows(); ++fi) {
    const auto &box = cache.triBoxes[fi];
    for (int j = 0; j < 3; ++j) {
      const ES::V3d v = obs.currentPositions().segment<3>(3 * F(fi, j));
      EXPECT_LE(box.lo.x(), v.x()); EXPECT_GE(box.hi.x(), v.x());
      EXPECT_LE(box.lo.y(), v.y()); EXPECT_GE(box.hi.y(), v.y());
      EXPECT_LE(box.lo.z(), v.z()); EXPECT_GE(box.hi.z(), v.z());
    }
  }

  // triHash query with a box covering triangle 0 must return 0.
  SpatialHashGrid::AABB queryBox;
  queryBox.init(obs.currentPositions().segment<3>(3 * F(0, 0)), 1e-3);
  queryBox.expand(obs.currentPositions().segment<3>(3 * F(0, 1)), 1e-3);
  queryBox.expand(obs.currentPositions().segment<3>(3 * F(0, 2)), 1e-3);
  std::vector<int> visited(F.rows(), 0);
  std::vector<int> candidates;
  cache.triHash.query(queryBox, -1, visited, 1, candidates);
  EXPECT_NE(std::find(candidates.begin(), candidates.end(), 0), candidates.end());

  // A second update(t) overwrites the cache cleanly (no stale entries).
  obs.update(1.0);
  const auto &cache2 = obs.cache();
  EXPECT_EQ(static_cast<int>(cache2.triBoxes.size()), F.rows());
  for (int fi = 0; fi < F.rows(); ++fi) {
    const auto &box = cache2.triBoxes[fi];
    const ES::V3d v0 = obs.currentPositions().segment<3>(3 * F(fi, 0));
    EXPECT_LE(box.lo.z(), v0.z());
    EXPECT_GE(box.hi.z(), v0.z());
  }
}

TEST(SurfaceIPCCoreGTest, ObstacleSurfaceKeepsOnlyFeatureEdgesForExternalEECache)
{
  using pgo::Contact::CIPC::ObstacleSurface;

  ES::MXd V(4, 3);
  V << 0.0, 0.0, 0.0,
       1.0, 0.0, 0.0,
       0.0, 1.0, 0.0,
       1.0, 1.0, 0.0;
  ES::MXi F(2, 3);
  F << 0, 1, 2,
       1, 3, 2;

  ES::VXd rest(V.rows() * 3);
  for (int vi = 0; vi < V.rows(); ++vi)
    rest.segment<3>(3 * vi) = V.row(vi).transpose();

  ObstacleSurface obs(V, F,
    pgo::Contact::CIPC::makeLinearTrajectorySampler(rest, ES::V3d::Zero()));
  obs.update(0.0);

  ASSERT_EQ(obs.uniqueEdges().rows(), 5);
  ASSERT_EQ(obs.contactEdges().rows(), 4);
  EXPECT_EQ(static_cast<int>(obs.cache().edgeBoxes.size()), obs.contactEdges().rows());
  EXPECT_EQ(static_cast<int>(obs.cache().edgeLengths.size()), obs.contactEdges().rows());

  for (int ei = 0; ei < obs.contactEdges().rows(); ++ei) {
    const int a = std::min(obs.contactEdges()(ei, 0), obs.contactEdges()(ei, 1));
    const int b = std::max(obs.contactEdges()(ei, 0), obs.contactEdges()(ei, 1));
    EXPECT_FALSE(a == 1 && b == 2);
  }
}

TEST(SurfaceIPCCoreGTest, ObstacleSurfaceKeepsSharpInteriorEdgesForExternalEECache)
{
  using pgo::Contact::CIPC::ObstacleSurface;

  ES::MXd V(4, 3);
  V << 0.0, 0.0, 0.0,
       1.0, 0.0, 0.0,
       0.0, 1.0, 0.0,
       1.0, 0.0, 1.0;
  ES::MXi F(2, 3);
  F << 0, 1, 2,
       1, 3, 2;

  ES::VXd rest(V.rows() * 3);
  for (int vi = 0; vi < V.rows(); ++vi)
    rest.segment<3>(3 * vi) = V.row(vi).transpose();

  ObstacleSurface obs(V, F,
    pgo::Contact::CIPC::makeLinearTrajectorySampler(rest, ES::V3d::Zero()));
  obs.update(0.0);

  ASSERT_EQ(obs.uniqueEdges().rows(), 5);
  ASSERT_EQ(obs.contactEdges().rows(), 5);

  bool sawSharedSharpEdge = false;
  for (int ei = 0; ei < obs.contactEdges().rows(); ++ei) {
    const int a = std::min(obs.contactEdges()(ei, 0), obs.contactEdges()(ei, 1));
    const int b = std::max(obs.contactEdges()(ei, 0), obs.contactEdges()(ei, 1));
    if (a == 1 && b == 2)
      sawSharedSharpEdge = true;
  }
  EXPECT_TRUE(sawSharedSharpEdge);
}
