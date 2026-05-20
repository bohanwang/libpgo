#include <gtest/gtest.h>

#include "ipc/broadPhase/surfaceIPCBroadPhase.h"
#include "ipc/core/surfaceIPCSelfBarrierAssembler.h"
#include "ipc/core/surfaceIPCExternalBarrierAssembler.h"
#include "ipc/core/surfaceIPCBarrierKernels.h"
#include "ipc/external/obstacleSurface.h"
#include "ipc/profiling/surfaceIPCProfiling.h"
#include "ipc/topology/surfaceIPCTopology.h"
#include "ipc/core/surfaceIPCCore.h"
#include "scopedProfileSection.h"

#include "testCIPCHelpers.h"

#include <algorithm>
#include <string_view>
#include <vector>

namespace
{
namespace ES = pgo::EigenSupport;
using pgo::Contact::CIPC::ExternalPairSet;
using pgo::Contact::CIPC::EEPair;
using pgo::Contact::CIPC::ObstacleSurface;
using pgo::Contact::CIPC::PTPair;
using pgo::Contact::CIPC::SelfPairSet;
using pgo::Contact::CIPC::SurfaceIPCCore;
using pgo::Contact::CIPC::SurfaceIPCTopology;
using pgo::Contact::CIPCTest::flattenPositions;
using pgo::Contact::CIPCTest::makeTwoTriangleMesh;
using pgo::Contact::CIPCTest::relativeError;
using pgo::Contact::CIPCTest::sparseToDense;
using pgo::Profiling::ProfileCounterStat;
using pgo::Profiling::ProfileStat;
namespace kernels = pgo::Contact::CIPC::barrier_kernels;

const ProfileStat *findStat(const std::vector<ProfileStat> &stats, std::string_view name)
{
  const auto it = std::find_if(stats.begin(), stats.end(),
    [name](const ProfileStat &stat) { return stat.name == name; });
  return it == stats.end() ? nullptr : &(*it);
}

const ProfileCounterStat *findCounterStat(const std::vector<ProfileCounterStat> &stats, std::string_view name)
{
  const auto it = std::find_if(stats.begin(), stats.end(),
    [name](const ProfileCounterStat &stat) { return stat.name == name; });
  return it == stats.end() ? nullptr : &(*it);
}

class SurfaceIPCBarrierAssemblerProfilingGTest : public ::testing::Test
{
protected:
  void SetUp() override
  {
    pgo::Profiling::setProfilingEnabled(false);
    pgo::Profiling::resetProfileStatistics();
  }

  void TearDown() override
  {
    pgo::Profiling::setProfilingEnabled(false);
    pgo::Profiling::resetProfileStatistics();
  }
};
}  // namespace

TEST(SurfaceIPCBarrierAssemblerGTest, HelperMatchesSurfaceIPCCoreAssembly)
{
  const auto [V, F] = makeTwoTriangleMesh();
  const ES::VXd x = flattenPositions(V);

  SurfaceIPCTopology topology;
  topology.setMesh(V, F);

  constexpr double dhat = 0.1;
  constexpr double kappa = 1.0;
  constexpr double epsEE = 0.0;
  SelfPairSet pairs;
  buildSelfPairs(topology, x, dhat, pairs);

  const double helperEnergy = computeSelfEnergy(x, pairs, topology.numVerts, dhat, kappa, epsEE);
  ES::VXd helperGradient(x.size());
  computeSelfGradient(x, pairs, topology.numVerts, dhat, kappa, epsEE, helperGradient);
  ES::SpMatD helperHessian;
  computeSelfHessian(x, pairs, topology.numVerts, dhat, kappa, epsEE, helperHessian);

  SurfaceIPCCore core;
  SurfaceIPCCore::Parameters params;
  params.dhat = dhat;
  params.kappa = kappa;
  params.eps_ee = epsEE;
  params.slackness = 0.9;
  core.setParameters(params);
  core.setMesh(V, F);

  ES::VXd coreGradient(x.size());
  core.computeGradient(x, coreGradient);
  ES::SpMatD coreHessian;
  core.computeHessian(x, coreHessian);

  EXPECT_NEAR(helperEnergy, core.computeEnergy(x), 1e-12);
  EXPECT_LT(relativeError(helperGradient, coreGradient), 1e-12);
  EXPECT_LT(relativeError(sparseToDense(helperHessian), sparseToDense(coreHessian)), 1e-12);
}

TEST(SurfaceIPCBarrierAssemblerGTest, KernelPTActiveAndInactive)
{
  // Vertex at origin, triangle 0.05 above — close enough to activate with dhat=0.1
  ES::V3d p(0.0, 0.0, 0.0);
  ES::V3d t0(0.0, 0.0, 0.05);
  ES::V3d t1(1.0, 0.0, 0.05);
  ES::V3d t2(0.0, 1.0, 0.05);

  const double dhat2 = 0.01;   // 0.1²
  const double kappa = 1.0;
  const double weight = 1.0;

  // Energy-only: should be active with positive energy
  auto kE = kernels::pointTriangle(p, t0, t1, t2, weight, dhat2, kappa, false, false);
  ASSERT_TRUE(kE.active);
  EXPECT_GT(kE.energy, 0.0);

  // Energy + gradient: gradient w.r.t. p should point toward triangle (+z)
  auto kG = kernels::pointTriangle(p, t0, t1, t2, weight, dhat2, kappa, true, false);
  ASSERT_TRUE(kG.active);
  EXPECT_NEAR(kG.energy, kE.energy, 1e-12);
  EXPECT_GT(kG.gradient[2], 0.0);  // p.z gradient positive

  // Energy + hessian (PSD-projected): p-p diagonal should be positive
  auto kH = kernels::pointTriangle(p, t0, t1, t2, weight, dhat2, kappa, false, true);
  ASSERT_TRUE(kH.active);
  EXPECT_GT(kH.hessian(2, 2), 0.0);

  // Too far: triangle 1.0 above — should be inactive
  ES::V3d far0(0.0, 0.0, 1.0);
  ES::V3d far1(1.0, 0.0, 1.0);
  ES::V3d far2(0.0, 1.0, 1.0);
  auto kFar = kernels::pointTriangle(p, far0, far1, far2, weight, dhat2, kappa, false, false);
  EXPECT_FALSE(kFar.active);

  // Coincident (d2=0) — should be inactive
  auto kZero = kernels::pointTriangle(p, p, t1, t2, weight, dhat2, kappa, false, false);
  EXPECT_FALSE(kZero.active);
}

TEST(SurfaceIPCBarrierAssemblerGTest, KernelEEActiveAndInactive)
{
  // Two edges 0.02 apart — close enough to activate with dhat=0.1
  ES::V3d ea0(0.0, -0.01, 0.0);
  ES::V3d ea1(1.0, -0.01, 0.0);
  ES::V3d eb0(0.0,  0.01, 0.0);
  ES::V3d eb1(1.0,  0.01, 0.0);

  const double dhat2 = 0.01;
  const double kappa = 1.0;
  const double weight = 1.0;

  // Without mollifier
  auto k = kernels::edgeEdge(ea0, ea1, eb0, eb1, weight, dhat2, kappa, 0.0, true, true);
  ASSERT_TRUE(k.active);
  EXPECT_GT(k.energy, 0.0);
  // Gradient should push edges apart (opposite y directions)
  EXPECT_GT(std::abs(k.gradient[1]), 0.0);
  EXPECT_GT(std::abs(k.gradient[7]), 0.0);
  EXPECT_LT(k.gradient[1] * k.gradient[7], 0.0);  // opposite signs
  // Hessian diagonal should be positive (PSD)
  EXPECT_GT(k.hessian(1, 1), 0.0);

  // With mollifier: energy should differ
  auto kM = kernels::edgeEdge(ea0, ea1, eb0, eb1, weight, dhat2, kappa, 1e-3, true, true);
  ASSERT_TRUE(kM.active);
  EXPECT_NE(kM.energy, k.energy);

  // Too far apart — inactive
  ES::V3d far0(0.0, -1.0, 0.0);
  ES::V3d far1(1.0, -1.0, 0.0);
  ES::V3d far2(0.0,  1.0, 0.0);
  ES::V3d far3(1.0,  1.0, 0.0);
  auto kFar = kernels::edgeEdge(far0, far1, far2, far3, weight, dhat2, kappa, 0.0, false, false);
  EXPECT_FALSE(kFar.active);
}

TEST(SurfaceIPCBarrierAssemblerGTest, ExternalDynamicPointKernelMatchesGenericSubBlock)
{
  const ES::V3d p(0.10, 0.20, 0.0);
  const ES::V3d t0(0.0, 0.0, 0.04);
  const ES::V3d t1(1.0, 0.0, 0.04);
  const ES::V3d t2(0.0, 1.0, 0.04);

  const double dhat2 = 0.01;
  const double kappa = 2.0;
  const double weight = 0.75;

  const auto generic = kernels::pointTriangle(p, t0, t1, t2, weight, dhat2, kappa, true, true);
  const auto external = kernels::pointStaticTriangle(p, t0, t1, t2, weight, dhat2, kappa, true, true);

  ASSERT_TRUE(generic.active);
  ASSERT_EQ(external.active, generic.active);
  EXPECT_NEAR(external.energy, generic.energy, 1e-12);
  EXPECT_LT((external.gradient - generic.gradient.head<3>()).norm(), 1e-12);
  EXPECT_LT((external.hessian - generic.hessian.block<3, 3>(0, 0)).norm(), 1e-12);
}

TEST(SurfaceIPCBarrierAssemblerGTest, ExternalDynamicTriangleKernelMatchesGenericSubBlock)
{
  const ES::V3d p(0.10, 0.20, 0.0);
  const ES::V3d t0(0.0, 0.0, 0.04);
  const ES::V3d t1(1.0, 0.0, 0.04);
  const ES::V3d t2(0.0, 1.0, 0.04);

  const double dhat2 = 0.01;
  const double kappa = 2.0;
  const double weight = 0.75;

  const auto generic = kernels::pointTriangle(p, t0, t1, t2, weight, dhat2, kappa, true, true);
  const auto external = kernels::staticPointTriangle(p, t0, t1, t2, weight, dhat2, kappa, true, true);

  ASSERT_TRUE(generic.active);
  ASSERT_EQ(external.active, generic.active);
  EXPECT_NEAR(external.energy, generic.energy, 1e-12);
  EXPECT_LT((external.gradient - generic.gradient.segment<9>(3)).norm(), 1e-12);
  EXPECT_LT((external.hessian - generic.hessian.block<9, 9>(3, 3)).norm(), 1e-12);
}

TEST(SurfaceIPCBarrierAssemblerGTest, ExternalDynamicEdgeKernelMatchesGenericSubBlock)
{
  const ES::V3d ea0(0.0, -0.01, 0.0);
  const ES::V3d ea1(1.0, -0.01, 0.0);
  const ES::V3d eb0(0.0,  0.01, 0.0);
  const ES::V3d eb1(1.0,  0.01, 0.0);

  const double dhat2 = 0.01;
  const double kappa = 1.5;
  const double weight = 0.8;
  const double epsEE = 1e-3;

  const auto generic = kernels::edgeEdge(ea0, ea1, eb0, eb1, weight, dhat2, kappa, epsEE, true, true);
  const auto external = kernels::edgeStaticEdge(ea0, ea1, eb0, eb1, weight, dhat2, kappa, epsEE, true, true);

  ASSERT_TRUE(generic.active);
  ASSERT_EQ(external.active, generic.active);
  EXPECT_NEAR(external.energy, generic.energy, 1e-12);
  EXPECT_LT((external.gradient - generic.gradient.head<6>()).norm(), 1e-12);
  EXPECT_LT((external.hessian - generic.hessian.block<6, 6>(0, 0)).norm(), 1e-12);
}

TEST(SurfaceIPCBarrierAssemblerGTest, SinglePairGradientScatterIsCorrect)
{
  // Create a minimal mesh with 2 vertices, build one PT pair manually
  // Then verify gradient is scattered to the correct DOFs
  ES::MXd V(2, 3);
  V << 0.0, 0.0, 0.0,
       0.0, 0.0, 0.1;
  ES::MXi F(0, 3);  // no triangles — we build the pair manually

  SurfaceIPCTopology topology;
  topology.setMesh(V, F);

  SelfPairSet pairs;
  PTPair ptPair = { 0, 1, 1, 1, 1.0 };  // vertex 0 vs (degenerate) "triangle" using vertex 1×3
  pairs.ptPairs.push_back(ptPair);

  const ES::VXd x = flattenPositions(V);
  const double dhat = 1.0;  // large enough to activate
  const double kappa = 1.0;
  const double epsEE = 0.0;

  ES::VXd grad = ES::VXd::Zero(x.size());
  computeSelfGradient(x, pairs, topology.numVerts, dhat, kappa, epsEE, grad);

  // The gradient should be non-zero (barrier active)
  double gradNorm = grad.norm();
  EXPECT_GT(gradNorm, 0.0);
}

TEST(SurfaceIPCBarrierAssemblerGTest, SinglePairHessianScatterProducesSymmetricMatrix)
{
  ES::MXd V(2, 3);
  V << 0.0, 0.0, 0.0,
       0.0, 0.0, 0.1;
  ES::MXi F(0, 3);

  SurfaceIPCTopology topology;
  topology.setMesh(V, F);

  SelfPairSet pairs;
  PTPair ptPair = { 0, 1, 1, 1, 1.0 };
  pairs.ptPairs.push_back(ptPair);

  const ES::VXd x = flattenPositions(V);
  const double dhat = 1.0;
  const double kappa = 1.0;
  const double epsEE = 0.0;

  ES::SpMatD hess;
  computeSelfHessian(x, pairs, topology.numVerts, dhat, kappa, epsEE, hess);

  // Hessian should be symmetric
  ES::MXd denseH = sparseToDense(hess);
  EXPECT_LT((denseH - denseH.transpose()).norm(), 1e-12);
  EXPECT_GT(hess.nonZeros(), 0);
}

TEST_F(SurfaceIPCBarrierAssemblerProfilingGTest, ComputeSelfAllRecordsBarrierBreakdownProfiling)
{
  ES::MXd V(8, 3);
  V << 0.0,  0.0,  0.0,
       0.0,  0.0,  0.05,
       1.0,  0.0,  0.05,
       0.0,  1.0,  0.05,
       0.0, -0.01, 0.0,
       1.0, -0.01, 0.0,
       0.0,  0.01, 0.0,
       1.0,  0.01, 0.0;
  const ES::VXd x = flattenPositions(V);

  SelfPairSet pairs;
  pairs.ptPairs.push_back({ 0, 1, 2, 3, 1.0 });
  pairs.eePairs.push_back({ 4, 5, 6, 7, 1.0 });

  pgo::Profiling::setProfilingEnabled(true);

  double energy = 0.0;
  ES::VXd grad;
  ES::SpMatD hess;
  computeSelfAll(x, pairs, static_cast<int>(V.rows()), 0.1, 1.0, 0.0, energy, grad, hess);

  const auto stats = pgo::Profiling::snapshotProfileStatistics();
  const auto counters = pgo::Profiling::snapshotProfileCounterStatistics();
  const ProfileStat *selfCombined = findStat(stats, pgo::Contact::SurfaceIPCProfileSections::kActiveSetSelfCombined);
  const ProfileStat *selfPT = findStat(stats, pgo::Contact::SurfaceIPCProfileSections::kActiveSetSelfPTCombined);
  const ProfileStat *selfEE = findStat(stats, pgo::Contact::SurfaceIPCProfileSections::kActiveSetSelfEECombined);
  const ProfileCounterStat *selfPTPairs = findCounterStat(counters, pgo::Contact::SurfaceIPCProfileSections::kActiveSetSelfPTPairCount);
  const ProfileCounterStat *selfEEPairs = findCounterStat(counters, pgo::Contact::SurfaceIPCProfileSections::kActiveSetSelfEEPairCount);

  ASSERT_NE(selfCombined, nullptr);
  ASSERT_NE(selfPT, nullptr);
  ASSERT_NE(selfEE, nullptr);
  ASSERT_NE(selfPTPairs, nullptr);
  ASSERT_NE(selfEEPairs, nullptr);
  EXPECT_EQ(selfCombined->callCount, 1u);
  EXPECT_EQ(selfPT->callCount, 1u);
  EXPECT_EQ(selfEE->callCount, 1u);
  EXPECT_EQ(selfPTPairs->sampleCount, 1u);
  EXPECT_EQ(selfPTPairs->total, pairs.ptPairs.size());
  EXPECT_EQ(selfPTPairs->max, pairs.ptPairs.size());
  EXPECT_EQ(selfEEPairs->sampleCount, 1u);
  EXPECT_EQ(selfEEPairs->total, pairs.eePairs.size());
  EXPECT_EQ(selfEEPairs->max, pairs.eePairs.size());
}

TEST_F(SurfaceIPCBarrierAssemblerProfilingGTest, ComputeExternalAllRecordsBarrierBreakdownProfiling)
{
  ES::MXd V(4, 3);
  V << 0.0, 0.0, 0.0,
       1.0, 0.0, 0.0,
       0.0, 0.0, 1.0,
       1.0, 0.0, 1.0;
  const ES::VXd x = flattenPositions(V);

  ES::MXd obsV(4, 3);
  obsV << 0.0, 0.05, 0.0,
          1.0, 0.05, 0.0,
          0.0, 0.05, 1.0,
          1.0, 0.05, 1.0;
  ES::MXi obsF(2, 3);
  obsF << 0, 1, 2,
          1, 3, 2;
  const ES::VXd obsRest = flattenPositions(obsV);
  ObstacleSurface obs(obsV, obsF,
    pgo::Contact::CIPC::makeLinearTrajectorySampler(obsRest, ES::V3d::Zero()));
  obs.setObjectId(0);
  obs.update(0.0);
  std::vector<ObstacleSurface> obstacles;
  obstacles.emplace_back(std::move(obs));

  ExternalPairSet pairs;
  pairs.ptPairs.push_back({ 0, 0, { 0, 1, 2 }, 1.0 });
  pairs.tpPairs.push_back({ 0, { 0, 1, 2 }, 0, 1.0 });
  pairs.eePairs.push_back({ 0, { 0, 1 }, { 0, 1 }, 1.0 });

  pgo::Profiling::setProfilingEnabled(true);

  double energy = 0.0;
  ES::VXd grad;
  ES::SpMatD hess;
  computeExternalAll(x, obstacles, pairs, static_cast<int>(V.rows()), 0.1, 1.0, 0.0, energy, grad, hess);

  const auto stats = pgo::Profiling::snapshotProfileStatistics();
  const auto counters = pgo::Profiling::snapshotProfileCounterStatistics();
  const ProfileStat *externalCombined = findStat(stats, pgo::Contact::SurfaceIPCProfileSections::kActiveSetExternalCombined);
  const ProfileStat *externalPT = findStat(stats, pgo::Contact::SurfaceIPCProfileSections::kActiveSetExternalPTCombined);
  const ProfileStat *externalTP = findStat(stats, pgo::Contact::SurfaceIPCProfileSections::kActiveSetExternalTPCombined);
  const ProfileStat *externalEE = findStat(stats, pgo::Contact::SurfaceIPCProfileSections::kActiveSetExternalEECombined);
  const ProfileCounterStat *externalPTPairs = findCounterStat(counters, pgo::Contact::SurfaceIPCProfileSections::kActiveSetExternalPTPairCount);
  const ProfileCounterStat *externalTPPairs = findCounterStat(counters, pgo::Contact::SurfaceIPCProfileSections::kActiveSetExternalTPPairCount);
  const ProfileCounterStat *externalEEPairs = findCounterStat(counters, pgo::Contact::SurfaceIPCProfileSections::kActiveSetExternalEEPairCount);

  ASSERT_NE(externalCombined, nullptr);
  ASSERT_NE(externalPT, nullptr);
  ASSERT_NE(externalTP, nullptr);
  ASSERT_NE(externalEE, nullptr);
  ASSERT_NE(externalPTPairs, nullptr);
  ASSERT_NE(externalTPPairs, nullptr);
  ASSERT_NE(externalEEPairs, nullptr);
  EXPECT_EQ(externalCombined->callCount, 1u);
  EXPECT_EQ(externalPT->callCount, 1u);
  EXPECT_EQ(externalTP->callCount, 1u);
  EXPECT_EQ(externalEE->callCount, 1u);
  EXPECT_EQ(externalPTPairs->sampleCount, 1u);
  EXPECT_EQ(externalPTPairs->total, pairs.ptPairs.size());
  EXPECT_EQ(externalPTPairs->max, pairs.ptPairs.size());
  EXPECT_EQ(externalTPPairs->sampleCount, 1u);
  EXPECT_EQ(externalTPPairs->total, pairs.tpPairs.size());
  EXPECT_EQ(externalTPPairs->max, pairs.tpPairs.size());
  EXPECT_EQ(externalEEPairs->sampleCount, 1u);
  EXPECT_EQ(externalEEPairs->total, pairs.eePairs.size());
  EXPECT_EQ(externalEEPairs->max, pairs.eePairs.size());
}
