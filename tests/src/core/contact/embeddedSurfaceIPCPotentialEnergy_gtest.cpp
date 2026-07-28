#include <gtest/gtest.h>

#include "embeddedSurfaceIPCPotentialEnergy.h"
#include "scopedProfileSection.h"
#include "ipc/core/surfaceIPCCore.h"
#include "ipc/profiling/surfaceIPCProfiling.h"
#include "testCIPCHelpers.h"

#include <algorithm>
#include <stdexcept>
#include <string_view>
#include <vector>

namespace
{
namespace ES = pgo::EigenSupport;
using pgo::Contact::CIPC::EmbeddedSurfaceIPCPotentialEnergy;
using pgo::Contact::CIPC::SurfaceIPCCore;
using pgo::Contact::CIPCTest::flattenPositions;
using pgo::Contact::CIPCTest::makeTwoTriangleMesh;
using pgo::Contact::CIPCTest::relativeError;
using pgo::Contact::CIPCTest::sparseToDense;
using pgo::Profiling::ProfileStat;

const ProfileStat *findStat(const std::vector<ProfileStat> &stats, std::string_view name)
{
  const auto it = std::find_if(stats.begin(), stats.end(),
    [name](const ProfileStat &stat) { return stat.name == name; });
  return it == stats.end() ? nullptr : &(*it);
}

SurfaceIPCCore::Parameters makeParams()
{
  SurfaceIPCCore::Parameters params;
  params.dhat = 0.1;
  params.kappa = 2.0;
  params.eps_ee = 0.0;
  params.slackness = 0.8;
  return params;
}

ES::SpMatD makeIdentityEmbedding(int n3)
{
  std::vector<ES::TripletD> triplets;
  triplets.reserve(n3);
  for (int i = 0; i < n3; ++i)
    triplets.emplace_back(i, i, 1.0);

  ES::SpMatD W(n3, n3);
  W.setFromTriplets(triplets.begin(), triplets.end());
  return W;
}
}  // namespace

TEST(EmbeddedSurfaceIPCPotentialEnergyGTest, IdentityEmbeddingMatchesSurfaceIPCCore)
{
  const auto [V, F] = makeTwoTriangleMesh();
  const ES::VXd rest = flattenPositions(V);
  ES::VXd u = ES::VXd::Zero(rest.size());
  for (int vi = 3; vi < 6; ++vi)
    u[3 * vi + 2] = 0.01;

  ES::VXd du = ES::VXd::Zero(rest.size());
  du[11] = -0.005;
  du[14] = -0.004;
  du[17] = -0.006;

  const auto params = makeParams();
  EmbeddedSurfaceIPCPotentialEnergy adapter(V, F, makeIdentityEmbedding(rest.size()), params);

  SurfaceIPCCore core(params);
  core.setMesh(V, F);
  const ES::VXd surfacePositions = rest + u;

  ES::VXd coreGradient(rest.size());
  core.computeGradient(surfacePositions, coreGradient);
  ES::SpMatD coreHessian;
  core.computeHessian(surfacePositions, coreHessian);

  ES::VXd adapterGradient(adapter.getNumDOFs());
  adapter.gradient(u, adapterGradient);
  ES::SpMatD adapterHessian;
  adapter.hessianDirect(u, adapterHessian);

  EXPECT_NEAR(adapter.func(u), core.computeEnergy(surfacePositions), 1e-10);
  EXPECT_LT(relativeError(adapterGradient, coreGradient), 1e-9);
  EXPECT_LT(relativeError(sparseToDense(adapterHessian), sparseToDense(coreHessian)), 1e-8);
  EXPECT_NEAR(adapter.computeMaxStepSize(u, du), core.computeMaxStepSize(surfacePositions, du), 1e-10);
}

TEST(EmbeddedSurfaceIPCPotentialEnergyGTest, SparseEmbeddingPullsBackGradientAndHessian)
{
  const auto [V, F] = makeTwoTriangleMesh();
  const ES::VXd surfaceRestPositions = flattenPositions(V);

  std::vector<ES::TripletD> triplets;
  triplets.emplace_back(0, 0, 1.0);
  triplets.emplace_back(1, 1, 1.0);
  triplets.emplace_back(2, 2, 1.0);
  triplets.emplace_back(3, 0, 0.25);
  triplets.emplace_back(3, 3, 0.75);
  triplets.emplace_back(4, 1, 0.25);
  triplets.emplace_back(4, 4, 0.75);
  triplets.emplace_back(5, 2, 0.25);
  triplets.emplace_back(5, 5, 0.75);
  triplets.emplace_back(6, 6, 1.0);
  triplets.emplace_back(7, 7, 1.0);
  triplets.emplace_back(8, 8, 1.0);
  triplets.emplace_back(9, 9, 1.0);
  triplets.emplace_back(10, 10, 1.0);
  triplets.emplace_back(11, 11, 1.0);
  triplets.emplace_back(12, 12, 1.0);
  triplets.emplace_back(13, 13, 1.0);
  triplets.emplace_back(14, 14, 1.0);
  triplets.emplace_back(15, 15, 1.0);
  triplets.emplace_back(16, 16, 1.0);
  triplets.emplace_back(17, 17, 1.0);

  ES::SpMatD W(surfaceRestPositions.size(), surfaceRestPositions.size());
  W.setFromTriplets(triplets.begin(), triplets.end());

  ES::VXd simulationDisplacements = ES::VXd::Zero(surfaceRestPositions.size());
  simulationDisplacements[11] = 0.01;
  simulationDisplacements[14] = 0.02;
  simulationDisplacements[17] = 0.015;

  EmbeddedSurfaceIPCPotentialEnergy adapter(V, F, W, makeParams());

  SurfaceIPCCore core(makeParams());
  core.setMesh(V, F);

  const ES::VXd surfacePositions = surfaceRestPositions + W * simulationDisplacements;

  ES::VXd surfaceGradient(surfacePositions.size());
  core.computeGradient(surfacePositions, surfaceGradient);
  ES::SpMatD surfaceHessian;
  core.computeHessian(surfacePositions, surfaceHessian);

  ES::VXd simulationGradient(adapter.getNumDOFs());
  adapter.gradient(simulationDisplacements, simulationGradient);
  ES::SpMatD simulationHessian;
  adapter.hessianDirect(simulationDisplacements, simulationHessian);

  EXPECT_LT(relativeError(simulationGradient, W.transpose() * surfaceGradient), 1e-12);
  EXPECT_LT(relativeError(sparseToDense(simulationHessian), sparseToDense(W.transpose() * surfaceHessian * W)), 1e-12);
}

TEST(EmbeddedSurfaceIPCPotentialEnergyGTest, ProfilingRecordsAdapterSections)
{
  const auto [V, F] = makeTwoTriangleMesh();
  const ES::VXd rest = flattenPositions(V);
  ES::VXd u = ES::VXd::Zero(rest.size());
  for (int vi = 3; vi < 6; ++vi)
    u[3 * vi + 2] = 0.01;

  EmbeddedSurfaceIPCPotentialEnergy adapter(V, F, makeIdentityEmbedding(rest.size()), makeParams());

  pgo::Profiling::setProfilingEnabled(true);
  pgo::Profiling::resetProfileStatistics();

  ES::VXd gradient(rest.size());
  adapter.gradient(u, gradient);
  ES::SpMatD hessian;
  adapter.hessianDirect(u, hessian);
  (void)adapter.func(u);

  const auto stats = pgo::Profiling::snapshotProfileStatistics();
  EXPECT_NE(findStat(stats, "contact.adapter.func"), nullptr);
  EXPECT_NE(findStat(stats, "contact.adapter.gradient"), nullptr);
  EXPECT_NE(findStat(stats, "contact.adapter.hessian_direct"), nullptr);
  EXPECT_NE(findStat(stats, "contact.adapter.map_to_surface"), nullptr);
  EXPECT_NE(findStat(stats, "contact.adapter.pullback_gradient"), nullptr);
  EXPECT_NE(findStat(stats, "contact.adapter.pullback_hessian"), nullptr);

  pgo::Profiling::setProfilingEnabled(false);
  pgo::Profiling::resetProfileStatistics();
}

TEST(EmbeddedSurfaceIPCPotentialEnergyGTest, ReusesPreparedPairsAcrossEnergyGradientHessianForSameState)
{
  const auto [V, F] = makeTwoTriangleMesh();
  const ES::VXd rest = flattenPositions(V);
  ES::VXd u = ES::VXd::Zero(rest.size());
  for (int vi = 3; vi < 6; ++vi)
    u[3 * vi + 2] = 0.01;

  EmbeddedSurfaceIPCPotentialEnergy adapter(V, F, makeIdentityEmbedding(rest.size()), makeParams());

  pgo::Profiling::setProfilingEnabled(true);
  pgo::Profiling::resetProfileStatistics();

  const double energy0 = adapter.func(u);
  ES::VXd gradient0 = ES::VXd::Zero(adapter.getNumDOFs());
  adapter.gradient(u, gradient0);
  ES::SpMatD hessian0;
  adapter.hessianDirect(u, hessian0);

  const double energy1 = adapter.func(u);
  ES::VXd gradient1 = ES::VXd::Zero(adapter.getNumDOFs());
  adapter.gradient(u, gradient1);
  ES::SpMatD hessian1;
  adapter.hessianDirect(u, hessian1);

  const auto stats = pgo::Profiling::snapshotProfileStatistics();
  const ProfileStat *pairBuild = findStat(stats, pgo::Contact::SurfaceIPCProfileSections::kPairBuildStatic);

  pgo::Profiling::setProfilingEnabled(false);
  pgo::Profiling::resetProfileStatistics();

  ASSERT_NE(pairBuild, nullptr);
  EXPECT_EQ(pairBuild->callCount, 1u);
  EXPECT_NEAR(energy1, energy0, 1e-12);
  EXPECT_LT(relativeError(gradient1, gradient0), 1e-12);
  EXPECT_LT(relativeError(sparseToDense(hessian1), sparseToDense(hessian0)), 1e-12);
}

TEST(EmbeddedSurfaceIPCPotentialEnergyGTest, InvalidEmbeddingRowsThrow)
{
  const auto [V, F] = makeTwoTriangleMesh();
  ES::SpMatD W(3 * V.rows() - 1, 3 * V.rows());

  EXPECT_THROW(
    EmbeddedSurfaceIPCPotentialEnergy(V, F, W, makeParams()),
    std::invalid_argument);
}
