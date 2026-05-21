#include <gtest/gtest.h>

#include "ipc/CIPC.h"
#include "ipc/core/surfaceIPCCore.h"
#include "scopedProfileSection.h"
#include "ipc/profiling/surfaceIPCProfiling.h"

#include "testCIPCHelpers.h"

#include <algorithm>
#include <string_view>
#include <vector>

namespace
{
namespace ES = pgo::EigenSupport;
using pgo::Contact::CIPC::CIPCPotentialEnergy;
using pgo::Contact::CIPC::SurfaceIPCCore;
using pgo::Contact::CIPCTest::computeFloorEnergy;
using pgo::Contact::CIPCTest::computeFloorGradient;
using pgo::Contact::CIPCTest::computeFloorHessian;
using pgo::Contact::CIPCTest::flattenPositions;
using pgo::Contact::CIPCTest::makeTwoTriangleMesh;
using pgo::Contact::CIPCTest::relativeError;
using pgo::Contact::CIPCTest::sparseToDense;
using pgo::NonlinearOptimization::PotentialEnergy;
using pgo::Profiling::ProfileStat;

const ProfileStat *findStat(const std::vector<ProfileStat> &stats, std::string_view name)
{
  const auto it = std::find_if(stats.begin(), stats.end(),
    [name](const ProfileStat &stat) { return stat.name == name; });
  return it == stats.end() ? nullptr : &(*it);
}

SurfaceIPCCore makeReferenceCore(const ES::MXd &V, const ES::MXi &F,
  double dhat, double kappa, double eps_ee, double slackness)
{
  SurfaceIPCCore core;
  SurfaceIPCCore::Parameters params;
  params.dhat = dhat;
  params.kappa = kappa;
  params.eps_ee = eps_ee;
  params.slackness = slackness;
  core.setParameters(params);
  core.setMesh(V, F);
  return core;
}
}  // namespace

TEST(CIPCPotentialEnergyGTest, AbsoluteInputWrapperMatchesCore)
{
  const auto [V, F] = makeTwoTriangleMesh();
  const ES::VXd x = flattenPositions(V);

  CIPCPotentialEnergy wrapper(0.1, 1.0, false);
  wrapper.setMesh(V, F);

  SurfaceIPCCore core = makeReferenceCore(V, F, wrapper.dhat, wrapper.kappa, wrapper.eps_ee, wrapper.slackness);

  ES::VXd wrapperGrad(x.size());
  wrapper.gradient(x, wrapperGrad);
  ES::SpMatD wrapperH;
  wrapper.hessianDirect(x, wrapperH);

  ES::VXd coreGrad(x.size());
  core.computeGradient(x, coreGrad);
  ES::SpMatD coreH;
  core.computeHessian(x, coreH);

  EXPECT_NEAR(wrapper.func(x), core.computeEnergy(x), 1e-10);
  EXPECT_LT(relativeError(wrapperGrad, coreGrad), 1e-12);
  EXPECT_LT(relativeError(sparseToDense(wrapperH), sparseToDense(coreH)), 1e-12);
}

TEST(CIPCPotentialEnergyGTest, DisplacementInputWrapperMatchesCore)
{
  const auto [V, F] = makeTwoTriangleMesh();
  ES::VXd x = ES::VXd::Zero(V.rows() * 3);
  for (int vi = 3; vi < 6; ++vi)
    x[3 * vi + 2] = 0.01;

  CIPCPotentialEnergy wrapper(0.1, 1.0, true);
  wrapper.setMesh(V, F);

  const ES::VXd xSurf = flattenPositions(V) + x;
  SurfaceIPCCore core = makeReferenceCore(V, F, wrapper.dhat, wrapper.kappa, wrapper.eps_ee, wrapper.slackness);

  ES::VXd wrapperGrad(x.size());
  wrapper.gradient(x, wrapperGrad);
  ES::SpMatD wrapperH;
  wrapper.hessianDirect(x, wrapperH);

  ES::VXd coreGrad(x.size());
  core.computeGradient(xSurf, coreGrad);
  ES::SpMatD coreH;
  core.computeHessian(xSurf, coreH);

  EXPECT_NEAR(wrapper.func(x), core.computeEnergy(xSurf), 1e-10);
  EXPECT_LT(relativeError(wrapperGrad, coreGrad), 1e-12);
  EXPECT_LT(relativeError(sparseToDense(wrapperH), sparseToDense(coreH)), 1e-12);
}

TEST(CIPCPotentialEnergyGTest, WrapperParameterChangesSyncIntoCore)
{
  const auto [V, F] = makeTwoTriangleMesh();
  const ES::VXd x = flattenPositions(V);

  CIPCPotentialEnergy wrapper(0.1, 1.0, false);
  wrapper.setMesh(V, F);
  wrapper.dhat = 0.12;
  wrapper.kappa = 2.0;
  wrapper.eps_ee = 0.0;
  wrapper.slackness = 0.8;

  SurfaceIPCCore core = makeReferenceCore(V, F, wrapper.dhat, wrapper.kappa, wrapper.eps_ee, wrapper.slackness);
  EXPECT_NEAR(wrapper.func(x), core.computeEnergy(x), 1e-10);
}

TEST(CIPCPotentialEnergyGTest, BarrierAndFloorActiveMatchesCorePlusFloorContribution)
{
  const auto [V, F] = makeTwoTriangleMesh();
  const ES::VXd x = flattenPositions(V);

  CIPCPotentialEnergy wrapper(0.1, 1.0, false, 0.0, true, 0.08, 3.0);
  wrapper.setMesh(V, F);

  SurfaceIPCCore core = makeReferenceCore(V, F, wrapper.dhat, wrapper.kappa, wrapper.eps_ee, wrapper.slackness);

  ES::VXd wrapperGrad(x.size());
  wrapper.gradient(x, wrapperGrad);
  ES::SpMatD wrapperH;
  wrapper.hessianDirect(x, wrapperH);

  ES::VXd coreGrad(x.size());
  core.computeGradient(x, coreGrad);
  ES::SpMatD coreH;
  core.computeHessian(x, coreH);

  const double floorE = computeFloorEnergy(x, wrapper.floorHeight, wrapper.floorKappa);
  const ES::VXd floorG = computeFloorGradient(x, wrapper.floorHeight, wrapper.floorKappa);
  const ES::MXd floorH = computeFloorHessian(x, wrapper.floorHeight, wrapper.floorKappa);

  EXPECT_NEAR(wrapper.func(x), core.computeEnergy(x) + floorE, 1e-10);
  EXPECT_LT(relativeError(wrapperGrad, coreGrad + floorG), 1e-12);
  EXPECT_LT(relativeError(sparseToDense(wrapperH), sparseToDense(coreH) + floorH), 1e-12);
}

TEST(CIPCPotentialEnergyGTest, BaseGradientHessianUsesHessianDirectDefault)
{
  const auto [V, F] = makeTwoTriangleMesh();
  const ES::VXd x = flattenPositions(V);

  CIPCPotentialEnergy wrapper(0.1, 1.0, false);
  wrapper.setMesh(V, F);

  const PotentialEnergy &base = wrapper;
  ES::VXd combinedGrad = ES::VXd::Zero(x.size());
  ES::SpMatD combinedH;
  EXPECT_NO_THROW(base.gradient_hessian(x, combinedGrad, combinedH));

  ES::VXd refGrad = ES::VXd::Zero(x.size());
  wrapper.gradient(x, refGrad);
  ES::SpMatD refH;
  wrapper.hessianDirect(x, refH);

  EXPECT_LT(relativeError(combinedGrad, refGrad), 1e-12);
  EXPECT_LT(relativeError(sparseToDense(combinedH), sparseToDense(refH)), 1e-12);
}

TEST(CIPCPotentialEnergyGTest, SeparateEvaluationsBuildIndependentActiveSetsForSameState)
{
  const auto [V, F] = makeTwoTriangleMesh();
  const ES::VXd x = flattenPositions(V);

  CIPCPotentialEnergy wrapper(0.1, 1.0, false);
  wrapper.setMesh(V, F);

  pgo::Profiling::setProfilingEnabled(true);
  pgo::Profiling::resetProfileStatistics();

  const double energy0 = wrapper.func(x);
  ES::VXd gradient0 = ES::VXd::Zero(x.size());
  wrapper.gradient(x, gradient0);
  ES::SpMatD hessian0;
  wrapper.hessianDirect(x, hessian0);

  const double energy1 = wrapper.func(x);
  ES::VXd gradient1 = ES::VXd::Zero(x.size());
  wrapper.gradient(x, gradient1);
  ES::SpMatD hessian1;
  wrapper.hessianDirect(x, hessian1);

  const auto stats = pgo::Profiling::snapshotProfileStatistics();
  const ProfileStat *activeSetBuild = findStat(stats, pgo::Contact::SurfaceIPCProfileSections::kBuildActiveSet);

  pgo::Profiling::setProfilingEnabled(false);
  pgo::Profiling::resetProfileStatistics();

  ASSERT_NE(activeSetBuild, nullptr);
  EXPECT_EQ(activeSetBuild->callCount, 6u);
  EXPECT_NEAR(energy1, energy0, 1e-12);
  EXPECT_LT(relativeError(gradient1, gradient0), 1e-12);
  EXPECT_LT(relativeError(sparseToDense(hessian1), sparseToDense(hessian0)), 1e-12);
}
