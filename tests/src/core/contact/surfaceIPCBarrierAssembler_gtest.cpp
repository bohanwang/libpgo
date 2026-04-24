#include <gtest/gtest.h>

#include "ipc/broadPhase/surfaceIPCSelfBroadPhase.h"
#include "ipc/core/surfaceIPCBarrierAssembler.h"
#include "ipc/topology/surfaceIPCTopology.h"
#include "ipc/core/surfaceIPCCore.h"

#include "testCIPCHelpers.h"

namespace
{
namespace ES = pgo::EigenSupport;
using pgo::Contact::CIPC::EEPair;
using pgo::Contact::CIPC::PTPair;
using pgo::Contact::CIPC::SurfaceIPCBarrierAssembler;
using pgo::Contact::CIPC::SurfaceIPCCore;
using pgo::Contact::CIPC::SurfaceIPCSelfBroadPhase;
using pgo::Contact::CIPC::SurfaceIPCTopology;
using pgo::Contact::CIPCTest::flattenPositions;
using pgo::Contact::CIPCTest::makeTwoTriangleMesh;
using pgo::Contact::CIPCTest::relativeError;
using pgo::Contact::CIPCTest::sparseToDense;
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
  std::vector<PTPair> ptPairs;
  std::vector<EEPair> eePairs;
  SurfaceIPCSelfBroadPhase().buildPairs(topology, x, dhat, ptPairs, eePairs);

  const SurfaceIPCBarrierAssembler assembler;
  const double helperEnergy = assembler.computeEnergy(x, ptPairs, eePairs, topology.numVerts, dhat, kappa, epsEE);
  ES::VXd helperGradient(x.size());
  assembler.computeGradient(x, ptPairs, eePairs, topology.numVerts, dhat, kappa, epsEE, helperGradient);
  ES::SpMatD helperHessian;
  assembler.computeHessian(x, ptPairs, eePairs, topology.numVerts, dhat, kappa, epsEE, helperHessian);

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
