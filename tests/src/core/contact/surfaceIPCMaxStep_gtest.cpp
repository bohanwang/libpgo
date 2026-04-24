#include <gtest/gtest.h>

#include "ipc/core/surfaceIPCMaxStep.h"
#include "ipc/topology/surfaceIPCTopology.h"
#include "ipc/core/surfaceIPCCore.h"

#include "testCIPCHelpers.h"

namespace
{
namespace ES = pgo::EigenSupport;
using pgo::Contact::CIPC::SurfaceIPCCore;
using pgo::Contact::CIPC::SurfaceIPCMaxStep;
using pgo::Contact::CIPC::SurfaceIPCTopology;
using pgo::Contact::CIPCTest::flattenPositions;
using pgo::Contact::CIPCTest::makeTwoTriangleMesh;
}  // namespace

TEST(SurfaceIPCMaxStepGTest, HelperMatchesSurfaceIPCCoreMaxStep)
{
  const auto [V, F] = makeTwoTriangleMesh();
  const ES::VXd x = flattenPositions(V);
  ES::VXd dx = ES::VXd::Zero(x.size());
  for (int vi = 3; vi < 6; ++vi)
    dx[3 * vi + 2] = -0.1;

  SurfaceIPCTopology topology;
  topology.setMesh(V, F);

  SurfaceIPCCore core;
  SurfaceIPCCore::Parameters params;
  params.dhat = 0.1;
  params.kappa = 1.0;
  params.eps_ee = 0.0;
  params.slackness = 0.9;
  core.setParameters(params);
  core.setMesh(V, F);

  const double helperAlpha = SurfaceIPCMaxStep().compute(topology, x, dx, params.dhat, params.slackness);
  const double coreAlpha = core.computeMaxStepSize(x, dx);

  EXPECT_GT(helperAlpha, 0.0);
  EXPECT_LT(helperAlpha, 1.0);
  EXPECT_NEAR(helperAlpha, coreAlpha, 1e-12);
}
