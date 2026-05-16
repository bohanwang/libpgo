#include "ipc/core/surfaceIPCMaxStep.h"
#include "ipc/core/surfaceIPCCore.h"
#include "ipc/external/obstacleSurface.h"
#include "ipc/topology/surfaceIPCTopology.h"

#include <gtest/gtest.h>

#include <memory>
#include <vector>

namespace ES = pgo::EigenSupport;
using pgo::Contact::CIPC::ObstacleSurface;
using pgo::Contact::CIPC::SurfaceIPCCore;
using pgo::Contact::CIPC::SurfaceIPCTopology;

TEST(SurfaceIPCExternalMaxStepGTest, HelperMatchesSurfaceIPCCoreExternalContribution)
{
  ES::MXd V(4, 3);
  V << 0.0, 0.0, 0.0,
       1.0, 0.0, 0.0,
       0.0, 0.0, 1.0,
       1.0, 0.0, 1.0;
  ES::MXi F(2, 3);
  F << 0, 1, 2,
       1, 3, 2;

  ES::MXd obsV(4, 3);
  obsV << -1.0, 0.5, -1.0,
           2.0, 0.5, -1.0,
          -1.0, 0.5,  2.0,
           2.0, 0.5,  2.0;
  ES::MXi obsF(2, 3);
  obsF << 0, 1, 2,
          1, 3, 2;

  ES::VXd obsRest(obsV.rows() * 3);
  for (int vi = 0; vi < obsV.rows(); ++vi)
    obsRest.segment<3>(3 * vi) = obsV.row(vi).transpose();

  auto obs = std::make_shared<ObstacleSurface>(
    obsV, obsF,
    pgo::Contact::CIPC::makeLinearTrajectorySampler(obsRest, ES::V3d(0.0, -1.0, 0.0)));
  obs->update(0.0, 1.0);

  SurfaceIPCCore::Parameters params;
  params.dhat_external = 0.5;
  params.slackness = 1.0;

  SurfaceIPCCore core(params);
  core.setMesh(V, F);
  core.addObstacleSurface(obs);

  ES::VXd x(V.rows() * 3);
  for (int vi = 0; vi < V.rows(); ++vi)
    x.segment<3>(3 * vi) = V.row(vi).transpose();
  const ES::VXd dx = ES::VXd::Zero(V.rows() * 3);

  SurfaceIPCTopology topology;
  topology.setMesh(V, F);
  const std::vector<std::shared_ptr<ObstacleSurface>> obstacles = { obs };

  const double selfAlpha = computeSelfMaxStep(topology, x, dx, params.dhat, params.slackness);
  ASSERT_NEAR(selfAlpha, 1.0, 1e-12);

  const double helperAlpha = computeExternalMaxStep(
    topology, x, dx, obstacles, params.dhat_external, params.slackness);
  const double coreAlpha = core.computeMaxStepLimit(x, dx).alpha;

  EXPECT_NEAR(helperAlpha, coreAlpha, 1e-12);
  EXPECT_LT(helperAlpha, 1.0);
  EXPECT_GE(helperAlpha, 0.0);
}
