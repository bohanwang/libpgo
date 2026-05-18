#include "ipc/core/surfaceIPCMaxStep.h"
#include "ipc/core/surfaceIPCCore.h"
#include "ipc/external/obstacleSurface.h"
#include "ipc/topology/surfaceIPCTopology.h"

#include <gtest/gtest.h>

#include <utility>
#include <vector>

namespace ES = pgo::EigenSupport;
using pgo::Contact::CIPC::ObstacleSurface;
using pgo::Contact::CIPC::SurfaceIPCCore;
using pgo::Contact::CIPC::SurfaceIPCTopology;

// Dynamic mesh moves into a fixed-pose obstacle: external max-step
// must clamp alpha < 1 and match SurfaceIPCCore's contribution exactly.
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

  // Obstacle is a plane at y = 0.5, 1 unit above the dynamic mesh.
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

  // Sampler is stationary; obstacle stays at its rest pose.
  auto makeObs = [&]() {
    ObstacleSurface o(obsV, obsF,
      pgo::Contact::CIPC::makeLinearTrajectorySampler(obsRest, ES::V3d::Zero()));
    o.update(0.0);
    return o;
  };

  SurfaceIPCCore::Parameters params;
  params.dhat_external = 0.5;
  params.slackness = 1.0;

  std::vector<ObstacleSurface> coreObstacles;
  coreObstacles.emplace_back(makeObs());
  SurfaceIPCCore core(params, std::move(coreObstacles));
  core.setMesh(V, F);

  ES::VXd x(V.rows() * 3);
  for (int vi = 0; vi < V.rows(); ++vi)
    x.segment<3>(3 * vi) = V.row(vi).transpose();
  // Dynamic mesh moves straight toward the obstacle (+y) by 1 unit.
  ES::VXd dx = ES::VXd::Zero(V.rows() * 3);
  for (int vi = 0; vi < V.rows(); ++vi)
    dx[3 * vi + 1] = 1.0;

  SurfaceIPCTopology topology;
  topology.setMesh(V, F);
  std::vector<ObstacleSurface> obstacles;
  obstacles.emplace_back(makeObs());
  obstacles.front().setObjectId(0);

  const double selfAlpha = computeSelfMaxStep(topology, x, dx, params.dhat, params.slackness);
  ASSERT_NEAR(selfAlpha, 1.0, 1e-12);

  const double helperAlpha = computeExternalMaxStep(
    topology, x, dx, obstacles, params.dhat_external, params.slackness);
  const double coreAlpha = core.computeMaxStepLimit(x, dx).alpha;

  EXPECT_NEAR(helperAlpha, coreAlpha, 1e-12);
  EXPECT_LT(helperAlpha, 1.0);
  EXPECT_GE(helperAlpha, 0.0);
}

// Regression intent: external max-step only sweeps the dynamic mesh.
// A moving-obstacle sampler combined with zero Newton displacement (dx == 0)
// must NOT clamp alpha. Intra-frame obstacle motion is intentionally NOT
// caught here — the barrier assembler enforces non-penetration on the next
// solve instead.
TEST(SurfaceIPCExternalMaxStepGTest, ZeroDynamicDisplacementReturnsUnitAlphaEvenWhenObstacleMoves)
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

  // Obstacle velocity drives it from y=0.5 toward y=-0.5 over a unit time;
  // the sampler is non-trivial but obstacle motion is invisible to max-step.
  ObstacleSurface obs(obsV, obsF,
    pgo::Contact::CIPC::makeLinearTrajectorySampler(obsRest, ES::V3d(0.0, -1.0, 0.0)));
  obs.setObjectId(0);
  obs.update(1.0);

  ES::VXd x(V.rows() * 3);
  for (int vi = 0; vi < V.rows(); ++vi)
    x.segment<3>(3 * vi) = V.row(vi).transpose();
  const ES::VXd dx = ES::VXd::Zero(V.rows() * 3);

  SurfaceIPCTopology topology;
  topology.setMesh(V, F);
  std::vector<ObstacleSurface> obstacles;
  obstacles.emplace_back(std::move(obs));

  const double alpha = computeExternalMaxStep(topology, x, dx, obstacles, 0.5, 1.0);
  EXPECT_DOUBLE_EQ(alpha, 1.0);
}
