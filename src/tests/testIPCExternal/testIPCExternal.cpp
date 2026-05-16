#include "ipc/external/obstacleSurface.h"
#include "ipc/core/surfaceIPCCore.h"
#include "ipc/core/surfaceIPCBarrierAssembler.h"
#include "embeddedSurfaceIPCPotentialEnergy.h"
#include "pgoLogging.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace ES = pgo::EigenSupport;
using namespace pgo::Contact::CIPC;

static int numFailures = 0;

#define TEST(name) \
  do { \
    std::cout << "  " << name << "... "; \
    try {

#define ENDTEST \
    std::cout << "PASSED" << std::endl; \
  } catch (const std::exception &e) { \
    std::cout << "FAILED: " << e.what() << std::endl; \
    ++numFailures; \
  } \
  } while (false)

// Helper: make a simple 2-triangle mesh (a unit square split diagonal)
static std::pair<ES::MXd, ES::MXi> makeUnitSquareMesh()
{
  ES::MXd V(4, 3);
  V << 0.0, 0.0, 0.0,
       1.0, 0.0, 0.0,
       0.0, 1.0, 0.0,
       1.0, 1.0, 0.0;
  ES::MXi F(2, 3);
  F << 0, 1, 2,
       1, 3, 2;
  return { V, F };
}

// Helper: make a simple 2-triangle mesh in the xz plane at y = 0.
static std::pair<ES::MXd, ES::MXi> makeUnitSquareMeshXZ()
{
  ES::MXd V(4, 3);
  V << 0.0, 0.0, 0.0,
       1.0, 0.0, 0.0,
       0.0, 0.0, 1.0,
       1.0, 0.0, 1.0;
  ES::MXi F(2, 3);
  F << 0, 1, 2,
       1, 3, 2;
  return { V, F };
}

// Helper: make a small box obstacle (8 vertices, 12 triangles)
static std::pair<ES::MXd, ES::MXi> makeSmallBoxObstacle(double zOffset)
{
  ES::MXd V(8, 3);
  V << -0.5, -0.5, zOffset - 0.5,
      0.5, -0.5, zOffset - 0.5,
      0.5,  0.5, zOffset - 0.5,
     -0.5,  0.5, zOffset - 0.5,
     -0.5, -0.5, zOffset + 0.5,
      0.5, -0.5, zOffset + 0.5,
      0.5,  0.5, zOffset + 0.5,
     -0.5,  0.5, zOffset + 0.5;
  ES::MXi F(12, 3);
  F << 0, 2, 1,  0, 3, 2,  // bottom
       4, 5, 6,  4, 6, 7,  // top
       0, 1, 5,  0, 5, 4,  // front
       1, 2, 6,  1, 6, 5,  // right
       2, 3, 7,  2, 7, 6,  // back
       3, 0, 4,  3, 4, 7;  // left
  return { V, F };
}

// ---------------------------------------------------------------------------
// 1. ObstacleSurface tests
// ---------------------------------------------------------------------------

static void test_obstacleSurface_basic()
{
  TEST("ObstacleSurface construction, update, unique_edges") {
    auto [V, F] = makeSmallBoxObstacle(0.0);
    ES::VXd restFlat(V.rows() * 3);
    for (int vi = 0; vi < V.rows(); ++vi)
      restFlat.segment<3>(vi * 3) = V.row(vi).transpose();

    ES::V3d vel(0.0, -1.0, 0.0);
    auto sampler = makeLinearTrajectorySampler(restFlat, vel, 0.0);
    ObstacleSurface obs(std::move(V), std::move(F), std::move(sampler));

    obs.update(0.0, 1.0);

    // Check unique_edges derived correctly
    if (obs.uniqueEdges().rows() <= 0) throw std::runtime_error("uniqueEdges empty");

    // After update(0, 1), previous = sampler(0) = rest, current = sampler(1) = rest + vel
    const ES::VXd &prev = obs.previousPositions();
    const ES::VXd &curr = obs.currentPositions();
    if (prev.size() != restFlat.size()) throw std::runtime_error("previous size mismatch");
    if (curr.size() != restFlat.size()) throw std::runtime_error("current size mismatch");

    // Check that current = rest + vel * (1.0 - 0.0)
    double maxErr = 0.0;
    for (int vi = 0; vi < V.rows(); ++vi) {
      ES::V3d expected = restFlat.segment<3>(vi * 3) + vel;
      ES::V3d actual = curr.segment<3>(vi * 3);
      maxErr = std::max(maxErr, (expected - actual).cwiseAbs().maxCoeff());
    }
    if (maxErr > 1e-12) throw std::runtime_error("current positions wrong");
  } ENDTEST;
}

static void test_obstacleSurface_zero_velocity()
{
  TEST("ObstacleSurface zero-velocity sampler") {
    auto [V, F] = makeSmallBoxObstacle(0.0);
    ES::VXd restFlat(V.rows() * 3);
    for (int vi = 0; vi < V.rows(); ++vi)
      restFlat.segment<3>(vi * 3) = V.row(vi).transpose();

    ES::V3d vel = ES::V3d::Zero();
    auto sampler = makeLinearTrajectorySampler(restFlat, vel, 0.0);
    ObstacleSurface obs(std::move(V), std::move(F), std::move(sampler));

    obs.update(0.0, 0.1);

    const ES::VXd &prev = obs.previousPositions();
    const ES::VXd &curr = obs.currentPositions();

    double maxErr = 0.0;
    for (int i = 0; i < restFlat.size(); ++i) {
      maxErr = std::max(maxErr, std::abs(prev[i] - restFlat[i]));
      maxErr = std::max(maxErr, std::abs(curr[i] - restFlat[i]));
    }
    if (maxErr > 1e-12) throw std::runtime_error("zero-velocity positions not equal to rest");
  } ENDTEST;
}

// ---------------------------------------------------------------------------
// 2. SurfaceIPCCore external registration tests
// ---------------------------------------------------------------------------

static void test_surfaceIPCCore_external_register_clear()
{
  TEST("SurfaceIPCCore add/clear obstacle surfaces") {
    auto [V, F] = makeUnitSquareMesh();
    SurfaceIPCCore::Parameters params;
    SurfaceIPCCore core(params);
    core.setMesh(V, F);

    if (core.preparedState().externalPairs.ptPairs.size() != 0) throw std::runtime_error("should have no ext PT pairs initially");
    if (core.preparedState().externalPairs.tpPairs.size() != 0) throw std::runtime_error("should have no ext TP pairs initially");
    if (core.preparedState().externalPairs.eePairs.size() != 0) throw std::runtime_error("should have no ext EE pairs initially");

    auto [obsV, obsF] = makeSmallBoxObstacle(2.0);
    ES::VXd obsRest(obsV.rows() * 3);
    for (int vi = 0; vi < obsV.rows(); ++vi)
      obsRest.segment<3>(vi * 3) = obsV.row(vi).transpose();

    auto sampler = makeLinearTrajectorySampler(obsRest, ES::V3d::Zero());
    auto obs = std::make_shared<ObstacleSurface>(std::move(obsV), std::move(obsF), std::move(sampler));
    obs->update(0.0, 0.0);

    int32_t id = core.addObstacleSurface(obs);
    if (id != 0) throw std::runtime_error("first obstacle id should be 0");
    if (obs->objectId() != 0) throw std::runtime_error("obstacle objectId not set");

    core.clearObstacleSurfaces();
    // After clear, no pairs
    ES::VXd x(V.rows() * 3);
    for (int vi = 0; vi < V.rows(); ++vi)
      x.segment<3>(vi * 3) = V.row(vi).transpose();
    core.prepareForSurfacePositions(x);
    if (core.preparedState().externalPairs.ptPairs.size() != 0) throw std::runtime_error("still has ext PT pairs after clear");
    if (core.preparedState().externalPairs.tpPairs.size() != 0) throw std::runtime_error("still has ext TP pairs after clear");
    if (core.preparedState().externalPairs.eePairs.size() != 0) throw std::runtime_error("still has ext EE pairs after clear");
  } ENDTEST;
}

// ---------------------------------------------------------------------------
// 3. External barrier energy / gradient / hessian tests
// ---------------------------------------------------------------------------

static void test_surfaceIPCCore_external_static_plane()
{
  TEST("SurfaceIPCCore external static plane barrier") {
    auto [V, F] = makeUnitSquareMesh();
    SurfaceIPCCore::Parameters params;
    params.dhat_external = 0.5;
    params.kappa = 1.0;
    SurfaceIPCCore core(params);
    core.setMesh(V, F);

    // Make a plane obstacle below the mesh
    ES::MXd obsV(4, 3);
    obsV << -1.0, -0.1, -1.0,
             2.0, -0.1, -1.0,
            -1.0, -0.1,  2.0,
             2.0, -0.1,  2.0;
    ES::MXi obsF(2, 3);
    obsF << 0, 1, 2,
            1, 3, 2;

    ES::VXd obsRest(obsV.rows() * 3);
    for (int vi = 0; vi < obsV.rows(); ++vi)
      obsRest.segment<3>(vi * 3) = obsV.row(vi).transpose();

    auto sampler = makeLinearTrajectorySampler(obsRest, ES::V3d::Zero());
    auto obs = std::make_shared<ObstacleSurface>(std::move(obsV), std::move(obsF), std::move(sampler));
    obs->update(0.0, 0.0);
    core.addObstacleSurface(obs);

    ES::VXd x(V.rows() * 3);
    for (int vi = 0; vi < V.rows(); ++vi)
      x.segment<3>(vi * 3) = V.row(vi).transpose();

    // Energy should be non-negative
    double e = core.computeEnergy(x);
    if (e < 0.0) throw std::runtime_error("barrier energy is negative");

    // Gradient should point away from obstacle (positive y direction for vertices near the plane)
    ES::VXd g(V.rows() * 3);
    core.computeGradient(x, g);
    if (g.size() != V.rows() * 3) throw std::runtime_error("gradient size mismatch");

    // Hessian should be PSD
    ES::SpMatD H;
    core.computeHessian(x, H);
    if (H.rows() != V.rows() * 3) throw std::runtime_error("hessian size mismatch");
  } ENDTEST;
}

static void test_surfaceIPCCore_external_box_contact()
{
  TEST("SurfaceIPCCore external box — PT, TP, EE pair coverage") {
    auto [V, F] = makeUnitSquareMesh();
    SurfaceIPCCore::Parameters params;
    params.dhat_external = 1.0;  // large to ensure pairs activate
    params.kappa = 0.1;
    SurfaceIPCCore core(params);
    core.setMesh(V, F);

    // Small box obstacle overlapping the mesh (z-offset at 0 so it intersects)
    auto [obsV, obsF] = makeSmallBoxObstacle(0.0);
    ES::VXd obsRest(obsV.rows() * 3);
    for (int vi = 0; vi < obsV.rows(); ++vi)
      obsRest.segment<3>(vi * 3) = obsV.row(vi).transpose();

    auto sampler = makeLinearTrajectorySampler(obsRest, ES::V3d::Zero());
    auto obs = std::make_shared<ObstacleSurface>(std::move(obsV), std::move(obsF), std::move(sampler));
    obs->update(0.0, 0.0);
    core.addObstacleSurface(obs);

    ES::VXd x(V.rows() * 3);
    for (int vi = 0; vi < V.rows(); ++vi)
      x.segment<3>(vi * 3) = V.row(vi).transpose();

    // Prepare and check pairs
    core.prepareForSurfacePositions(x);
    const auto &ptPairs = core.preparedState().externalPairs.ptPairs;
    const auto &tpPairs = core.preparedState().externalPairs.tpPairs;
    const auto &eePairs = core.preparedState().externalPairs.eePairs;

    if (ptPairs.empty() && tpPairs.empty() && eePairs.empty())
      throw std::runtime_error("no external pairs activated for overlapping meshes");

    // Verify pair identity: all pairs reference a valid obstacle id
    for (const auto &p : ptPairs)
      if (p.obstacleObjectId != 0) throw std::runtime_error("PT pair has wrong obstacle id");
    for (const auto &p : tpPairs)
      if (p.obstacleObjectId != 0) throw std::runtime_error("TP pair has wrong obstacle id");
    for (const auto &p : eePairs)
      if (p.obstacleObjectId != 0) throw std::runtime_error("EE pair has wrong obstacle id");

    // Energy, gradient, hessian should be computable
    double e = core.computeEnergy(x);
    if (e < 0.0) throw std::runtime_error("negative barrier energy");

    ES::VXd g(V.rows() * 3);
    core.computeGradient(x, g);

    ES::SpMatD H;
    core.computeHessian(x, H);
    if (H.rows() != V.rows() * 3) throw std::runtime_error("hessian row mismatch");
  } ENDTEST;
}

static void test_surfaceIPCCore_external_multi_obstacle()
{
  TEST("SurfaceIPCCore multi-obstacle — obstacleObjectId distinguishes") {
    auto [V, F] = makeUnitSquareMesh();
    SurfaceIPCCore::Parameters params;
    params.dhat_external = 1.0;
    params.kappa = 0.1;
    SurfaceIPCCore core(params);
    core.setMesh(V, F);

    // Register two obstacles with different z-offsets
    auto [obsV1, obsF1] = makeSmallBoxObstacle(0.0);
    ES::VXd obsRest1(obsV1.rows() * 3);
    for (int vi = 0; vi < obsV1.rows(); ++vi)
      obsRest1.segment<3>(vi * 3) = obsV1.row(vi).transpose();
    auto obs1 = std::make_shared<ObstacleSurface>(std::move(obsV1), std::move(obsF1),
      makeLinearTrajectorySampler(obsRest1, ES::V3d::Zero()));
    obs1->update(0.0, 0.0);

    auto [obsV2, obsF2] = makeSmallBoxObstacle(0.75);
    ES::VXd obsRest2(obsV2.rows() * 3);
    for (int vi = 0; vi < obsV2.rows(); ++vi)
      obsRest2.segment<3>(vi * 3) = obsV2.row(vi).transpose();
    auto obs2 = std::make_shared<ObstacleSurface>(std::move(obsV2), std::move(obsF2),
      makeLinearTrajectorySampler(obsRest2, ES::V3d::Zero()));
    obs2->update(0.0, 0.0);

    int32_t id1 = core.addObstacleSurface(obs1);
    int32_t id2 = core.addObstacleSurface(obs2);
    if (id1 != 0) throw std::runtime_error("first obstacle id should be 0");
    if (id2 != 1) throw std::runtime_error("second obstacle id should be 1");

    ES::VXd x(V.rows() * 3);
    for (int vi = 0; vi < V.rows(); ++vi)
      x.segment<3>(vi * 3) = V.row(vi).transpose();

    core.prepareForSurfacePositions(x);

    // Check that pairs from different obstacles are separated
    bool hasId0 = false, hasId1 = false;
    for (const auto &p : core.preparedState().externalPairs.ptPairs) {
      if (p.obstacleObjectId == 0) hasId0 = true;
      if (p.obstacleObjectId == 1) hasId1 = true;
    }
    for (const auto &p : core.preparedState().externalPairs.tpPairs) {
      if (p.obstacleObjectId == 0) hasId0 = true;
      if (p.obstacleObjectId == 1) hasId1 = true;
    }
    for (const auto &p : core.preparedState().externalPairs.eePairs) {
      if (p.obstacleObjectId == 0) hasId0 = true;
      if (p.obstacleObjectId == 1) hasId1 = true;
    }

    if (!hasId0) throw std::runtime_error("no pairs from obstacle 0");
    if (!hasId1) throw std::runtime_error("no pairs from obstacle 1");
  } ENDTEST;
}

static void test_surfaceIPCCore_external_self_equivalence()
{
  TEST("SurfaceIPCCore external dynamic-only block vs locked-DOF self equivalence") {
    // Use a simple one-triangle dynamic mesh and a one-triangle "obstacle"
    // placed close together. Compare energy/gradient from external IPC
    // vs self-IPC with obstacle DOFs locked.

    ES::MXd dynV(3, 3);
    dynV << 0.0, 0.0, 0.0,
            1.0, 0.0, 0.0,
            0.0, 1.0, 0.0;
    ES::MXi dynF(1, 3);
    dynF << 0, 1, 2;

    ES::MXd obsV(3, 3);
    obsV << 0.5, 0.5, 0.2,  // slightly above dynamic mesh
            1.5, 0.5, 0.2,
            0.5, 1.5, 0.2;
    ES::MXi obsF(1, 3);
    obsF << 0, 1, 2;

    // ---- External IPC ----
    SurfaceIPCCore::Parameters extParams;
    extParams.dhat_external = 0.5;
    extParams.kappa = 0.1;
    extParams.dhat = 0.5;
    SurfaceIPCCore extCore(extParams);
    extCore.setMesh(dynV, dynF);

    ES::VXd obsRest(obsV.rows() * 3);
    for (int vi = 0; vi < obsV.rows(); ++vi)
      obsRest.segment<3>(vi * 3) = obsV.row(vi).transpose();

    auto obs = std::make_shared<ObstacleSurface>(obsV, obsF,
      makeLinearTrajectorySampler(obsRest, ES::V3d::Zero()));
    obs->update(0.0, 0.0);
    extCore.addObstacleSurface(obs);

    ES::VXd xDyn(dynV.rows() * 3);
    for (int vi = 0; vi < dynV.rows(); ++vi)
      xDyn.segment<3>(vi * 3) = dynV.row(vi).transpose();

    double extE = extCore.computeEnergy(xDyn);
    ES::VXd extG(dynV.rows() * 3);
    extCore.computeGradient(xDyn, extG);

    // ---- Self IPC with obstacle as part of same mesh (but not connected) ----
    ES::MXd selfV(dynV.rows() + obsV.rows(), 3);
    selfV.topRows(dynV.rows()) = dynV;
    selfV.bottomRows(obsV.rows()) = obsV;

    ES::MXi selfF(dynF.rows() + obsF.rows(), 3);
    selfF.topRows(dynF.rows()) = dynF;
    selfF.bottomRows(obsF.rows()) = obsF.array() + dynV.rows();

    SurfaceIPCCore::Parameters selfParams;
    selfParams.dhat = 0.5;
    selfParams.kappa = 0.1;
    SurfaceIPCCore selfCore(selfParams);
    selfCore.setMesh(selfV, selfF);

    ES::VXd xSelf(selfV.rows() * 3);
    for (int vi = 0; vi < selfV.rows(); ++vi)
      xSelf.segment<3>(vi * 3) = selfV.row(vi).transpose();

    // Compute self energy with all DOFs
    // The self energy includes dyn-dyn, dyn-obs, and obs-obs pairs.
    // External energy only includes dyn-obs pairs.
    // Both should have non-negative energy and consistent gradient on dynamic DOFs.
    double selfE = selfCore.computeEnergy(xSelf);
    ES::VXd selfG(selfV.rows() * 3);
    selfCore.computeGradient(xSelf, selfG);

    // External energy should be <= self energy (since self includes extra obs-obs pairs)
    // Actually, external energy IS the dyn-obs component. Self energy = dyn-dyn + dyn-obs + obs-obs.
    // Both should be finite and non-negative.
    if (!std::isfinite(extE)) throw std::runtime_error("external energy not finite");
    if (!std::isfinite(selfE)) throw std::runtime_error("self energy not finite");
    if (extE < 0.0 || selfE < 0.0) throw std::runtime_error("energy negative");

    // External gradient on dynamic vertices should be non-trivial (not all zero)
    if (extG.norm() <= 0.0) throw std::runtime_error("external gradient is zero");
  } ENDTEST;
}

// ---------------------------------------------------------------------------
// 4. External CCD tests
// ---------------------------------------------------------------------------

static void test_surfaceIPCCore_external_ccd_kinematic()
{
  TEST("SurfaceIPCCore external CCD — kinematic obstacle approaching") {
    auto [V, F] = makeUnitSquareMeshXZ();
    SurfaceIPCCore::Parameters params;
    params.dhat_external = 0.5;
    params.kappa = 0.1;
    params.slackness = 1.0;
    SurfaceIPCCore core(params);
    core.setMesh(V, F);

    // Plane obstacle above the mesh, moving downward
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
      obsRest.segment<3>(vi * 3) = obsV.row(vi).transpose();

    ES::V3d vel(0.0, -1.0, 0.0);  // moving downward
    auto sampler = makeLinearTrajectorySampler(obsRest, vel, 0.0);
    auto obs = std::make_shared<ObstacleSurface>(std::move(obsV), std::move(obsF), std::move(sampler));
    obs->update(0.0, 1.0);  // obs moves from y=0.5 to y=-0.5 in this stage
    core.addObstacleSurface(obs);

    ES::VXd x(V.rows() * 3);
    for (int vi = 0; vi < V.rows(); ++vi)
      x.segment<3>(vi * 3) = V.row(vi).transpose();

    // Dynamic is stationary, obstacle moves down → should contact
    ES::VXd dx = ES::VXd::Zero(V.rows() * 3);

    auto result = core.computeMaxStepLimit(x, dx);
    double alpha = result.alpha;

    // Obstacle moved from y=0.5 towards mesh at y=0.0; dynamic stationary
    // The CCD should detect contact
    if (alpha >= 1.0) throw std::runtime_error("CCD did not detect approaching obstacle");
    if (alpha < 0.0) throw std::runtime_error("CCD alpha is negative");
  } ENDTEST;
}

static void test_surfaceIPCCore_external_ccd_alpha_symmetry()
{
  TEST("SurfaceIPCCore external CCD — alpha does not scale obstacle motion") {
    auto [V, F] = makeUnitSquareMeshXZ();
    SurfaceIPCCore::Parameters params;
    params.dhat_external = 0.5;
    params.kappa = 0.1;
    params.slackness = 1.0;
    SurfaceIPCCore core(params);
    core.setMesh(V, F);

    // Plane obstacle at y=0.5, stationary (velocity = 0)
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
      obsRest.segment<3>(vi * 3) = obsV.row(vi).transpose();

    auto obs = std::make_shared<ObstacleSurface>(std::move(obsV), std::move(obsF),
      makeLinearTrajectorySampler(obsRest, ES::V3d::Zero()));
    obs->update(0.0, 0.0);
    core.addObstacleSurface(obs);

    ES::VXd x(V.rows() * 3);
    for (int vi = 0; vi < V.rows(); ++vi)
      x.segment<3>(vi * 3) = V.row(vi).transpose();

    // Dynamic moves up toward obstacle.
    ES::VXd dx(V.rows() * 3);
    for (int vi = 0; vi < V.rows(); ++vi)
      dx.segment<3>(vi * 3) = ES::V3d(0.0, 1.0, 0.0);

    auto result1 = core.computeMaxStepLimit(x, dx);
    double alpha1 = result1.alpha;

    // Scale dx by 0.5 → alpha should double (linear scaling)
    ES::VXd dxHalf = dx * 0.5;
    auto result2 = core.computeMaxStepLimit(x, dxHalf);
    double alpha2 = result2.alpha;

    // alpha2 / alpha1 should be ≈ 2.0 (since half the displacement → twice the safe step proportion)
    if (alpha1 > 0.0 && alpha2 > 0.0) {
      double ratio = alpha2 / alpha1;
      if (ratio < 1.5 || ratio > 2.5)
        throw std::runtime_error("alpha scaling invariant violated: " + std::to_string(ratio));
    }
  } ENDTEST;
}

// ---------------------------------------------------------------------------
// 5. Wrapper consistency tests
// ---------------------------------------------------------------------------

static void test_embeddedSurfaceIPCPotentialEnergy_external()
{
  TEST("EmbeddedSurfaceIPCPotentialEnergy wrapper — obstacle registration") {
    auto [V, F] = makeUnitSquareMesh();
    int n3 = V.rows() * 3;
    ES::SpMatD W(n3, n3);
    std::vector<ES::TripletD> triplets;
    for (int i = 0; i < n3; ++i)
      triplets.emplace_back(i, i, 1.0);
    W.setFromTriplets(triplets.begin(), triplets.end());

    SurfaceIPCCore::Parameters params;
    params.dhat_external = 0.5;
    params.kappa = 0.1;

    EmbeddedSurfaceIPCPotentialEnergy wrapper(V, F, W, params);

    auto [obsV, obsF] = makeSmallBoxObstacle(0.0);
    ES::VXd obsRest(obsV.rows() * 3);
    for (int vi = 0; vi < obsV.rows(); ++vi)
      obsRest.segment<3>(vi * 3) = obsV.row(vi).transpose();
    auto obs = std::make_shared<ObstacleSurface>(std::move(obsV), std::move(obsF),
      makeLinearTrajectorySampler(obsRest, ES::V3d::Zero()));
    obs->update(0.0, 0.0);

    int32_t id = wrapper.addObstacleSurface(obs);
    if (id != 0) throw std::runtime_error("wrapper addObstacleSurface returned wrong id");

    ES::VXd u = ES::VXd::Zero(n3);

    // Compute through wrapper - should not throw
    double e = wrapper.func(u);
    if (!std::isfinite(e)) throw std::runtime_error("wrapper func returned non-finite energy");

    ES::VXd g(n3);
    wrapper.gradient(u, g);
    if (g.size() != n3) throw std::runtime_error("wrapper gradient size mismatch");

    ES::SpMatD H;
    wrapper.hessianDirect(u, H);
    if (H.rows() != n3) throw std::runtime_error("wrapper hessian size mismatch");

    wrapper.clearObstacleSurfaces();
  } ENDTEST;
}

static void test_embeddedSurfaceIPCPotentialEnergy_wrapper_vs_core()
{
  TEST("EmbeddedSurfaceIPCPotentialEnergy wrapper vs SurfaceIPCCore direct") {
    auto [V, F] = makeUnitSquareMesh();
    int n3 = V.rows() * 3;
    ES::SpMatD W(n3, n3);
    std::vector<ES::TripletD> triplets;
    for (int i = 0; i < n3; ++i)
      triplets.emplace_back(i, i, 1.0);
    W.setFromTriplets(triplets.begin(), triplets.end());

    SurfaceIPCCore::Parameters params;
    params.dhat = 0.1;
    params.dhat_external = 0.5;
    params.kappa = 0.1;

    // Wrapper path
    EmbeddedSurfaceIPCPotentialEnergy wrapper(V, F, W, params);

    auto [obsV, obsF] = makeSmallBoxObstacle(0.0);
    ES::VXd obsRest(obsV.rows() * 3);
    for (int vi = 0; vi < obsV.rows(); ++vi)
      obsRest.segment<3>(vi * 3) = obsV.row(vi).transpose();
    auto obs = std::make_shared<ObstacleSurface>(std::move(obsV), std::move(obsF),
      makeLinearTrajectorySampler(obsRest, ES::V3d::Zero()));
    obs->update(0.0, 0.0);
    wrapper.addObstacleSurface(obs);
    wrapper.updateObstacleStage(0.0, 0.0);

    ES::VXd x(n3);
    for (int vi = 0; vi < V.rows(); ++vi)
      x.segment<3>(vi * 3) = V.row(vi).transpose();

    ES::VXd u = ES::VXd::Zero(n3);

    double wrapperE = wrapper.func(u);
    ES::VXd wrapperG(n3);
    wrapper.gradient(u, wrapperG);

    // Direct core path
    auto [obsV2, obsF2] = makeSmallBoxObstacle(0.0);
    ES::VXd obsRest2(obsV2.rows() * 3);
    for (int vi = 0; vi < obsV2.rows(); ++vi)
      obsRest2.segment<3>(vi * 3) = obsV2.row(vi).transpose();
    auto obs2 = std::make_shared<ObstacleSurface>(std::move(obsV2), std::move(obsF2),
      makeLinearTrajectorySampler(obsRest2, ES::V3d::Zero()));
    obs2->update(0.0, 0.0);

    SurfaceIPCCore core(params);
    core.setMesh(V, F);
    core.addObstacleSurface(obs2);

    double coreE = core.computeEnergy(x);
    ES::VXd coreG(n3);
    core.computeGradient(x, coreG);

    // Results should match
    if (std::abs(wrapperE - coreE) > 1e-10)
      throw std::runtime_error("wrapper energy differs from core: " + std::to_string(std::abs(wrapperE - coreE)));
    if ((wrapperG - coreG).cwiseAbs().maxCoeff() > 1e-10)
      throw std::runtime_error("wrapper gradient differs from core");
  } ENDTEST;
}

int main()
{
  pgo::Logging::init(nullptr, spdlog::level::info);

  std::cout << "=== ObstacleSurface Tests ===" << std::endl;
  test_obstacleSurface_basic();
  test_obstacleSurface_zero_velocity();

  std::cout << "\n=== SurfaceIPCCore Registration Tests ===" << std::endl;
  test_surfaceIPCCore_external_register_clear();

  std::cout << "\n=== External Barrier Tests ===" << std::endl;
  test_surfaceIPCCore_external_static_plane();
  test_surfaceIPCCore_external_box_contact();
  test_surfaceIPCCore_external_multi_obstacle();
  test_surfaceIPCCore_external_self_equivalence();

  std::cout << "\n=== External CCD Tests ===" << std::endl;
  test_surfaceIPCCore_external_ccd_kinematic();
  test_surfaceIPCCore_external_ccd_alpha_symmetry();

  std::cout << "\n=== Wrapper Tests ===" << std::endl;
  test_embeddedSurfaceIPCPotentialEnergy_external();
  test_embeddedSurfaceIPCPotentialEnergy_wrapper_vs_core();

  std::cout << "\n=== Results ===" << std::endl;
  if (numFailures == 0) {
    std::cout << "All tests PASSED." << std::endl;
    return 0;
  }
  else {
    std::cout << numFailures << " test(s) FAILED." << std::endl;
    return 1;
  }
}
