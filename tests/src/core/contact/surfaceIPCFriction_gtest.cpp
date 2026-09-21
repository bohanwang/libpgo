#include <gtest/gtest.h>

#include "ipc/core/surfaceIPCFriction.h"
#include "ipc/core/surfaceIPCCore.h"
#include "embeddedSurfaceIPCPotentialEnergy.h"
#include "ipc/geometry/ipcBarrier.h"
#include "testCIPCHelpers.h"

#include <Eigen/Eigenvalues>
#include <limits>

namespace
{
namespace ES = pgo::EigenSupport;
using namespace pgo::Contact::CIPC;
using namespace pgo::Contact::CIPCTest;

ES::VXd pointTriangle()
{
  ES::MXd V(4, 3);
  V << .2, .3, .05, 0, 0, 0, 1, 0, 0, 0, 1, 0;
  return flattenPositions(V);
}

ES::VXd gradient(const SurfaceIPCFriction &friction, const ES::VXd &x)
{
  ES::VXd g = ES::VXd::Zero(x.size());
  friction.addGradient(x, g);
  return g;
}

void checkDerivatives(const SurfaceIPCFriction &friction, const ES::VXd &x, double h = 1e-7)
{
  const ES::VXd g = gradient(friction, x);
  const ES::VXd fdg = finiteDifferenceGradient([&](const ES::VXd &q) { return friction.computeEnergy(q); }, x, h);
  EXPECT_LT(relativeError(g, fdg), 2e-7);
  ES::SpMatD H;
  friction.computeHessian(x, H);
  const ES::MXd dense = H;
  const ES::MXd fdh = finiteDifferenceHessian([&](const ES::VXd &q) { return gradient(friction, q); }, x, h);
  EXPECT_TRUE(dense.allFinite());
  EXPECT_LT(relativeError(dense, fdh), 2e-5);
  EXPECT_LT((dense - dense.transpose()).norm(), 1e-8 * std::max(1.0, dense.norm()));
  EXPECT_GE(Eigen::SelfAdjointEigenSolver<ES::MXd>(dense).eigenvalues().minCoeff(), -1e-9 * std::max(1.0, dense.norm()));
  ES::V3d total = ES::V3d::Zero();
  for (int i = 0; i < x.size() / 3; ++i)
    total += g.segment<3>(3 * i);
  EXPECT_LT(total.norm(), 1e-10 * std::max(1.0, g.norm()));
}
}  // namespace

TEST(SurfaceIPCFriction, PointTriangleDerivativesAcrossClosestFeaturesAndSlipRegimes)
{
  for (const ES::V3d p : { ES::V3d(.2, .3, .05), ES::V3d(.4, -.02, .05), ES::V3d(-.02, -.03, .05) }) {
    ES::VXd reference = pointTriangle();
    reference.head<3>() = p;
    SurfaceIPCFriction friction;
    friction.update(reference, reference, { { 0, 1, 2, 3, 2.0 } }, {}, .2, 3, 0, .3, .1, .1);
    ASSERT_EQ(friction.numPairs(), 1);
    for (double displacement : { 0.0, .002, .03 }) {
      ES::VXd x = reference;
      x[0] += displacement;
      x[1] += .4 * displacement;
      checkDerivatives(friction, x);
    }
  }
}

TEST(SurfaceIPCFriction, EdgeEdgeDerivativesIncludeParallelAndBoundaryCases)
{
  for (double offset : { 0.3, 1.03 }) {
    for (double angle : { 0.0, .03, 1.0 }) {
      ES::MXd V(4, 3);
      V << 0, 0, 0, 1, 0, 0, offset, -.2 * angle, .05, offset + 1 - angle, .8 * angle, .05;
      const ES::VXd reference = flattenPositions(V);
      SurfaceIPCFriction friction;
      friction.update(reference, reference, {}, { { 0, 1, 2, 3, .7 } }, .2, 3, 0, .3, .1, .1);
      ASSERT_EQ(friction.numPairs(), 1);
      for (double displacement : { 0.0, .002, .03 }) {
        ES::VXd x = reference;
        x[0] += displacement;
        x[4] += displacement * .3;
        checkDerivatives(friction, x);
      }
    }
  }
}

TEST(SurfaceIPCFriction, CoulombForceBarycentricDistributionAndZeroSlipLimit)
{
  const ES::VXd reference = pointTriangle();
  SurfaceIPCFriction friction;
  friction.update(reference, reference, { { 0, 1, 2, 3, 2.0 } }, {}, .2, 3, 0, .3, .1, .1);
  const double muLambda = .3 * (-2 * .05 * 2 * 3 * barrier::dbds(.05 * .05, .2 * .2));
  ES::VXd x = reference;
  x[0] += .02;
  const ES::VXd g = gradient(friction, x);
  EXPECT_NEAR(g[0], muLambda, 1e-10);
  EXPECT_NEAR(g[3], -.5 * muLambda, 1e-10);
  EXPECT_NEAR(g[6], -.2 * muLambda, 1e-10);
  EXPECT_NEAR(g[9], -.3 * muLambda, 1e-10);
  EXPECT_NEAR(g[2], 0, 1e-12);
  ES::SpMatD H;
  friction.computeHessian(reference, H);
  EXPECT_NEAR(H.coeff(0, 0), 2 * muLambda / .01, 1e-8);
  EXPECT_NEAR(H.coeff(2, 2), 0, 1e-12);
  EXPECT_EQ(gradient(friction, reference).norm(), 0);
  for (double t : { .1, .5, .999999, 1.0, 1.000001 }) {
    x = reference;
    x[0] += t * .01;
    EXPECT_LE(gradient(friction, x).head<3>().norm(), muLambda * (1 + 1e-12));
    checkDerivatives(friction, x, 1e-8);
  }
}

TEST(SurfaceIPCFriction, InvarianceAndReferenceIsSeparateFromLaggedGeometry)
{
  const ES::VXd reference = pointTriangle();
  ES::VXd lagged = reference;
  lagged[0] += .03;
  SurfaceIPCFriction friction;
  friction.update(reference, lagged, { { 0, 1, 2, 3, 2.0 } }, {}, .2, 3, 0, .3, .1, .1);
  EXPECT_GT(gradient(friction, lagged).norm(), 1.0);
  const double energy = friction.computeEnergy(lagged);
  ES::VXd translated = lagged;
  for (int i = 0; i < 4; ++i)
    translated.segment<3>(3 * i) += ES::V3d(.1, -.2, .3);
  EXPECT_NEAR(friction.computeEnergy(translated), energy, 1e-11);
  translated = lagged;
  translated[2] += .07;
  EXPECT_NEAR(friction.computeEnergy(translated), energy, 1e-11);
  checkDerivatives(friction, translated);
  friction.update(lagged, lagged, { { 0, 1, 2, 3, 2.0 } }, {}, .2, 3, 0, .3, .1, .1);
  EXPECT_EQ(gradient(friction, lagged).norm(), 0);
}

TEST(SurfaceIPCFriction, EdgeMollifierScalesNormalLoad)
{
  ES::MXd V(4, 3);
  V << 0, 0, 0, 1, 0, 0, .3, -.01, .05, .8, .02, .05;
  const ES::VXd reference = flattenPositions(V);
  SurfaceIPCFriction plain, mollified;
  plain.update(reference, reference, {}, { { 0, 1, 2, 3, .7 } }, .2, 3, 0, .3, .1, .1);
  mollified.update(reference, reference, {}, { { 0, 1, 2, 3, .7 } }, .2, 3, .01, .3, .1, .1);
  const double m = distance::eeMollifier(V.row(0), V.row(1), V.row(2), V.row(3), .01);
  EXPECT_GT(m, 0);
  EXPECT_LT(m, 1);
  ES::VXd x = reference;
  x[0] += .02;
  EXPECT_NEAR(mollified.computeEnergy(x), m * plain.computeEnergy(x), 1e-12);
  EXPECT_LT((gradient(mollified, x) - m * gradient(plain, x)).norm(), 1e-10);
}

TEST(SurfaceIPCFriction, CrossingEdgesDistributeCoulombForceAtTheirMidpoints)
{
  ES::MXd V(4, 3);
  V << -1, 0, 0, 1, 0, 0, 0, -1, .05, 0, 1, .05;
  const ES::VXd reference = flattenPositions(V);
  SurfaceIPCFriction friction;
  friction.update(reference, reference, {}, { { 0, 1, 2, 3, .7 } }, .2, 3, 0, .3, .1, .1);
  const double muLambda = .3 * (-2 * .05 * .7 * 3 * barrier::dbds(.05 * .05, .2 * .2));
  ES::VXd x = reference;
  x[0] += .02;
  x[3] += .02;
  const ES::VXd g = gradient(friction, x);
  for (int i = 0; i < 4; ++i) {
    EXPECT_NEAR(g[3 * i], (i < 2 ? .5 : -.5) * muLambda, 1e-10);
    EXPECT_NEAR(g[3 * i + 1], 0, 1e-12);
    EXPECT_NEAR(g[3 * i + 2], 0, 1e-12);
  }
  ES::SpMatD H;
  friction.computeHessian(reference, H);
  EXPECT_NEAR(H.coeff(0, 0), .25 * 2 * muLambda / .01, 1e-8);
}

TEST(SurfaceIPCFriction, CoreFrozenStateCopyInvalidationAndZeroCoefficient)
{
  const auto [V, F] = makeTwoTriangleMesh();
  const ES::VXd reference = flattenPositions(V);
  SurfaceIPCCore::Parameters params;
  params.dhat = .2;
  params.frictionCoeff = .3;
  params.projectHessianToPSD = false;
  SurfaceIPCCore core(params);
  core.setMesh(V, F);
  EXPECT_THROW(core.computeEnergy(reference), std::logic_error);
  core.updateFriction(reference, reference, .01);
  EXPECT_GT(core.getNumFrictionPairs(), 0);
  ES::VXd x = reference;
  x[9] += .001;
  const double energy = core.computeEnergy(x);
  core.invalidatePreparedState();
  EXPECT_DOUBLE_EQ(core.computeEnergy(x), energy);
  SurfaceIPCCore copy(core), assigned;
  assigned = core;
  EXPECT_DOUBLE_EQ(copy.computeEnergy(x), energy);
  EXPECT_DOUBLE_EQ(assigned.computeEnergy(x), energy);
  ES::VXd g;
  ES::SpMatD H, H2;
  double combined;
  core.computeAll(x, combined, g, H);
  EXPECT_DOUBLE_EQ(combined, energy);
  const ES::VXd fdg = finiteDifferenceGradient([&](const ES::VXd &q) { return core.computeEnergy(q); }, x, 1e-7);
  EXPECT_LT(relativeError(g, fdg), 1e-6);
  const ES::MXd fdh = finiteDifferenceHessian([&](const ES::VXd &q) {
    ES::VXd out(q.size());
    core.computeGradient(q, out);
    return out;
  },
    x, 1e-7);
  EXPECT_LT(relativeError(ES::MXd(H), fdh), 1e-5);
  core.computeHessian(x, H2);
  EXPECT_LT((H - H2).norm(), 1e-10);
  core.setParameters(params);
  EXPECT_THROW(core.computeEnergy(x), std::logic_error);
  core.updateFriction(reference, reference, .01);
  core.setMesh(V, F);
  EXPECT_THROW(core.computeEnergy(x), std::logic_error);
  params.frictionCoeff = 0;
  core.setParameters(params);
  SurfaceIPCCore baseline;
  params.frictionCoeff = 0;
  baseline.setParameters(params);
  baseline.setMesh(V, F);
  EXPECT_DOUBLE_EQ(core.computeEnergy(x), baseline.computeEnergy(x));
  EXPECT_DOUBLE_EQ(core.computeMaxStepSize(x, reference - x), baseline.computeMaxStepSize(x, reference - x));
}

TEST(SurfaceIPCFriction, EmbeddedExternalSurfaceAndNonidentityMapDerivatives)
{
  const auto [V, F] = makeTwoTriangleMesh();
  ES::SpMatD W(18, 9);
  std::vector<ES::TripletD> triplets;
  for (int i = 0; i < 9; ++i) {
    triplets.emplace_back(i, i, .8);
    triplets.emplace_back(i, (i + 3) % 9, .2);
  }
  W.setFromTriplets(triplets.begin(), triplets.end());
  SurfaceIPCCore::Parameters params;
  params.dhat = .2;
  params.frictionCoeff = .3;
  params.projectHessianToPSD = false;
  EmbeddedSurfaceIPCPotentialEnergy energy(V, F, W, { 1, 1, 1, 0, 0, 0 }, params);
  const ES::VXd reference = ES::VXd::Zero(9);
  ES::VXd x = reference;
  x[0] = .001;
  energy.updateFriction(reference, reference, .01);
  EXPECT_GT(energy.getNumFrictionPairs(), 0);
  ES::VXd g(9);
  energy.gradient(x, g);
  EXPECT_LT(relativeError(g, finiteDifferenceGradient([&](const ES::VXd &q) { return energy.func(q); }, x, 1e-7)), 1e-6);
  ES::SpMatD H;
  energy.hessianDirect(x, H);
  const ES::MXd fd = finiteDifferenceHessian([&](const ES::VXd &q) { ES::VXd out(9); energy.gradient(q, out); return out; }, x, 1e-7);
  EXPECT_LT(relativeError(ES::MXd(H), fd), 1e-5);
}

TEST(SurfaceIPCFriction, RejectsInvalidParametersAndSeparations)
{
  SurfaceIPCCore::Parameters params;
  params.frictionCoeff = -.1;
  EXPECT_THROW(SurfaceIPCCore{ params }, std::invalid_argument);
  params.frictionCoeff = std::numeric_limits<double>::infinity();
  EXPECT_THROW(SurfaceIPCCore{ params }, std::invalid_argument);
  params.frictionCoeff = .3;
  params.frictionEpsV = 0;
  EXPECT_THROW(SurfaceIPCCore{ params }, std::invalid_argument);
  const ES::VXd reference = pointTriangle();
  SurfaceIPCFriction friction;
  EXPECT_THROW(friction.update(reference, reference, {}, {}, .2, 3, 0, .3, .1, 0), std::invalid_argument);
  ES::VXd intersecting = reference;
  intersecting[2] = 0;
  EXPECT_THROW(friction.update(reference, intersecting, { { 0, 1, 2, 3, 1 } }, {}, .2, 3, 0, .3, .1, .1), std::runtime_error);
}
