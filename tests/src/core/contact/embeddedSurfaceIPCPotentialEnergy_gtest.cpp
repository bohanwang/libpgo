#include <gtest/gtest.h>

#include "embeddedSurfaceIPCPotentialEnergy.h"
#include "ipc/core/surfaceIPCCore.h"
#include "potentialEnergies.h"
#include "testCIPCHelpers.h"

#include <cmath>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace
{
namespace ES = pgo::EigenSupport;
using pgo::Contact::CIPC::EmbeddedSurfaceIPCPotentialEnergy;
using pgo::Contact::CIPC::SurfaceIPCCore;
using pgo::Contact::CIPCTest::flattenPositions;
using pgo::Contact::CIPCTest::finiteDifferenceHessian;
using pgo::Contact::CIPCTest::makeTwoTriangleMesh;
using pgo::Contact::CIPCTest::relativeError;
using pgo::Contact::CIPCTest::sparseToDense;
using pgo::NonlinearOptimization::PotentialEnergies;

class FixedDiagonalEnergy final : public pgo::NonlinearOptimization::PotentialEnergy
{
public:
  FixedDiagonalEnergy(int numDOFs, double diagonal):
    numDOFs_(numDOFs), diagonal_(diagonal), dofs_(numDOFs)
  {
    std::iota(dofs_.begin(), dofs_.end(), 0);
  }

  double func(ES::ConstRefVecXd) const override { return 0.0; }
  void gradient(ES::ConstRefVecXd, ES::RefVecXd gradient) const override { gradient.setZero(); }
  void hessian(ES::ConstRefVecXd, ES::SpMatD &hessian) const override
  {
    for (int i = 0; i < numDOFs_; ++i)
      hessian.coeffRef(i, i) = diagonal_;
  }
  void createHessian(ES::SpMatD &hessian) const override
  {
    hessian.resize(numDOFs_, numDOFs_);
    hessian.reserve(numDOFs_);
    for (int i = 0; i < numDOFs_; ++i)
      hessian.insert(i, i) = 0.0;
    hessian.makeCompressed();
  }
  void getDOFs(std::vector<int> &dofs) const override { dofs = dofs_; }
  int getNumDOFs() const override { return numDOFs_; }

private:
  int numDOFs_;
  double diagonal_;
  std::vector<int> dofs_;
};

class DynamicDiagonalEnergy final : public pgo::NonlinearOptimization::PotentialEnergy
{
public:
  explicit DynamicDiagonalEnergy(int numDOFs):
    numDOFs_(numDOFs), dofs_(numDOFs)
  {
    std::iota(dofs_.begin(), dofs_.end(), 0);
  }

  double func(ES::ConstRefVecXd) const override { return 0.0; }
  void gradient(ES::ConstRefVecXd, ES::RefVecXd gradient) const override { gradient.setZero(); }
  void hessian(ES::ConstRefVecXd, ES::SpMatD &) const override
  {
    throw std::runtime_error("DynamicDiagonalEnergy::hessian() must not be called.");
  }
  void createHessian(ES::SpMatD &) const override
  {
    throw std::runtime_error("DynamicDiagonalEnergy::createHessian() must not be called.");
  }
  void hessianDirect(ES::ConstRefVecXd x, ES::SpMatD &hessian) const override
  {
    ++hessianDirectCalls_;
    hessian.resize(numDOFs_, numDOFs_);
    hessian.insert(x[0] < 0.0 ? 0 : numDOFs_ - 1, x[0] < 0.0 ? 0 : numDOFs_ - 1) = 7.0;
    hessian.makeCompressed();
  }
  void getDOFs(std::vector<int> &dofs) const override { dofs = dofs_; }
  int getNumDOFs() const override { return numDOFs_; }
  int isHessianTopologyFixed() const override { return 0; }

  int hessianDirectCalls() const { return hessianDirectCalls_; }

private:
  int numDOFs_;
  std::vector<int> dofs_;
  mutable int hessianDirectCalls_ = 0;
};

class StepLimitedEnergy final : public pgo::NonlinearOptimization::PotentialEnergy
{
public:
  explicit StepLimitedEnergy(double maxStep):
    maxStep_(maxStep)
  {
  }

  double func(ES::ConstRefVecXd) const override { return 0.0; }
  void gradient(ES::ConstRefVecXd, ES::RefVecXd gradient) const override { gradient.setZero(); }
  void hessian(ES::ConstRefVecXd, ES::SpMatD &) const override {}
  void createHessian(ES::SpMatD &hessian) const override { hessian = ES::SpMatD(1, 1); }
  void getDOFs(std::vector<int> &dofs) const override { dofs = { 0 }; }
  int getNumDOFs() const override { return 1; }
  double computeMaxStepSize(ES::ConstRefVecXd, ES::ConstRefVecXd) const override { return maxStep_; }

private:
  double maxStep_;
};
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

ES::SpMatD makeCoupledEmbedding(int n3)
{
  std::vector<ES::TripletD> triplets;
  triplets.reserve(n3 + 1);
  for (int i = 0; i < n3; ++i) {
    if (i != 3)
      triplets.emplace_back(i, i, 1.0);
  }
  triplets.emplace_back(3, 0, 0.25);
  triplets.emplace_back(3, 3, 0.75);

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

TEST(EmbeddedSurfaceIPCPotentialEnergyGTest, GradientMatchesEnergyDirectionalFiniteDifference)
{
  const auto [V, F] = makeTwoTriangleMesh();
  const ES::VXd rest = flattenPositions(V);
  ES::VXd u = ES::VXd::Zero(rest.size());
  for (int vi = 3; vi < 6; ++vi)
    u[3 * vi + 2] = 0.01;

  ES::VXd direction(rest.size());
  for (Eigen::Index i = 0; i < direction.size(); ++i)
    direction[i] = 0.1 + 0.03 * static_cast<double>((5 * i) % 11);
  direction.normalize();

  EmbeddedSurfaceIPCPotentialEnergy adapter(V, F, makeIdentityEmbedding(rest.size()), makeParams());
  ES::VXd gradient = ES::VXd::Zero(adapter.getNumDOFs());
  adapter.gradient(u, gradient);

  constexpr double kEpsilon = 1e-6;
  const double fdDirectionalGradient =
    (adapter.func(u + kEpsilon * direction) - adapter.func(u - kEpsilon * direction)) / (2.0 * kEpsilon);
  EXPECT_NEAR(gradient.dot(direction), fdDirectionalGradient,
    2e-5 * std::max(1.0, std::abs(fdDirectionalGradient)));
}

TEST(EmbeddedSurfaceIPCPotentialEnergyGTest, UnprojectedHessianMatchesGradientFiniteDifference)
{
  const auto [V, F] = makeTwoTriangleMesh();
  const ES::VXd rest = flattenPositions(V);
  ES::VXd u = ES::VXd::Zero(rest.size());
  for (int vi = 3; vi < 6; ++vi)
    u[3 * vi + 2] = 0.01;

  auto params = makeParams();
  params.projectHessianToPSD = false;
  EmbeddedSurfaceIPCPotentialEnergy adapter(V, F, makeIdentityEmbedding(rest.size()), params);

  ES::SpMatD hessian;
  adapter.hessianDirect(u, hessian);
  const ES::MXd analyticHessian = sparseToDense(hessian);
  const ES::MXd fdHessian = finiteDifferenceHessian(
    [&adapter](const ES::VXd &state) {
      ES::VXd gradient = ES::VXd::Zero(adapter.getNumDOFs());
      adapter.gradient(state, gradient);
      return gradient;
    },
    u,
    1e-5);

  EXPECT_LT(relativeError(analyticHessian, analyticHessian.transpose()), 1e-12);
  EXPECT_LT(relativeError(analyticHessian, fdHessian), 2e-5);
}

TEST(EmbeddedSurfaceIPCPotentialEnergyGTest, UnprojectedHessianWithCoupledEmbeddingMatchesGradientFiniteDifference)
{
  const auto [V, F] = makeTwoTriangleMesh();
  const ES::VXd rest = flattenPositions(V);
  ES::VXd u = ES::VXd::Zero(rest.size());
  u[11] = 0.01;
  u[14] = 0.02;
  u[17] = 0.015;

  auto params = makeParams();
  params.projectHessianToPSD = false;
  EmbeddedSurfaceIPCPotentialEnergy adapter(V, F, makeCoupledEmbedding(rest.size()), params);

  ES::SpMatD hessian;
  adapter.hessianDirect(u, hessian);
  const ES::MXd analyticHessian = sparseToDense(hessian);
  const ES::MXd fdHessian = finiteDifferenceHessian(
    [&adapter](const ES::VXd &state) {
      ES::VXd gradient = ES::VXd::Zero(adapter.getNumDOFs());
      adapter.gradient(state, gradient);
      return gradient;
    },
    u,
    1e-5);

  EXPECT_LT(relativeError(analyticHessian, analyticHessian.transpose()), 1e-12);
  EXPECT_LT(relativeError(analyticHessian, fdHessian), 2e-5);
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

TEST(EmbeddedSurfaceIPCPotentialEnergyGTest, ReusesPreparedPairsAcrossEnergyGradientHessianForSameState)
{
  const auto [V, F] = makeTwoTriangleMesh();
  const ES::VXd rest = flattenPositions(V);
  ES::VXd u = ES::VXd::Zero(rest.size());
  for (int vi = 3; vi < 6; ++vi)
    u[3 * vi + 2] = 0.01;

  EmbeddedSurfaceIPCPotentialEnergy adapter(V, F, makeIdentityEmbedding(rest.size()), makeParams());

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

  EXPECT_NEAR(energy1, energy0, 1e-12);
  EXPECT_LT(relativeError(gradient1, gradient0), 1e-12);
  EXPECT_LT(relativeError(sparseToDense(hessian1), sparseToDense(hessian0)), 1e-12);
}

TEST(EmbeddedSurfaceIPCPotentialEnergyGTest, DynamicHessianCreateRemainsInvalidButAggregateUsesDirectPath)
{
  const auto [V, F] = makeTwoTriangleMesh();
  const ES::VXd rest = flattenPositions(V);
  ES::VXd u = ES::VXd::Zero(rest.size());
  for (int vi = 3; vi < 6; ++vi)
    u[3 * vi + 2] = 0.01;

  auto adapter = std::make_shared<EmbeddedSurfaceIPCPotentialEnergy>(V, F, makeIdentityEmbedding(rest.size()), makeParams());
  ES::SpMatD unusedPattern;
  EXPECT_THROW(adapter->createHessian(unusedPattern), std::runtime_error);

  PotentialEnergies aggregate(adapter->getNumDOFs());
  aggregate.addPotentialEnergy(adapter);
  EXPECT_NO_THROW(aggregate.init());

  ES::SpMatD aggregateHessian;
  EXPECT_NO_THROW(aggregate.hessianDirect(u, aggregateHessian));

  ES::SpMatD adapterHessian;
  adapter->hessianDirect(u, adapterHessian);
  EXPECT_LT(relativeError(sparseToDense(aggregateHessian), sparseToDense(adapterHessian)), 1e-12);
}

TEST(EmbeddedSurfaceIPCPotentialEnergyGTest, AggregateAddsFixedAndDynamicHessiansOnce)
{
  auto fixed = std::make_shared<FixedDiagonalEnergy>(2, 4.0);
  auto dynamic = std::make_shared<DynamicDiagonalEnergy>(2);

  PotentialEnergies aggregate(2);
  aggregate.addPotentialEnergy(fixed, 2.0);
  aggregate.addPotentialEnergy(dynamic, 3.0);
  EXPECT_NO_THROW(aggregate.init());

  ES::VXd x(2);
  x << -1.0, 0.0;
  ES::SpMatD hessian;
  aggregate.hessianDirect(x, hessian);
  EXPECT_EQ(dynamic->hessianDirectCalls(), 1);
  EXPECT_EQ(hessian.rows(), 2);
  EXPECT_EQ(hessian.cols(), 2);
  EXPECT_TRUE(sparseToDense(hessian).isApprox(sparseToDense(hessian).transpose(), 1e-12));
  EXPECT_TRUE(sparseToDense(hessian).allFinite());
  EXPECT_NEAR(hessian.coeff(0, 0), 29.0, 1e-12);
  EXPECT_NEAR(hessian.coeff(1, 1), 8.0, 1e-12);

  x[0] = 1.0;
  aggregate.hessianDirect(x, hessian);
  EXPECT_EQ(dynamic->hessianDirectCalls(), 2);
  EXPECT_NEAR(hessian.coeff(0, 0), 8.0, 1e-12);
  EXPECT_NEAR(hessian.coeff(1, 1), 29.0, 1e-12);
}

TEST(EmbeddedSurfaceIPCPotentialEnergyGTest, ZeroCoefficientSkipsDynamicHessianEvaluation)
{
  auto fixed = std::make_shared<FixedDiagonalEnergy>(2, 4.0);
  auto dynamic = std::make_shared<DynamicDiagonalEnergy>(2);

  PotentialEnergies aggregate(2);
  aggregate.addPotentialEnergy(fixed);
  aggregate.addPotentialEnergy(dynamic, 0.0);
  EXPECT_NO_THROW(aggregate.init());

  const ES::VXd x = ES::VXd::Zero(2);
  ES::SpMatD hessian;
  aggregate.hessianDirect(x, hessian);

  EXPECT_EQ(dynamic->hessianDirectCalls(), 0);
  EXPECT_NEAR(hessian.coeff(0, 0), 4.0, 1e-12);
  EXPECT_NEAR(hessian.coeff(1, 1), 4.0, 1e-12);
}

TEST(EmbeddedSurfaceIPCPotentialEnergyGTest, AggregateMaxStepUsesMinimumAndSkipsZeroCoefficient)
{
  auto unlimited = std::make_shared<FixedDiagonalEnergy>(1, 1.0);
  auto limited = std::make_shared<StepLimitedEnergy>(0.25);
  PotentialEnergies aggregate(1);
  aggregate.addPotentialEnergy(unlimited);
  aggregate.addPotentialEnergy(limited);
  aggregate.init();

  ES::VXd x(1);
  ES::VXd dx(1);
  x << 0.0;
  dx << 1.0;
  EXPECT_DOUBLE_EQ(aggregate.computeMaxStepSize(x, dx), 0.25);

  aggregate.setEnergyCoeffs(1, 0.0);
  EXPECT_DOUBLE_EQ(aggregate.computeMaxStepSize(x, dx), 1.0);
}

TEST(EmbeddedSurfaceIPCPotentialEnergyGTest, InvalidEmbeddingRowsThrow)
{
  const auto [V, F] = makeTwoTriangleMesh();
  ES::SpMatD W(3 * V.rows() - 1, 3 * V.rows());

  EXPECT_THROW(
    EmbeddedSurfaceIPCPotentialEnergy(V, F, W, makeParams()),
    std::invalid_argument);
}

TEST(EmbeddedSurfaceIPCPotentialEnergyGTest, NonFiniteRestPositionsAndEmbeddingThrow)
{
  const auto [V, F] = makeTwoTriangleMesh();

  ES::MXd nonFiniteVertices = V;
  nonFiniteVertices(0, 0) = std::numeric_limits<double>::infinity();
  EXPECT_THROW(
    EmbeddedSurfaceIPCPotentialEnergy(nonFiniteVertices, F, makeIdentityEmbedding(3 * V.rows()), makeParams()),
    std::invalid_argument);

  ES::SpMatD nonFiniteEmbedding = makeIdentityEmbedding(3 * V.rows());
  nonFiniteEmbedding.coeffRef(0, 0) = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(
    EmbeddedSurfaceIPCPotentialEnergy(V, F, nonFiniteEmbedding, makeParams()),
    std::invalid_argument);

  ES::MXi invalidTriangles = F;
  invalidTriangles(0, 0) = V.rows();
  EXPECT_THROW(
    EmbeddedSurfaceIPCPotentialEnergy(V, invalidTriangles, makeIdentityEmbedding(3 * V.rows()), makeParams()),
    std::invalid_argument);
}
