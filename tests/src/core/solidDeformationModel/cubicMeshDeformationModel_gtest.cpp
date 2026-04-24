#include <gtest/gtest.h>

#include "cubicMeshDeformationModel.h"

#include "elasticModelHillTypeMaterial.h"
#include "elasticModelStableNeoHookeanMaterial.h"
#include "plasticModel3D6DOF.h"
#include "plasticModel3DConstant.h"

#include "EigenSupport.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <iostream>
#include <memory>
#include <vector>

namespace
{
namespace ES = pgo::EigenSupport;
using pgo::SolidDeformationModel::CubicMeshDeformationModel;
using pgo::SolidDeformationModel::ElasticModelHillTypeMaterial;
using pgo::SolidDeformationModel::ElasticModelStableNeoHookeanMaterial;
using pgo::SolidDeformationModel::PlasticModel3D6DOF;
using pgo::SolidDeformationModel::PlasticModel3DConstant;

constexpr int kNumVertices = 8;
constexpr int kNumDofs = 24;
constexpr int kNumGaussPoints = 8;

std::array<double, kNumDofs> makeUnitCubeRestPositions()
{
  return {
    0.0, 0.0, 0.0,
    1.0, 0.0, 0.0,
    1.0, 1.0, 0.0,
    0.0, 1.0, 0.0,
    0.0, 0.0, 1.0,
    1.0, 0.0, 1.0,
    1.0, 1.0, 1.0,
    0.0, 1.0, 1.0,
  };
}

ES::V24d makePerturbedState(const std::array<double, kNumDofs> &rest)
{
  ES::V24d x = Eigen::Map<const ES::V24d>(rest.data());
  for (int vi = 0; vi < kNumVertices; vi++) {
    x[vi * 3 + 0] += 0.015 * ((vi % 3) - 1);
    x[vi * 3 + 1] += 0.010 * ((vi % 5) - 2);
    x[vi * 3 + 2] += 0.012 * ((vi % 7) - 3);
  }
  return x;
}

template<typename Model>
double computeEnergy(Model &model, const double *x, const double *plasticParam, const double *materialParam)
{
  std::unique_ptr<pgo::SolidDeformationModel::DeformationModel::CacheData, std::function<void(pgo::SolidDeformationModel::DeformationModel::CacheData *)>> cache(
    model.allocateCacheData(),
    [&](pgo::SolidDeformationModel::DeformationModel::CacheData *ptr) { model.freeCacheData(ptr); });

  model.prepareData(x, plasticParam, materialParam, cache.get());
  return model.computeEnergy(cache.get());
}

template<typename Model>
ES::V24d computeGradient(Model &model, const double *x, const double *plasticParam, const double *materialParam)
{
  std::unique_ptr<pgo::SolidDeformationModel::DeformationModel::CacheData, std::function<void(pgo::SolidDeformationModel::DeformationModel::CacheData *)>> cache(
    model.allocateCacheData(),
    [&](pgo::SolidDeformationModel::DeformationModel::CacheData *ptr) { model.freeCacheData(ptr); });

  model.prepareData(x, plasticParam, materialParam, cache.get());
  ES::V24d grad = ES::V24d::Zero();
  model.compute_dE_dx(cache.get(), grad.data());
  return grad;
}

template<typename Model>
ES::M24d computeHessian(Model &model, const double *x, const double *plasticParam, const double *materialParam)
{
  std::unique_ptr<pgo::SolidDeformationModel::DeformationModel::CacheData, std::function<void(pgo::SolidDeformationModel::DeformationModel::CacheData *)>> cache(
    model.allocateCacheData(),
    [&](pgo::SolidDeformationModel::DeformationModel::CacheData *ptr) { model.freeCacheData(ptr); });

  model.prepareData(x, plasticParam, materialParam, cache.get());
  ES::M24d hess = ES::M24d::Zero();
  model.compute_d2E_dx2(cache.get(), hess.data());
  return hess;
}

template<typename Model>
Eigen::Matrix<double, kNumDofs, 6> computeMixedDxDa(Model &model, const double *x, const double *plasticParam, const double *materialParam)
{
  std::unique_ptr<pgo::SolidDeformationModel::DeformationModel::CacheData, std::function<void(pgo::SolidDeformationModel::DeformationModel::CacheData *)>> cache(
    model.allocateCacheData(),
    [&](pgo::SolidDeformationModel::DeformationModel::CacheData *ptr) { model.freeCacheData(ptr); });

  model.prepareData(x, plasticParam, materialParam, cache.get());
  Eigen::Matrix<double, kNumDofs, 6> mixed = Eigen::Matrix<double, kNumDofs, 6>::Zero();
  model.compute_d2E_dxda(cache.get(), mixed.data());
  return mixed;
}

template<typename Model>
Eigen::Matrix<double, kNumDofs, 1> computeMixedDxDb(Model &model, const double *x, const double *plasticParam, const double *materialParam)
{
  std::unique_ptr<pgo::SolidDeformationModel::DeformationModel::CacheData, std::function<void(pgo::SolidDeformationModel::DeformationModel::CacheData *)>> cache(
    model.allocateCacheData(),
    [&](pgo::SolidDeformationModel::DeformationModel::CacheData *ptr) { model.freeCacheData(ptr); });

  model.prepareData(x, plasticParam, materialParam, cache.get());
  Eigen::Matrix<double, kNumDofs, 1> mixed = Eigen::Matrix<double, kNumDofs, 1>::Zero();
  model.compute_d2E_dxdb(cache.get(), mixed.data());
  return mixed;
}

struct ErrorSummary
{
  double maxAbs = 0.0;
  double maxRel = 0.0;
};

template<typename DerivedA, typename DerivedB>
ErrorSummary summarizeError(const Eigen::MatrixBase<DerivedA> &actual, const Eigen::MatrixBase<DerivedB> &expected)
{
  ErrorSummary err;
  for (Eigen::Index i = 0; i < actual.size(); i++) {
    const double absErr = std::abs(actual.derived().data()[i] - expected.derived().data()[i]);
    const double relErr = absErr / std::max(1e-9, std::abs(expected.derived().data()[i]));
    err.maxAbs = std::max(err.maxAbs, absErr);
    err.maxRel = std::max(err.maxRel, relErr);
  }
  return err;
}

void expectFinite(const ES::V24d &v)
{
  for (int i = 0; i < v.size(); i++) {
    EXPECT_TRUE(std::isfinite(v[i])) << "Non-finite value at " << i;
  }
}

void expectFinite(const ES::M24d &m)
{
  for (int i = 0; i < m.size(); i++) {
    EXPECT_TRUE(std::isfinite(m.data()[i])) << "Non-finite value at " << i;
  }
}

void expectFinite(const Eigen::VectorXd &v)
{
  for (int i = 0; i < v.size(); i++) {
    EXPECT_TRUE(std::isfinite(v[i])) << "Non-finite value at " << i;
  }
}
}  // namespace

TEST(CubicMeshDeformationModelGTest, RestStateGeometryAndStressAreFinite)
{
  const auto rest = makeUnitCubeRestPositions();
  ElasticModelStableNeoHookeanMaterial elastic(1200.0, 1800.0);
  PlasticModel3D6DOF plastic;
  CubicMeshDeformationModel model(rest.data(), &elastic, &plastic);

  const ES::M3d identity = ES::M3d::Identity();
  double plasticParam[6];
  plastic.toParam(identity.data(), plasticParam);

  std::unique_ptr<pgo::SolidDeformationModel::DeformationModel::CacheData, std::function<void(pgo::SolidDeformationModel::DeformationModel::CacheData *)>> cache(
    model.allocateCacheData(),
    [&](pgo::SolidDeformationModel::DeformationModel::CacheData *ptr) { model.freeCacheData(ptr); });
  model.prepareData(rest.data(), plasticParam, nullptr, cache.get());

  double totalWeight = 0.0;
  for (int q = 0; q < kNumGaussPoints; q++) {
    ES::M3d F = ES::M3d::Zero();
    model.computeFe(cache.get(), q, F.data());
    EXPECT_NEAR((F - ES::M3d::Identity()).norm(), 0.0, 1e-12);

    const double weightDetJ = model.getWeightDetJ(q);
    EXPECT_GT(weightDetJ, 0.0);
    totalWeight += weightDetJ;
  }
  EXPECT_NEAR(totalWeight, 1.0, 1e-12);

  const double restEnergy = model.computeEnergy(cache.get());
  EXPECT_TRUE(std::isfinite(restEnergy));
  EXPECT_GE(restEnergy, -1e-10);
  EXPECT_NEAR(restEnergy, 0.0, 1e-10);

  int nPt = 0;
  double stresses[kNumGaussPoints];
  double strains[kNumGaussPoints];
  model.vonMisesStress(cache.get(), nPt, stresses);
  EXPECT_EQ(nPt, kNumGaussPoints);
  model.maxStrain(cache.get(), nPt, strains);
  EXPECT_EQ(nPt, kNumGaussPoints);
  for (int q = 0; q < kNumGaussPoints; q++) {
    EXPECT_TRUE(std::isfinite(stresses[q]));
    EXPECT_TRUE(std::isfinite(strains[q]));
  }

  ES::V24d x = makePerturbedState(rest);
  model.prepareData(x.data(), plasticParam, nullptr, cache.get());
  const double perturbedEnergy = model.computeEnergy(cache.get());
  EXPECT_TRUE(std::isfinite(perturbedEnergy));

  model.vonMisesStress(cache.get(), nPt, stresses);
  model.maxStrain(cache.get(), nPt, strains);
  EXPECT_EQ(nPt, kNumGaussPoints);
  for (int q = 0; q < kNumGaussPoints; q++) {
    EXPECT_TRUE(std::isfinite(stresses[q]));
    EXPECT_TRUE(std::isfinite(strains[q]));
  }
}

TEST(CubicMeshDeformationModelGTest, PositionAndPlasticDerivativesMatchFiniteDifferences)
{
  const auto rest = makeUnitCubeRestPositions();
  ElasticModelStableNeoHookeanMaterial elastic(1200.0, 1800.0);
  PlasticModel3D6DOF plastic;
  CubicMeshDeformationModel model(rest.data(), &elastic, &plastic);

  ES::V24d x = makePerturbedState(rest);
  Eigen::Matrix<double, 6, 1> plasticParam;
  plasticParam << 1.05, 0.02, -0.01, 0.96, 0.015, 1.02;

  const double epsX = 1e-6;
  const double epsA = 1e-6;

  const ES::V24d analyticGrad = computeGradient(model, x.data(), plasticParam.data(), nullptr);
  const ES::M24d analyticHess = computeHessian(model, x.data(), plasticParam.data(), nullptr);
  const Eigen::Matrix<double, kNumDofs, 6> analyticDxDa = computeMixedDxDa(model, x.data(), plasticParam.data(), nullptr);

  ES::V24d fdGrad = ES::V24d::Zero();
  for (int i = 0; i < kNumDofs; i++) {
    ES::V24d xPlus = x;
    ES::V24d xMinus = x;
    xPlus[i] += epsX;
    xMinus[i] -= epsX;
    fdGrad[i] = (computeEnergy(model, xPlus.data(), plasticParam.data(), nullptr) -
      computeEnergy(model, xMinus.data(), plasticParam.data(), nullptr)) / (2.0 * epsX);
  }

  ES::M24d fdHess = ES::M24d::Zero();
  for (int i = 0; i < kNumDofs; i++) {
    ES::V24d xPlus = x;
    ES::V24d xMinus = x;
    xPlus[i] += epsX;
    xMinus[i] -= epsX;
    fdHess.col(i) = (computeGradient(model, xPlus.data(), plasticParam.data(), nullptr) -
      computeGradient(model, xMinus.data(), plasticParam.data(), nullptr)) / (2.0 * epsX);
  }

  Eigen::Matrix<double, kNumDofs, 6> fdDxDa = Eigen::Matrix<double, kNumDofs, 6>::Zero();
  for (int i = 0; i < 6; i++) {
    Eigen::Matrix<double, 6, 1> aPlus = plasticParam;
    Eigen::Matrix<double, 6, 1> aMinus = plasticParam;
    aPlus[i] += epsA;
    aMinus[i] -= epsA;
    fdDxDa.col(i) = (computeGradient(model, x.data(), aPlus.data(), nullptr) -
      computeGradient(model, x.data(), aMinus.data(), nullptr)) / (2.0 * epsA);
  }

  const ErrorSummary gradErr = summarizeError(analyticGrad, fdGrad);
  const ErrorSummary hessErr = summarizeError(analyticHess, fdHess);
  const ErrorSummary mixedErr = summarizeError(analyticDxDa, fdDxDa);
  std::cout << "dx grad abs=" << gradErr.maxAbs << " rel=" << gradErr.maxRel
            << " hess abs=" << hessErr.maxAbs << " rel=" << hessErr.maxRel
            << " dxda abs=" << mixedErr.maxAbs << " rel=" << mixedErr.maxRel << std::endl;

  expectFinite(analyticGrad);
  expectFinite(analyticHess);
  EXPECT_LT((analyticHess - analyticHess.transpose()).norm(), 1e-7);
  EXPECT_LT(gradErr.maxAbs, 2e-5);
  EXPECT_LT(gradErr.maxRel, 2e-4);
  EXPECT_LT(hessErr.maxAbs, 3e-4);
  EXPECT_LT(hessErr.maxRel, 3e-3);
  EXPECT_LT(mixedErr.maxAbs, 3e-4);
  EXPECT_LT(mixedErr.maxRel, 5e-3);
}

TEST(CubicMeshDeformationModelGTest, MaterialDerivativeMatchesFiniteDifferences)
{
  const auto rest = makeUnitCubeRestPositions();
  const ES::V24d x = makePerturbedState(rest);
  const double fiberDirection[3] = { 1.0, 0.0, 0.0 };

  ElasticModelHillTypeMaterial elastic(0.35, 2.0, 1.0, fiberDirection);
  const ES::M3d identity = ES::M3d::Identity();
  PlasticModel3DConstant plastic(identity.data());
  CubicMeshDeformationModel model(rest.data(), &elastic, &plastic);

  const double materialParam[1] = { 0.75 };
  const double epsB = 1e-6;

  const Eigen::Matrix<double, kNumDofs, 1> analyticDxDb = computeMixedDxDb(model, x.data(), nullptr, materialParam);

  double bPlus[1] = { materialParam[0] + epsB };
  double bMinus[1] = { materialParam[0] - epsB };
  const Eigen::Matrix<double, kNumDofs, 1> fdDxDb =
    (computeGradient(model, x.data(), nullptr, bPlus) - computeGradient(model, x.data(), nullptr, bMinus)) / (2.0 * epsB);

  const ErrorSummary mixedErr = summarizeError(analyticDxDb, fdDxDb);
  std::cout << "dxdb abs=" << mixedErr.maxAbs << " rel=" << mixedErr.maxRel << std::endl;

  EXPECT_LT(mixedErr.maxAbs, 2e-5);
  EXPECT_LT(mixedErr.maxRel, 2e-4);
}

TEST(CubicMeshDeformationModelGTest, MultipleCachesRemainStableAcrossRepeatedPrepareDataCalls)
{
  const auto rest = makeUnitCubeRestPositions();
  const double fiberDirection[3] = { 1.0, 0.0, 0.0 };

  ElasticModelHillTypeMaterial elastic(0.35, 2.0, 1.0, fiberDirection);
  PlasticModel3D6DOF plastic;
  CubicMeshDeformationModel model(rest.data(), &elastic, &plastic);

  ES::V24d x1 = makePerturbedState(rest);
  ES::V24d x2 = x1;
  x2[1] += 0.008;
  x2[7] -= 0.006;
  x2[18] += 0.004;

  Eigen::Matrix<double, 6, 1> plasticParam1;
  plasticParam1 << 1.02, 0.01, -0.02, 0.97, 0.015, 1.01;
  Eigen::Matrix<double, 6, 1> plasticParam2;
  plasticParam2 << 1.06, -0.015, 0.025, 0.95, -0.01, 1.03;

  const double materialParam1[1] = { 0.70 };
  const double materialParam2[1] = { 0.92 };

  std::unique_ptr<pgo::SolidDeformationModel::DeformationModel::CacheData, std::function<void(pgo::SolidDeformationModel::DeformationModel::CacheData *)>> cacheA(
    model.allocateCacheData(),
    [&](pgo::SolidDeformationModel::DeformationModel::CacheData *ptr) { model.freeCacheData(ptr); });
  std::unique_ptr<pgo::SolidDeformationModel::DeformationModel::CacheData, std::function<void(pgo::SolidDeformationModel::DeformationModel::CacheData *)>> cacheB(
    model.allocateCacheData(),
    [&](pgo::SolidDeformationModel::DeformationModel::CacheData *ptr) { model.freeCacheData(ptr); });

  model.prepareData(x1.data(), plasticParam1.data(), materialParam1, cacheA.get());
  const double energyA1 = model.computeEnergy(cacheA.get());
  ES::V24d gradA1 = ES::V24d::Zero();
  model.compute_dE_dx(cacheA.get(), gradA1.data());

  model.prepareData(x2.data(), plasticParam2.data(), materialParam2, cacheB.get());
  const double energyB2 = model.computeEnergy(cacheB.get());
  ES::V24d gradB2 = ES::V24d::Zero();
  model.compute_dE_dx(cacheB.get(), gradB2.data());

  model.prepareData(x2.data(), plasticParam2.data(), materialParam2, cacheA.get());
  const double energyA2 = model.computeEnergy(cacheA.get());
  ES::V24d gradA2 = ES::V24d::Zero();
  model.compute_dE_dx(cacheA.get(), gradA2.data());

  model.prepareData(x1.data(), plasticParam1.data(), materialParam1, cacheB.get());
  const double energyB1 = model.computeEnergy(cacheB.get());
  ES::V24d gradB1 = ES::V24d::Zero();
  model.compute_dE_dx(cacheB.get(), gradB1.data());

  expectFinite(gradA1);
  expectFinite(gradA2);
  expectFinite(gradB1);
  expectFinite(gradB2);
  EXPECT_TRUE(std::isfinite(energyA1));
  EXPECT_TRUE(std::isfinite(energyA2));
  EXPECT_TRUE(std::isfinite(energyB1));
  EXPECT_TRUE(std::isfinite(energyB2));

  EXPECT_NEAR(energyA2, energyB2, 1e-12);
  EXPECT_NEAR(energyA1, energyB1, 1e-12);
  EXPECT_LT((gradA2 - gradB2).norm(), 1e-12);
  EXPECT_LT((gradA1 - gradB1).norm(), 1e-12);
}
