#include <gtest/gtest.h>

#include "tetMeshDeformationModel.h"

#include "elasticModelHillTypeMaterial.h"
#include "elasticModelStableNeoHookeanMaterial.h"
#include "plasticModel3D6DOF.h"
#include "plasticModel3DConstant.h"

#include "EigenSupport.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <memory>

namespace
{
namespace ES = pgo::EigenSupport;
using pgo::SolidDeformationModel::ElasticModelHillTypeMaterial;
using pgo::SolidDeformationModel::ElasticModelStableNeoHookeanMaterial;
using pgo::SolidDeformationModel::PlasticModel3D6DOF;
using pgo::SolidDeformationModel::PlasticModel3DConstant;
using pgo::SolidDeformationModel::TetMeshDeformationModel;

constexpr int kNumVertices = 4;
constexpr int kNumDOFs = 12;

std::array<double, kNumDOFs> makeRestPositions()
{
  return {
    0.0, 0.0, 0.0,
    1.0, 0.0, 0.0,
    0.0, 1.0, 0.0,
    0.0, 0.0, 1.0,
  };
}

ES::V12d makePerturbedState(const std::array<double, kNumDOFs> &rest)
{
  ES::V12d x = Eigen::Map<const ES::V12d>(rest.data());
  for (int vi = 0; vi < kNumVertices; ++vi) {
    x[3 * vi + 0] += 0.012 * ((vi % 3) - 1);
    x[3 * vi + 1] += 0.009 * ((vi % 5) - 2);
    x[3 * vi + 2] += 0.011 * ((vi % 7) - 3);
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
ES::V12d computeGradient(Model &model, const double *x, const double *plasticParam, const double *materialParam)
{
  std::unique_ptr<pgo::SolidDeformationModel::DeformationModel::CacheData, std::function<void(pgo::SolidDeformationModel::DeformationModel::CacheData *)>> cache(
    model.allocateCacheData(),
    [&](pgo::SolidDeformationModel::DeformationModel::CacheData *ptr) { model.freeCacheData(ptr); });
  model.prepareData(x, plasticParam, materialParam, cache.get());
  ES::V12d gradient = ES::V12d::Zero();
  model.compute_dE_dx(cache.get(), gradient.data());
  return gradient;
}

template<typename Model>
ES::M12d computeHessian(Model &model, const double *x, const double *plasticParam, const double *materialParam)
{
  std::unique_ptr<pgo::SolidDeformationModel::DeformationModel::CacheData, std::function<void(pgo::SolidDeformationModel::DeformationModel::CacheData *)>> cache(
    model.allocateCacheData(),
    [&](pgo::SolidDeformationModel::DeformationModel::CacheData *ptr) { model.freeCacheData(ptr); });
  model.prepareData(x, plasticParam, materialParam, cache.get());
  ES::M12d hessian = ES::M12d::Zero();
  model.compute_d2E_dx2(cache.get(), hessian.data());
  return hessian;
}

template<typename Model>
Eigen::Matrix<double, kNumDOFs, 6> computeMixedDxDa(Model &model, const double *x, const double *plasticParam, const double *materialParam)
{
  std::unique_ptr<pgo::SolidDeformationModel::DeformationModel::CacheData, std::function<void(pgo::SolidDeformationModel::DeformationModel::CacheData *)>> cache(
    model.allocateCacheData(),
    [&](pgo::SolidDeformationModel::DeformationModel::CacheData *ptr) { model.freeCacheData(ptr); });
  model.prepareData(x, plasticParam, materialParam, cache.get());
  Eigen::Matrix<double, kNumDOFs, 6> mixed = Eigen::Matrix<double, kNumDOFs, 6>::Zero();
  model.compute_d2E_dxda(cache.get(), mixed.data());
  return mixed;
}

template<typename Model>
Eigen::Matrix<double, kNumDOFs, 1> computeMixedDxDb(Model &model, const double *x, const double *plasticParam, const double *materialParam)
{
  std::unique_ptr<pgo::SolidDeformationModel::DeformationModel::CacheData, std::function<void(pgo::SolidDeformationModel::DeformationModel::CacheData *)>> cache(
    model.allocateCacheData(),
    [&](pgo::SolidDeformationModel::DeformationModel::CacheData *ptr) { model.freeCacheData(ptr); });
  model.prepareData(x, plasticParam, materialParam, cache.get());
  Eigen::Matrix<double, kNumDOFs, 1> mixed = Eigen::Matrix<double, kNumDOFs, 1>::Zero();
  model.compute_d2E_dxdb(cache.get(), mixed.data());
  return mixed;
}

template<typename DerivedA, typename DerivedB>
void expectClose(const Eigen::MatrixBase<DerivedA> &actual, const Eigen::MatrixBase<DerivedB> &expected,
  double absoluteTolerance, double relativeTolerance)
{
  for (Eigen::Index i = 0; i < actual.size(); ++i) {
    const double error = std::abs(actual.derived().coeff(i) - expected.derived().coeff(i));
    const double scale = std::max(1.0, std::abs(expected.derived().coeff(i)));
    EXPECT_LE(error, absoluteTolerance + relativeTolerance * scale)
      << "entry " << i << ", analytic=" << actual.derived().coeff(i)
      << ", finite difference=" << expected.derived().coeff(i);
  }
}
}  // namespace

TEST(TetMeshDeformationModelGTest, PositionAndPlasticDerivativesMatchFiniteDifferences)
{
  const auto rest = makeRestPositions();
  ElasticModelStableNeoHookeanMaterial elastic(1200.0, 1800.0);
  PlasticModel3D6DOF plastic;
  TetMeshDeformationModel model(rest.data(), rest.data() + 3, rest.data() + 6, rest.data() + 9, &elastic, &plastic);

  const ES::V12d x = makePerturbedState(rest);
  Eigen::Matrix<double, 6, 1> plasticParam;
  plasticParam << 1.05, 0.02, -0.01, 0.96, 0.015, 1.02;

  constexpr double kEpsX = 1e-6;
  constexpr double kEpsA = 1e-6;
  const ES::V12d analyticGradient = computeGradient(model, x.data(), plasticParam.data(), nullptr);
  const ES::M12d analyticHessian = computeHessian(model, x.data(), plasticParam.data(), nullptr);
  const auto analyticDxDa = computeMixedDxDa(model, x.data(), plasticParam.data(), nullptr);

  ES::V12d fdGradient = ES::V12d::Zero();
  ES::M12d fdHessian = ES::M12d::Zero();
  for (int i = 0; i < kNumDOFs; ++i) {
    ES::V12d xPlus = x;
    ES::V12d xMinus = x;
    xPlus[i] += kEpsX;
    xMinus[i] -= kEpsX;
    fdGradient[i] = (computeEnergy(model, xPlus.data(), plasticParam.data(), nullptr) -
      computeEnergy(model, xMinus.data(), plasticParam.data(), nullptr)) / (2.0 * kEpsX);
    fdHessian.col(i) = (computeGradient(model, xPlus.data(), plasticParam.data(), nullptr) -
      computeGradient(model, xMinus.data(), plasticParam.data(), nullptr)) / (2.0 * kEpsX);
  }

  Eigen::Matrix<double, kNumDOFs, 6> fdDxDa = Eigen::Matrix<double, kNumDOFs, 6>::Zero();
  for (int i = 0; i < plasticParam.size(); ++i) {
    auto aPlus = plasticParam;
    auto aMinus = plasticParam;
    aPlus[i] += kEpsA;
    aMinus[i] -= kEpsA;
    fdDxDa.col(i) = (computeGradient(model, x.data(), aPlus.data(), nullptr) -
      computeGradient(model, x.data(), aMinus.data(), nullptr)) / (2.0 * kEpsA);
  }

  expectClose(analyticGradient, fdGradient, 2e-5, 2e-4);
  expectClose(analyticHessian, fdHessian, 3e-4, 3e-3);
  expectClose(analyticDxDa, fdDxDa, 3e-4, 5e-3);
}

TEST(TetMeshDeformationModelGTest, MaterialDerivativeMatchesFiniteDifferences)
{
  const auto rest = makeRestPositions();
  const double fiberDirection[3] = { 1.0, 0.0, 0.0 };
  const ES::M3d identity = ES::M3d::Identity();
  ElasticModelHillTypeMaterial elastic(0.35, 2.0, 1.0, fiberDirection);
  PlasticModel3DConstant plastic(identity.data());
  TetMeshDeformationModel model(rest.data(), rest.data() + 3, rest.data() + 6, rest.data() + 9, &elastic, &plastic);

  const ES::V12d x = makePerturbedState(rest);
  double materialParam = 0.75;
  constexpr double kEpsB = 1e-6;

  const auto analyticDxDb = computeMixedDxDb(model, x.data(), nullptr, &materialParam);
  const double bPlus = materialParam + kEpsB;
  const double bMinus = materialParam - kEpsB;
  const ES::V12d fdDxDb = (computeGradient(model, x.data(), nullptr, &bPlus) -
    computeGradient(model, x.data(), nullptr, &bMinus)) / (2.0 * kEpsB);

  expectClose(analyticDxDb, fdDxDb, 2e-5, 2e-4);
}
