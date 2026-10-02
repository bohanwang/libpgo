#include "cudaKoiterDeformationModelEnergy.h"
#include "cudaDeviceBuffer.h"
#include "cudaFrozenHessianEnergy.h"
#include "deformationModelAssembler.h"
#include "deformationModelEnergy.h"
#include "deformationModelManager.h"
#include "simulationMesh.h"
#include "triMeshGeo.h"
#include "pgoLogging.h"
#include "cudaElasticEquivalence.h"

#include <gtest/gtest.h>
#include <Eigen/Geometry>

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <vector>

namespace
{
namespace ES = pgo::EigenSupport;
using namespace pgo::SolidDeformationModel;
using pgo::CudaUtilities::DeviceBuffer;

class CudaKoiterDeformationModelEnergyGTest : public ::testing::Test
{
protected:
  void SetUp() override
  {
    int count = 0;
    if (cudaGetDeviceCount(&count) != cudaSuccess || count == 0)
      GTEST_SKIP() << "A CUDA device is required.";
    static const bool loggingInitialized = [] {
      pgo::Logging::init();
      return true;
    }();
    (void)loggingInitialized;
  }
};

struct Fixture
{
  std::shared_ptr<SimulationMesh> mesh;
  std::shared_ptr<DeformationModelManager> manager;
  std::shared_ptr<DeformationModelAssembler> assembler;
  std::unique_ptr<DeformationModelEnergy> cpu;
  std::shared_ptr<CudaKoiterDeformationModelEnergy> gpu;
  ES::VXd rest, parameters;
};

Fixture makeFixture(int type, bool zeroRestCurvature = false, double scale = 1)
{
  std::vector<double> vertices;
  std::vector<int> triangles;
  if (type == 0) {
    vertices = { 0, 0, 0, 1.1, 0.1, 0, 0.2, 0.9, 0.03 };
    triangles = { 0, 1, 2 };
  }
  else if (type == 1) {
    // Curved rest hinge with missing neighbors on its four boundary edges.
    vertices = { 0, 0, 0, 1.1, 0.1, 0, 0.2, 0.9, 0.03, 1, 1.1, 0.45 };
    triangles = { 0, 1, 2, 1, 3, 2 };
  }
  else if (type == 2) {
    // On a closed tetrahedral surface each patch repeats its fourth vertex
    // in all three neighbor slots, exercising global diagonal assembly.
    vertices = { 0, 0, 0, 1.1, 0.1, 0, 0.2, 0.9, 0.03, 0.15, 0.2, 0.95 };
    triangles = { 0, 2, 1, 0, 1, 3, 0, 3, 2, 1, 2, 3 };
  }
  else {
    vertices = { 0, 0, 0, 1.1, 0.1, 0, 0.2, 0.9, 0.03,
      0.5, -0.7, 0.2, 1.1, 1.1, -0.15, -0.6, 0.4, 0.3 };
    triangles = { 0, 1, 2, 1, 0, 3, 2, 1, 4, 0, 2, 5 };
  }
  for (double &value : vertices)
    value *= scale;
  pgo::Mesh::TriMeshGeo surface(vertices.size() / 3, vertices.data(), triangles.size() / 3, triangles.data());
  SimulationMeshENuhMaterial material(10000, 0.4, 0.001);
  Fixture fixture;
  fixture.mesh.reset(loadShellMesh(surface, &material));
  fixture.rest = Eigen::Map<const ES::VXd>(vertices.data(), vertices.size());
  fixture.parameters.resize(5 * fixture.mesh->getNumElements());
  for (int e = 0; e < fixture.mesh->getNumElements(); ++e)
    fixture.parameters.segment<5>(5 * e) << 1300 + 90 * e, 0.27 + 0.025 * e, 8700 - 200 * e, 0.18 + 0.015 * e, 0.08 * scale;
  fixture.manager = std::make_shared<DeformationModelManager>();
  fixture.manager->setMesh(fixture.mesh.get());
  fixture.manager->init(DeformationModelPlasticMaterial::SHELL_FF_DOF0, DeformationModelElasticMaterial::KOITER_STVK);
  fixture.manager->setEnforceSPD(1);
  if (zeroRestCurvature)
    fixture.manager->zeroShellRestCurvature();
  fixture.assembler = std::make_shared<DeformationModelAssembler>(fixture.manager);
  fixture.assembler->setDeterministic(true);
  fixture.cpu = std::make_unique<DeformationModelEnergy>(fixture.assembler, &fixture.rest);
  fixture.cpu->setElasticParams(fixture.parameters);
  fixture.gpu = std::make_shared<CudaKoiterDeformationModelEnergy>(*fixture.mesh, fixture.parameters, zeroRestCurvature);
  return fixture;
}

ES::VXd direction(int n)
{
  ES::VXd result(n);
  for (int i = 0; i < n; ++i)
    result[i] = std::sin(0.57 * i + 0.13) + 0.23 * std::cos(0.19 * i + 0.2);
  return result;
}

ES::VXd deformed(const Fixture &fixture, double amplitude = 0.1)
{
  return amplitude * fixture.rest.norm() * direction(fixture.rest.size());
}

void expectVectorNear(const ES::VXd &actual, const ES::VXd &expected, double tolerance)
{
  EXPECT_LE((actual - expected).norm(), tolerance * std::max(1e-16, expected.norm()));
}

void expectParity(Fixture &fixture, const ES::VXd &x)
{
  const int n = x.size();
  DeviceBuffer<double> deviceX, deviceV, deviceOut;
  deviceX.upload(x.data(), n);
  const double expectedEnergy = fixture.cpu->func(x);
  EXPECT_NEAR(fixture.gpu->func(deviceX.data()), expectedEnergy, 3e-10 * std::max(1e-18, std::abs(expectedEnergy)));
  ES::VXd gradient(n), result(n), baseline = 0.013 * direction(n), v = direction(n).reverse().eval();
  fixture.cpu->gradient(x, gradient);
  deviceOut.upload(baseline.data(), n);
  fixture.gpu->addGradient(deviceX.data(), deviceOut.data());
  deviceOut.download(result.data(), n);
  expectVectorNear(result, baseline + gradient, 4e-10);

  ES::SpMatD H;
  fixture.cpu->createHessian(H);
  fixture.cpu->hessian(x, H);
  fixture.gpu->prepareHessian(deviceX.data());
  deviceV.upload(v.data(), n);
  deviceOut.upload(baseline.data(), n);
  fixture.gpu->addHessianVector(deviceV.data(), deviceOut.data());
  deviceOut.download(result.data(), n);
  expectVectorNear(result, baseline + H * v, 2e-8);
  EXPECT_GE(v.dot(result - baseline), -1e-9 * (H * v).norm());
  deviceOut.upload(baseline.data(), n);
  fixture.gpu->addDiagonal(deviceOut.data());
  deviceOut.download(result.data(), n);
  expectVectorNear(result, baseline + ES::VXd(H.diagonal()), 2e-8);
  EXPECT_EQ(fixture.gpu->computeMaxStepSize(deviceX.data(), deviceV.data()), fixture.cpu->computeMaxStepSize(x, v));

  const double h = 2e-7 * fixture.rest.norm();
  const ES::VXd plus = x + h * v, minus = x - h * v;
  deviceX.upload(plus.data(), n);
  const double ep = fixture.gpu->func(deviceX.data());
  deviceX.upload(minus.data(), n);
  const double em = fixture.gpu->func(deviceX.data());
  EXPECT_NEAR((ep - em) / (2 * h), gradient.dot(v), 3e-7 * std::max(1e-14, std::abs(gradient.dot(v))));
}
}  // namespace

TEST_F(CudaKoiterDeformationModelEnergyGTest, BoundaryTriangleMatchesCpuMembraneAndProjection)
{
  auto fixture = makeFixture(0);
  expectParity(fixture, deformed(fixture));
  ES::VXd x = ES::VXd::Zero(fixture.rest.size());
  for (int v = 0; v < x.size() / 3; ++v)
    x.segment<3>(3 * v) = -0.6 * fixture.rest.segment<3>(3 * v);
  expectParity(fixture, x);
}

TEST_F(CudaKoiterDeformationModelEnergyGTest, QuantitativeFullHessianColumnsMatchCpu)
{
  nlohmann::json reports = nlohmann::json::array();
  for (int type : { 0, 1, 2, 3 }) {
    for (bool zeroRestCurvature : { false, true }) {
      auto fixture = makeFixture(type, zeroRestCurvature);
      const ES::VXd x = deformed(fixture);
      auto report = pgo::Tests::ElasticEquivalence::compare(*fixture.cpu, *fixture.manager, *fixture.gpu, x, true);
      report["element_type"] = "koiter_shell";
      report["material"] = "koiter-stvk";
      report["fixture_type"] = type;
      report["zero_rest_curvature"] = zeroRestCurvature;
      report["num_elements"] = fixture.mesh->getNumElements();
      report["projection_scope"] = "full 18x18 six-vertex element Hessian, matching CPU; repeated global vertex IDs accumulated";
      EXPECT_TRUE(report.at("passed").get<bool>()) << report.dump(2);
      reports.push_back(std::move(report));
    }
  }
  RecordProperty("elastic_equivalence_measurements", reports.dump());
}

TEST_F(CudaKoiterDeformationModelEnergyGTest, CurvedHingeMatchesCpuWithSeparateMaterialParameters)
{
  auto fixture = makeFixture(1);
  expectParity(fixture, deformed(fixture));
  expectParity(fixture, deformed(fixture, -0.17));
}

TEST_F(CudaKoiterDeformationModelEnergyGTest, ClosedPatchesAccumulateRepeatedNeighborDiagonal)
{
  auto fixture = makeFixture(2);
  expectParity(fixture, deformed(fixture));
}

TEST_F(CudaKoiterDeformationModelEnergyGTest, InteriorPatchWithSixDistinctVerticesMatchesCpu)
{
  auto fixture = makeFixture(3);
  expectParity(fixture, deformed(fixture, 0.04));
}

TEST_F(CudaKoiterDeformationModelEnergyGTest, ZeroRestCurvatureMatchesCpu)
{
  for (int type : { 1, 2, 3 }) {
    auto fixture = makeFixture(type, true);
    expectParity(fixture, ES::VXd::Zero(fixture.rest.size()));
    expectParity(fixture, deformed(fixture));
  }
}

TEST_F(CudaKoiterDeformationModelEnergyGTest, SmallScaleRotatedTranslatedShellMatchesCpu)
{
  auto fixture = makeFixture(1, false, 1e-5);
  ES::VXd x = deformed(fixture);
  const ES::M3d rotation = Eigen::AngleAxisd(0.71, ES::V3d(1, 2, -1).normalized()).toRotationMatrix();
  for (int v = 0; v < x.size() / 3; ++v)
    x.segment<3>(3 * v) = rotation * (fixture.rest.segment<3>(3 * v) + x.segment<3>(3 * v)) -
      fixture.rest.segment<3>(3 * v) + ES::V3d(2e-4, -1e-4, 3e-4);
  expectParity(fixture, x);
}

TEST_F(CudaKoiterDeformationModelEnergyGTest, ExactHessianMatchesCpuAndGradientDirectionalDerivative)
{
  for (int type : { 1, 2, 3 }) {
    auto fixture = makeFixture(type, true);
    fixture.manager->setEnforceSPD(0);
    const ES::VXd x = deformed(fixture), v = direction(x.size()).reverse().eval();
    const int n = x.size();
    DeviceBuffer<double> deviceX, deviceV, deviceOut(n);
    deviceX.upload(x.data(), n);
    deviceV.upload(v.data(), n);
    fixture.gpu->prepareExactHessian(deviceX.data());
    deviceOut.zero();
    fixture.gpu->addHessianVector(deviceV.data(), deviceOut.data());
    ES::VXd product(n), plusGradient(n), minusGradient(n);
    deviceOut.download(product.data(), n);
    ES::SpMatD H;
    fixture.cpu->createHessian(H);
    fixture.cpu->hessian(x, H);
    expectVectorNear(product, H * v, 3e-10);
    const double step = 1e-6 * fixture.rest.norm();
    const ES::VXd plus = x + step * v, minus = x - step * v;
    deviceX.upload(plus.data(), n);
    deviceOut.zero();
    fixture.gpu->addGradient(deviceX.data(), deviceOut.data());
    deviceOut.download(plusGradient.data(), n);
    deviceX.upload(minus.data(), n);
    deviceOut.zero();
    fixture.gpu->addGradient(deviceX.data(), deviceOut.data());
    deviceOut.download(minusGradient.data(), n);
    expectVectorNear(product, (plusGradient - minusGradient) / (2 * step), 2e-8);
  }
}

TEST_F(CudaKoiterDeformationModelEnergyGTest, PairedEnergyChangeRetainsSmallDecreases)
{
  for (bool zeroRest : { false, true }) {
    auto fixture = makeFixture(1, zeroRest);
    const ES::VXd x = deformed(fixture, 0.08), v = direction(x.size()).reverse().eval();
    const int n = x.size();
    DeviceBuffer<double> base, trial, out(n);
    base.upload(x.data(), n);
    out.zero();
    fixture.gpu->addGradient(base.data(), out.data());
    ES::VXd gradient(n);
    out.download(gradient.data(), n);
    for (double step : { 1e-5, 1e-12, -1e-12 }) {
      const ES::VXd y = x + step * v;
      trial.upload(y.data(), n);
      const auto change = fixture.gpu->energyChange(base.data(), trial.data());
      const auto reverse = fixture.gpu->energyChange(trial.data(), base.data());
      EXPECT_TRUE(std::isfinite(change.value));
      EXPECT_GE(change.roundoffScale, std::abs(change.value) * (1 - 1e-12));
      EXPECT_NEAR(change.value + reverse.value, 0, 2e-12 * (change.roundoffScale + reverse.roundoffScale));
      if (std::abs(step) < 1e-8) {
        const double linear = gradient.dot(y - x);
        EXPECT_NEAR(change.value, linear, 2e-8 * std::abs(linear));
      }
      else {
        const double reference = fixture.cpu->func(y) - fixture.cpu->func(x);
        EXPECT_NEAR(change.value, reference, 1e-8 * std::abs(reference));
      }
    }
    const auto zero = fixture.gpu->energyChange(base.data(), base.data());
    EXPECT_EQ(zero.value, 0);
    EXPECT_EQ(zero.roundoffScale, 0);
  }
}

TEST_F(CudaKoiterDeformationModelEnergyGTest, RestEnergyAndUniformTranslationAreZero)
{
  auto fixture = makeFixture(2);
  ES::VXd x(fixture.rest.size());
  for (int v = 0; v < x.size() / 3; ++v)
    x.segment<3>(3 * v) = ES::V3d(8, -16, 32);
  DeviceBuffer<double> deviceX, out(x.size());
  deviceX.upload(x.data(), x.size());
  EXPECT_EQ(fixture.gpu->func(deviceX.data()), 0);
  out.zero();
  fixture.gpu->addGradient(deviceX.data(), out.data());
  ES::VXd gradient(x.size());
  out.download(gradient.data(), gradient.size());
  EXPECT_EQ(gradient.norm(), 0);
}

TEST_F(CudaKoiterDeformationModelEnergyGTest, FrozenStiffnessDampingUsesCpuProjectedShellHessian)
{
  auto fixture = makeFixture(1);
  const ES::VXd reference = deformed(fixture), x = deformed(fixture, 0.12), v = direction(x.size());
  const int n = x.size();
  ES::SpMatD H;
  fixture.cpu->createHessian(H);
  fixture.cpu->hessian(reference, H);
  pgo::NonlinearOptimization::CudaFrozenHessianEnergy damping(fixture.gpu);
  constexpr double scale = 2.7;
  damping.setReference(reference, scale);
  DeviceBuffer<double> deviceX, deviceV, out(n);
  deviceX.upload(x.data(), n);
  deviceV.upload(v.data(), n);
  EXPECT_NEAR(damping.func(deviceX.data()), 0.5 * scale * (x - reference).dot(H * (x - reference)), 1e-9);
  damping.prepareHessian(deviceX.data());
  out.zero();
  damping.addHessianVector(deviceV.data(), out.data());
  ES::VXd actual(n);
  out.download(actual.data(), n);
  expectVectorNear(actual, scale * H * v, 2e-8);
  out.zero();
  damping.addGradient(deviceX.data(), out.data());
  out.download(actual.data(), n);
  expectVectorNear(actual, scale * H * (x - reference), 2e-8);
}

TEST_F(CudaKoiterDeformationModelEnergyGTest, RejectsUnsupportedParametersAndUnpreparedOperator)
{
  auto fixture = makeFixture(1);
  EXPECT_THROW(CudaKoiterDeformationModelEnergy(*fixture.mesh, fixture.parameters.head(5)), std::invalid_argument);
  for (int parameter : { 0, 1, 2, 3, 4 }) {
    ES::VXd invalid = fixture.parameters;
    invalid[parameter] = parameter == 1 || parameter == 3 ? 0.5 : 0;
    EXPECT_THROW(CudaKoiterDeformationModelEnergy(*fixture.mesh, invalid), std::invalid_argument);
  }
  ES::VXd invalid = fixture.parameters;
  invalid[0] = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(CudaKoiterDeformationModelEnergy(*fixture.mesh, invalid), std::invalid_argument);
  EXPECT_THROW(fixture.gpu->addHessianVector(nullptr, nullptr), std::logic_error);
  EXPECT_THROW(fixture.gpu->addDiagonal(nullptr), std::logic_error);
}
