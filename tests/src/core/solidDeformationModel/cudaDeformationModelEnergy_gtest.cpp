#include "cudaDeformationModelEnergy.h"
#include "cudaDeviceBuffer.h"
#include "deformationModelAssembler.h"
#include "deformationModelEnergy.h"
#include "deformationModelManager.h"
#include "simulationMesh.h"
#include "pgoLogging.h"
#include "hexJacobianValidation.h"
#include "cudaElasticEquivalence.h"

#include <gtest/gtest.h>
#include <Eigen/Geometry>

#include <algorithm>
#include <cmath>
#include <memory>
#include <limits>
#include <stdexcept>
#include <vector>

namespace
{
namespace ES = pgo::EigenSupport;
using namespace pgo::SolidDeformationModel;
using pgo::CudaUtilities::DeviceBuffer;
using CubicValidity = CudaDeformationModelEnergy::CubicValidity;

class CudaDeformationModelEnergyGTest : public ::testing::Test
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
  std::unique_ptr<CudaDeformationModelEnergy> gpu;
  ES::VXd rest;
  bool projectedTangentIsPSD = true;
  CubicValidity cubicValidity = CubicValidity::Gauss;
};

Fixture makeFixture(double scale = 1.0, bool mirrored = false, bool multiple = true, bool maxStep = true,
  DeformationModelElasticMaterial model = DeformationModelElasticMaterial::STABLE_NEO, double compression = 10000)
{
  std::vector<double> vertices{
    0.0, 0.0, 0.0,
    1.0, 0.0, 0.0,
    0.0, 1.0, 0.0,
    0.0, 0.0, 1.0,
    1.0, 1.0, 1.0
  };
  for (std::size_t i = 0; i < vertices.size(); ++i)
    vertices[i] *= scale * ((mirrored && i % 3 == 0) ? -1.0 : 1.0);
  std::vector<int> indices{ 0, 1, 2, 3 };
  if (multiple)
    indices.insert(indices.end(), { 1, 2, 3, 4 });
  std::vector<int> materialIndices(indices.size() / 4, 0);
  SimulationMeshENuMaterial material(1200.0, 0.35, compression);
  const SimulationMeshMaterial *materials[]{ &material };
  Fixture fixture;
  fixture.mesh = std::make_shared<SimulationMesh>(static_cast<int>(vertices.size() / 3), vertices.data(),
    static_cast<int>(indices.size() / 4), 4, indices.data(), materialIndices.data(), 1, materials, SimulationMeshType::TET);
  fixture.rest = Eigen::Map<const ES::VXd>(vertices.data(), vertices.size());
  fixture.manager = std::make_shared<DeformationModelManager>();
  fixture.manager->setMesh(fixture.mesh.get());
  fixture.manager->init(DeformationModelPlasticMaterial::VOLUMETRIC_DOF0, model);
  fixture.manager->setEnforceSPD(1);
  fixture.assembler = std::make_shared<DeformationModelAssembler>(fixture.manager);
  fixture.assembler->setDeterministic(true);
  fixture.cpu = std::make_unique<DeformationModelEnergy>(fixture.assembler, &fixture.rest);
  fixture.cpu->setEnableMaterialMaxStep(maxStep);
  fixture.gpu = std::make_unique<CudaDeformationModelEnergy>(*fixture.mesh, model, maxStep);
  fixture.projectedTangentIsPSD = model != DeformationModelElasticMaterial::STVK_VOL;
  return fixture;
}

Fixture makeCubicFixture(double scale = 1.0, bool mirrored = false, bool distorted = true, bool maxStep = true,
  DeformationModelElasticMaterial model = DeformationModelElasticMaterial::STABLE_NEO, double compression = 10000,
  CubicValidity cubicValidity = CubicValidity::Gauss)
{
  std::vector<double> vertices;
  for (int z = 0; z < 2; ++z) {
    for (int y = 0; y < 2; ++y) {
      for (int x = 0; x < 3; ++x) {
        const int vertex = x + 3 * (y + 2 * z);
        const double position[3] = { double(x), double(y), double(z) };
        for (int c = 0; c < 3; ++c) {
          const double warp = distorted ? 0.07 * std::sin(0.81 * vertex + 1.3 * c) : 0.0;
          vertices.push_back(scale * (position[c] + warp) * (mirrored && c == 0 ? -1.0 : 1.0));
        }
      }
    }
  }
  const int indices[] = { 0, 1, 4, 3, 6, 7, 10, 9, 1, 2, 5, 4, 7, 8, 11, 10 };
  const int materialIndices[] = { 0, 1 };
  SimulationMeshENuMaterial material(1200.0, 0.35, compression), secondMaterial(850.0, 0.2, 0.7 * compression);
  const SimulationMeshMaterial *materials[]{ &material, &secondMaterial };
  Fixture fixture;
  fixture.mesh = std::make_shared<SimulationMesh>(12, vertices.data(), 2, 8, indices,
    materialIndices, 2, materials, SimulationMeshType::CUBIC);
  fixture.rest = Eigen::Map<const ES::VXd>(vertices.data(), vertices.size());
  fixture.manager = std::make_shared<DeformationModelManager>();
  fixture.manager->setMesh(fixture.mesh.get());
  fixture.manager->init(DeformationModelPlasticMaterial::VOLUMETRIC_DOF0, model);
  fixture.manager->setEnforceSPD(1);
  fixture.assembler = std::make_shared<DeformationModelAssembler>(fixture.manager);
  fixture.assembler->setDeterministic(true);
  fixture.cpu = std::make_unique<DeformationModelEnergy>(fixture.assembler, &fixture.rest);
  fixture.cpu->setEnableMaterialMaxStep(maxStep);
  fixture.gpu = std::make_unique<CudaDeformationModelEnergy>(*fixture.mesh, model, maxStep, cubicValidity);
  fixture.projectedTangentIsPSD = model != DeformationModelElasticMaterial::STVK_VOL;
  fixture.cubicValidity = cubicValidity;
  return fixture;
}

ES::VXd affineDisplacement(const Fixture &fixture, const ES::M3d &F)
{
  ES::VXd x(fixture.rest.size());
  for (int i = 0; i < x.size() / 3; ++i)
    x.segment<3>(3 * i) = (F - ES::M3d::Identity()) * fixture.rest.segment<3>(3 * i);
  return x;
}

void expectVectorNear(const ES::VXd &actual, const ES::VXd &expected, double tolerance)
{
  EXPECT_LE((actual - expected).norm(), tolerance * std::max(1e-14, expected.norm()));
}

void expectDerivatives(Fixture &fixture, const ES::VXd &x)
{
  const int n = fixture.gpu->getNumDOFs();
  DeviceBuffer<double> deviceX, deviceV, deviceOut;
  deviceX.upload(x.data(), x.size());
  const double referenceEnergy = fixture.cpu->func(x);
  EXPECT_NEAR(fixture.gpu->func(deviceX.data()), referenceEnergy, 2e-10 * std::max(1e-10, std::abs(referenceEnergy)));

  ES::VXd referenceGradient(n), result(n), baseline(n), direction(n);
  for (int i = 0; i < n; ++i) {
    baseline[i] = 0.13 * std::sin(0.7 * i);
    direction[i] = std::cos(0.31 * i + 0.4);
  }
  fixture.cpu->gradient(x, referenceGradient);
  deviceOut.upload(baseline.data(), n);
  fixture.gpu->addGradient(deviceX.data(), deviceOut.data());
  deviceOut.download(result.data(), n);
  expectVectorNear(result, referenceGradient + baseline, 2e-10);

  ES::SpMatD hessian;
  fixture.cpu->createHessian(hessian);
  fixture.cpu->hessian(x, hessian);
  fixture.gpu->prepareHessian(deviceX.data());
  deviceV.upload(direction.data(), n);
  deviceOut.upload(baseline.data(), n);
  fixture.gpu->addHessianVector(deviceV.data(), deviceOut.data());
  deviceOut.download(result.data(), n);
  const ES::VXd referenceProduct = hessian * direction;
  expectVectorNear(result, referenceProduct + baseline, 3e-8);
  if (fixture.projectedTangentIsPSD) {
    EXPECT_GE(direction.dot(result - baseline), -1e-8 * std::max(1.0, referenceProduct.norm()));
  }

  deviceOut.upload(baseline.data(), n);
  fixture.gpu->addDiagonal(deviceOut.data());
  deviceOut.download(result.data(), n);
  const ES::VXd referenceDiagonal = hessian.diagonal();
  expectVectorNear(result, referenceDiagonal + baseline, 3e-8);

  // Gradient checks concern the unmodified energy. The projected Hessian is
  // compared above with the CPU projection, not with finite differences.
  const double h = 1e-6 * std::max(1e-5, fixture.rest.norm());
  const ES::VXd plus = x + h * direction, minus = x - h * direction;
  deviceX.upload(plus.data(), n);
  const double ePlus = fixture.gpu->func(deviceX.data());
  deviceX.upload(minus.data(), n);
  const double eMinus = fixture.gpu->func(deviceX.data());
  const double derivative = referenceGradient.dot(direction);
  EXPECT_NEAR((ePlus - eMinus) / (2.0 * h), derivative, 2e-6 * std::max(1e-12, std::abs(derivative)));
}

double compareMaxStep(Fixture &fixture, const ES::VXd &x, const ES::VXd &dx)
{
  DeviceBuffer<double> deviceX, deviceDx;
  deviceX.upload(x.data(), x.size());
  deviceDx.upload(dx.data(), dx.size());
  const double expected = fixture.cpu->computeMaxStepSize(x, dx);
  const double actual = fixture.gpu->computeMaxStepSize(deviceX.data(), deviceDx.data());
  if (fixture.cubicValidity == CubicValidity::WholeCell || dx.isZero()) {
    // Whole-cell cubic validity is stricter than the CPU material's Gauss-only
    // bound. Zero directions now also reject invalid initial configurations.
    EXPECT_GE(actual, 0);
    EXPECT_LE(actual, expected + 2e-11 * std::max(1e-12, expected));
  }
  else
    EXPECT_NEAR(actual, expected, 2e-11 * std::max(1e-12, expected));
  return actual;
}

double gpuMaxStep(const CudaDeformationModelEnergy &energy, const ES::VXd &x, const ES::VXd &dx)
{
  DeviceBuffer<double> deviceX, deviceDx;
  deviceX.upload(x.data(), x.size());
  deviceDx.upload(dx.data(), dx.size());
  return energy.computeMaxStepSize(deviceX.data(), deviceDx.data());
}

namespace HV = pgo::RunSim::HexValidation;
HV::Positions hexPositions(const SimulationMesh &mesh, int element, const ES::VXd &x)
{
  HV::Positions positions;
  for (int local = 0; local < 8; ++local) {
    const int vertex = mesh.getVertexIndex(element, local);
    double rest[3];
    mesh.getVertex(vertex, rest);
    for (int c = 0; c < 3; ++c)
      positions[local][c] = HV::Interval(rest[c] + x[3 * vertex + c]);
  }
  return positions;
}

void expectIndependentHexCertificate(const SimulationMesh &mesh, const ES::VXd &from, const ES::VXd &to)
{
  const ES::VXd zero = ES::VXd::Zero(from.size());
  for (int e = 0; e < mesh.getNumElements(); ++e) {
    const auto rest = hexPositions(mesh, e, zero);
    const double sign = HV::jacobian(rest, { HV::Interval(0.5), HV::Interval(0.5), HV::Interval(0.5) }).lo < 0 ? -1 : 1;
    const auto a = hexPositions(mesh, e, from), b = hexPositions(mesh, e, to);
    EXPECT_TRUE(HV::certify(a, &b, sign).certified) << "element=" << e;
  }
}
}  // namespace

TEST_F(CudaDeformationModelEnergyGTest, MatchesCpuEnergyGradientAndProjectedTangent)
{
  auto fixture = makeFixture();
  ES::M3d F;
  F << 1.15, 0.12, -0.08,
    0.06, 0.84, 0.09,
    0.02, -0.04, 1.08;
  ES::VXd x = affineDisplacement(fixture, F);
  for (int i = 0; i < x.size(); ++i)
    x[i] += 0.02 * std::sin(0.91 * i);
  expectDerivatives(fixture, x);

  F.setIdentity();
  F(0, 0) = 0.05;
  F(1, 1) = 0.2;
  expectDerivatives(fixture, affineDisplacement(fixture, F));
  F(0, 0) = -0.25;
  expectDerivatives(fixture, affineDisplacement(fixture, F));
}

TEST_F(CudaDeformationModelEnergyGTest, QuantitativeFullHessianColumnsMatchCpu)
{
  nlohmann::json reports = nlohmann::json::array();
  for (auto material : { DeformationModelElasticMaterial::STABLE_NEO, DeformationModelElasticMaterial::STVK_VOL }) {
    for (bool cubic : { false, true }) {
      auto fixture = cubic ? makeCubicFixture(1, false, true, true, material) : makeFixture(1, false, true, true, material);
      for (bool inverted : { false, true }) {
        ES::M3d F;
        F << (inverted ? -0.25 : 0.65), 0.12, -0.04,
          0.03, 0.78, 0.06, -0.02, 0.08, 1.03;
        ES::VXd x = affineDisplacement(fixture, F);
        for (int i = 0; i < x.size(); ++i)
          x[i] += 0.012 * std::sin(0.83 * i);
        auto report = pgo::Tests::ElasticEquivalence::compare(*fixture.cpu, *fixture.manager, *fixture.gpu, x, true);
        report["element_type"] = cubic ? "cubic" : "tet";
        report["material"] = material == DeformationModelElasticMaterial::STABLE_NEO ? "stable-neo" : "stvk-vol";
        report["state"] = inverted ? "inverted_affine_with_nonaffine_perturbation" : "compressed_affine_with_nonaffine_perturbation";
        report["num_elements"] = fixture.mesh->getNumElements();
        report["projection_scope"] = material == DeformationModelElasticMaterial::STABLE_NEO ?
          "9x9 material tangent at each quadrature point" : "StVK child 9x9 tangent projection plus unprojected volumetric penalty, matching CPU; not necessarily PSD";
        EXPECT_TRUE(report.at("passed").get<bool>()) << report.dump(2);
        reports.push_back(std::move(report));
      }
    }
  }
  RecordProperty("elastic_equivalence_measurements", reports.dump());
}

TEST_F(CudaDeformationModelEnergyGTest, TinyAndMirroredRestTetsMatchCpu)
{
  for (bool mirrored : { false, true }) {
    auto fixture = makeFixture(1e-5, mirrored);
    ES::M3d F = ES::M3d::Identity();
    F(0, 0) = 0.73;
    F(0, 1) = 0.18;
    expectDerivatives(fixture, affineDisplacement(fixture, F));
  }
}

TEST_F(CudaDeformationModelEnergyGTest, ExactTangentMatchesCpuAndGradientDerivative)
{
  for (bool cubic : { false, true }) {
    auto fixture = cubic ? makeCubicFixture() : makeFixture();
    fixture.manager->setEnforceSPD(0);
    const int n = fixture.gpu->getNumDOFs();
    ES::VXd direction(n), baseline(n), result(n), plusGradient(n), minusGradient(n);
    for (int i = 0; i < n; ++i) {
      direction[i] = std::sin(0.63 * i + 0.2);
      baseline[i] = 0.07 * std::cos(0.47 * i);
    }
    DeviceBuffer<double> deviceX, deviceV, deviceOut;
    deviceV.upload(direction.data(), n);
    for (bool inverted : { false, true }) {
      ES::M3d F;
      F << (inverted ? -0.25 : 0.55), 0.12, -0.04,
        0.03, 0.71, 0.06,
        -0.02, 0.08, 1.03;
      ES::VXd x = affineDisplacement(fixture, F);
      for (int i = 0; i < n; ++i)
        x[i] += 0.015 * std::sin(0.83 * i);
      ES::SpMatD hessian;
      fixture.cpu->createHessian(hessian);
      fixture.cpu->hessian(x, hessian);
      deviceX.upload(x.data(), n);
      fixture.gpu->prepareExactHessian(deviceX.data());
      deviceOut.upload(baseline.data(), n);
      fixture.gpu->addHessianVector(deviceV.data(), deviceOut.data());
      deviceOut.download(result.data(), n);
      const ES::VXd exactProduct = result - baseline;
      expectVectorNear(exactProduct, hessian * direction, 3e-10);

      deviceOut.upload(baseline.data(), n);
      fixture.gpu->addDiagonal(deviceOut.data());
      deviceOut.download(result.data(), n);
      const ES::VXd referenceDiagonal = hessian.diagonal();
      expectVectorNear(result, referenceDiagonal + baseline, 3e-10);

      // The unprojected operator must differentiate the actual gradient, even
      // where compression or inversion makes the material tangent indefinite.
      const double step = 1e-6;
      const ES::VXd plus = x + step * direction, minus = x - step * direction;
      deviceX.upload(plus.data(), n);
      deviceOut.zero();
      fixture.gpu->addGradient(deviceX.data(), deviceOut.data());
      deviceOut.download(plusGradient.data(), n);
      deviceX.upload(minus.data(), n);
      deviceOut.zero();
      fixture.gpu->addGradient(deviceX.data(), deviceOut.data());
      deviceOut.download(minusGradient.data(), n);
      expectVectorNear(exactProduct, (plusGradient - minusGradient) / (2.0 * step), 3e-8);

      // Re-preparation must replace the exact cache with the usual projected
      // operator, as required after an indefinite Newton system is rejected.
      fixture.manager->setEnforceSPD(1);
      fixture.cpu->hessian(x, hessian);
      deviceX.upload(x.data(), n);
      fixture.gpu->prepareHessian(deviceX.data());
      deviceOut.zero();
      fixture.gpu->addHessianVector(deviceV.data(), deviceOut.data());
      deviceOut.download(result.data(), n);
      expectVectorNear(result, hessian * direction, 3e-8);
      EXPECT_GT((result - exactProduct).norm(), 1e-3);
      fixture.manager->setEnforceSPD(0);
    }
  }
}

TEST_F(CudaDeformationModelEnergyGTest, MaterialMaxStepFindsFirstBoundaryAndMatchesCpu)
{
  for (double scale : { 1.0, 1e-5 }) {
    for (bool mirrored : { false, true }) {
      auto fixture = makeFixture(scale, mirrored);
      const ES::VXd zero = ES::VXd::Zero(fixture.rest.size());
      EXPECT_DOUBLE_EQ(compareMaxStep(fixture, zero, zero), 1.0);
      ES::VXd translation = zero;
      for (int i = 0; i < translation.size() / 3; ++i)
        translation.segment<3>(3 * i) << 0.1, -0.3, 0.2;
      EXPECT_DOUBLE_EQ(compareMaxStep(fixture, zero, translation), 1.0);

      ES::M3d endF = ES::M3d::Identity();
      endF(0, 0) = -2.0;
      endF(1, 1) = -1.0;
      // Both endpoint determinants are positive; the first interior crossing
      // must still limit this step.
      const double alpha = compareMaxStep(fixture, zero, affineDisplacement(fixture, endF));
      EXPECT_GT(alpha, 0.32);
      EXPECT_LT(alpha, 1.0 / 3.0);

      ES::M3d illegalF = ES::M3d::Identity();
      illegalF(0, 0) = -0.1;
      const ES::VXd illegal = affineDisplacement(fixture, illegalF);
      EXPECT_NEAR(compareMaxStep(fixture, illegal, translation), 1e-12, 1e-24);
      EXPECT_DOUBLE_EQ(compareMaxStep(fixture, illegal, zero), 0.0);
    }
  }
}

TEST_F(CudaDeformationModelEnergyGTest, MaterialMaxStepCanBeDisabled)
{
  auto fixture = makeFixture(1.0, false, true, false);
  const ES::VXd zero = ES::VXd::Zero(fixture.rest.size());
  ES::M3d endF = ES::M3d::Identity();
  endF(0, 0) = -2.0;
  EXPECT_DOUBLE_EQ(compareMaxStep(fixture, zero, affineDisplacement(fixture, endF)), 1.0);
}

TEST_F(CudaDeformationModelEnergyGTest, NearRestEnergyRetainsQuadraticScaling)
{
  auto fixture = makeFixture(1.0, false, false);
  constexpr double strain = 1e-9;
  const ES::VXd x = strain * fixture.rest;
  DeviceBuffer<double> deviceX, deviceGradient(x.size());
  deviceX.upload(x.data(), x.size());
  const auto *material = dynamic_cast<const SimulationMeshENuMaterial *>(fixture.mesh->getElementMaterial(0, 0));
  ASSERT_NE(material, nullptr);
  const double mu = material->getMuLame(), lambda = material->getLambdaLame() + mu;
  // Leading term of isotropic strain energy: volume * .5 * (9 lambda - 3 mu) s^2.
  // The relative higher-order correction is O(s); constant subtraction in the
  // original expression destroys this small positive energy entirely.
  const double expected = (1.0 / 6.0) * 0.5 * (9.0 * lambda - 3.0 * mu) * strain * strain;
  const double energy = fixture.gpu->func(deviceX.data());
  EXPECT_GT(energy, 0.0);
  EXPECT_NEAR(energy, expected, 1e-7 * expected);
  deviceGradient.zero();
  fixture.gpu->addGradient(deviceX.data(), deviceGradient.data());
  ES::VXd gradient(x.size());
  deviceGradient.download(gradient.data(), gradient.size());
  EXPECT_NEAR(gradient.dot(fixture.rest), 2.0 * expected / strain, 1e-7 * expected / strain);
}

TEST_F(CudaDeformationModelEnergyGTest, CubicMatchesCpuQuadratureAndProjectedTangent)
{
  auto fixture = makeCubicFixture();
  ES::M3d F;
  F << 1.12, 0.15, -0.06,
    -0.02, 0.82, 0.11,
    0.04, 0.07, 0.93;
  ES::VXd displacement = affineDisplacement(fixture, F);
  for (int i = 0; i < displacement.size(); ++i)
    displacement[i] += 0.08 * std::sin(1.21 * i + 0.3);
  expectDerivatives(fixture, displacement);
  F.setIdentity();
  F(0, 0) = 0.05;
  F(1, 1) = 0.2;
  expectDerivatives(fixture, affineDisplacement(fixture, F));
  F(0, 0) = -0.25;
  expectDerivatives(fixture, affineDisplacement(fixture, F));
}

TEST_F(CudaDeformationModelEnergyGTest, TinyAndMirroredCubicRestElementsMatchCpu)
{
  for (bool mirrored : { false, true }) {
    auto fixture = makeCubicFixture(1e-5, mirrored);
    ES::M3d F = ES::M3d::Identity();
    F(0, 0) = 0.71;
    F(0, 1) = 0.17;
    expectDerivatives(fixture, affineDisplacement(fixture, F));
  }
}

TEST_F(CudaDeformationModelEnergyGTest, CubicMaterialMaxStepMatchesCpuGaussPointBound)
{
  for (double scale : { 1.0, 1e-5 }) {
    for (bool mirrored : { false, true }) {
      auto fixture = makeCubicFixture(scale, mirrored);
      const ES::VXd zero = ES::VXd::Zero(fixture.rest.size());
      EXPECT_DOUBLE_EQ(compareMaxStep(fixture, zero, zero), 1.0);
      ES::VXd translation = zero;
      for (int i = 0; i < translation.size() / 3; ++i)
        translation.segment<3>(3 * i) << 0.1, -0.3, 0.2;
      EXPECT_DOUBLE_EQ(compareMaxStep(fixture, zero, translation), 1.0);
      ES::M3d endF = ES::M3d::Identity();
      endF(0, 0) = -2.0;
      endF(1, 1) = -1.0;
      EXPECT_LT(compareMaxStep(fixture, zero, affineDisplacement(fixture, endF)), 1.0 / 3.0);
      ES::VXd cornerMotion = zero;
      cornerMotion[3 * 10 + 1] = -5.0 * scale;
      const double bound = compareMaxStep(fixture, zero, cornerMotion);
      EXPECT_GT(bound, 0.0);
      EXPECT_LT(bound, 1.0);
      ES::VXd randomMotion = zero;
      for (int i = 0; i < randomMotion.size(); ++i)
        randomMotion[i] = scale * 1.5 * std::sin(0.713 * i + 0.24);
      compareMaxStep(fixture, zero, randomMotion);
      endF.setIdentity();
      endF(0, 0) = -0.1;
      const ES::VXd illegal = affineDisplacement(fixture, endF);
      EXPECT_NEAR(compareMaxStep(fixture, illegal, translation), 1e-12, 1e-24);
      EXPECT_DOUBLE_EQ(compareMaxStep(fixture, illegal, zero), 0.0);
    }
  }
  auto disabled = makeCubicFixture(1.0, false, true, false);
  const ES::VXd zero = ES::VXd::Zero(disabled.rest.size());
  EXPECT_DOUBLE_EQ(compareMaxStep(disabled, zero, -2.0 * disabled.rest), 1.0);
}

TEST_F(CudaDeformationModelEnergyGTest, WholeCubicCellRejectsCornerInversionHiddenFromQuadrature)
{
  for (double scale : { 1.0, 1e-5 }) {
    for (bool mirrored : { false, true }) {
      auto fixture = makeCubicFixture(scale, mirrored, false, true,
        DeformationModelElasticMaterial::STABLE_NEO, 10000, CubicValidity::WholeCell);
      const ES::VXd zero = ES::VXd::Zero(fixture.rest.size());
      ES::VXd target = zero;
      target[3 * 10 + 2] = -1.2 * scale;
      // Each affected regular cell has det(F)=1-1.2*u*v: positive at
      // all eight Gauss points, but -0.2 at its upper corner.
      EXPECT_DOUBLE_EQ(fixture.cpu->computeMaxStepSize(zero, target), 1.0);
      const double bound = gpuMaxStep(*fixture.gpu, zero, target);
      EXPECT_GT(bound, 0.7);
      EXPECT_LT(bound, 1.0 / 1.2);
      expectIndependentHexCertificate(*fixture.mesh, zero, bound * target);
      EXPECT_DOUBLE_EQ(gpuMaxStep(*fixture.gpu, target, zero), 0.0);
      target[3 * 10 + 2] = -(1 - 0.5e-8) * scale;
      EXPECT_DOUBLE_EQ(gpuMaxStep(*fixture.gpu, target, zero), 0.0) << "Existing 1e-8 determinant floor applies to the full cell.";
    }
  }
}

TEST_F(CudaDeformationModelEnergyGTest, WholeCubicCellChecksSweptMotionAndRespectsDisabledGuard)
{
  auto fixture = makeCubicFixture(1, false, false, true,
    DeformationModelElasticMaterial::STABLE_NEO, 10000, CubicValidity::WholeCell);
  const ES::VXd zero = ES::VXd::Zero(fixture.rest.size());
  ES::M3d F = ES::M3d::Identity();
  F(0, 0) = -1;
  F(1, 1) = -2;
  const auto target = affineDisplacement(fixture, F);
  EXPECT_DOUBLE_EQ(gpuMaxStep(*fixture.gpu, target, zero), 1.0);  // Positive endpoint determinant.
  const double bound = gpuMaxStep(*fixture.gpu, zero, target);
  EXPECT_GT(bound, 0);
  EXPECT_LT(bound, 1.0 / 3.0);
  expectIndependentHexCertificate(*fixture.mesh, zero, bound * target);
  auto disabled = makeCubicFixture(1, false, false, false,
    DeformationModelElasticMaterial::STABLE_NEO, 10000, CubicValidity::WholeCell);
  ES::VXd inverted = zero;
  inverted[3 * 10 + 2] = -1.2;
  EXPECT_DOUBLE_EQ(gpuMaxStep(*disabled.gpu, inverted, zero), 1.0);
  EXPECT_DOUBLE_EQ(gpuMaxStep(*disabled.gpu, zero, inverted), 1.0);
}

TEST_F(CudaDeformationModelEnergyGTest, WholeCubicCellRejectsInvalidReferenceAndNonfiniteMotion)
{
  const int indices[8]{ 0, 1, 2, 3, 4, 5, 6, 7 }, materialIndex = 0;
  double vertices[24]{ 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0,
    0, 0, 1, 1, 0, 1, 1, 1, -0.2, 0, 1, 1 };
  SimulationMeshENuMaterial material(1200, 0.35, 10000);
  const SimulationMeshMaterial *materials[]{ &material };
  SimulationMesh invalid(8, vertices, 1, 8, indices, &materialIndex, 1, materials, SimulationMeshType::CUBIC);
  EXPECT_THROW((void)CudaDeformationModelEnergy(invalid, true, CubicValidity::WholeCell), std::invalid_argument);
  EXPECT_THROW((void)CudaDeformationModelEnergy(invalid, false, CubicValidity::WholeCell), std::invalid_argument);
  // The same reference is invertible at every CPU Gauss location. Gauss mode
  // must not impose an additional whole-cell reference-shape restriction.
  CudaDeformationModelEnergy gauss(invalid);
  const ES::VXd referenceZero = ES::VXd::Zero(24);
  EXPECT_DOUBLE_EQ(gpuMaxStep(gauss, referenceZero, referenceZero), 1.0);
  for (auto validity : { CubicValidity::Gauss, CubicValidity::WholeCell }) {
    auto fixture = makeCubicFixture(1, false, false, true,
      DeformationModelElasticMaterial::STABLE_NEO, 10000, validity);
    const ES::VXd zero = ES::VXd::Zero(fixture.rest.size());
    for (double value : { std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity(), std::numeric_limits<double>::max() }) {
      ES::VXd direction = zero;
      direction[0] = value;
      EXPECT_DOUBLE_EQ(gpuMaxStep(*fixture.gpu, zero, direction), 0.0);
    }
  }
}

TEST_F(CudaDeformationModelEnergyGTest, CubicGaussValidityIsDefaultAndPreservesCpuMaterial)
{
  for (auto material : { DeformationModelElasticMaterial::STABLE_NEO, DeformationModelElasticMaterial::STVK_VOL }) {
    auto fixture = makeCubicFixture(1, false, false, true, material);
    const ES::VXd zero = ES::VXd::Zero(fixture.rest.size());
    ES::VXd target = zero;
    target[3 * 10 + 2] = -1.2;
    CudaDeformationModelEnergy defaultMode(*fixture.mesh, material);
    CudaDeformationModelEnergy gauss(*fixture.mesh, material, true, CubicValidity::Gauss);
    CudaDeformationModelEnergy whole(*fixture.mesh, material, true, CubicValidity::WholeCell);
    EXPECT_DOUBLE_EQ(compareMaxStep(fixture, zero, target), 1.0);
    EXPECT_DOUBLE_EQ(gpuMaxStep(defaultMode, zero, target), 1.0);
    EXPECT_DOUBLE_EQ(gpuMaxStep(gauss, target, zero), 1.0);
    EXPECT_LT(gpuMaxStep(whole, zero, target), 1.0 / 1.2);
    EXPECT_DOUBLE_EQ(gpuMaxStep(whole, target, zero), 0.0);
    // Feasibility scope does not change the original finite material energy,
    // gradient, projected Hessian or diagonal, including at an inverted corner.
    expectDerivatives(fixture, target);
    DeviceBuffer<double> deviceX;
    deviceX.upload(target.data(), target.size());
    EXPECT_NEAR(whole.func(deviceX.data()), gauss.func(deviceX.data()), 1e-10);
  }
  auto fixture = makeCubicFixture(1, false, false);
  const ES::VXd zero = ES::VXd::Zero(fixture.rest.size());
  ES::VXd target = zero;
  target[3 * 10 + 2] = -1.2;
  CudaDeformationModelEnergy defaultStableNeo(*fixture.mesh);
  EXPECT_DOUBLE_EQ(gpuMaxStep(defaultStableNeo, zero, target), 1.0);
}

TEST_F(CudaDeformationModelEnergyGTest, RotatedNearRestEnergyRetainsQuadraticScaling)
{
  constexpr double strain = 1e-9;
  for (bool cubic : { false, true }) {
    auto fixture = cubic ? makeCubicFixture(1.0, false, false) : makeFixture(1.0, false, false);
    double coefficient = 0.0;
    for (int e = 0; e < fixture.mesh->getNumElements(); ++e) {
      const auto *material = dynamic_cast<const SimulationMeshENuMaterial *>(fixture.mesh->getElementMaterial(e, 0));
      ASSERT_NE(material, nullptr);
      const double mu = material->getMuLame(), lambda = material->getLambdaLame() + mu;
      coefficient += (cubic ? 1.0 : 1.0 / 6.0) * 0.5 * (9.0 * lambda - 3.0 * mu);
    }
    const double expected = coefficient * strain * strain;
    for (double angle : { 0.7, 2.2, 3.141592653589793 }) {
      const ES::M3d rotation = Eigen::AngleAxisd(angle, ES::V3d(1.0, 2.0, -0.5).normalized()).toRotationMatrix();
      ES::VXd displacement = affineDisplacement(fixture, (1.0 + strain) * rotation);
      ES::VXd direction(displacement.size());
      for (int i = 0; i < displacement.size() / 3; ++i) {
        displacement.segment<3>(3 * i) += ES::V3d(3.0, -4.0, 5.0);
        direction.segment<3>(3 * i) = rotation * fixture.rest.segment<3>(3 * i);
      }
      DeviceBuffer<double> deviceX, deviceGradient(displacement.size());
      deviceX.upload(displacement.data(), displacement.size());
      const double energy = fixture.gpu->func(deviceX.data());
      EXPECT_GT(energy, 0.0);
      EXPECT_NEAR(energy, expected, 2e-5 * expected);
      deviceGradient.zero();
      fixture.gpu->addGradient(deviceX.data(), deviceGradient.data());
      ES::VXd gradient(displacement.size());
      deviceGradient.download(gradient.data(), gradient.size());
      EXPECT_NEAR(gradient.dot(direction), 2.0 * expected / strain, 2e-5 * expected / strain);
    }
  }
}

TEST_F(CudaDeformationModelEnergyGTest, RotatedModerateStrainRetainsSmallEnergyDifferences)
{
  auto fixture = makeFixture(1.0, false, false);
  const auto *material = dynamic_cast<const SimulationMeshENuMaterial *>(fixture.mesh->getElementMaterial(0, 0));
  ASSERT_NE(material, nullptr);
  const double mu = material->getMuLame(), lambda = material->getLambdaLame() + mu;
  constexpr double stretch = 0.9375;
  const double delta = std::ldexp(1.0, -34);
  // Uniform stretch has an analytic derivative, independently of its rotation.
  const double derivative = (3.0 * mu * stretch + 3.0 * stretch * stretch * (lambda * (stretch * stretch * stretch - 1.0) - mu)) / 6.0;
  DeviceBuffer<double> deviceX;
  for (int axis = 0; axis < 3; ++axis) {
    // Exactly representable quarter turns and binary stretch increments isolate
    // energy evaluation error from trigonometric/input-rounding error.
    ES::M3d rotation = ES::M3d::Zero();
    rotation(axis, axis) = 1.0;
    rotation((axis + 1) % 3, (axis + 2) % 3) = -1.0;
    rotation((axis + 2) % 3, (axis + 1) % 3) = 1.0;
    const ES::VXd plus = affineDisplacement(fixture, (stretch + delta) * rotation);
    const ES::VXd minus = affineDisplacement(fixture, (stretch - delta) * rotation);
    deviceX.upload(plus.data(), plus.size());
    const double energyPlus = fixture.gpu->func(deviceX.data());
    deviceX.upload(minus.data(), minus.size());
    const double energyMinus = fixture.gpu->func(deviceX.data());
    EXPECT_NEAR((energyPlus - energyMinus) / (2.0 * delta), derivative, 3e-7 * std::abs(derivative));
  }
}

TEST_F(CudaDeformationModelEnergyGTest, PairedEnergyChangesResolveSmallTrialDifferences)
{
  constexpr double stretch = 0.9375;
  const double delta = std::ldexp(1.0, -36);
  ES::M3d rotation;
  rotation << 0, -1, 0, 1, 0, 0, 0, 0, 1;
  for (bool cubic : { false, true }) {
    auto fixture = cubic ? makeCubicFixture(1.0, false, false) : makeFixture(1.0, false, false);
    const ES::VXd base = affineDisplacement(fixture, stretch * rotation);
    const ES::VXd trial = affineDisplacement(fixture, (stretch + delta) * rotation);
    DeviceBuffer<double> deviceBase, deviceTrial;
    deviceBase.upload(base.data(), base.size());
    deviceTrial.upload(trial.data(), trial.size());
    long double expected = 0.0L;
    const long double s = stretch, d = delta;
    const long double deltaJ = 3 * s * s * d + 3 * s * d * d + d * d * d;
    for (int e = 0; e < fixture.mesh->getNumElements(); ++e) {
      const auto *material = dynamic_cast<const SimulationMeshENuMaterial *>(fixture.mesh->getElementMaterial(e, 0));
      ASSERT_NE(material, nullptr);
      const long double mu = material->getMuLame(), lambda = material->getLambdaLame() + material->getMuLame();
      const long double densityChange = 3 * mu * s * d + 1.5L * mu * d * d +
        (lambda * (s * s * s - 1) - mu) * deltaJ + 0.5L * lambda * deltaJ * deltaJ;
      expected += (cubic ? 1.0L : 1.0L / 6.0L) * densityChange;
    }
    const auto change = fixture.gpu->energyChange(deviceBase.data(), deviceTrial.data());
    EXPECT_NEAR(change.value, double(expected), 2e-11 * std::abs(double(expected)));
    EXPECT_GE(change.roundoffScale, std::abs(change.value));
    const auto reverse = fixture.gpu->energyChange(deviceTrial.data(), deviceBase.data());
    EXPECT_NEAR(reverse.value, -change.value, 2e-11 * std::abs(change.value));
    const auto zero = fixture.gpu->energyChange(deviceBase.data(), deviceBase.data());
    EXPECT_DOUBLE_EQ(zero.value, 0.0);
    EXPECT_DOUBLE_EQ(zero.roundoffScale, 0.0);

    // Exercise the same exact increment at a large non-affine change, including
    // inversion. No near-rest approximation may be used by energyChange().
    ES::VXd large = -1.3 * fixture.rest;
    for (int i = 0; i < large.size(); ++i)
      large[i] += 0.04 * std::sin(0.72 * i);
    deviceTrial.upload(large.data(), large.size());
    const auto largeChange = fixture.gpu->energyChange(deviceBase.data(), deviceTrial.data());
    const double reference = fixture.cpu->func(large) - fixture.cpu->func(base);
    EXPECT_NEAR(largeChange.value, reference, 2e-10 * std::max(1.0, std::abs(reference)));
  }
}

TEST_F(CudaDeformationModelEnergyGTest, StvkVolumeTetAndCubicMatchCpuCombinedMaterial)
{
  for (bool cubic : { false, true }) {
    for (double compression : { 0.0, 10000.0 }) {
      auto fixture = cubic ? makeCubicFixture(1, false, true, true, DeformationModelElasticMaterial::STVK_VOL, compression) :
                             makeFixture(1, false, true, true, DeformationModelElasticMaterial::STVK_VOL, compression);
      for (double stretch : { 1.2, 0.23, -0.25 }) {
        ES::M3d F;
        F << stretch, 0.12, -0.08,
          0.03, 0.91, 0.06,
          -0.02, 0.04, 1.07;
        ES::VXd x = affineDisplacement(fixture, F);
        for (int i = 0; i < x.size(); ++i)
          x[i] += 0.006 * std::sin(0.79 * i);
        expectDerivatives(fixture, x);
      }
    }
  }
}

TEST_F(CudaDeformationModelEnergyGTest, StvkVolumeExactTangentAndPairedEnergyChange)
{
  for (bool cubic : { false, true }) {
    auto fixture = cubic ? makeCubicFixture(1, false, true, true, DeformationModelElasticMaterial::STVK_VOL) :
                           makeFixture(1, false, true, true, DeformationModelElasticMaterial::STVK_VOL);
    fixture.manager->setEnforceSPD(0);
    const int n = fixture.rest.size();
    ES::VXd direction(n), result(n), gradient(n);
    for (int i = 0; i < n; ++i)
      direction[i] = std::sin(0.57 * i + 0.41);
    DeviceBuffer<double> deviceX, deviceTrial, deviceV, deviceOut(n);
    deviceV.upload(direction.data(), n);
    for (double stretch : { 0.65, -0.65 }) {
      ES::M3d F = ES::M3d::Identity();
      F(0, 0) = stretch;
      F(0, 1) = 0.12;
      const ES::VXd x = affineDisplacement(fixture, F);
      deviceX.upload(x.data(), n);
      fixture.gpu->prepareExactHessian(deviceX.data());
      deviceOut.zero();
      fixture.gpu->addHessianVector(deviceV.data(), deviceOut.data());
      deviceOut.download(result.data(), n);
      ES::SpMatD hessian;
      fixture.cpu->createHessian(hessian);
      fixture.cpu->hessian(x, hessian);
      expectVectorNear(result, hessian * direction, 3e-8);
      const double h = 1e-6;
      ES::VXd plusGradient(n), minusGradient(n);
      fixture.cpu->gradient(x + h * direction, plusGradient);
      fixture.cpu->gradient(x - h * direction, minusGradient);
      expectVectorNear(result, (plusGradient - minusGradient) / (2 * h), 2e-7);
      fixture.cpu->gradient(x, gradient);
      for (double step : { 1e-3, 1e-10 }) {
        const ES::VXd trial = x + step * direction;
        deviceTrial.upload(trial.data(), n);
        const auto change = fixture.gpu->energyChange(deviceX.data(), deviceTrial.data());
        const double expected = step > 1e-6 ? fixture.cpu->func(trial) - fixture.cpu->func(x) : gradient.dot(trial - x);
        EXPECT_NEAR(change.value, expected, 2e-7 * std::max(1e-12, std::abs(expected)));
        EXPECT_GE(change.roundoffScale, std::abs(change.value) * (1 - 1e-12));
      }
    }
  }
}

TEST_F(CudaDeformationModelEnergyGTest, StvkVolumePreservesMaterialMaxStepAndTinyMeshParity)
{
  for (bool cubic : { false, true }) {
    auto fixture = cubic ? makeCubicFixture(1e-5, true, true, true, DeformationModelElasticMaterial::STVK_VOL) :
                           makeFixture(1e-5, true, true, true, DeformationModelElasticMaterial::STVK_VOL);
    ES::M3d F = ES::M3d::Identity();
    F(0, 0) = 0.72;
    F(0, 1) = 0.14;
    const ES::VXd x = affineDisplacement(fixture, F);
    expectDerivatives(fixture, x);
    F(0, 0) = -1;
    const ES::VXd trial = affineDisplacement(fixture, F);
    EXPECT_LT(compareMaxStep(fixture, x, trial - x), 1.0);
    EXPECT_DOUBLE_EQ(compareMaxStep(fixture, x, ES::VXd::Zero(x.size())), 1.0);
  }
}

TEST_F(CudaDeformationModelEnergyGTest, RejectsUnsupportedMaterialAndNonfiniteCompression)
{
  const auto fixture = makeFixture();
  EXPECT_THROW(CudaDeformationModelEnergy(*fixture.mesh, DeformationModelElasticMaterial::MOONEY_RIVLIN), std::invalid_argument);
  EXPECT_THROW(CudaDeformationModelEnergy(*fixture.mesh, true, static_cast<CubicValidity>(-1)), std::invalid_argument);
  EXPECT_THROW(makeFixture(1, false, true, true, DeformationModelElasticMaterial::STVK_VOL,
                 std::numeric_limits<double>::quiet_NaN()),
    std::invalid_argument);
}

namespace
{
struct DragonHexFixture
{
  std::unique_ptr<SimulationMesh> mesh;
  ES::VXd from, inverted, candidate;
};
DragonHexFixture makeDragonHexFixture()
{
  // Canonical dragon element3972, frames575->576. The old Gauss-only bound
  // returned1: min detF at Gauss points is .4895 while corner detF is -.0211.
  const double rest[24]{
    0.0813736696624768, 0.148012127311545, 0.155870242315661,
    0.103742651022811, 0.148012127311545, 0.155870242315661,
    0.103742651022811, 0.170381108671879, 0.155870242315661,
    0.0813736696624768, 0.170381108671879, 0.155870242315661,
    0.0813736696624768, 0.148012127311545, 0.178239223675995,
    0.103742651022811, 0.148012127311545, 0.178239223675995,
    0.103742651022811, 0.170381108671879, 0.178239223675995,
    0.0813736696624768, 0.170381108671879, 0.178239223675995
  };
  const double fromValues[24]{
    -1.0185909227048155, -0.3618984178498372, -0.11587397771659622,
    -1.0184956098394886, -0.3548022943942232, -0.11759147388298245,
    -1.0290257154699949, -0.3552472323204903, -0.1295604148332937,
    -1.0274536757881638, -0.3661603081455091, -0.12426328969031009,
    -1.0232443954816934, -0.3571139085908797, -0.11641333042806232,
    -1.021659649530752, -0.35285589274508256, -0.11724822719380806,
    -1.0382190112716956, -0.363054863522392, -0.12875549760137756,
    -1.0333905934099086, -0.3586207106288724, -0.12406921182925805
  };
  const double toValues[24]{
    -1.0177105499906298, -0.36230014815420514, -0.11563167519882342,
    -1.0177047052972295, -0.35480989648039163, -0.11711667711193667,
    -1.028618262324014, -0.35560782124872975, -0.12937633492200137,
    -1.0268316835601308, -0.36667757440554755, -0.12403644888328544,
    -1.0225486131216899, -0.35763248359342775, -0.116153287370494,
    -1.0210652453865183, -0.3529111425704744, -0.11670780603864973,
    -1.0374911765642507, -0.36369441405863967, -0.12846811925388604,
    -1.0328047472396806, -0.35913730505710334, -0.12377895996255678
  };
  const double candidateValues[24]{
    -1.0182152468291574, -0.36207684947456165, -0.11576640599804874,
    -1.0181740070407583, -0.35484641250147519, -0.11739697497061835,
    -1.0287885213503098, -0.35540015934703606, -0.1294707918960851,
    -1.0271565646107628, -0.36636706887926895, -0.12416446384536993,
    -1.0229313404766798, -0.35732415057475952, -0.11629120573582032,
    -1.0214210574924658, -0.35290918587456616, -0.11704918032321743,
    -1.0379665903670994, -0.36323368501231318, -0.12864739096459618,
    -1.0331294838232785, -0.35884014376778972, -0.12394548245235952
  };
  const int indices[8]{ 0, 1, 2, 3, 4, 5, 6, 7 }, materialIndex = 0;
  SimulationMeshENuMaterial material(1200, 0.35, 10000);
  const SimulationMeshMaterial *materials[]{ &material };
  auto mesh = std::make_unique<SimulationMesh>(8, rest, 1, 8, indices, &materialIndex, 1, materials, SimulationMeshType::CUBIC);
  return { std::move(mesh), Eigen::Map<const ES::VXd>(fromValues, 24),
    Eigen::Map<const ES::VXd>(toValues, 24), Eigen::Map<const ES::VXd>(candidateValues, 24) };
}
}  // namespace

TEST_F(CudaDeformationModelEnergyGTest, CubicGaussValidityMatchesCpuOnRecordedDragonState)
{
  auto fixture = makeDragonHexFixture();
  const auto &mesh = *fixture.mesh;
  ES::VXd rest(24);
  for (int vertex = 0; vertex < 8; ++vertex)
    mesh.getVertex(vertex, rest.data() + 3 * vertex);
  const ES::VXd zero = ES::VXd::Zero(24);
  for (auto material : { DeformationModelElasticMaterial::STABLE_NEO, DeformationModelElasticMaterial::STVK_VOL }) {
    auto manager = std::make_shared<DeformationModelManager>();
    manager->setMesh(&mesh);
    manager->init(DeformationModelPlasticMaterial::VOLUMETRIC_DOF0, material);
    auto assembler = std::make_shared<DeformationModelAssembler>(manager);
    DeformationModelEnergy cpu(assembler, &rest);
    cpu.setEnableMaterialMaxStep(true);
    CudaDeformationModelEnergy gpu(mesh, material);
    EXPECT_DOUBLE_EQ(cpu.computeMaxStepSize(fixture.from, fixture.inverted - fixture.from), 1.0);
    EXPECT_DOUBLE_EQ(gpuMaxStep(gpu, fixture.from, fixture.inverted - fixture.from), 1.0);
    EXPECT_DOUBLE_EQ(gpuMaxStep(gpu, fixture.inverted, zero), 1.0);
    DeviceBuffer<double> deviceX;
    deviceX.upload(fixture.inverted.data(), 24);
    const double reference = cpu.func(fixture.inverted);
    EXPECT_NEAR(gpu.func(deviceX.data()), reference, 2e-10 * std::max(1e-10, std::abs(reference)));
  }
}

TEST_F(CudaDeformationModelEnergyGTest, WholeCubicCellRejectsDragonCornerInversionFromAcceptedGaussStates)
{
  auto fixture = makeDragonHexFixture();
  const auto &mesh = *fixture.mesh;
  const auto &from = fixture.from, &to = fixture.inverted;
  const ES::VXd zero = ES::VXd::Zero(24);
  for (auto materialType : { DeformationModelElasticMaterial::STABLE_NEO, DeformationModelElasticMaterial::STVK_VOL }) {
    CudaDeformationModelEnergy energy(mesh, materialType, true, CubicValidity::WholeCell);
    EXPECT_DOUBLE_EQ(gpuMaxStep(energy, from, zero), 1.0) << "Last independently valid checkpoint must remain usable.";
    EXPECT_DOUBLE_EQ(gpuMaxStep(energy, to, zero), 0.0);
    const double bound = gpuMaxStep(energy, from, to - from);
    EXPECT_GT(bound, 0);
    EXPECT_LT(bound, 0.3);
    expectIndependentHexCertificate(mesh, from, from + bound * (to - from));
  }
}

TEST_F(CudaDeformationModelEnergyGTest, WholeCubicCellCertifiesPositiveDragonChordsNearMaterialFloor)
{
  auto fixture = makeDragonHexFixture();
  const auto &mesh = *fixture.mesh;
  const ES::VXd zero = ES::VXd::Zero(24);
  const auto reference = HV::certify(hexPositions(mesh, 0, zero));
  ASSERT_TRUE(reference.certified);
  const auto fromPositions = hexPositions(mesh, 0, fixture.from);
  std::vector<ES::VXd> targets;
  // The depth24 guard accepted this candidate at an artificial interval
  // boundary: extending the same CPU-certified chord by one part per million
  // failed, even though its minimum det(F) exceeded .0019.
  for (double scale : { 1.0, 1.000001, 1.00001, 1.001, 1.01 })
    targets.push_back(fixture.from + scale * (fixture.candidate - fixture.from));

  const auto cornerDeterminant = [&](const ES::VXd &x) {
    Eigen::Vector3d positions[8];
    for (int v = 0; v < 8; ++v) {
      mesh.getVertex(v, positions[v].data());
      positions[v] += x.segment<3>(3 * v);
    }
    ES::M3d jacobian;
    jacobian.col(0) = positions[6] - positions[7];
    jacobian.col(1) = positions[6] - positions[5];
    jacobian.col(2) = positions[6] - positions[2];
    return jacobian.determinant();
  };
  const double restDeterminant = cornerDeterminant(zero);
  for (double desired : { 1e-4, 1e-6, 1e-7, 2e-8 }) {
    double lower = 0, upper = 1;
    for (int iteration = 0; iteration < 60; ++iteration) {
      const double middle = 0.5 * (lower + upper);
      const ES::VXd x = fixture.from + middle * (fixture.inverted - fixture.from);
      if (cornerDeterminant(x) / restDeterminant > desired)
        lower = middle;
      else
        upper = middle;
    }
    targets.push_back(fixture.from + lower * (fixture.inverted - fixture.from));
  }
  // Independent CPU Bernstein certificates establish the full original
  // chords, including a positive margin above the unchanged1e-8 floor.
  for (const auto &target : targets) {
    const auto toPositions = hexPositions(mesh, 0, target);
    const auto certificate = HV::certify(fromPositions, &toPositions);
    ASSERT_TRUE(certificate.certified);
    ASSERT_GT(certificate.lower / reference.upper, 1e-8);
  }
  for (auto material : { DeformationModelElasticMaterial::STABLE_NEO, DeformationModelElasticMaterial::STVK_VOL }) {
    CudaDeformationModelEnergy energy(mesh, material, true, CubicValidity::WholeCell);
    for (std::size_t i = 0; i < targets.size(); ++i) {
      SCOPED_TRACE(i);
      EXPECT_DOUBLE_EQ(gpuMaxStep(energy, fixture.from, targets[i] - fixture.from), 1.0);
      EXPECT_DOUBLE_EQ(gpuMaxStep(energy, targets[i], zero), 1.0);
    }
  }
}
