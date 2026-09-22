#include <gtest/gtest.h>

#include "EigenSupport.h"
#include "TRBDF2TimeIntegratorHelper.h"
#include "TRBDF2TimeIntegrator.h"
#include "deformationModelAssembler.h"
#include "deformationModelEnergy.h"
#include "deformationModelManager.h"
#include "implicitBackwardEulerTimeIntegratorHelper.h"
#include "implicitBackwardEulerTimeIntegrator.h"
#include "pgoLogging.h"
#include "potentialEnergy.h"
#include "simulationMesh.h"
#include "tetMeshDeformationModel.h"
#include "triMeshGeo.h"
#include "cubicMeshDeformationModel.h"
#include "materialMaxStepPolynomialUtils.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
namespace ES = pgo::EigenSupport;
using pgo::NonlinearOptimization::PotentialEnergy;
using pgo::Simulation::ImplicitBackwardEulerTimeIntegrator;
using pgo::Simulation::TRBDF2TimeIntegrator;
using pgo::SolidDeformationModel::CubicMeshDeformationModel;
using pgo::SolidDeformationModel::DeformationModelAssembler;
using pgo::SolidDeformationModel::DeformationModelElasticMaterial;
using pgo::SolidDeformationModel::DeformationModelEnergy;
using pgo::SolidDeformationModel::DeformationModelManager;
using pgo::SolidDeformationModel::DeformationModelPlasticMaterial;
using pgo::SolidDeformationModel::SimulationMesh;
using pgo::SolidDeformationModel::SimulationMeshENuMaterial;
using pgo::SolidDeformationModel::SimulationMeshENuhMaterial;
using pgo::SolidDeformationModel::SimulationMeshMaterial;
using pgo::SolidDeformationModel::SimulationMeshType;
using pgo::SolidDeformationModel::TetMeshDeformationModel;

constexpr const char *kShellObjPath = LIBPGO_TEST_SHELL_OBJ;

void initializeLogging()
{
  static const bool initialized = []() {
    pgo::Logging::init();
    return true;
  }();
  (void)initialized;
}

struct EnergyFixture
{
  std::shared_ptr<SimulationMesh> mesh;
  std::shared_ptr<DeformationModelManager> manager;
  std::shared_ptr<DeformationModelAssembler> assembler;
  std::shared_ptr<DeformationModelEnergy> energy;
  ES::VXd restPositions;
};

ES::VXd gatherRestPositions(const SimulationMesh &mesh)
{
  ES::VXd rest(mesh.getNumVertices() * 3);
  for (int vi = 0; vi < mesh.getNumVertices(); vi++) {
    double p[3];
    mesh.getVertex(vi, p);
    rest.segment<3>(vi * 3) << p[0], p[1], p[2];
  }
  return rest;
}

EnergyFixture makeTetFixture(const std::vector<double> &vertices, const std::vector<int> &elementVertices)
{
  initializeLogging();

  std::vector<int> elementMaterialIndices(elementVertices.size() / 4, 0);
  SimulationMeshENuMaterial baseMaterial(1200.0, 0.45);
  const SimulationMeshMaterial *materials[] = { &baseMaterial };

  EnergyFixture fixture;
  fixture.mesh = std::shared_ptr<SimulationMesh>(new SimulationMesh(
    static_cast<int>(vertices.size() / 3), vertices.data(),
    static_cast<int>(elementVertices.size() / 4), 4, elementVertices.data(),
    elementMaterialIndices.data(), 1, materials,
    SimulationMeshType::TET));

  fixture.manager = std::make_shared<DeformationModelManager>();
  fixture.manager->setMesh(fixture.mesh.get());
  fixture.manager->init(DeformationModelPlasticMaterial::VOLUMETRIC_DOF6, DeformationModelElasticMaterial::STABLE_NEO);

  fixture.assembler = std::make_shared<DeformationModelAssembler>(fixture.manager, nullptr);
  fixture.restPositions = gatherRestPositions(*fixture.mesh);
  fixture.energy = std::make_shared<DeformationModelEnergy>(fixture.assembler, &fixture.restPositions, 0);
  return fixture;
}

EnergyFixture makeSingleTetFixture()
{
  const std::vector<double> vertices = {
    0.0, 0.0, 0.0,
    1.0, 0.0, 0.0,
    0.0, 1.0, 0.0,
    0.0, 0.0, 1.0,
  };
  const std::vector<int> elementVertices = { 0, 1, 2, 3 };
  return makeTetFixture(vertices, elementVertices);
}

EnergyFixture makeCubicFixture(const std::vector<double> &vertices, const std::vector<int> &elementVertices)
{
  initializeLogging();

  std::vector<int> elementMaterialIndices(elementVertices.size() / 8, 0);
  SimulationMeshENuMaterial baseMaterial(1200.0, 0.45);
  const SimulationMeshMaterial *materials[] = { &baseMaterial };

  EnergyFixture fixture;
  fixture.mesh = std::shared_ptr<SimulationMesh>(new SimulationMesh(
    static_cast<int>(vertices.size() / 3), vertices.data(),
    static_cast<int>(elementVertices.size() / 8), 8, elementVertices.data(),
    elementMaterialIndices.data(), 1, materials,
    SimulationMeshType::CUBIC));

  fixture.manager = std::make_shared<DeformationModelManager>();
  fixture.manager->setMesh(fixture.mesh.get());
  fixture.manager->init(DeformationModelPlasticMaterial::VOLUMETRIC_DOF6, DeformationModelElasticMaterial::STABLE_NEO);

  fixture.assembler = std::make_shared<DeformationModelAssembler>(fixture.manager, nullptr);
  fixture.restPositions = gatherRestPositions(*fixture.mesh);
  fixture.energy = std::make_shared<DeformationModelEnergy>(fixture.assembler, &fixture.restPositions, 0);
  return fixture;
}

EnergyFixture makeSingleCubicFixture()
{
  const std::vector<double> vertices = {
    0.0, 0.0, 0.0,
    1.0, 0.0, 0.0,
    1.0, 1.0, 0.0,
    0.0, 1.0, 0.0,
    0.0, 0.0, 1.0,
    1.0, 0.0, 1.0,
    1.0, 1.0, 1.0,
    0.0, 1.0, 1.0,
  };
  const std::vector<int> elementVertices = { 0, 1, 2, 3, 4, 5, 6, 7 };
  return makeCubicFixture(vertices, elementVertices);
}

EnergyFixture makeShellFixture()
{
  initializeLogging();

  pgo::Mesh::TriMeshGeo surfaceMesh;
  if (!surfaceMesh.load(kShellObjPath))
    throw std::runtime_error("Failed to load shell regression mesh.");

  SimulationMeshENuhMaterial shellMaterial(1000.0, 0.45, 1e-3);

  EnergyFixture fixture;
  fixture.mesh = std::shared_ptr<SimulationMesh>(pgo::SolidDeformationModel::loadShellMesh(surfaceMesh, &shellMaterial));
  fixture.manager = std::make_shared<DeformationModelManager>();
  fixture.manager->setMesh(fixture.mesh.get());
  fixture.manager->init(DeformationModelPlasticMaterial::SHELL_FF_DOF1, DeformationModelElasticMaterial::KOITER_STVK);

  fixture.assembler = std::make_shared<DeformationModelAssembler>(fixture.manager, nullptr);
  fixture.restPositions = gatherRestPositions(*fixture.mesh);
  fixture.energy = std::make_shared<DeformationModelEnergy>(fixture.assembler, &fixture.restPositions, 0);
  return fixture;
}

ES::VXd makeTetFlipDirection(int numVertices, int topVertex, double dz)
{
  ES::VXd dx = ES::VXd::Zero(numVertices * 3);
  dx[topVertex * 3 + 2] = dz;
  return dx;
}

void applyUniformTranslation(ES::VXd &dx, double tx, double ty, double tz)
{
  for (Eigen::Index vi = 0; vi < dx.size() / 3; vi++) {
    dx[vi * 3 + 0] = tx;
    dx[vi * 3 + 1] = ty;
    dx[vi * 3 + 2] = tz;
  }
}

ES::VXd makeCubicTopFaceDirection(int numVertices, int baseVertex, double dz)
{
  ES::VXd dx = ES::VXd::Zero(numVertices * 3);
  for (int localVertex = 4; localVertex < 8; localVertex++) {
    dx[(baseVertex + localVertex) * 3 + 2] = dz;
  }
  return dx;
}

double tetDeterminant(const SimulationMesh &mesh, int ele, const ES::VXd &absolutePositions)
{
  std::array<double, 12> localPositions{};
  for (int j = 0; j < 4; j++) {
    const int vi = mesh.getVertexIndex(ele, j);
    localPositions[j * 3 + 0] = absolutePositions[vi * 3 + 0];
    localPositions[j * 3 + 1] = absolutePositions[vi * 3 + 1];
    localPositions[j * 3 + 2] = absolutePositions[vi * 3 + 2];
  }

  std::array<double, 9> Ds{};
  TetMeshDeformationModel::computeDs(localPositions.data(), Ds.data());
  return Eigen::Map<const ES::M3d>(Ds.data()).determinant();
}

double minCubicDeterminant(const SimulationMesh &mesh, const DeformationModelManager &manager,
  int ele, const ES::VXd &absolutePositions)
{
  const auto *model = dynamic_cast<const CubicMeshDeformationModel *>(manager.getDeformationModel(ele));
  if (model == nullptr)
    return -std::numeric_limits<double>::infinity();

  std::array<double, 24> localPositions{};
  for (int j = 0; j < 8; j++) {
    const int vi = mesh.getVertexIndex(ele, j);
    localPositions[j * 3 + 0] = absolutePositions[vi * 3 + 0];
    localPositions[j * 3 + 1] = absolutePositions[vi * 3 + 1];
    localPositions[j * 3 + 2] = absolutePositions[vi * 3 + 2];
  }

  double minDet = std::numeric_limits<double>::infinity();
  for (int q = 0; q < model->getNumMaterialLocations(); q++) {
    std::array<double, 9> F{};
    model->computeF(localPositions.data(), q, F.data());
    minDet = std::min(minDet, Eigen::Map<const ES::M3d>(F.data()).determinant());
  }

  return minDet;
}

class FixedMaxStepEnergy : public PotentialEnergy
{
public:
  FixedMaxStepEnergy(int numDOFs, double maxStep):
    numDOFs_(numDOFs), maxStep_(maxStep)
  {
    dofs_.resize(numDOFs_);
    std::iota(dofs_.begin(), dofs_.end(), 0);
  }

  double func(ES::ConstRefVecXd) const override { return 0.0; }
  void gradient(ES::ConstRefVecXd, ES::RefVecXd grad) const override { grad.setZero(); }
  void hessian(ES::ConstRefVecXd, ES::SpMatD &) const override {}
  void createHessian(ES::SpMatD &hess) const override { hess = ES::SpMatD(numDOFs_, numDOFs_); }
  void getDOFs(std::vector<int> &dofs) const override { dofs = dofs_; }
  int getNumDOFs() const override { return numDOFs_; }
  double computeMaxStepSize(ES::ConstRefVecXd, ES::ConstRefVecXd) const override { return maxStep_; }

private:
  int numDOFs_;
  double maxStep_;
  std::vector<int> dofs_;
};

class TestImplicitBackwardEulerTimeIntegrator : public ImplicitBackwardEulerTimeIntegrator
{
public:
  using ImplicitBackwardEulerTimeIntegrator::ImplicitBackwardEulerTimeIntegrator;
  using pgo::Simulation::TimeIntegrator::assembleImplicitModels;
  void prepareStepEnergy()
  {
    assembleImplicitModels();
    updateD();
    updateA();
    updateb();
  }
};

class ZeroElasticEnergy : public FixedMaxStepEnergy
{
public:
  explicit ZeroElasticEnergy(int n): FixedMaxStepEnergy(n, 1.0) {}
  void createHessian(ES::SpMatD &hessian) const override
  {
    hessian.resize(getNumDOFs(), getNumDOFs());
    hessian.setIdentity();
    hessian.coeffs().setZero();
  }
};

class TestTRBDF2TimeIntegrator : public TRBDF2TimeIntegrator
{
public:
  using TRBDF2TimeIntegrator::TRBDF2TimeIntegrator;
  using pgo::Simulation::TimeIntegrator::assembleImplicitModels;

  std::shared_ptr<const PotentialEnergy> getTRStageEnergy() const
  {
    return std::static_pointer_cast<const PotentialEnergy>(trEnergy);
  }
};
}  // namespace

TEST(DeformationModelEnergyMaxStepGTest, ZeroDirectionReturnsOne)
{
  EnergyFixture fixture = makeSingleTetFixture();
  const ES::VXd x = ES::VXd::Zero(fixture.restPositions.size());
  const ES::VXd dx = ES::VXd::Zero(fixture.restPositions.size());

  EXPECT_DOUBLE_EQ(fixture.energy->computeMaxStepSize(x, dx), 1.0);
}

TEST(DeformationModelEnergyMaxStepGTest, TetPureTranslationReturnsOne)
{
  EnergyFixture fixture = makeSingleTetFixture();
  const ES::VXd x = ES::VXd::Zero(fixture.restPositions.size());
  ES::VXd dx = ES::VXd::Zero(fixture.restPositions.size());
  applyUniformTranslation(dx, 1.0, 2.0, 3.0);

  EXPECT_DOUBLE_EQ(fixture.energy->computeMaxStepSize(x, dx), 1.0);
}

TEST(DeformationModelEnergyMaxStepGTest, TetShrinksBeforeInversion)
{
  EnergyFixture fixture = makeSingleTetFixture();
  const ES::VXd x = ES::VXd::Zero(fixture.restPositions.size());
  const ES::VXd dx = makeTetFlipDirection(fixture.mesh->getNumVertices(), 3, -2.0);

  const double alpha = fixture.energy->computeMaxStepSize(x, dx);
  ASSERT_LT(alpha, 1.0);
  ASSERT_GT(alpha, 0.0);

  const ES::VXd updatedPositions = fixture.restPositions + alpha * dx;
  EXPECT_LT(tetDeterminant(*fixture.mesh, 0, fixture.restPositions + dx), 0.0);
  EXPECT_GT(tetDeterminant(*fixture.mesh, 0, updatedPositions), 0.0);
}

TEST(DeformationModelEnergyMaxStepGTest, DisabledMaterialMaxStepSkipsTetClamp)
{
  EnergyFixture fixture = makeSingleTetFixture();
  fixture.energy->setEnableMaterialMaxStep(false);

  const ES::VXd x = ES::VXd::Zero(fixture.restPositions.size());
  const ES::VXd dx = makeTetFlipDirection(fixture.mesh->getNumVertices(), 3, -2.0);

  EXPECT_DOUBLE_EQ(fixture.energy->computeMaxStepSize(x, dx), 1.0);
}

TEST(DeformationModelEnergyMaxStepGTest, TetIllegalInitialStateWarnsEachCallAndClamps)
{
  EnergyFixture fixture = makeSingleTetFixture();
  ES::VXd x = ES::VXd::Zero(fixture.restPositions.size());
  ES::VXd dx = ES::VXd::Zero(fixture.restPositions.size());
  x[3 * 3 + 2] = -2.2;
  dx[0] = 0.1;

  testing::internal::CaptureStdout();
  const double alpha1 = fixture.energy->computeMaxStepSize(x, dx);
  const double alpha2 = fixture.energy->computeMaxStepSize(x, dx);
  const std::string logOutput = testing::internal::GetCapturedStdout();

  EXPECT_GT(alpha1, 0.0);
  EXPECT_LT(alpha1, 1e-9);
  EXPECT_DOUBLE_EQ(alpha1, alpha2);
  EXPECT_NE(logOutput.find("Phase 1.5 material max step encountered illegal initial state"), std::string::npos);
}

TEST(DeformationModelEnergyMaxStepGTest, TetSmallAlphaRemainsFeasible)
{
  EnergyFixture fixture = makeSingleTetFixture();
  const ES::VXd x = ES::VXd::Zero(fixture.restPositions.size());
  const ES::VXd dx = makeTetFlipDirection(fixture.mesh->getNumVertices(), 3, -200.0);

  const double alpha = fixture.energy->computeMaxStepSize(x, dx);

  EXPECT_GT(alpha, 0.0);
  EXPECT_LT(alpha, 0.01);
}

TEST(DeformationModelEnergyMaxStepGTest, TetMultipleElementsReturnEarliestClamp)
{
  const std::vector<double> vertices = {
    0.0, 0.0, 0.0,
    1.0, 0.0, 0.0,
    0.0, 1.0, 0.0,
    0.0, 0.0, 1.0,
    3.0, 0.0, 0.0,
    4.0, 0.0, 0.0,
    3.0, 1.0, 0.0,
    3.0, 0.0, 1.0,
  };
  const std::vector<int> elementVertices = {
    0, 1, 2, 3,
    4, 5, 6, 7,
  };

  EnergyFixture multiFixture = makeTetFixture(vertices, elementVertices);
  const ES::VXd x = ES::VXd::Zero(multiFixture.restPositions.size());
  ES::VXd dx = ES::VXd::Zero(multiFixture.restPositions.size());
  dx[3 * 3 + 2] = -2.0;
  dx[7 * 3 + 2] = -1.2;

  const double alpha = multiFixture.energy->computeMaxStepSize(x, dx);

  EnergyFixture singleFixture = makeSingleTetFixture();
  const ES::VXd singleX = ES::VXd::Zero(singleFixture.restPositions.size());
  const double alphaA = singleFixture.energy->computeMaxStepSize(singleX, makeTetFlipDirection(singleFixture.mesh->getNumVertices(), 3, -2.0));
  const double alphaB = singleFixture.energy->computeMaxStepSize(singleX, makeTetFlipDirection(singleFixture.mesh->getNumVertices(), 3, -1.2));

  EXPECT_NEAR(alpha, std::min(alphaA, alphaB), 1e-12);
}

TEST(DeformationModelEnergyMaxStepGTest, CubicShrinksBeforeInversion)
{
  EnergyFixture fixture = makeSingleCubicFixture();
  const ES::VXd x = ES::VXd::Zero(fixture.restPositions.size());
  const ES::VXd dx = makeCubicTopFaceDirection(fixture.mesh->getNumVertices(), 0, -2.0);

  const double alpha = fixture.energy->computeMaxStepSize(x, dx);
  ASSERT_LT(alpha, 1.0);
  ASSERT_GT(alpha, 0.0);

  const ES::VXd updatedPositions = fixture.restPositions + alpha * dx;
  EXPECT_LT(minCubicDeterminant(*fixture.mesh, *fixture.manager, 0, fixture.restPositions + dx), 0.0);
  EXPECT_GT(minCubicDeterminant(*fixture.mesh, *fixture.manager, 0, updatedPositions), 0.0);
}

TEST(DeformationModelEnergyMaxStepGTest, CubicFeasibleDirectionReturnsOne)
{
  EnergyFixture fixture = makeSingleCubicFixture();
  const ES::VXd x = ES::VXd::Zero(fixture.restPositions.size());
  const ES::VXd dx = makeCubicTopFaceDirection(fixture.mesh->getNumVertices(), 0, -0.2);

  EXPECT_DOUBLE_EQ(fixture.energy->computeMaxStepSize(x, dx), 1.0);
}

TEST(DeformationModelEnergyMaxStepGTest, CubicMultipleElementsReturnEarliestClamp)
{
  const std::vector<double> vertices = {
    0.0, 0.0, 0.0,
    1.0, 0.0, 0.0,
    1.0, 1.0, 0.0,
    0.0, 1.0, 0.0,
    0.0, 0.0, 1.0,
    1.0, 0.0, 1.0,
    1.0, 1.0, 1.0,
    0.0, 1.0, 1.0,
    3.0, 0.0, 0.0,
    4.0, 0.0, 0.0,
    4.0, 1.0, 0.0,
    3.0, 1.0, 0.0,
    3.0, 0.0, 1.0,
    4.0, 0.0, 1.0,
    4.0, 1.0, 1.0,
    3.0, 1.0, 1.0,
  };
  const std::vector<int> elementVertices = {
    0, 1, 2, 3, 4, 5, 6, 7,
    8, 9, 10, 11, 12, 13, 14, 15,
  };

  EnergyFixture multiFixture = makeCubicFixture(vertices, elementVertices);
  const ES::VXd x = ES::VXd::Zero(multiFixture.restPositions.size());
  ES::VXd dx = ES::VXd::Zero(multiFixture.restPositions.size());
  dx += makeCubicTopFaceDirection(multiFixture.mesh->getNumVertices(), 0, -2.0);
  dx += makeCubicTopFaceDirection(multiFixture.mesh->getNumVertices(), 8, -1.2);

  const double alpha = multiFixture.energy->computeMaxStepSize(x, dx);

  EnergyFixture singleFixture = makeSingleCubicFixture();
  const ES::VXd singleX = ES::VXd::Zero(singleFixture.restPositions.size());
  const double alphaA = singleFixture.energy->computeMaxStepSize(singleX, makeCubicTopFaceDirection(singleFixture.mesh->getNumVertices(), 0, -2.0));
  const double alphaB = singleFixture.energy->computeMaxStepSize(singleX, makeCubicTopFaceDirection(singleFixture.mesh->getNumVertices(), 0, -1.2));

  EXPECT_NEAR(alpha, std::min(alphaA, alphaB), 1e-12);
}

TEST(DeformationModelEnergyMaxStepGTest, ShellKeepsUnitStep)
{
  EnergyFixture fixture = makeShellFixture();
  const ES::VXd x = ES::VXd::Zero(fixture.restPositions.size());
  ES::VXd dx = ES::VXd::Zero(fixture.restPositions.size());
  applyUniformTranslation(dx, 0.1, -0.05, 0.2);

  EXPECT_DOUBLE_EQ(fixture.energy->computeMaxStepSize(x, dx), 1.0);
}

TEST(DeformationModelEnergyMaxStepGTest, BackwardEulerResolvesSmallStepsAtLargeDisplacement)
{
  initializeLogging();
  ES::SpMatD mass(3, 3);
  mass.setIdentity();
  mass.coeffRef(0, 0) = 2;
  mass.coeffRef(1, 1) = 3;
  mass.coeffRef(2, 2) = 4;
  constexpr double h = 1e-3, damping = .4;
  TestImplicitBackwardEulerTimeIntegrator integrator(mass, std::make_shared<ZeroElasticEnergy>(3), damping, 0, h, 10, 1e-8);
  const ES::VXd q = ES::VXd::Constant(3, -16384);
  const ES::V3d velocity(.2, -.1, .5), force(.4, -.7, .2);
  integrator.setqState(q, velocity, ES::V3d::Zero());
  integrator.setExternalForce(force.data());
  integrator.prepareStepEnergy();
  const auto energy = integrator.getInternalEnergy();
  const ES::V3d du(std::ldexp(1., -12), -std::ldexp(1., -13), std::ldexp(1., -14));
  const ES::VXd x = q + du;
  const ES::M3d A = ES::M3d(mass) * (1 / (h * h) + damping / h);
  const ES::V3d b = mass * velocity / h + force;
  const double expected = .5 * du.dot(A * du) - du.dot(b);
  EXPECT_DOUBLE_EQ(energy->func(q), 0);
  EXPECT_NEAR(energy->func(x), expected, 1e-13);
  ES::VXd g(3);
  energy->gradient(x, g);
  EXPECT_LT((g - (A * du - b)).norm(), 1e-10);
  ES::SpMatD H;
  energy->hessianDirect(x, H);
  EXPECT_LT((ES::M3d(H) - A).norm(), 1e-9);
  for (int i = 0; i < 3; ++i) {
    constexpr double delta = 1.0 / 65536;
    ES::VXd xp = x, xm = x;
    xp[i] += delta;
    xm[i] -= delta;
    EXPECT_NEAR((energy->func(xp) - energy->func(xm)) / (2 * delta), g[i], 1e-9);
  }
  integrator.setSolution(x);
  integrator.proceedTimestep();
  EXPECT_NEAR(energy->func(x), expected, 1e-13);
  integrator.prepareStepEnergy();
  EXPECT_DOUBLE_EQ(energy->func(x), 0);
}

TEST(DeformationModelEnergyMaxStepGTest, ImplicitBackwardEulerTakesMinWithOtherEnergy)
{
  EnergyFixture fixture = makeSingleTetFixture();
  const ES::VXd x = ES::VXd::Zero(fixture.restPositions.size());
  const ES::VXd dx = makeTetFlipDirection(fixture.mesh->getNumVertices(), 3, -2.0);
  const double materialAlpha = fixture.energy->computeMaxStepSize(x, dx);
  ASSERT_LT(materialAlpha, 0.95);

  ES::SpMatD mass(fixture.restPositions.size(), fixture.restPositions.size());
  mass.setIdentity();

  TestImplicitBackwardEulerTimeIntegrator integrator(mass, fixture.energy, 0.0, 0.0, 0.01, 0, 1e-6);

  integrator.addGeneralImplicitForceModel(std::make_shared<FixedMaxStepEnergy>(fixture.restPositions.size(), 0.95));
  integrator.assembleImplicitModels();
  EXPECT_NEAR(integrator.getInternalEnergy()->computeMaxStepSize(x, dx), materialAlpha, 1e-12);

  integrator.clearGeneralImplicitForceModel();
  integrator.addGeneralImplicitForceModel(std::make_shared<FixedMaxStepEnergy>(fixture.restPositions.size(), 0.25));
  integrator.assembleImplicitModels();
  EXPECT_DOUBLE_EQ(integrator.getInternalEnergy()->computeMaxStepSize(x, dx), 0.25);
}

TEST(DeformationModelEnergyMaxStepGTest, TRBDF2TakesMinWithOtherEnergy)
{
  EnergyFixture fixture = makeSingleTetFixture();
  const ES::VXd x = ES::VXd::Zero(fixture.restPositions.size());
  const ES::VXd dx = makeTetFlipDirection(fixture.mesh->getNumVertices(), 3, -2.0);
  const double materialAlpha = fixture.energy->computeMaxStepSize(x, dx);
  ASSERT_LT(materialAlpha, 0.9);

  ES::SpMatD mass(fixture.restPositions.size(), fixture.restPositions.size());
  mass.setIdentity();

  TestTRBDF2TimeIntegrator integrator(mass, fixture.energy, 0.0, 0.0, 0.5, 0.01, 0, 1e-6);

  integrator.addGeneralImplicitForceModel(std::make_shared<FixedMaxStepEnergy>(fixture.restPositions.size(), 0.9));
  integrator.assembleImplicitModels();
  EXPECT_NEAR(integrator.getTRStageEnergy()->computeMaxStepSize(x, dx), materialAlpha, 1e-12);

  integrator.clearGeneralImplicitForceModel();
  integrator.addGeneralImplicitForceModel(std::make_shared<FixedMaxStepEnergy>(fixture.restPositions.size(), 0.2));
  integrator.assembleImplicitModels();
  EXPECT_DOUBLE_EQ(integrator.getTRStageEnergy()->computeMaxStepSize(x, dx), 0.2);
}

TEST(DeformationModelEnergyMaxStepGTest, SafetyClampNeverPassesTheFeasibleBound)
{
  using namespace pgo::SolidDeformationModel;
  // Linear "determinant" 1e-6 - 2e6 alpha: feasible up to alpha = 5e-13, which
  // is below the 1e-12 minimum clamp. The clamp must not round past the root.
  pgo::BasicAlgorithms::CubicPolynomial poly;
  poly.c0 = 1e-6;
  poly.c1 = -2e6;
  const ConservativeFeasibleAlphaResult result = findConservativeFeasibleAlpha(poly, 0.0);
  EXPECT_FALSE(result.illegalInitialState);
  EXPECT_GT(result.alpha, 0.0);
  EXPECT_LE(result.alpha, 5e-13);
  EXPECT_GE(poly.eval(result.alpha), 0.0);

  // ordinary roots keep the interior safety factor
  pgo::BasicAlgorithms::CubicPolynomial half;
  half.c0 = 0.5;
  half.c1 = -1.0;
  EXPECT_NEAR(findConservativeFeasibleAlpha(half, 0.0).alpha, 0.5 * kMaterialMaxStepInteriorSafety, 1e-12);

  EXPECT_LE(applyMaterialMaxStepSafetyClamp(5e-13), 5e-13);
  EXPECT_DOUBLE_EQ(applyMaterialMaxStepSafetyClamp(2e-12), 2e-12 * kMaterialMaxStepInteriorSafety);
  EXPECT_DOUBLE_EQ(applyMaterialMaxStepSafetyClamp(1.5), 1.0);
  EXPECT_DOUBLE_EQ(applyMaterialMaxStepSafetyClamp(0.0), 0.0);
}

TEST(DeformationModelEnergyMaxStepGTest, NegativelyOrientedRestTetIsLegalAndLimitedLikeItsMirror)
{
  const std::vector<double> vertices = {
    0.0, 0.0, 0.0,
    1.0, 0.0, 0.0,
    0.0, 1.0, 0.0,
    0.0, 0.0, 1.0,
  };
  // Same geometry, vertices 1 and 2 swapped: det(Dm) < 0 while F = I at rest.
  EnergyFixture reversed = makeTetFixture(vertices, { 0, 2, 1, 3 });
  EnergyFixture forward = makeTetFixture(vertices, { 0, 1, 2, 3 });
  const ES::VXd zeroU = ES::VXd::Zero(12);

  const auto atRest = reversed.assembler->computeMaxStepObservation(reversed.restPositions.data(), zeroU.data());
  EXPECT_FALSE(atRest.hasIllegalInitialState);
  EXPECT_DOUBLE_EQ(atRest.alpha, 1.0);

  // pushing the apex through the base plane inverts either tet at alpha = 0.5
  ES::VXd du = ES::VXd::Zero(12);
  du.segment<3>(9) = ES::V3d(0.0, 0.0, -2.0);
  const double reversedAlpha = reversed.energy->computeMaxStepSize(zeroU, du);
  const double forwardAlpha = forward.energy->computeMaxStepSize(zeroU, du);
  EXPECT_GT(reversedAlpha, 0.0);
  EXPECT_LT(reversedAlpha, 1.0);
  EXPECT_NEAR(reversedAlpha, forwardAlpha, 1e-12);

  // a harmless stretch keeps the full step for both orientations
  ES::VXd stretch = ES::VXd::Zero(12);
  stretch.segment<3>(9) = ES::V3d(0.0, 0.0, 0.5);
  EXPECT_DOUBLE_EQ(reversed.energy->computeMaxStepSize(zeroU, stretch), 1.0);
  EXPECT_DOUBLE_EQ(forward.energy->computeMaxStepSize(zeroU, stretch), 1.0);
}

TEST(DeformationModelEnergyMaxStepGTest, TinyTetIsLimitedExactlyLikeItsUnitSizedCopy)
{
  // A 0.1 mm tet has det(Dm) = 1e-12, so in det(Ds) units its whole volume sat
  // inside the limiter's absolute 1e-12 boundary tolerance and a step ending
  // on an inverted element passed as "on the boundary" with alpha = 1. The
  // limiter now works in det(F) units and must be scale invariant.
  const double scale = 1e-4;
  const std::vector<double> unitVertices = {
    0.0, 0.0, 0.0,
    1.0, 0.0, 0.0,
    0.0, 1.0, 0.0,
    0.0, 0.0, 1.0,
  };
  std::vector<double> tinyVertices = unitVertices;
  for (double &v : tinyVertices)
    v *= scale;

  EnergyFixture unit = makeTetFixture(unitVertices, { 0, 1, 2, 3 });
  EnergyFixture tiny = makeTetFixture(tinyVertices, { 0, 1, 2, 3 });
  const ES::VXd zeroU = ES::VXd::Zero(12);

  // Push the apex 1.5 heights downward: det(Ds(1)) = -0.5 det(Dm), which for
  // the tiny tet is -5e-13, inside the old absolute tolerance band.
  const ES::VXd unitDx = makeTetFlipDirection(4, 3, -1.5);
  const ES::VXd tinyDx = makeTetFlipDirection(4, 3, -1.5 * scale);
  EXPECT_LT(tetDeterminant(*tiny.mesh, 0, tiny.restPositions + tinyDx), 0.0);
  EXPECT_GT(tetDeterminant(*tiny.mesh, 0, tiny.restPositions + tinyDx), -1e-12);

  const double unitAlpha = unit.energy->computeMaxStepSize(zeroU, unitDx);
  const double tinyAlpha = tiny.energy->computeMaxStepSize(zeroU, tinyDx);
  EXPECT_LT(unitAlpha, 1.0);
  EXPECT_GT(unitAlpha, 0.0);
  EXPECT_NEAR(tinyAlpha, unitAlpha, 1e-9);
  EXPECT_GT(tetDeterminant(*tiny.mesh, 0, tiny.restPositions + tinyAlpha * tinyDx), 0.0);

  // The recovery path is scale invariant too: an already inverted tiny tet is
  // reported as illegal, not accepted because its determinant is "small".
  const ES::VXd invertedU = makeTetFlipDirection(4, 3, -2.0 * scale);
  const auto observation = tiny.assembler->computeMaxStepObservation((tiny.restPositions + invertedU).eval().data(), tinyDx.data());
  EXPECT_TRUE(observation.hasIllegalInitialState);
}
