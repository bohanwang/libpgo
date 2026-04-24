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

std::size_t countOccurrences(const std::string &haystack, const std::string &needle)
{
  std::size_t count = 0;
  std::size_t pos = 0;
  while ((pos = haystack.find(needle, pos)) != std::string::npos) {
    count++;
    pos += needle.size();
  }
  return count;
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

TEST(DeformationModelEnergyMaxStepGTest, ZeroDirectionReturnsOneAndDoesNotClamp)
{
  EnergyFixture fixture = makeSingleTetFixture();
  const ES::VXd x = ES::VXd::Zero(fixture.restPositions.size());
  const ES::VXd dx = ES::VXd::Zero(fixture.restPositions.size());

  EXPECT_DOUBLE_EQ(fixture.energy->computeMaxStepSize(x, dx), 1.0);
  EXPECT_EQ(fixture.energy->getMaterialClampCount(), 0);
}

TEST(DeformationModelEnergyMaxStepGTest, TetPureTranslationReturnsOne)
{
  EnergyFixture fixture = makeSingleTetFixture();
  const ES::VXd x = ES::VXd::Zero(fixture.restPositions.size());
  ES::VXd dx = ES::VXd::Zero(fixture.restPositions.size());
  applyUniformTranslation(dx, 1.0, 2.0, 3.0);

  EXPECT_DOUBLE_EQ(fixture.energy->computeMaxStepSize(x, dx), 1.0);
  EXPECT_EQ(fixture.energy->getMaterialClampCount(), 0);
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
  EXPECT_EQ(fixture.energy->getMaterialClampCount(), 1);
}

TEST(DeformationModelEnergyMaxStepGTest, DisabledMaterialMaxStepSkipsTetClamp)
{
  EnergyFixture fixture = makeSingleTetFixture();
  fixture.energy->setEnableMaterialMaxStep(false);

  const ES::VXd x = ES::VXd::Zero(fixture.restPositions.size());
  const ES::VXd dx = makeTetFlipDirection(fixture.mesh->getNumVertices(), 3, -2.0);

  EXPECT_DOUBLE_EQ(fixture.energy->computeMaxStepSize(x, dx), 1.0);
  EXPECT_EQ(fixture.energy->getMaterialClampCount(), 0);
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
  EXPECT_EQ(fixture.energy->getMaterialClampCount(), 2);
  EXPECT_EQ(countOccurrences(logOutput, "Phase 1.5 material max step encountered illegal initial state"), 2u);
}

TEST(DeformationModelEnergyMaxStepGTest, TetSmallAlphaWarnsAndTracksSolveMinimumAlpha)
{
  EnergyFixture fixture = makeSingleTetFixture();
  const ES::VXd x = ES::VXd::Zero(fixture.restPositions.size());
  const ES::VXd dx = makeTetFlipDirection(fixture.mesh->getNumVertices(), 3, -200.0);

  testing::internal::CaptureStdout();
  const double alpha = fixture.energy->computeMaxStepSize(x, dx);
  const std::string logOutput = testing::internal::GetCapturedStdout();

  EXPECT_GT(alpha, 0.0);
  EXPECT_LT(alpha, 0.01);
  EXPECT_EQ(fixture.energy->getMaterialClampCount(), 1);
  EXPECT_DOUBLE_EQ(fixture.energy->getMinMaterialFeasibleAlphaThisSolve(), alpha);
  EXPECT_NE(logOutput.find("materialFeasibleAlpha"), std::string::npos);
}

TEST(DeformationModelEnergyMaxStepGTest, ResetMaterialMaxStepStatsClearsCountAndMinimumAlpha)
{
  EnergyFixture fixture = makeSingleTetFixture();
  const ES::VXd x = ES::VXd::Zero(fixture.restPositions.size());
  const ES::VXd dx = makeTetFlipDirection(fixture.mesh->getNumVertices(), 3, -2.0);

  const double alpha = fixture.energy->computeMaxStepSize(x, dx);
  ASSERT_LT(alpha, 1.0);
  ASSERT_EQ(fixture.energy->getMaterialClampCount(), 1);

  fixture.energy->resetMaterialMaxStepStats();
  EXPECT_EQ(fixture.energy->getMaterialClampCount(), 0);
  EXPECT_DOUBLE_EQ(fixture.energy->getMinMaterialFeasibleAlphaThisSolve(), 1.0);
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
  EXPECT_EQ(multiFixture.energy->getMaterialClampCount(), 1);
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
  EXPECT_EQ(fixture.energy->getMaterialClampCount(), 1);
}

TEST(DeformationModelEnergyMaxStepGTest, CubicFeasibleDirectionReturnsOne)
{
  EnergyFixture fixture = makeSingleCubicFixture();
  const ES::VXd x = ES::VXd::Zero(fixture.restPositions.size());
  const ES::VXd dx = makeCubicTopFaceDirection(fixture.mesh->getNumVertices(), 0, -0.2);

  EXPECT_DOUBLE_EQ(fixture.energy->computeMaxStepSize(x, dx), 1.0);
  EXPECT_EQ(fixture.energy->getMaterialClampCount(), 0);
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
  EXPECT_EQ(multiFixture.energy->getMaterialClampCount(), 1);
}

TEST(DeformationModelEnergyMaxStepGTest, ShellKeepsUnitStep)
{
  EnergyFixture fixture = makeShellFixture();
  const ES::VXd x = ES::VXd::Zero(fixture.restPositions.size());
  ES::VXd dx = ES::VXd::Zero(fixture.restPositions.size());
  applyUniformTranslation(dx, 0.1, -0.05, 0.2);

  EXPECT_DOUBLE_EQ(fixture.energy->computeMaxStepSize(x, dx), 1.0);
  EXPECT_EQ(fixture.energy->getMaterialClampCount(), 0);
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
  EXPECT_NEAR(std::static_pointer_cast<const pgo::Simulation::ImplicitBackwardEulerEnergy>(
    integrator.getInternalEnergy())->getMinFeasibleAlphaThisSolve(), materialAlpha, 1e-12);

  integrator.clearGeneralImplicitForceModel();
  integrator.addGeneralImplicitForceModel(std::make_shared<FixedMaxStepEnergy>(fixture.restPositions.size(), 0.25));
  integrator.assembleImplicitModels();
  EXPECT_DOUBLE_EQ(integrator.getInternalEnergy()->computeMaxStepSize(x, dx), 0.25);
  EXPECT_DOUBLE_EQ(std::static_pointer_cast<const pgo::Simulation::ImplicitBackwardEulerEnergy>(
    integrator.getInternalEnergy())->getMinFeasibleAlphaThisSolve(), 0.25);
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
  EXPECT_NEAR(std::static_pointer_cast<const pgo::Simulation::TRBDF2TimeIntegratorEnergy>(
    integrator.getTRStageEnergy())->getMinFeasibleAlphaThisSolve(), materialAlpha, 1e-12);

  integrator.clearGeneralImplicitForceModel();
  integrator.addGeneralImplicitForceModel(std::make_shared<FixedMaxStepEnergy>(fixture.restPositions.size(), 0.2));
  integrator.assembleImplicitModels();
  EXPECT_DOUBLE_EQ(integrator.getTRStageEnergy()->computeMaxStepSize(x, dx), 0.2);
  EXPECT_DOUBLE_EQ(std::static_pointer_cast<const pgo::Simulation::TRBDF2TimeIntegratorEnergy>(
    integrator.getTRStageEnergy())->getMinFeasibleAlphaThisSolve(), 0.2);
}
