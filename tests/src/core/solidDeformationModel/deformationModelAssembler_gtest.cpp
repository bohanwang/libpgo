#include <gtest/gtest.h>

#include "deformationModelAssembler.h"
#include "deformationModel.h"
#include "deformationModelManager.h"
#include "pgoLogging.h"
#include "simulationMesh.h"
#include "plasticModel3DDeformationGradient.h"
#include "cubicMesh.h"
#include "tetMesh.h"
#include "triMeshGeo.h"

#include <cmath>
#include <memory>

namespace
{
namespace ES = pgo::EigenSupport;
using pgo::SolidDeformationModel::DeformationModelAssembler;
using pgo::SolidDeformationModel::DeformationModelElasticMaterial;
using pgo::SolidDeformationModel::DeformationModelManager;
using pgo::SolidDeformationModel::DeformationModelPlasticMaterial;
using pgo::SolidDeformationModel::PlasticModel3DDeformationGradient;
using pgo::SolidDeformationModel::SimulationMesh;
using pgo::SolidDeformationModel::SimulationMeshENuhMaterial;
using pgo::SolidDeformationModel::SimulationMeshENuMaterial;
using pgo::SolidDeformationModel::SimulationMeshType;

constexpr const char *kTorusVegPath = LIBPGO_TEST_TORUS_VEG;
constexpr const char *kShellObjPath = LIBPGO_TEST_SHELL_OBJ;
constexpr const char *kCubicBoxVegPath = LIBPGO_TEST_CUBIC_BOX_VEG;

const double *dataOrNull(const ES::VXd &v)
{
  return v.size() ? v.data() : nullptr;
}

ES::VXd makePerturbedRestPositions(const SimulationMesh &mesh)
{
  ES::VXd x(mesh.getNumVertices() * 3);
  for (int vi = 0; vi < mesh.getNumVertices(); vi++) {
    double p[3];
    mesh.getVertex(vi, p);
    x[vi * 3 + 0] = p[0] + 1e-3 * ((vi % 3) - 1);
    x[vi * 3 + 1] = p[1] + 5e-4 * ((vi % 5) - 2);
    x[vi * 3 + 2] = p[2] + 7.5e-4 * ((vi % 7) - 3);
  }
  return x;
}

void expectAllFinite(const ES::VXd &v)
{
  for (Eigen::Index i = 0; i < v.size(); i++) {
    EXPECT_TRUE(std::isfinite(v[i])) << "Non-finite vector entry at " << i;
  }
}

void expectAllFinite(const ES::SpMatD &m)
{
  for (Eigen::Index i = 0; i < m.nonZeros(); i++) {
    EXPECT_TRUE(std::isfinite(m.valuePtr()[i])) << "Non-finite sparse entry at " << i;
  }
}

std::shared_ptr<SimulationMesh> makeSingleElementCubicSimulationMesh()
{
  const double vertices[] = {
    0.0, 0.0, 0.0,
    1.0, 0.0, 0.0,
    1.0, 1.0, 0.0,
    0.0, 1.0, 0.0,
    0.0, 0.0, 1.0,
    1.0, 0.0, 1.0,
    1.0, 1.0, 1.0,
    0.0, 1.0, 1.0,
  };
  const int elementVertices[] = { 0, 1, 2, 3, 4, 5, 6, 7 };
  const int elementMaterialIndices[] = { 0 };

  SimulationMeshENuMaterial baseMaterial(1200.0, 0.45);
  const pgo::SolidDeformationModel::SimulationMeshMaterial *materials[] = { &baseMaterial };

  return std::shared_ptr<SimulationMesh>(new SimulationMesh(
    8, vertices,
    1, 8, elementVertices,
    elementMaterialIndices, 1, materials,
    SimulationMeshType::CUBIC));
}

std::shared_ptr<SimulationMesh> makeSingleElementTetSimulationMesh()
{
  const double vertices[] = {
    0.0, 0.0, 0.0,
    1.0, 0.0, 0.0,
    0.0, 1.0, 0.0,
    0.0, 0.0, 1.0,
  };
  const int elementVertices[] = { 0, 1, 2, 3 };
  const int elementMaterialIndices[] = { 0 };
  SimulationMeshENuMaterial baseMaterial(1200.0, 0.45);
  const pgo::SolidDeformationModel::SimulationMeshMaterial *materials[] = { &baseMaterial };

  return std::shared_ptr<SimulationMesh>(new SimulationMesh(
    4, vertices,
    1, 4, elementVertices,
    elementMaterialIndices, 1, materials,
    SimulationMeshType::TET));
}

void expectAssemblerDirectionalDerivatives(const std::shared_ptr<SimulationMesh> &mesh)
{
  auto dmm = std::make_shared<DeformationModelManager>();
  dmm->setMesh(mesh.get());
  dmm->init(DeformationModelPlasticMaterial::VOLUMETRIC_DOF6, DeformationModelElasticMaterial::STABLE_NEO);
  DeformationModelAssembler assembler(dmm, nullptr);

  ES::VXd x = makePerturbedRestPositions(*mesh);
  ES::VXd plasticParams(dmm->getNumPlasticParameters() * mesh->getNumElements());
  const auto *plasticModel = dynamic_cast<const PlasticModel3DDeformationGradient *>(
    dmm->getDeformationModel(0)->getPlasticModel());
  ASSERT_NE(plasticModel, nullptr);
  const ES::M3d identity = ES::M3d::Identity();
  for (int ei = 0; ei < mesh->getNumElements(); ++ei)
    plasticModel->toParam(identity.data(), plasticParams.data() + ei * dmm->getNumPlasticParameters());
  ES::VXd elasticParams = ES::VXd::Zero(dmm->getNumElasticParameters() * mesh->getNumElements());

  ES::VXd direction(x.size());
  for (Eigen::Index i = 0; i < direction.size(); ++i)
    direction[i] = 0.15 + 0.07 * static_cast<double>((3 * i) % 7);
  direction.normalize();

  ES::VXd gradient = ES::VXd::Zero(assembler.getNumDOFs());
  assembler.computeGradient(x.data(), dataOrNull(plasticParams), dataOrNull(elasticParams), gradient.data());
  ES::SpMatD hessian = assembler.getHessianTemplate();
  assembler.computeHessian(x.data(), dataOrNull(plasticParams), dataOrNull(elasticParams), hessian);

  constexpr double kEpsilon = 1e-6;
  const ES::VXd xPlus = x + kEpsilon * direction;
  const ES::VXd xMinus = x - kEpsilon * direction;
  const double fdDirectionalGradient =
    (assembler.computeEnergy(xPlus.data(), dataOrNull(plasticParams), dataOrNull(elasticParams)) -
      assembler.computeEnergy(xMinus.data(), dataOrNull(plasticParams), dataOrNull(elasticParams))) /
    (2.0 * kEpsilon);
  EXPECT_NEAR(gradient.dot(direction), fdDirectionalGradient,
    2e-5 * std::max(1.0, std::abs(fdDirectionalGradient)));

  ES::VXd gradientPlus = ES::VXd::Zero(assembler.getNumDOFs());
  ES::VXd gradientMinus = ES::VXd::Zero(assembler.getNumDOFs());
  assembler.computeGradient(xPlus.data(), dataOrNull(plasticParams), dataOrNull(elasticParams), gradientPlus.data());
  assembler.computeGradient(xMinus.data(), dataOrNull(plasticParams), dataOrNull(elasticParams), gradientMinus.data());
  const ES::VXd fdHessianVector = (gradientPlus - gradientMinus) / (2.0 * kEpsilon);
  const ES::VXd analyticHessianVector = hessian * direction;
  EXPECT_LT((analyticHessianVector - fdHessianVector).norm(),
    3e-4 * std::max(1.0, fdHessianVector.norm()));
}
}

TEST(DeformationModelAssemblerGTest, TetAssemblerRegression)
{
  pgo::Logging::init();

  pgo::VolumetricMeshes::TetMesh tetMesh(kTorusVegPath);
  std::shared_ptr<SimulationMesh> mesh(pgo::SolidDeformationModel::loadTetMesh(&tetMesh));
  ASSERT_NE(mesh, nullptr);

  auto dmm = std::make_shared<DeformationModelManager>();
  dmm->setMesh(mesh.get());
  dmm->init(DeformationModelPlasticMaterial::VOLUMETRIC_DOF6, DeformationModelElasticMaterial::STABLE_NEO);

  auto assembler = std::make_shared<DeformationModelAssembler>(dmm, nullptr);

  ES::VXd x = makePerturbedRestPositions(*mesh);
  ES::VXd plasticParams(dmm->getNumPlasticParameters() * mesh->getNumElements());
  ES::VXd elasticParams(dmm->getNumElasticParameters() * mesh->getNumElements());
  elasticParams.setZero();

  const auto *plasticModel = dynamic_cast<const PlasticModel3DDeformationGradient *>(
    dmm->getDeformationModel(0)->getPlasticModel());
  ASSERT_NE(plasticModel, nullptr);

  ES::M3d identity = ES::M3d::Identity();
  for (int ei = 0; ei < mesh->getNumElements(); ei++) {
    plasticModel->toParam(identity.data(), plasticParams.data() + ei * dmm->getNumPlasticParameters());
  }

  ES::VXd grad = ES::VXd::Zero(assembler->getNumDOFs());
  assembler->computeGradient(x.data(), dataOrNull(plasticParams), dataOrNull(elasticParams), grad.data());
  EXPECT_EQ(grad.size(), assembler->getNumDOFs());
  expectAllFinite(grad);

  ES::SpMatD hess = assembler->getHessianTemplate();
  assembler->computeHessian(x.data(), dataOrNull(plasticParams), dataOrNull(elasticParams), hess);
  EXPECT_EQ(hess.rows(), assembler->getNumDOFs());
  EXPECT_EQ(hess.cols(), assembler->getNumDOFs());
  expectAllFinite(hess);

  ES::SpMatD dfda = assembler->get_dfda_Template();
  assembler->compute_df_da(x.data(), dataOrNull(plasticParams), dataOrNull(elasticParams), dfda);
  EXPECT_EQ(dfda.rows(), assembler->getNumDOFs());
  EXPECT_EQ(dfda.cols(), mesh->getNumElements() * dmm->getNumPlasticParameters());
  expectAllFinite(dfda);
}

TEST(DeformationModelAssemblerGTest, ShellAssemblerRegression)
{
  pgo::Logging::init();

  pgo::Mesh::TriMeshGeo surfaceMesh;
  ASSERT_TRUE(surfaceMesh.load(kShellObjPath));

  SimulationMeshENuhMaterial mat(1000.0, 0.45, 1e-3);
  std::shared_ptr<SimulationMesh> mesh(pgo::SolidDeformationModel::loadShellMesh(surfaceMesh, &mat));
  ASSERT_NE(mesh, nullptr);

  auto dmm = std::make_shared<DeformationModelManager>();
  dmm->setMesh(mesh.get());
  dmm->init(DeformationModelPlasticMaterial::SHELL_FF_DOF1, DeformationModelElasticMaterial::KOITER_STVK);

  auto assembler = std::make_shared<DeformationModelAssembler>(dmm, nullptr);

  ES::VXd x = makePerturbedRestPositions(*mesh);
  ES::VXd plasticParams = ES::VXd::Constant(dmm->getNumPlasticParameters() * mesh->getNumElements(), 1.2);
  ES::VXd elasticParams(dmm->getNumElasticParameters() * mesh->getNumElements());

  ASSERT_EQ(dmm->getNumElasticParameters(), 5);
  for (int ei = 0; ei < mesh->getNumElements(); ei++) {
    elasticParams.segment<5>(ei * 5) << 20000.0, 0.45, 10000.0, 0.3, 1e-3;
  }

  ES::VXd grad = ES::VXd::Zero(assembler->getNumDOFs());
  assembler->computeGradient(x.data(), dataOrNull(plasticParams), dataOrNull(elasticParams), grad.data());
  EXPECT_EQ(grad.size(), assembler->getNumDOFs());
  expectAllFinite(grad);

  ES::SpMatD hess = assembler->getHessianTemplate();
  assembler->computeHessian(x.data(), dataOrNull(plasticParams), dataOrNull(elasticParams), hess);
  EXPECT_EQ(hess.rows(), assembler->getNumDOFs());
  EXPECT_EQ(hess.cols(), assembler->getNumDOFs());
  expectAllFinite(hess);

  ES::SpMatD dfda = assembler->get_dfda_Template();
  assembler->compute_df_da(x.data(), dataOrNull(plasticParams), dataOrNull(elasticParams), dfda);
  EXPECT_EQ(dfda.rows(), assembler->getNumDOFs());
  EXPECT_EQ(dfda.cols(), mesh->getNumElements() * dmm->getNumPlasticParameters());
  expectAllFinite(dfda);

  if (dmm->getNumElasticParameters() > 0) {
    ES::SpMatD dfdb = assembler->get_dfdb_Template();
    assembler->compute_df_db(x.data(), dataOrNull(plasticParams), dataOrNull(elasticParams), dfdb);
    EXPECT_EQ(dfdb.rows(), assembler->getNumDOFs());
    EXPECT_EQ(dfdb.cols(), mesh->getNumElements() * dmm->getNumElasticParameters());
    expectAllFinite(dfdb);
  }
}

TEST(DeformationModelAssemblerGTest, CubicAssemblerSmokeRegression)
{
  pgo::Logging::init();

  pgo::VolumetricMeshes::CubicMesh cubicMesh(kCubicBoxVegPath);
  std::shared_ptr<SimulationMesh> mesh(pgo::SolidDeformationModel::loadCubicMesh(&cubicMesh));
  ASSERT_NE(mesh, nullptr);

  auto dmm = std::make_shared<DeformationModelManager>();
  dmm->setMesh(mesh.get());
  dmm->init(DeformationModelPlasticMaterial::VOLUMETRIC_DOF6, DeformationModelElasticMaterial::STABLE_NEO);

  auto assembler = std::make_shared<DeformationModelAssembler>(dmm, nullptr);

  ES::VXd x = makePerturbedRestPositions(*mesh);
  ES::VXd plasticParams(dmm->getNumPlasticParameters() * mesh->getNumElements());
  ES::VXd elasticParams(dmm->getNumElasticParameters() * mesh->getNumElements());
  elasticParams.setZero();

  const auto *plasticModel = dynamic_cast<const PlasticModel3DDeformationGradient *>(
    dmm->getDeformationModel(0)->getPlasticModel());
  ASSERT_NE(plasticModel, nullptr);

  ES::M3d identity = ES::M3d::Identity();
  for (int ei = 0; ei < mesh->getNumElements(); ei++) {
    plasticModel->toParam(identity.data(), plasticParams.data() + ei * dmm->getNumPlasticParameters());
  }

  ES::VXd grad = ES::VXd::Zero(assembler->getNumDOFs());
  assembler->computeGradient(x.data(), dataOrNull(plasticParams), dataOrNull(elasticParams), grad.data());
  EXPECT_EQ(grad.size(), assembler->getNumDOFs());
  expectAllFinite(grad);

  ES::SpMatD hess = assembler->getHessianTemplate();
  assembler->computeHessian(x.data(), dataOrNull(plasticParams), dataOrNull(elasticParams), hess);
  EXPECT_EQ(hess.rows(), assembler->getNumDOFs());
  EXPECT_EQ(hess.cols(), assembler->getNumDOFs());
  expectAllFinite(hess);

  ES::SpMatD dfda = assembler->get_dfda_Template();
  assembler->compute_df_da(x.data(), dataOrNull(plasticParams), dataOrNull(elasticParams), dfda);
  EXPECT_EQ(dfda.rows(), assembler->getNumDOFs());
  EXPECT_EQ(dfda.cols(), mesh->getNumElements() * dmm->getNumPlasticParameters());
  expectAllFinite(dfda);
}

TEST(DeformationModelAssemblerGTest, CubicAssemblerMaterialParamRegression)
{
  pgo::Logging::init();

  std::shared_ptr<SimulationMesh> mesh = makeSingleElementCubicSimulationMesh();
  ASSERT_NE(mesh, nullptr);

  pgo::SolidDeformationModel::SimulationMeshHillMaterial hillMaterial(2500.0, 0.35, 1.0);
  mesh->appendMaterialToAllElements(&hillMaterial);

  ES::VXd elementFiberDirections = ES::VXd::Zero(mesh->getNumElements() * 3);
  for (int ei = 0; ei < mesh->getNumElements(); ei++) {
    elementFiberDirections.segment<3>(ei * 3) << 1.0, 0.0, 0.0;
  }

  ES::VXd vertexFiberDirections = ES::VXd::Zero(mesh->getNumVertices() * 3);
  for (int vi = 0; vi < mesh->getNumVertices(); vi++) {
    vertexFiberDirections.segment<3>(vi * 3) << 1.0, 0.0, 0.0;
  }

  auto dmm = std::make_shared<DeformationModelManager>();
  dmm->setMesh(mesh.get(), elementFiberDirections.data(), vertexFiberDirections.data());
  dmm->init(DeformationModelPlasticMaterial::VOLUMETRIC_DOF6, DeformationModelElasticMaterial::HILL_STABLE_NEO);

  ASSERT_EQ(dmm->getNumElasticParameters(), 1);

  auto assembler = std::make_shared<DeformationModelAssembler>(dmm, nullptr);

  ES::VXd x = makePerturbedRestPositions(*mesh);
  ES::VXd plasticParams(dmm->getNumPlasticParameters() * mesh->getNumElements());
  ES::VXd elasticParams = ES::VXd::Constant(dmm->getNumElasticParameters() * mesh->getNumElements(), 0.75);

  const auto *plasticModel = dynamic_cast<const PlasticModel3DDeformationGradient *>(
    dmm->getDeformationModel(0)->getPlasticModel());
  ASSERT_NE(plasticModel, nullptr);

  ES::M3d identity = ES::M3d::Identity();
  for (int ei = 0; ei < mesh->getNumElements(); ei++) {
    plasticModel->toParam(identity.data(), plasticParams.data() + ei * dmm->getNumPlasticParameters());
  }

  ES::SpMatD dfdb = assembler->get_dfdb_Template();
  assembler->compute_df_db(x.data(), dataOrNull(plasticParams), dataOrNull(elasticParams), dfdb);
  EXPECT_EQ(dfdb.rows(), assembler->getNumDOFs());
  EXPECT_EQ(dfdb.cols(), mesh->getNumElements() * dmm->getNumElasticParameters());
  expectAllFinite(dfdb);
  EXPECT_GT(dfdb.norm(), 0.0);
}

TEST(DeformationModelAssemblerGTest, TetAndCubicDirectionalDerivativesMatchFiniteDifferences)
{
  pgo::Logging::init();

  expectAssemblerDirectionalDerivatives(makeSingleElementTetSimulationMesh());
  expectAssemblerDirectionalDerivatives(makeSingleElementCubicSimulationMesh());
}
