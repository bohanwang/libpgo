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

ES::VXd makeRestPositions(const SimulationMesh &mesh)
{
  ES::VXd x(mesh.getNumVertices() * 3);
  for (int vi = 0; vi < mesh.getNumVertices(); vi++) {
    double p[3];
    mesh.getVertex(vi, p);
    x.segment<3>(vi * 3) = ES::V3d(p[0], p[1], p[2]);
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

std::unique_ptr<SimulationMesh> makeSingleElementCubicSimulationMesh()
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

  return std::unique_ptr<SimulationMesh>(new SimulationMesh(
    8, vertices,
    1, 8, elementVertices,
    elementMaterialIndices, 1, materials,
    SimulationMeshType::CUBIC));
}
}

TEST(DeformationModelAssemblerGTest, TetAssemblerRegression)
{
  pgo::Logging::init();

  pgo::VolumetricMeshes::TetMesh tetMesh(kTorusVegPath);
  std::unique_ptr<SimulationMesh> mesh = pgo::SolidDeformationModel::loadTetMesh(&tetMesh);
  ASSERT_NE(mesh, nullptr);

  const int nele = mesh->getNumElements();
  const int nvtx = mesh->getNumVertices();
  const int n3 = nvtx * 3;

  auto dmm = std::make_unique<DeformationModelManager>();
  dmm->setMesh(std::move(mesh), nullptr, nullptr);
  dmm->init(DeformationModelPlasticMaterial::VOLUMETRIC_DOF6, DeformationModelElasticMaterial::STABLE_NEO);

  const int numPlasticParams = dmm->getNumPlasticParameters();
  const int numElasticParams = dmm->getNumElasticParameters();
  const auto *plasticModel = dynamic_cast<const PlasticModel3DDeformationGradient *>(dmm->getDeformationModel(0)->getPlasticModel());

  auto assembler = std::make_unique<DeformationModelAssembler>(std::move(dmm), nullptr);

  ES::VXd x = makePerturbedRestPositions(*assembler->getDeformationModelManager().getMesh());
  ES::VXd plasticParams(numPlasticParams * nele);
  ES::VXd elasticParams(numElasticParams * nele);
  elasticParams.setZero();

  ASSERT_NE(plasticModel, nullptr);

  ES::M3d identity = ES::M3d::Identity();
  for (int ei = 0; ei < nele; ei++) {
    plasticModel->toParam(identity.data(), plasticParams.data() + ei * numPlasticParams);
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
  EXPECT_EQ(dfda.cols(), nele * numPlasticParams);
  expectAllFinite(dfda);
}

TEST(DeformationModelAssemblerGTest, TetVonMisesStressIsZeroAtRestAndNonzeroUnderStretch)
{
  pgo::Logging::init();

  pgo::VolumetricMeshes::TetMesh tetMesh(kTorusVegPath);
  std::unique_ptr<SimulationMesh> mesh = pgo::SolidDeformationModel::loadTetMesh(&tetMesh);
  ASSERT_NE(mesh, nullptr);

  const int nele = mesh->getNumElements();
  const int nvtx = mesh->getNumVertices();

  auto dmm = std::make_unique<DeformationModelManager>();
  dmm->setMesh(std::move(mesh), nullptr, nullptr);
  dmm->init(DeformationModelPlasticMaterial::VOLUMETRIC_DOF6, DeformationModelElasticMaterial::STABLE_NEO);

  const int numPlasticParams = dmm->getNumPlasticParameters();
  const int numElasticParams = dmm->getNumElasticParameters();
  const auto *plasticModel = dynamic_cast<const PlasticModel3DDeformationGradient *>(dmm->getDeformationModel(0)->getPlasticModel());

  auto assembler = std::make_unique<DeformationModelAssembler>(std::move(dmm), nullptr);

  ES::VXd plasticParams(numPlasticParams * nele);
  ES::VXd elasticParams(numElasticParams * nele);
  elasticParams.setZero();

  ASSERT_NE(plasticModel, nullptr);

  ES::M3d identity = ES::M3d::Identity();
  for (int ei = 0; ei < nele; ei++) {
    plasticModel->toParam(identity.data(), plasticParams.data() + ei * numPlasticParams);
  }

  const auto &meshPtr = *assembler->getDeformationModelManager().getMesh();
  ES::VXd rest = makeRestPositions(meshPtr);
  ES::VXd stresses = ES::VXd::Constant(nele, -1.0);
  assembler->computeVonMisesStresses(rest.data(), dataOrNull(plasticParams), dataOrNull(elasticParams), stresses.data());
  expectAllFinite(stresses);
  EXPECT_LE(stresses.cwiseAbs().maxCoeff(), 1e-8);

  ES::VXd stretched = rest;
  for (int vi = 0; vi < nvtx; vi++) {
    stretched[vi * 3] *= 1.01;
  }
  stresses.setConstant(-1.0);
  assembler->computeVonMisesStresses(stretched.data(), dataOrNull(plasticParams), dataOrNull(elasticParams), stresses.data());
  expectAllFinite(stresses);
  EXPECT_GE(stresses.minCoeff(), 0.0);
  EXPECT_GT(stresses.maxCoeff(), 1e-8);
}

TEST(DeformationModelAssemblerGTest, ShellAssemblerRegression)
{
  pgo::Logging::init();

  pgo::Mesh::TriMeshGeo surfaceMesh;
  ASSERT_TRUE(surfaceMesh.load(kShellObjPath));

  SimulationMeshENuhMaterial mat(1000.0, 0.45, 1e-3);
  std::unique_ptr<SimulationMesh> mesh = pgo::SolidDeformationModel::loadShellMesh(surfaceMesh, &mat);
  ASSERT_NE(mesh, nullptr);

  const int nele = mesh->getNumElements();

  auto dmm = std::make_unique<DeformationModelManager>();
  dmm->setMesh(std::move(mesh), nullptr, nullptr);
  dmm->init(DeformationModelPlasticMaterial::SHELL_FF_DOF1, DeformationModelElasticMaterial::KOITER_STVK);

  const int numPlasticParams = dmm->getNumPlasticParameters();
  const int numElasticParams = dmm->getNumElasticParameters();

  auto assembler = std::make_unique<DeformationModelAssembler>(std::move(dmm), nullptr);

  const auto &meshPtr = *assembler->getDeformationModelManager().getMesh();
  ES::VXd x = makePerturbedRestPositions(meshPtr);
  ES::VXd plasticParams = ES::VXd::Constant(numPlasticParams * nele, 1.2);
  ES::VXd elasticParams(numElasticParams * nele);

  ASSERT_EQ(numElasticParams, 5);
  for (int ei = 0; ei < nele; ei++) {
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
  EXPECT_EQ(dfda.cols(), nele * numPlasticParams);
  expectAllFinite(dfda);

  if (numElasticParams > 0) {
    ES::SpMatD dfdb = assembler->get_dfdb_Template();
    assembler->compute_df_db(x.data(), dataOrNull(plasticParams), dataOrNull(elasticParams), dfdb);
    EXPECT_EQ(dfdb.rows(), assembler->getNumDOFs());
    EXPECT_EQ(dfdb.cols(), nele * numElasticParams);
    expectAllFinite(dfdb);
  }
}

TEST(DeformationModelAssemblerGTest, CubicAssemblerSmokeRegression)
{
  pgo::Logging::init();

  pgo::VolumetricMeshes::CubicMesh cubicMesh(kCubicBoxVegPath);
  std::unique_ptr<SimulationMesh> mesh = pgo::SolidDeformationModel::loadCubicMesh(&cubicMesh);
  ASSERT_NE(mesh, nullptr);

  const int nele = mesh->getNumElements();

  auto dmm = std::make_unique<DeformationModelManager>();
  dmm->setMesh(std::move(mesh), nullptr, nullptr);
  dmm->init(DeformationModelPlasticMaterial::VOLUMETRIC_DOF6, DeformationModelElasticMaterial::STABLE_NEO);

  const int numPlasticParams = dmm->getNumPlasticParameters();
  const int numElasticParams = dmm->getNumElasticParameters();
  const auto *plasticModel = dynamic_cast<const PlasticModel3DDeformationGradient *>(dmm->getDeformationModel(0)->getPlasticModel());

  auto assembler = std::make_unique<DeformationModelAssembler>(std::move(dmm), nullptr);

  const auto &meshPtr = *assembler->getDeformationModelManager().getMesh();
  ES::VXd x = makePerturbedRestPositions(meshPtr);
  ES::VXd plasticParams(numPlasticParams * nele);
  ES::VXd elasticParams(numElasticParams * nele);
  elasticParams.setZero();

  ASSERT_NE(plasticModel, nullptr);

  ES::M3d identity = ES::M3d::Identity();
  for (int ei = 0; ei < nele; ei++) {
    plasticModel->toParam(identity.data(), plasticParams.data() + ei * numPlasticParams);
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
  EXPECT_EQ(dfda.cols(), nele * numPlasticParams);
  expectAllFinite(dfda);
}

TEST(DeformationModelAssemblerGTest, CubicAssemblerMaterialParamRegression)
{
  pgo::Logging::init();

  std::unique_ptr<SimulationMesh> mesh = makeSingleElementCubicSimulationMesh();
  ASSERT_NE(mesh, nullptr);

  pgo::SolidDeformationModel::SimulationMeshHillMaterial hillMaterial(2500.0, 0.35, 1.0);
  mesh->appendMaterialToAllElements(&hillMaterial);

  const int nele = mesh->getNumElements();
  const int nvtx = mesh->getNumVertices();

  ES::VXd elementFiberDirections = ES::VXd::Zero(nele * 3);
  for (int ei = 0; ei < nele; ei++) {
    elementFiberDirections.segment<3>(ei * 3) << 1.0, 0.0, 0.0;
  }

  ES::VXd vertexFiberDirections = ES::VXd::Zero(nvtx * 3);
  for (int vi = 0; vi < nvtx; vi++) {
    vertexFiberDirections.segment<3>(vi * 3) << 1.0, 0.0, 0.0;
  }

  auto dmm = std::make_unique<DeformationModelManager>();
  dmm->setMesh(std::move(mesh), elementFiberDirections.data(), vertexFiberDirections.data());
  dmm->init(DeformationModelPlasticMaterial::VOLUMETRIC_DOF6, DeformationModelElasticMaterial::HILL_STABLE_NEO);

  const int numPlasticParams = dmm->getNumPlasticParameters();
  const int numElasticParams = dmm->getNumElasticParameters();
  ASSERT_EQ(numElasticParams, 1);

  auto assembler = std::make_unique<DeformationModelAssembler>(std::move(dmm), nullptr);

  const auto &meshPtr = *assembler->getDeformationModelManager().getMesh();
  ES::VXd x = makePerturbedRestPositions(meshPtr);
  ES::VXd plasticParams(numPlasticParams * nele);
  ES::VXd elasticParams = ES::VXd::Constant(numElasticParams * nele, 0.75);

  const auto *plasticModel = dynamic_cast<const PlasticModel3DDeformationGradient *>(
    assembler->getDeformationModelManager().getDeformationModel(0)->getPlasticModel());
  ASSERT_NE(plasticModel, nullptr);

  ES::M3d identity = ES::M3d::Identity();
  for (int ei = 0; ei < nele; ei++) {
    plasticModel->toParam(identity.data(), plasticParams.data() + ei * numPlasticParams);
  }

  ES::SpMatD dfdb = assembler->get_dfdb_Template();
  assembler->compute_df_db(x.data(), dataOrNull(plasticParams), dataOrNull(elasticParams), dfdb);
  EXPECT_EQ(dfdb.rows(), assembler->getNumDOFs());
  EXPECT_EQ(dfdb.cols(), nele * numElasticParams);
  expectAllFinite(dfdb);
  EXPECT_GT(dfdb.norm(), 0.0);
}
