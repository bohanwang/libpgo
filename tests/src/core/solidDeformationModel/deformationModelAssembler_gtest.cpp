#include <gtest/gtest.h>

#include "deformationModelAssembler.h"
#include "deformationModel.h"
#include "deformationModelManager.h"
#include "pgoLogging.h"
#include "simulationMesh.h"
#include "plasticModel3DDeformationGradient.h"
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

constexpr const char *kTorusVegPath = LIBPGO_TEST_TORUS_VEG;
constexpr const char *kShellObjPath = LIBPGO_TEST_SHELL_OBJ;

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
