#include <gtest/gtest.h>

#include "deformationModelAssembler.h"
#include "deformationModelEnergy.h"
#include "potentialEnergies.h"
#include "deformationModel.h"
#include "deformationModelManager.h"
#include "pgoLogging.h"
#include "simulationMesh.h"
#include "plasticModel3DDeformationGradient.h"
#include "cubicMesh.h"
#include "tetMesh.h"
#include "triMeshGeo.h"

#include <algorithm>
#include <cmath>
#include <cstring>
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

struct AssemblerCase
{
  std::shared_ptr<SimulationMesh> mesh;
  std::shared_ptr<DeformationModelManager> dmm;
  std::shared_ptr<DeformationModelAssembler> assembler;
  ES::VXd x, plasticParams, elasticParams;
};

void fillIdentityPlasticity(AssemblerCase &c)
{
  c.plasticParams.resize(c.dmm->getNumPlasticParameters() * c.mesh->getNumElements());
  const auto *plasticModel = dynamic_cast<const PlasticModel3DDeformationGradient *>(
    c.dmm->getDeformationModel(0)->getPlasticModel());
  ASSERT_NE(plasticModel, nullptr);
  const ES::M3d identity = ES::M3d::Identity();
  for (int ei = 0; ei < c.mesh->getNumElements(); ei++)
    plasticModel->toParam(identity.data(), c.plasticParams.data() + ei * c.dmm->getNumPlasticParameters());
}

AssemblerCase makeVolumetricCase(std::shared_ptr<SimulationMesh> mesh)
{
  AssemblerCase c;
  c.mesh = std::move(mesh);
  c.dmm = std::make_shared<DeformationModelManager>();
  c.dmm->setMesh(c.mesh.get());
  c.dmm->init(DeformationModelPlasticMaterial::VOLUMETRIC_DOF6, DeformationModelElasticMaterial::STABLE_NEO);
  c.assembler = std::make_shared<DeformationModelAssembler>(c.dmm, nullptr);
  c.x = makePerturbedRestPositions(*c.mesh);
  fillIdentityPlasticity(c);
  c.elasticParams = ES::VXd::Zero(c.dmm->getNumElasticParameters() * c.mesh->getNumElements());
  return c;
}

AssemblerCase makeShellCase()
{
  AssemblerCase c;
  pgo::Mesh::TriMeshGeo surfaceMesh;
  EXPECT_TRUE(surfaceMesh.load(kShellObjPath));
  SimulationMeshENuhMaterial mat(1000.0, 0.45, 1e-3);
  c.mesh.reset(pgo::SolidDeformationModel::loadShellMesh(surfaceMesh, &mat));
  c.dmm = std::make_shared<DeformationModelManager>();
  c.dmm->setMesh(c.mesh.get());
  c.dmm->init(DeformationModelPlasticMaterial::SHELL_FF_DOF1, DeformationModelElasticMaterial::KOITER_STVK);
  c.assembler = std::make_shared<DeformationModelAssembler>(c.dmm, nullptr);
  c.x = makePerturbedRestPositions(*c.mesh);
  c.plasticParams = ES::VXd::Constant(c.dmm->getNumPlasticParameters() * c.mesh->getNumElements(), 1.2);
  c.elasticParams.resize(c.dmm->getNumElasticParameters() * c.mesh->getNumElements());
  EXPECT_EQ(c.dmm->getNumElasticParameters(), 5);
  for (int ei = 0; ei < c.mesh->getNumElements(); ei++)
    c.elasticParams.segment<5>(ei * 5) << 20000.0, 0.45, 10000.0, 0.3, 1e-3;
  return c;
}

struct AssemblyResult
{
  double energy = 0.0;
  ES::VXd gradient;
  ES::SpMatD hessian;
};

AssemblyResult assemble(const AssemblerCase &c)
{
  AssemblyResult r;
  r.energy = c.assembler->computeEnergy(c.x.data(), dataOrNull(c.plasticParams), dataOrNull(c.elasticParams));
  r.gradient = ES::VXd::Zero(c.assembler->getNumDOFs());
  c.assembler->computeGradient(c.x.data(), dataOrNull(c.plasticParams), dataOrNull(c.elasticParams), r.gradient.data());
  r.hessian = c.assembler->getHessianTemplate();
  c.assembler->computeHessian(c.x.data(), dataOrNull(c.plasticParams), dataOrNull(c.elasticParams), r.hessian);
  return r;
}

bool bitwiseEqual(const ES::VXd &a, const ES::VXd &b)
{
  return a.size() == b.size() && std::memcmp(a.data(), b.data(), sizeof(double) * a.size()) == 0;
}

bool bitwiseEqual(const ES::SpMatD &a, const ES::SpMatD &b)
{
  return a.nonZeros() == b.nonZeros() && std::memcmp(a.valuePtr(), b.valuePtr(), sizeof(double) * a.nonZeros()) == 0;
}

// The deterministic mode only changes the summation order, so it must agree
// with the default mode to roundoff and reproduce itself exactly.
void expectDeterministicAssemblyMatches(const AssemblerCase &c)
{
  c.assembler->setDeterministic(false);
  EXPECT_FALSE(c.assembler->isDeterministic());
  const AssemblyResult base = assemble(c);
  expectAllFinite(base.gradient);
  expectAllFinite(base.hessian);
  EXPECT_GT(base.gradient.norm(), 0.0);

  c.assembler->setDeterministic(true);
  EXPECT_TRUE(c.assembler->isDeterministic());
  const AssemblyResult det = assemble(c);
  const AssemblyResult again = assemble(c);

  EXPECT_NEAR(det.energy, base.energy, 1e-10 * std::max(1.0, std::abs(base.energy)));
  EXPECT_LT((det.gradient - base.gradient).norm(), 1e-10 * std::max(1.0, base.gradient.norm()));
  ASSERT_EQ(det.hessian.nonZeros(), base.hessian.nonZeros());
  const Eigen::Map<const ES::VXd> detValues(det.hessian.valuePtr(), det.hessian.nonZeros());
  const Eigen::Map<const ES::VXd> baseValues(base.hessian.valuePtr(), base.hessian.nonZeros());
  EXPECT_LT((detValues - baseValues).norm(), 1e-10 * std::max(1.0, baseValues.norm()));

  EXPECT_EQ(det.energy, again.energy);
  EXPECT_TRUE(bitwiseEqual(det.gradient, again.gradient));
  EXPECT_TRUE(bitwiseEqual(det.hessian, again.hessian));

  // switching back restores the default path
  c.assembler->setDeterministic(false);
  const AssemblyResult back = assemble(c);
  EXPECT_NEAR(back.energy, base.energy, 1e-10 * std::max(1.0, std::abs(base.energy)));
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

TEST(DeformationModelAssemblerGTest, DeterministicAssemblyMatchesDefaultAndReproducesBitwise)
{
  pgo::Logging::init();

  pgo::VolumetricMeshes::TetMesh tetMesh(kTorusVegPath);
  std::shared_ptr<SimulationMesh> tet(pgo::SolidDeformationModel::loadTetMesh(&tetMesh));
  ASSERT_NE(tet, nullptr);
  expectDeterministicAssemblyMatches(makeVolumetricCase(tet));

  pgo::VolumetricMeshes::CubicMesh cubicMesh(kCubicBoxVegPath);
  std::shared_ptr<SimulationMesh> cubic(pgo::SolidDeformationModel::loadCubicMesh(&cubicMesh));
  ASSERT_NE(cubic, nullptr);
  expectDeterministicAssemblyMatches(makeVolumetricCase(cubic));

  const AssemblerCase shell = makeShellCase();
  ASSERT_NE(shell.mesh, nullptr);
  expectDeterministicAssemblyMatches(shell);
}

TEST(DeformationModelAssemblerGTest, EnergyWithAggregateOffsetIndexesItsOwnDOFs)
{
  pgo::Logging::init();
  const AssemblerCase c = makeVolumetricCase(makeSingleElementTetSimulationMesh());
  ES::VXd rest(c.assembler->getNumDOFs());
  for (int vi = 0; vi < c.mesh->getNumVertices(); vi++) {
    double p[3];
    c.mesh->getVertex(vi, p);
    rest.segment<3>(vi * 3) = ES::V3d(p[0], p[1], p[2]);
  }
  const int n = c.assembler->getNumDOFs();
  const ES::VXd u = c.x - rest;  // displacement of the perturbed state

  auto makeEnergy = [&](int offset) {
    auto energy = std::make_shared<pgo::SolidDeformationModel::DeformationModelEnergy>(c.assembler, &rest, offset);
    energy->setPlasticParams(c.plasticParams);
    energy->setElasticParams(c.elasticParams);
    return energy;
  };

  const auto reference = makeEnergy(0);
  const double energy0 = reference->func(u);
  ES::VXd gradient0 = ES::VXd::Zero(n);
  reference->gradient(u, gradient0);
  ES::SpMatD hessian0;
  reference->createHessian(hessian0);
  reference->hessian(u, hessian0);
  EXPECT_GT(gradient0.norm(), 0.0);

  // Placed after 5 unrelated DOFs inside an aggregate, the energy must see the
  // same local vector and scatter into DOFs 5..5+n.
  constexpr int offset = 5;
  pgo::NonlinearOptimization::PotentialEnergies aggregate(offset + n);
  aggregate.addPotentialEnergy(makeEnergy(offset));
  aggregate.init();

  ES::VXd xAll = ES::VXd::Constant(offset + n, 7.0);  // garbage in the unrelated DOFs
  xAll.tail(n) = u;
  EXPECT_NEAR(aggregate.func(xAll), energy0, 1e-12 * std::max(1.0, std::abs(energy0)));

  ES::VXd gradientAll = ES::VXd::Zero(offset + n);
  aggregate.gradient(xAll, gradientAll);
  EXPECT_EQ(gradientAll.head(offset).norm(), 0.0);
  EXPECT_LT((gradientAll.tail(n) - gradient0).norm(), 1e-12 * std::max(1.0, gradient0.norm()));

  ES::SpMatD hessianAll;
  aggregate.createHessian(hessianAll);
  aggregate.hessian(xAll, hessianAll);
  const ES::MXd block = ES::MXd(hessianAll).block(offset, offset, n, n);
  EXPECT_LT((block - ES::MXd(hessian0)).norm(), 1e-12 * std::max(1.0, ES::MXd(hessian0).norm()));
  EXPECT_EQ(ES::MXd(hessianAll).topLeftCorner(offset, offset).norm(), 0.0);

  // the material step limit is evaluated on the same local slice
  ES::VXd du = ES::VXd::Zero(n);
  du.segment<3>(9) = ES::V3d(0.0, 0.0, -2.0);  // pushes the apex through the base plane
  ES::VXd duAll = ES::VXd::Zero(offset + n);
  duAll.tail(n) = du;
  const double alpha0 = reference->computeMaxStepSize(u, du);
  EXPECT_LT(alpha0, 1.0);
  EXPECT_DOUBLE_EQ(aggregate.computeMaxStepSize(xAll, duAll), alpha0);
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
