#include <gtest/gtest.h>

#include "NewtonSolver.h"
#include "deformationModelAssembler.h"
#include "deformationModel.h"
#include "deformationModelEnergy.h"
#include "deformationModelManager.h"
#include "embeddedSurfaceIPCPotentialEnergy.h"
#include "plasticModel3DDeformationGradient.h"
#include "pgoLogging.h"
#include "potentialEnergies.h"
#include "simulationMesh.h"

#include <cmath>
#include <memory>
#include <stdexcept>
#include <vector>

namespace
{
namespace ES = pgo::EigenSupport;
namespace NO = pgo::NonlinearOptimization;
namespace SDM = pgo::SolidDeformationModel;

class FixedQuarticEnergy final : public NO::PotentialEnergy
{
public:
  double func(ES::ConstRefVecXd x) const override
  {
    return 0.25 * x[0] * x[0] * x[0] * x[0];
  }

  void gradient(ES::ConstRefVecXd x, ES::RefVecXd gradient) const override
  {
    gradient[0] = x[0] * x[0] * x[0];
  }

  void hessian(ES::ConstRefVecXd x, ES::SpMatD &hessian) const override
  {
    hessian.coeffRef(0, 0) = 3.0 * x[0] * x[0];
  }

  void createHessian(ES::SpMatD &hessian) const override
  {
    hessian.resize(1, 1);
    hessian.insert(0, 0) = 0.0;
    hessian.makeCompressed();
  }

  void getDOFs(std::vector<int> &dofs) const override { dofs = { 0 }; }
  int getNumDOFs() const override { return 1; }
  double computeMaxStepSize(ES::ConstRefVecXd, ES::ConstRefVecXd) const override { return 1.0; }
};

class InspectableNewtonSolver final : public NO::NewtonSolver
{
public:
  using NewtonSolver::NewtonSolver;

  const void *solverAddress() const { return solver.get(); }
};

class DynamicQuarticEnergy final : public NO::PotentialEnergy
{
public:
  double func(ES::ConstRefVecXd x) const override
  {
    return 0.25 * (x[0] * x[0] * x[0] * x[0] + x[1] * x[1] * x[1] * x[1]);
  }

  void gradient(ES::ConstRefVecXd x, ES::RefVecXd gradient) const override
  {
    gradient[0] = x[0] * x[0] * x[0];
    gradient[1] = x[1] * x[1] * x[1];
  }

  void hessian(ES::ConstRefVecXd, ES::SpMatD &) const override
  {
    throw std::runtime_error("DynamicQuarticEnergy::hessian() must not be called.");
  }

  void createHessian(ES::SpMatD &) const override
  {
    throw std::runtime_error("DynamicQuarticEnergy::createHessian() must not be called.");
  }

  void hessianDirect(ES::ConstRefVecXd x, ES::SpMatD &hessian) const override
  {
    ++hessianDirectCalls_;
    hessian.resize(2, 2);
    hessian.insert(0, 0) = 3.0 * x[0] * x[0];
    hessian.insert(1, 1) = 3.0 * x[1] * x[1];
    if (x[0] > 1.0) {
      hessian.insert(0, 1) = 0.0;
      hessian.insert(1, 0) = 0.0;
    }
    hessian.makeCompressed();
    hessianNonzeroCounts_.push_back(hessian.nonZeros());
  }

  void getDOFs(std::vector<int> &dofs) const override { dofs = { 0, 1 }; }
  int getNumDOFs() const override { return 2; }
  int isHessianTopologyFixed() const override { return 0; }
  double computeMaxStepSize(ES::ConstRefVecXd, ES::ConstRefVecXd) const override { return 1.0; }

  int hessianDirectCalls() const { return hessianDirectCalls_; }
  const std::vector<Eigen::Index> &hessianNonzeroCounts() const { return hessianNonzeroCounts_; }

private:
  mutable int hessianDirectCalls_ = 0;
  mutable std::vector<Eigen::Index> hessianNonzeroCounts_;
};
}  // namespace

TEST(NewtonSolverGTest, FixedTopologyKeepsSymbolicFactorizationAcrossIterations)
{
  pgo::Logging::init();
  const auto energy = std::make_shared<FixedQuarticEnergy>();
  NO::NewtonSolver::SolverParam params;
  params.sst = NO::NewtonSolver::SST_SUBITERATION_ONE;
  params.stopAfterIncrease = 0;

  double x = 3.0;
  InspectableNewtonSolver solver(&x, params, energy, {});
  const void *initialSolver = solver.solverAddress();
  ASSERT_NE(initialSolver, nullptr);

  std::vector<const void *> solverAddresses;
  solver.setStepFunc([&](const ES::VXd &, int) {
    solverAddresses.push_back(solver.solverAddress());
  });
  ASSERT_EQ(solver.solve(&x, 3, 0.0, 0), 0);

  ASSERT_EQ(solverAddresses.size(), 3u);
  for (const void *address : solverAddresses)
    EXPECT_EQ(address, initialSolver);
  EXPECT_NEAR(x, 3.0 * std::pow(2.0 / 3.0, 3.0), 1e-12);
}

TEST(NewtonSolverGTest, DynamicTopologyRebuildsChangingHessianPatterns)
{
  pgo::Logging::init();
  const auto energy = std::make_shared<DynamicQuarticEnergy>();
  NO::NewtonSolver::SolverParam params;
  params.sst = NO::NewtonSolver::SST_SUBITERATION_ONE;
  params.stopAfterIncrease = 0;

  ES::VXd x(2);
  x << 3.0, 2.0;
  NO::NewtonSolver solver(x.data(), params, energy, {});
  ASSERT_EQ(solver.solve(x.data(), 4, 0.0, 0), 0);

  ASSERT_EQ(energy->hessianDirectCalls(), 4);
  ASSERT_EQ(energy->hessianNonzeroCounts().size(), 4u);
  EXPECT_EQ(energy->hessianNonzeroCounts().front(), 4);
  EXPECT_EQ(energy->hessianNonzeroCounts().back(), 2);
  EXPECT_TRUE(x.allFinite());
}

TEST(NewtonSolverGTest, TinyTetWithIpcStaticSolveProducesFiniteResult)
{
  pgo::Logging::init();

  const double vertices[] = {
    0.0, 0.0, 0.0,
    1.0, 0.0, 0.0,
    0.0, 1.0, 0.0,
    0.0, 0.0, 1.0,
  };
  const int elementVertices[] = { 0, 1, 2, 3 };
  const int elementMaterialIndices[] = { 0 };
  SDM::SimulationMeshENuMaterial material(1200.0, 0.45);
  const SDM::SimulationMeshMaterial *materials[] = { &material };
  auto mesh = std::make_shared<SDM::SimulationMesh>(
    4, vertices, 1, 4, elementVertices, elementMaterialIndices, 1, materials, SDM::SimulationMeshType::TET);

  auto manager = std::make_shared<SDM::DeformationModelManager>();
  manager->setMesh(mesh.get());
  manager->init(SDM::DeformationModelPlasticMaterial::VOLUMETRIC_DOF6, SDM::DeformationModelElasticMaterial::STABLE_NEO);
  auto assembler = std::make_shared<SDM::DeformationModelAssembler>(manager, nullptr);

  ES::VXd restPositions(12);
  for (int vi = 0; vi < 4; ++vi)
    restPositions.segment<3>(3 * vi) = ES::Mp<const ES::V3d>(vertices + 3 * vi);
  auto deformationEnergy = std::make_shared<SDM::DeformationModelEnergy>(assembler, &restPositions, 0);
  const auto *plasticModel = dynamic_cast<const SDM::PlasticModel3DDeformationGradient *>(
    manager->getDeformationModel(0)->getPlasticModel());
  ASSERT_NE(plasticModel, nullptr);
  ES::VXd plasticParams(manager->getNumPlasticParameters() * mesh->getNumElements());
  const ES::M3d identity = ES::M3d::Identity();
  plasticModel->toParam(identity.data(), plasticParams.data());
  deformationEnergy->setPlasticParams(plasticParams);
  deformationEnergy->setElasticParams(ES::VXd::Zero(manager->getNumElasticParameters() * mesh->getNumElements()));

  ES::MXd surfaceVertices(4, 3);
  for (int vi = 0; vi < 4; ++vi)
    surfaceVertices.row(vi) = ES::Mp<const ES::V3d>(vertices + 3 * vi).transpose();
  ES::MXi surfaceTriangles(4, 3);
  surfaceTriangles <<
    0, 2, 1,
    0, 1, 3,
    0, 3, 2,
    1, 2, 3;
  ES::SpMatD identityEmbedding(12, 12);
  identityEmbedding.setIdentity();
  pgo::Contact::CIPC::SurfaceIPCCore::Parameters ipcParams;
  ipcParams.dhat = 0.1;
  ipcParams.kappa = 10.0;
  auto ipcEnergy = std::make_shared<pgo::Contact::CIPC::EmbeddedSurfaceIPCPotentialEnergy>(
    surfaceVertices, surfaceTriangles, identityEmbedding, ipcParams);

  auto aggregate = std::make_shared<NO::PotentialEnergies>(12);
  aggregate->addPotentialEnergy(deformationEnergy);
  aggregate->addPotentialEnergy(ipcEnergy);
  ASSERT_EQ(aggregate->isHessianTopologyFixed(), 0);
  ASSERT_NO_THROW(aggregate->init());

  NO::NewtonSolver::SolverParam params;
  params.sst = NO::NewtonSolver::SST_SUBITERATION_ONE;
  params.stopAfterIncrease = 0;
  ES::VXd displacements = ES::VXd::Zero(12);
  const std::vector<int> fixedDOFs = { 0, 1, 2, 3, 4, 5, 6, 7, 8 };
  const std::vector<double> fixedValues(fixedDOFs.size(), 0.0);
  NO::NewtonSolver solver(displacements.data(), params, aggregate, fixedDOFs, fixedValues.data());

  ASSERT_EQ(solver.solve(displacements.data(), 1, 0.0, 0), 0);
  EXPECT_TRUE(displacements.allFinite());
  EXPECT_TRUE(std::isfinite(aggregate->func(displacements)));
  ES::VXd gradient = ES::VXd::Zero(12);
  aggregate->gradient(displacements, gradient);
  EXPECT_TRUE(gradient.allFinite());
}
