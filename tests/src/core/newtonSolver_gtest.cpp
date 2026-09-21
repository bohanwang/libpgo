#include <gtest/gtest.h>

#include "NewtonSolver.h"
#include "deformationModelAssembler.h"
#include "deformationModel.h"
#include "deformationModelEnergy.h"
#include "deformationModelManager.h"
#include "embeddedSurfaceIPCPotentialEnergy.h"
#include "linearPotentialEnergy.h"
#include "plasticModel3DDeformationGradient.h"
#include "pgoLogging.h"
#include "potentialEnergies.h"
#include "simulationMesh.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
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
};

// -0.5 x^2: the Newton direction always increases the energy.
class ConcaveEnergy final : public NO::PotentialEnergy
{
public:
  double func(ES::ConstRefVecXd x) const override { return -0.5 * x[0] * x[0]; }
  void gradient(ES::ConstRefVecXd x, ES::RefVecXd gradient) const override { gradient[0] = -x[0]; }
  void hessian(ES::ConstRefVecXd, ES::SpMatD &hessian) const override { hessian.coeffRef(0, 0) = -1.0; }
  void createHessian(ES::SpMatD &hessian) const override
  {
    hessian.resize(1, 1);
    hessian.insert(0, 0) = 0.0;
    hessian.makeCompressed();
  }
  void getDOFs(std::vector<int> &dofs) const override { dofs = { 0 }; }
  int getNumDOFs() const override { return 1; }
};

// 0.5 (x + 1)^2 for x >= 3 and +inf below: Newton from x = 3 always steps into
// the wall, so every trial step has infinite energy and the line search fails.
class WallEnergy final : public NO::PotentialEnergy
{
public:
  double func(ES::ConstRefVecXd x) const override
  {
    return x[0] >= 3.0 ? 0.5 * (x[0] + 1.0) * (x[0] + 1.0) : std::numeric_limits<double>::infinity();
  }
  void gradient(ES::ConstRefVecXd x, ES::RefVecXd gradient) const override { gradient[0] = x[0] + 1.0; }
  void hessian(ES::ConstRefVecXd, ES::SpMatD &hessian) const override { hessian.coeffRef(0, 0) = 1.0; }
  void createHessian(ES::SpMatD &hessian) const override
  {
    hessian.resize(1, 1);
    hessian.insert(0, 0) = 0.0;
    hessian.makeCompressed();
  }
  void getDOFs(std::vector<int> &dofs) const override { dofs = { 0 }; }
  int getNumDOFs() const override { return 1; }
};

class NonFiniteGradientEnergy final : public NO::PotentialEnergy
{
public:
  double func(ES::ConstRefVecXd x) const override { return 0.5 * x[0] * x[0]; }
  void gradient(ES::ConstRefVecXd, ES::RefVecXd gradient) const override
  {
    gradient[0] = std::numeric_limits<double>::quiet_NaN();
  }
  void hessian(ES::ConstRefVecXd, ES::SpMatD &hessian) const override { hessian.coeffRef(0, 0) = 1.0; }
  void createHessian(ES::SpMatD &hessian) const override
  {
    hessian.resize(1, 1);
    hessian.insert(0, 0) = 0.0;
    hessian.makeCompressed();
  }
  void getDOFs(std::vector<int> &dofs) const override { dofs = { 0 }; }
  int getNumDOFs() const override { return 1; }
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

class RecordingPotentialEnergy final : public NO::PotentialEnergy
{
public:
  explicit RecordingPotentialEnergy(NO::PotentialEnergy_p energy): energy_(std::move(energy)) {}

  double func(ES::ConstRefVecXd x) const override { return energy_->func(x); }
  void gradient(ES::ConstRefVecXd x, ES::RefVecXd gradient) const override { energy_->gradient(x, gradient); }
  void hessian(ES::ConstRefVecXd x, ES::SpMatD &hessian) const override { energy_->hessian(x, hessian); }
  void hessianDirect(ES::ConstRefVecXd x, ES::SpMatD &hessian) const override { energy_->hessianDirect(x, hessian); }
  void createHessian(ES::SpMatD &hessian) const override { energy_->createHessian(hessian); }
  void getDOFs(std::vector<int> &dofs) const override { energy_->getDOFs(dofs); }
  int getNumDOFs() const override { return energy_->getNumDOFs(); }
  int isQuadratic() const override { return energy_->isQuadratic(); }
  int hasHessianVector() const override { return energy_->hasHessianVector(); }
  int hasHessian() const override { return energy_->hasHessian(); }
  int isHessianTopologyFixed() const override { return energy_->isHessianTopologyFixed(); }

  double computeMaxStepSize(ES::ConstRefVecXd x, ES::ConstRefVecXd dx) const override
  {
    const double alpha = energy_->computeMaxStepSize(x, dx);
    ++maxStepCalls_;
    minimumReturnedAlpha_ = std::min(minimumReturnedAlpha_, alpha);
    return alpha;
  }

  int maxStepCalls() const { return maxStepCalls_; }
  double minimumReturnedAlpha() const { return minimumReturnedAlpha_; }

private:
  NO::PotentialEnergy_p energy_;
  mutable int maxStepCalls_ = 0;
  mutable double minimumReturnedAlpha_ = 1.0;
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
  // A zero tolerance can never be met, so the solve ends at the iteration limit.
  ASSERT_EQ(solver.solve(&x, 3, 0.0, 0), NO::NewtonSolver::SOLVE_NOT_CONVERGED);

  ASSERT_EQ(solverAddresses.size(), 3u);
  for (const void *address : solverAddresses)
    EXPECT_EQ(address, initialSolver);
  EXPECT_NEAR(x, 3.0 * std::pow(2.0 / 3.0, 3.0), 1e-12);
}

TEST(NewtonSolverGTest, PotentialEnergyDefaultMaxStepSizeIsOne)
{
  FixedQuarticEnergy energy;
  ES::VXd x(1);
  ES::VXd dx(1);
  x << 0.0;
  dx << 1.0;
  EXPECT_DOUBLE_EQ(energy.computeMaxStepSize(x, dx), 1.0);
}

TEST(NewtonSolverGTest, ReportsGradientAtReturnedStateAndDistinguishesIterationLimit)
{
  pgo::Logging::init();
  const auto energy = std::make_shared<FixedQuarticEnergy>();
  NO::NewtonSolver::SolverParam params;
  double x = 3.0;
  NO::NewtonSolver solver(&x, params, energy, {});

  testing::internal::CaptureStdout();
  const int limitedReturn = solver.solve(&x, 1, 1e-4, 1);
  const std::string limitedLog = testing::internal::GetCapturedStdout();
  EXPECT_EQ(limitedReturn, NO::NewtonSolver::SOLVE_NOT_CONVERGED);
  EXPECT_STREQ(solver.getLastStopReason(), "iteration_limit");
  EXPECT_DOUBLE_EQ(solver.getLastGradientNorm(), 8.0);
  EXPECT_NEAR(x, 2.0, 1e-12);
  EXPECT_NE(limitedLog.find("stop=iteration_limit, iteration=1, gradient_inf=8,"), std::string::npos);
  EXPECT_NE(limitedLog.find("converged=false, status=not_converged"), std::string::npos);

  testing::internal::CaptureStdout();
  const int convergedReturn = solver.solve(&x, 50, 1e-4, 1);
  const std::string convergedLog = testing::internal::GetCapturedStdout();
  EXPECT_EQ(convergedReturn, NO::NewtonSolver::SOLVE_CONVERGED);
  EXPECT_STREQ(solver.getLastStopReason(), "gradient_tolerance");
  EXPECT_LT(std::abs(x * x * x), 1e-4);
  EXPECT_NE(convergedLog.find("stop=gradient_tolerance"), std::string::npos);
  EXPECT_NE(convergedLog.find("converged=true, status=converged"), std::string::npos);
}

TEST(NewtonSolverGTest, StopsAfterFailedLineSearchRegardlessOfVerbosity)
{
  pgo::Logging::init();
  // Every trial step hits the wall, so the backtracking line search returns
  // its full, infinitely expensive step and the solve must stop at the first
  // iteration without applying it, for any verbosity level.
  for (int verbose : { 0, 1, 2 }) {
    const auto energy = std::make_shared<WallEnergy>();
    NO::NewtonSolver::SolverParam params;
    params.lsm = NO::NewtonSolver::LSM_BACKTRACK;
    double x = 3.0;
    NO::NewtonSolver solver(&x, params, energy, {});
    int iterationsRun = 0;
    solver.setStepFunc([&](const ES::VXd &, int) { ++iterationsRun; });

    testing::internal::CaptureStdout();
    const int status = solver.solve(&x, 20, 1e-8, verbose);
    testing::internal::GetCapturedStdout();

    EXPECT_EQ(status, NO::NewtonSolver::SOLVE_NOT_CONVERGED) << "verbose=" << verbose;
    EXPECT_STREQ(solver.getLastStopReason(), "line_search_failed") << "verbose=" << verbose;
    EXPECT_EQ(iterationsRun, 0) << "verbose=" << verbose;
    EXPECT_DOUBLE_EQ(x, 3.0) << "verbose=" << verbose;
    EXPECT_DOUBLE_EQ(solver.getLastGradientNorm(), 4.0) << "verbose=" << verbose;
  }

  // On a concave energy the backtracking halves the step until its energy
  // change drops below roundoff, so the solve ends through either the failed
  // line search or the tiny-step exit. Both must stop at the first iteration.
  for (int verbose : { 0, 1, 2 }) {
    const auto energy = std::make_shared<ConcaveEnergy>();
    NO::NewtonSolver::SolverParam params;
    double x = 3.0;
    NO::NewtonSolver solver(&x, params, energy, {});
    int iterationsRun = 0;
    solver.setStepFunc([&](const ES::VXd &, int) { ++iterationsRun; });

    testing::internal::CaptureStdout();
    const int status = solver.solve(&x, 20, 1e-8, verbose);
    testing::internal::GetCapturedStdout();

    const std::string reason = solver.getLastStopReason();
    EXPECT_EQ(status, NO::NewtonSolver::SOLVE_NOT_CONVERGED) << "verbose=" << verbose;
    EXPECT_TRUE(reason == "line_search_failed" || reason == "tiny_step") << "verbose=" << verbose << " reason=" << reason;
    EXPECT_EQ(iterationsRun, 0) << "verbose=" << verbose;
    EXPECT_DOUBLE_EQ(x, 3.0) << "verbose=" << verbose;
  }
}

TEST(NewtonSolverGTest, ReportsNumericalFailureInsteadOfAborting)
{
  pgo::Logging::init();
  const auto energy = std::make_shared<NonFiniteGradientEnergy>();
  NO::NewtonSolver::SolverParam params;
  double x = 3.0;
  NO::NewtonSolver solver(&x, params, energy, {});

  testing::internal::CaptureStdout();
  const int status = solver.solve(&x, 5, 1e-8, 0);
  testing::internal::GetCapturedStdout();

  EXPECT_EQ(status, NO::NewtonSolver::SOLVE_NUMERICAL_FAILURE);
  EXPECT_STREQ(solver.getLastStopReason(), "non_finite_energy");
  EXPECT_STREQ(NO::NewtonSolver::solveStatusName(status), "numerical_failure");
  EXPECT_DOUBLE_EQ(x, 3.0);
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
  ASSERT_EQ(solver.solve(x.data(), 4, 0.0, 0), NO::NewtonSolver::SOLVE_NOT_CONVERGED);

  ASSERT_EQ(energy->hessianDirectCalls(), 4);
  ASSERT_EQ(energy->hessianNonzeroCounts().size(), 4u);
  EXPECT_EQ(energy->hessianNonzeroCounts().front(), 4);
  EXPECT_EQ(energy->hessianNonzeroCounts().back(), 2);
  EXPECT_TRUE(x.allFinite());
}

TEST(NewtonSolverGTest, StaticIpcSelfContactUsesCcdMaxStepAndProducesFiniteResult)
{
  pgo::Logging::init();

  const double vertices[] = {
    -1.0, -1.0, 0.0,
     1.0, -1.0, 0.0,
     0.0,  1.0, 0.0,
     0.0,  0.0, -1.0,
    -1.0, -1.0, 1.0,
     1.0, -1.0, 1.0,
     0.0,  1.0, 1.0,
     0.0,  0.0, 0.2,
  };
  const int elementVertices[] = {
    0, 2, 1, 3,
    4, 6, 5, 7,
  };
  const int elementMaterialIndices[] = { 0, 0 };
  SDM::SimulationMeshENuMaterial material(1200.0, 0.45);
  const SDM::SimulationMeshMaterial *materials[] = { &material };
  auto mesh = std::make_shared<SDM::SimulationMesh>(
    8, vertices, 2, 4, elementVertices, elementMaterialIndices, 1, materials, SDM::SimulationMeshType::TET);

  auto manager = std::make_shared<SDM::DeformationModelManager>();
  manager->setMesh(mesh.get());
  manager->init(SDM::DeformationModelPlasticMaterial::VOLUMETRIC_DOF6, SDM::DeformationModelElasticMaterial::STABLE_NEO);
  manager->setEnforceSPD(1);
  auto assembler = std::make_shared<SDM::DeformationModelAssembler>(manager, nullptr);

  ES::VXd restPositions(24);
  for (int vi = 0; vi < 8; ++vi)
    restPositions.segment<3>(3 * vi) = ES::Mp<const ES::V3d>(vertices + 3 * vi);
  auto deformationEnergy = std::make_shared<SDM::DeformationModelEnergy>(assembler, &restPositions, 0);
  const auto *plasticModel = dynamic_cast<const SDM::PlasticModel3DDeformationGradient *>(
    manager->getDeformationModel(0)->getPlasticModel());
  ASSERT_NE(plasticModel, nullptr);
  ES::VXd plasticParams = ES::VXd::Zero(manager->getNumPlasticParameters() * mesh->getNumElements());
  const ES::M3d identity = ES::M3d::Identity();
  for (int element = 0; element < mesh->getNumElements(); ++element)
    plasticModel->toParam(identity.data(), plasticParams.data() + element * manager->getNumPlasticParameters());
  deformationEnergy->setPlasticParams(plasticParams);

  ES::MXd surfaceVertices(8, 3);
  for (int vi = 0; vi < 8; ++vi)
    surfaceVertices.row(vi) = ES::Mp<const ES::V3d>(vertices + 3 * vi).transpose();
  ES::MXi surfaceTriangles(8, 3);
  surfaceTriangles <<
    0, 2, 1,
    0, 1, 3,
    0, 3, 2,
    1, 2, 3,
    4, 6, 5,
    4, 5, 7,
    4, 7, 6,
    5, 6, 7;
  ES::SpMatD identityEmbedding(24, 24);
  identityEmbedding.setIdentity();
  pgo::Contact::CIPC::SurfaceIPCCore::Parameters ipcParams;
  ipcParams.dhat = 0.3;
  ipcParams.kappa = 1.0e4;
  auto embeddedIpcEnergy = std::make_shared<pgo::Contact::CIPC::EmbeddedSurfaceIPCPotentialEnergy>(
    surfaceVertices, surfaceTriangles, identityEmbedding, ipcParams);
  auto ipcEnergy = std::make_shared<RecordingPotentialEnergy>(embeddedIpcEnergy);

  ES::VXd downwardLoad = ES::VXd::Zero(24);
  downwardLoad[3 * 7 + 2] = -50.0;
  auto loadEnergy = std::make_shared<pgo::PredefinedPotentialEnergies::LinearPotentialEnergy>(downwardLoad);

  auto aggregate = std::make_shared<NO::PotentialEnergies>(24);
  aggregate->addPotentialEnergy(deformationEnergy);
  aggregate->addPotentialEnergy(loadEnergy, -1.0);
  aggregate->addPotentialEnergy(ipcEnergy);
  ASSERT_EQ(aggregate->isHessianTopologyFixed(), 0);
  ASSERT_NO_THROW(aggregate->init());

  ES::VXd displacements = ES::VXd::Zero(24);
  ES::VXd crossingDirection = ES::VXd::Zero(24);
  crossingDirection[3 * 7 + 2] = -1.0;
  const double ccdAlpha = aggregate->computeMaxStepSize(displacements, crossingDirection);
  EXPECT_GT(ccdAlpha, 0.0);
  EXPECT_LT(ccdAlpha, 1.0);

  const double initialEnergy = aggregate->func(displacements);
  EXPECT_GT(initialEnergy, 0.0);
  ES::VXd initialGradient = ES::VXd::Zero(24);
  aggregate->gradient(displacements, initialGradient);

  NO::NewtonSolver::SolverParam params;
  params.sst = NO::NewtonSolver::SST_SUBITERATION_LINE_SEARCH;
  params.lsm = NO::NewtonSolver::LSM_SIMPLE;
  params.stopAfterIncrease = 0;
  const std::vector<int> fixedDOFs = {
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,
    12, 13, 14, 15, 16, 17, 18, 19, 20,
  };
  const std::vector<double> fixedValues(fixedDOFs.size(), 0.0);
  NO::NewtonSolver solver(displacements.data(), params, aggregate, fixedDOFs, fixedValues.data());

  ASSERT_NE(solver.solve(displacements.data(), 12, 1e-8, 0), NO::NewtonSolver::SOLVE_NUMERICAL_FAILURE);
  EXPECT_TRUE(displacements.allFinite());
  const double finalEnergy = aggregate->func(displacements);
  EXPECT_TRUE(std::isfinite(finalEnergy));
  EXPECT_LT(finalEnergy, initialEnergy);
  ES::VXd gradient = ES::VXd::Zero(24);
  aggregate->gradient(displacements, gradient);
  EXPECT_TRUE(gradient.allFinite());
  for (int dof : fixedDOFs)
    EXPECT_DOUBLE_EQ(displacements[dof], 0.0);
  EXPECT_LT(gradient.tail<3>().norm(), initialGradient.tail<3>().norm());
  EXPECT_GT(ipcEnergy->maxStepCalls(), 1);
  EXPECT_GT(ipcEnergy->minimumReturnedAlpha(), 0.0);
  EXPECT_LT(ipcEnergy->minimumReturnedAlpha(), 1.0);
}
