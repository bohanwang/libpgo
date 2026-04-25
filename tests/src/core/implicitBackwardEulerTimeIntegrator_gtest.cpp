#include <gtest/gtest.h>

#include "implicitBackwardEulerTimeIntegrator.h"
#include "NewtonSolver.h"
#include "pgoLogging.h"

#include <numeric>
#include <stdexcept>

namespace
{
namespace ES = pgo::EigenSupport;
using pgo::NonlinearOptimization::NewtonSolver;
using pgo::NonlinearOptimization::PotentialEnergy;
using pgo::Simulation::ImplicitBackwardEulerTimeIntegrator;

void initializeLogging()
{
  static const bool initialized = []() {
    pgo::Logging::init();
    return true;
  }();
  (void)initialized;
}

class TestQuadraticEnergy : public PotentialEnergy
{
public:
  explicit TestQuadraticEnergy(int n, double maxStep = 1.0): n(n), maxStep(maxStep) {}

  double func(ES::ConstRefVecXd x) const override
  {
    return 0.5 * x.squaredNorm();
  }

  void gradient(ES::ConstRefVecXd x, ES::RefVecXd grad) const override
  {
    grad = x;
  }

  void hessian(ES::ConstRefVecXd, ES::SpMatD &hess) const override
  {
    hess.setIdentity();
  }

  void createHessian(ES::SpMatD &hess) const override
  {
    hess.resize(n, n);
    hess.setIdentity();
  }

  void getDOFs(std::vector<int> &dofs) const override
  {
    dofs.resize(n);
    std::iota(dofs.begin(), dofs.end(), 0);
  }

  int getNumDOFs() const override { return n; }
  double computeMaxStepSize(ES::ConstRefVecXd, ES::ConstRefVecXd) const override { return maxStep; }

private:
  int n;
  double maxStep;
};

ES::SpMatD identityMass(int n)
{
  ES::SpMatD mass(n, n);
  mass.setIdentity();
  return mass;
}
}  // namespace

TEST(ImplicitBackwardEulerTimeIntegratorGTest, TryTimestepRejectsFailedNewtonSolveWithoutAdvancingState)
{
  initializeLogging();

  auto energy = std::make_shared<TestQuadraticEnergy>(1);
  ImplicitBackwardEulerTimeIntegrator integrator(identityMass(1), energy, 0.0, 0.0, 0.1, 0, 1e-8);
  integrator.setTimestepID(7);

  const double force[1] = { 1.0 };
  integrator.setExternalForce(force);

  ES::VXd qBefore(1), qvelBefore(1), qaccBefore(1);
  integrator.getq(qBefore);
  integrator.getqvel(qvelBefore);
  integrator.getqacc(qaccBefore);

  testing::internal::CaptureStdout();
  const int ret = integrator.tryTimestep(1, 1, 0);
  const std::string output = testing::internal::GetCapturedStdout();

  EXPECT_EQ(ret, static_cast<int>(NewtonSolver::SolveStatus::MaxIterations));
  EXPECT_EQ(integrator.getTimestepID(), 7u);

  ES::VXd qAfter(1), qvelAfter(1), qaccAfter(1);
  integrator.getq(qAfter);
  integrator.getqvel(qvelAfter);
  integrator.getqacc(qaccAfter);

  EXPECT_TRUE(qAfter.isApprox(qBefore));
  EXPECT_TRUE(qvelAfter.isApprox(qvelBefore));
  EXPECT_TRUE(qaccAfter.isApprox(qaccBefore));
  EXPECT_NE(output.find("ImplicitBackwardEuler timestep begin: T7"), std::string::npos);
  EXPECT_NE(output.find("solverRet=MaxIterations"), std::string::npos);
  EXPECT_NE(output.find("accepted=false"), std::string::npos);
}

TEST(ImplicitBackwardEulerTimeIntegratorGTest, DoTimestepThrowsOnFailedNewtonSolve)
{
  initializeLogging();

  auto energy = std::make_shared<TestQuadraticEnergy>(1);
  ImplicitBackwardEulerTimeIntegrator integrator(identityMass(1), energy, 0.0, 0.0, 0.1, 0, 1e-8);

  const double force[1] = { 1.0 };
  integrator.setExternalForce(force);

  EXPECT_THROW(integrator.doTimestep(1, 0, 0), std::runtime_error);
}

TEST(ImplicitBackwardEulerTimeIntegratorGTest, TryTimestepAcceptsNearConvergedStepTooSmallWithWarning)
{
  initializeLogging();

  auto energy = std::make_shared<TestQuadraticEnergy>(1, 0.0);
  ImplicitBackwardEulerTimeIntegrator integrator(identityMass(1), energy, 0.0, 0.0, 0.1, 8, 1e-8);
  integrator.setTimestepID(11);

  const double force[1] = { 1.5e-8 };
  integrator.setExternalForce(force);

  testing::internal::CaptureStdout();
  const int ret = integrator.tryTimestep(1, 1, 0);
  const std::string output = testing::internal::GetCapturedStdout();

  EXPECT_EQ(ret, 0);
  EXPECT_EQ(integrator.getSolverReturn(), static_cast<int>(NewtonSolver::SolveStatus::StepTooSmall));
  EXPECT_EQ(integrator.getTimestepID(), 12u);
  EXPECT_NE(output.find("solverRet=StepTooSmall"), std::string::npos);
  EXPECT_NE(output.find("accepted=true"), std::string::npos);
  EXPECT_NE(output.find("near-converged stalled solve accepted"), std::string::npos);
}
