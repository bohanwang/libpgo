#include <gtest/gtest.h>

#include "implicitBackwardEulerTimeIntegrator.h"
#include "pgoLogging.h"
#include "solverResult.h"
#include "TRBDF2TimeIntegrator.h"

#include <numeric>
#include <stdexcept>

namespace
{
namespace ES = pgo::EigenSupport;
using pgo::NonlinearOptimization::MaxStepResult;
using pgo::NonlinearOptimization::PotentialEnergy;
using pgo::NonlinearOptimization::SolverResult;
using pgo::NonlinearOptimization::SolveStatus;
using pgo::Simulation::ImplicitBackwardEulerTimeIntegrator;
using pgo::Simulation::TRBDF2TimeIntegrator;

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
    funcCalls++;
    return 0.5 * x.squaredNorm();
  }

  void gradient(ES::ConstRefVecXd x, ES::RefVecXd grad) const override
  {
    gradientCalls++;
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
  MaxStepResult computeMaxStepLimit(ES::ConstRefVecXd, ES::ConstRefVecXd) const override { return MaxStepResult::material(maxStep); }

  mutable int funcCalls = 0;
  mutable int gradientCalls = 0;

private:
  int n;
  double maxStep;
};

class TestNonFixedQuadraticEnergy : public PotentialEnergy
{
public:
  explicit TestNonFixedQuadraticEnergy(int n): n(n) {}

  double func(ES::ConstRefVecXd x) const override
  {
    return 0.5 * x.squaredNorm();
  }

  void gradient(ES::ConstRefVecXd x, ES::RefVecXd grad) const override
  {
    gradientCalls++;
    grad = x;
  }

  void hessian(ES::ConstRefVecXd, ES::SpMatD &hess) const override
  {
    hessianCalls++;
    hess.resize(n, n);
    hess.setIdentity();
  }

  void gradient_hessian(ES::ConstRefVecXd x, ES::RefVecXd grad, ES::SpMatD &hess) const override
  {
    gradientHessianCalls++;
    grad = x;
    hess.resize(n, n);
    hess.setIdentity();
  }

  double func_grad_hessian(ES::ConstRefVecXd x, ES::RefVecXd grad, ES::SpMatD &hess) const override
  {
    gradient_hessian(x, grad, hess);
    return func(x);
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
  int isHessianTopologyFixed() const override { return 0; }
  MaxStepResult computeMaxStepLimit(ES::ConstRefVecXd, ES::ConstRefVecXd) const override { return MaxStepResult::unconstrained(); }

  void resetCounts() const
  {
    gradientCalls = 0;
    hessianCalls = 0;
    gradientHessianCalls = 0;
  }

  mutable int gradientCalls = 0;
  mutable int hessianCalls = 0;
  mutable int gradientHessianCalls = 0;

private:
  int n;
};

ES::SpMatD identityMass(int n)
{
  ES::SpMatD mass(n, n);
  mass.setIdentity();
  return mass;
}
}  // namespace

TEST(ImplicitBackwardEulerTimeIntegratorGTest, TryTimestepAcceptsMaxIterationsAndAdvancesTimestep)
{
  initializeLogging();

  auto energy = std::make_shared<TestQuadraticEnergy>(1);
  ImplicitBackwardEulerTimeIntegrator integrator(identityMass(1), energy, 0.0, 0.0, 0.1, 0, 1e-8);
  integrator.setTimestepID(7);

  const double force[1] = { 1.0 };
  integrator.setExternalForce(force);

  testing::internal::CaptureStdout();
  const SolverResult result = integrator.tryTimestep(1, 1, 0);
  const std::string output = testing::internal::GetCapturedStdout();

  EXPECT_EQ(result.status, SolveStatus::MaxIterations);
  EXPECT_EQ(integrator.getSolverStatus(), SolveStatus::MaxIterations);
  EXPECT_EQ(integrator.getLastSolveResult().rawStatusCode, static_cast<int>(SolveStatus::MaxIterations));
  EXPECT_EQ(integrator.getTimestepID(), 8u);

  EXPECT_NE(output.find("ImplicitBackwardEuler timestep begin: T7"), std::string::npos);
  EXPECT_NE(output.find("status=MaxIterations"), std::string::npos);
  EXPECT_NE(output.find("accepted=true"), std::string::npos);
}

TEST(ImplicitBackwardEulerTimeIntegratorGTest, DoTimestepDoesNotThrowOnAcceptedMaxIterations)
{
  initializeLogging();

  auto energy = std::make_shared<TestQuadraticEnergy>(1);
  ImplicitBackwardEulerTimeIntegrator integrator(identityMass(1), energy, 0.0, 0.0, 0.1, 0, 1e-8);

  const double force[1] = { 1.0 };
  integrator.setExternalForce(force);

  EXPECT_NO_THROW(integrator.doTimestep(1, 0, 0));
  EXPECT_EQ(integrator.getSolverStatus(), SolveStatus::MaxIterations);
  EXPECT_EQ(integrator.getTimestepID(), 1u);
}

TEST(ImplicitBackwardEulerTimeIntegratorGTest, ResidualPrintReusesNewtonFinalGradient)
{
  initializeLogging();

  auto energy = std::make_shared<TestQuadraticEnergy>(1);
  ImplicitBackwardEulerTimeIntegrator integrator(identityMass(1), energy, 0.0, 0.0, 0.1, 10, 1e-8);

  testing::internal::CaptureStdout();
  const SolverResult result = integrator.tryTimestep(0, 1, 1);
  const std::string output = testing::internal::GetCapturedStdout();

  EXPECT_EQ(result.status, SolveStatus::Converged);
  EXPECT_EQ(integrator.getSolverStatus(), SolveStatus::Converged);
  EXPECT_NE(output.find("residual=0"), std::string::npos);
  EXPECT_NE(output.find("sub 0: 0"), std::string::npos);
  EXPECT_EQ(energy->funcCalls, 1);
  EXPECT_EQ(energy->gradientCalls, 1);
}

TEST(ImplicitBackwardEulerTimeIntegratorGTest, WrapperGradientHessianDispatchesToNonFixedModel)
{
  initializeLogging();

  auto elastic = std::make_shared<TestQuadraticEnergy>(1);
  auto nonFixed = std::make_shared<TestNonFixedQuadraticEnergy>(1);
  ImplicitBackwardEulerTimeIntegrator integrator(identityMass(1), elastic, 0.0, 0.0, 0.1, 0, 1e-8);
  integrator.addGeneralImplicitForceModel(nonFixed, 0.0, 0.0);

  testing::internal::CaptureStdout();
  (void)integrator.tryTimestep(0, 0, 0);
  (void)testing::internal::GetCapturedStdout();
  nonFixed->resetCounts();

  const auto energy = integrator.getInternalEnergy();
  ES::VXd x(1);
  x[0] = 0.25;
  ES::VXd combinedGrad = ES::VXd::Zero(1);
  ES::SpMatD combinedH;
  energy->gradient_hessian(x, combinedGrad, combinedH);

  EXPECT_EQ(nonFixed->gradientHessianCalls, 1);
  EXPECT_EQ(nonFixed->gradientCalls, 0);
  EXPECT_EQ(nonFixed->hessianCalls, 0);

  nonFixed->resetCounts();
  ES::VXd refGrad = ES::VXd::Zero(1);
  energy->gradient(x, refGrad);
  ES::SpMatD refH;
  energy->hessianDirect(x, refH);

  EXPECT_LT((combinedGrad - refGrad).norm(), 1e-12);
  EXPECT_LT((ES::MXd(combinedH) - ES::MXd(refH)).norm(), 1e-12);
}

TEST(TRBDF2TimeIntegratorGTest, NewtonSolveDispatchesWrapperGradientHessianToNonFixedModel)
{
  initializeLogging();

  auto elastic = std::make_shared<TestQuadraticEnergy>(1);
  auto nonFixed = std::make_shared<TestNonFixedQuadraticEnergy>(1);
  TRBDF2TimeIntegrator integrator(identityMass(1), elastic, 0.0, 0.0, 0.5, 0.1, 1, 1e-8);
  integrator.addGeneralImplicitForceModel(nonFixed, 0.0, 0.0);

  const double force[1] = { 1.0 };
  integrator.setExternalForce(force);

  testing::internal::CaptureStdout();
  integrator.doTimestep(0, 0, 0);
  (void)testing::internal::GetCapturedStdout();

  EXPECT_GT(nonFixed->gradientHessianCalls, 0);
}
