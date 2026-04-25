#include <gtest/gtest.h>

#include "NewtonSolver.h"
#include "pgoLogging.h"

#include <numeric>

namespace
{
namespace ES = pgo::EigenSupport;
using pgo::NonlinearOptimization::NewtonSolver;
using pgo::NonlinearOptimization::PotentialEnergy;

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
}  // namespace

TEST(NewtonSolverGTest, ConvergedSolveReturnsConvergedStatus)
{
  initializeLogging();

  auto energy = std::make_shared<TestQuadraticEnergy>(2);
  ES::VXd x(2);
  x[0] = 2.0;
  x[1] = 0.0;

  NewtonSolver::SolverParam solverParam;
  const std::vector<int> fixedDOFs = { 1 };
  const double fixedValues[1] = { 0.0 };
  NewtonSolver solver(x.data(), solverParam, energy, fixedDOFs, fixedValues);

  const int ret = solver.solve(x.data(), 8, 1e-10, 0);

  EXPECT_EQ(ret, static_cast<int>(NewtonSolver::SolveStatus::Converged));
  EXPECT_NEAR(x[0], 0.0, 1e-10);
  EXPECT_NEAR(x[1], 0.0, 1e-10);
}

TEST(NewtonSolverGTest, ZeroFeasibleStepWithLargeResidualReturnsStepTooSmall)
{
  initializeLogging();

  auto energy = std::make_shared<TestQuadraticEnergy>(2, 0.0);
  ES::VXd x(2);
  x[0] = 2.0;
  x[1] = 0.0;

  NewtonSolver::SolverParam solverParam;
  const std::vector<int> fixedDOFs = { 1 };
  const double fixedValues[1] = { 0.0 };
  NewtonSolver solver(x.data(), solverParam, energy, fixedDOFs, fixedValues);

  testing::internal::CaptureStdout();
  const int ret = solver.solve(x.data(), 8, 1e-10, 1);
  const std::string output = testing::internal::GetCapturedStdout();

  EXPECT_EQ(ret, static_cast<int>(NewtonSolver::SolveStatus::StepTooSmall));
  EXPECT_NE(output.find("status=StepTooSmall"), std::string::npos);
  EXPECT_EQ(output.find("T2330"), std::string::npos);
}

TEST(NewtonSolverGTest, SolveStatusToStringReturnsStableNames)
{
  EXPECT_STREQ(NewtonSolver::solveStatusToString(static_cast<int>(NewtonSolver::SolveStatus::Converged)), "Converged");
  EXPECT_STREQ(NewtonSolver::solveStatusToString(static_cast<int>(NewtonSolver::SolveStatus::MaxIterations)), "MaxIterations");
  EXPECT_STREQ(NewtonSolver::solveStatusToString(static_cast<int>(NewtonSolver::SolveStatus::LineSearchFailed)), "LineSearchFailed");
  EXPECT_STREQ(NewtonSolver::solveStatusToString(static_cast<int>(NewtonSolver::SolveStatus::StepTooSmall)), "StepTooSmall");
  EXPECT_STREQ(NewtonSolver::solveStatusToString(static_cast<int>(NewtonSolver::SolveStatus::NonFinite)), "NonFinite");
  EXPECT_STREQ(NewtonSolver::solveStatusToString(static_cast<int>(NewtonSolver::SolveStatus::LinearSolveFailed)), "LinearSolveFailed");
  EXPECT_STREQ(NewtonSolver::solveStatusToString(999), "Unknown");
}
