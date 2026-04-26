#include <gtest/gtest.h>

#include "NewtonSolver.h"
#include "pgoLogging.h"
#include "solveDiagnostics.h"

#include <numeric>

namespace
{
namespace ES = pgo::EigenSupport;
using pgo::NonlinearOptimization::NewtonSolver;
using pgo::NonlinearOptimization::MaxStepResult;
using pgo::NonlinearOptimization::PotentialEnergy;
using pgo::NonlinearOptimization::SolveDiagnostics;

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
  explicit TestQuadraticEnergy(int n, MaxStepResult maxStep = MaxStepResult::unconstrained()): n(n), maxStep(maxStep) {}

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
  MaxStepResult computeMaxStepLimit(ES::ConstRefVecXd, ES::ConstRefVecXd) const override { return maxStep; }

private:
  int n;
  MaxStepResult maxStep;
};
}  // namespace

TEST(SolveDiagnosticsGTest, RecordsAndResetsMaxStepAndLineSearch)
{
  SolveDiagnostics diagnostics;

  diagnostics.recordMaxStep(MaxStepResult::material(0.4));
  diagnostics.recordMaxStep(MaxStepResult::contact(0.25));
  diagnostics.recordLineSearch(0.25, 0.5, 0.125);

  EXPECT_EQ(diagnostics.materialClampCount, 1);
  EXPECT_EQ(diagnostics.contactClampCount, 1);
  EXPECT_DOUBLE_EQ(diagnostics.minFeasibleAlpha, 0.25);
  EXPECT_DOUBLE_EQ(diagnostics.minMaterialFeasibleAlpha, 0.4);
  EXPECT_DOUBLE_EQ(diagnostics.minContactFeasibleAlpha, 0.25);
  EXPECT_DOUBLE_EQ(diagnostics.minLineSearchAlpha, 0.5);
  EXPECT_DOUBLE_EQ(diagnostics.minEffectiveAlpha, 0.125);
  EXPECT_DOUBLE_EQ(diagnostics.currentMaterialAlpha, 1.0);
  EXPECT_DOUBLE_EQ(diagnostics.currentContactAlpha, 0.25);

  diagnostics.reset();

  EXPECT_EQ(diagnostics.materialClampCount, 0);
  EXPECT_EQ(diagnostics.contactClampCount, 0);
  EXPECT_DOUBLE_EQ(diagnostics.minFeasibleAlpha, 1.0);
  EXPECT_DOUBLE_EQ(diagnostics.minMaterialFeasibleAlpha, 1.0);
  EXPECT_DOUBLE_EQ(diagnostics.minContactFeasibleAlpha, 1.0);
  EXPECT_DOUBLE_EQ(diagnostics.minLineSearchAlpha, 1.0);
  EXPECT_DOUBLE_EQ(diagnostics.minEffectiveAlpha, 1.0);
  EXPECT_DOUBLE_EQ(diagnostics.currentMaterialAlpha, 1.0);
  EXPECT_DOUBLE_EQ(diagnostics.currentContactAlpha, 1.0);
}

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

  auto energy = std::make_shared<TestQuadraticEnergy>(2, MaxStepResult::contact(0.0));
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

TEST(NewtonSolverGTest, SolveDiagnosticsRecordsMaxStepBreakdown)
{
  initializeLogging();

  auto energy = std::make_shared<TestQuadraticEnergy>(2, MaxStepResult::material(0.25));
  ES::VXd x(2);
  x[0] = 2.0;
  x[1] = 0.0;

  NewtonSolver::SolverParam solverParam;
  const std::vector<int> fixedDOFs = { 1 };
  const double fixedValues[1] = { 0.0 };
  NewtonSolver solver(x.data(), solverParam, energy, fixedDOFs, fixedValues);

  testing::internal::CaptureStdout();
  const int ret = solver.solve(x.data(), 1, 1e-10, 2);
  const std::string output = testing::internal::GetCapturedStdout();

  const SolveDiagnostics &diagnostics = solver.getSolveDiagnostics();
  EXPECT_EQ(ret, static_cast<int>(NewtonSolver::SolveStatus::MaxIterations));
  EXPECT_EQ(diagnostics.materialClampCount, 1);
  EXPECT_EQ(diagnostics.contactClampCount, 0);
  EXPECT_DOUBLE_EQ(diagnostics.minFeasibleAlpha, 0.25);
  EXPECT_DOUBLE_EQ(diagnostics.minMaterialFeasibleAlpha, 0.25);
  EXPECT_DOUBLE_EQ(diagnostics.minContactFeasibleAlpha, 1.0);
  EXPECT_NE(output.find("feasible alpha clamped: material:0.25 contact:1"), std::string::npos);
}

TEST(NewtonSolverGTest, VerboseIterationLogLabelsMaxGradientNorm)
{
  initializeLogging();

  auto energy = std::make_shared<TestQuadraticEnergy>(2);
  ES::VXd x(2);
  x[0] = 3.0;
  x[1] = 4.0;

  NewtonSolver::SolverParam solverParam;
  const std::vector<int> fixedDOFs;
  NewtonSolver solver(x.data(), solverParam, energy, fixedDOFs);

  testing::internal::CaptureStdout();
  solver.solve(x.data(), 1, 0.0, 2);
  const std::string output = testing::internal::GetCapturedStdout();

  EXPECT_NE(output.find("||grad||_max=4"), std::string::npos);
  EXPECT_EQ(output.find("; ||grad||="), std::string::npos);
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
