#pragma once

#include "solveDiagnostics.h"

#include <string>

namespace pgo::NonlinearOptimization
{

enum class SolveStatus : int
{
  Converged = 0,
  MaxIterations = 1,
  LineSearchFailed = 2,
  StepTooSmall = 3,
  NonFinite = 4,
  LinearSolveFailed = 5,
  ExternalSolverFailure = 100,
  UnsupportedBackend = 101
};

struct SolverResult
{
  SolveStatus status = SolveStatus::MaxIterations;
  int iterations = 0;
  int rawStatusCode = static_cast<int>(SolveStatus::MaxIterations);
  double finalGradientNorm = 0.0;
  double finalGradientMaxNorm = 0.0;
  bool hasFinalGradientStats = false;
  SolveDiagnostics diagnostics;

  bool converged() const { return status == SolveStatus::Converged; }
};

const char *solveStatusToString(SolveStatus status);
const char *solveStatusToString(int status);
bool isSolveStatusCode(int status);
std::string formatSolverResultSummary(const SolverResult &result);

}  // namespace pgo::NonlinearOptimization

namespace pgo
{
using NonlinearOptimization::SolveStatus;
using NonlinearOptimization::SolverResult;
}  // namespace pgo
