#include "solverResult.h"

#include <sstream>

namespace pgo::NonlinearOptimization
{

const char *solveStatusToString(SolveStatus status)
{
  switch (status) {
    case SolveStatus::Converged:
      return "Converged";
    case SolveStatus::MaxIterations:
      return "MaxIterations";
    case SolveStatus::LineSearchFailed:
      return "LineSearchFailed";
    case SolveStatus::StepTooSmall:
      return "StepTooSmall";
    case SolveStatus::NonFinite:
      return "NonFinite";
    case SolveStatus::LinearSolveFailed:
      return "LinearSolveFailed";
    case SolveStatus::ExternalSolverFailure:
      return "ExternalSolverFailure";
    case SolveStatus::UnsupportedBackend:
      return "UnsupportedBackend";
  }

  return "Unknown";
}

const char *solveStatusToString(int status)
{
  return isSolveStatusCode(status) ? solveStatusToString(static_cast<SolveStatus>(status)) : "Unknown";
}

bool isSolveStatusCode(int status)
{
  switch (static_cast<SolveStatus>(status)) {
    case SolveStatus::Converged:
    case SolveStatus::MaxIterations:
    case SolveStatus::LineSearchFailed:
    case SolveStatus::StepTooSmall:
    case SolveStatus::NonFinite:
    case SolveStatus::LinearSolveFailed:
    case SolveStatus::ExternalSolverFailure:
    case SolveStatus::UnsupportedBackend:
      return true;
  }

  return false;
}

SolverResult makeSolverResult(SolveStatus status, int rawStatusCode, int iterations)
{
  SolverResult result;
  result.status = status;
  result.rawStatusCode = rawStatusCode;
  result.iterations = iterations;
  result.diagnostics.reset();
  return result;
}

SolverResult makeIpoptSolverResult(int rawStatusCode)
{
  switch (rawStatusCode) {
    case 0:   // Solve_Succeeded
    case 1:   // Solved_To_Acceptable_Level
    case 6:   // Feasible_Point_Found
      return makeSolverResult(SolveStatus::Converged, rawStatusCode);
    case -1:  // Maximum_Iterations_Exceeded
      return makeSolverResult(SolveStatus::MaxIterations, rawStatusCode);
    case 3:   // Search_Direction_Becomes_Too_Small
      return makeSolverResult(SolveStatus::StepTooSmall, rawStatusCode);
    case -13: // Invalid_Number_Detected
      return makeSolverResult(SolveStatus::NonFinite, rawStatusCode);
    case -3:  // Error_In_Step_Computation
      return makeSolverResult(SolveStatus::LinearSolveFailed, rawStatusCode);
    default:
      return makeSolverResult(SolveStatus::ExternalSolverFailure, rawStatusCode);
  }
}

SolverResult makeKnitroSolverResult(int rawStatusCode)
{
  if (rawStatusCode == 0 || rawStatusCode == -100 || rawStatusCode == -101 || rawStatusCode == -102) {
    return makeSolverResult(SolveStatus::Converged, rawStatusCode);
  }

  if (rawStatusCode == -400 || rawStatusCode == -401 || rawStatusCode == -402) {
    return makeSolverResult(SolveStatus::MaxIterations, rawStatusCode);
  }

  if (rawStatusCode <= -500 && rawStatusCode > -600) {
    return makeSolverResult(SolveStatus::LinearSolveFailed, rawStatusCode);
  }

  return makeSolverResult(SolveStatus::ExternalSolverFailure, rawStatusCode);
}

bool acceptsDynamicSolveStatus(SolveStatus status)
{
  return status == SolveStatus::Converged ||
    status == SolveStatus::MaxIterations ||
    status == SolveStatus::StepTooSmall;
}

bool acceptsStrictSolveStatus(SolveStatus status)
{
  return status == SolveStatus::Converged;
}

std::string formatSolverResultSummary(const SolverResult &result)
{
  std::ostringstream ss;
  ss << "status=" << solveStatusToString(result.status)
     << " iterations=" << result.iterations
     << " rawStatusCode=" << result.rawStatusCode;
  if (result.hasFinalGradientStats) {
    ss << " finalGradient=" << result.finalGradientNorm
       << " finalGradientMax=" << result.finalGradientMaxNorm;
  }
  return ss.str();
}

}  // namespace pgo::NonlinearOptimization
