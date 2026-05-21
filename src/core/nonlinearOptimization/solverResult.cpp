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
