#pragma once

#include <algorithm>
#include <cstdint>

namespace pgo::NonlinearOptimization
{

struct MaxStepResult
{
  double alpha = 1.0;
  double materialAlpha = 1.0;
  double contactAlpha = 1.0;
  bool materialClamped = false;
  bool contactClamped = false;

  static MaxStepResult unconstrained()
  {
    return {};
  }

  static MaxStepResult material(double alpha)
  {
    MaxStepResult result;
    result.alpha = alpha;
    result.materialAlpha = alpha;
    result.materialClamped = alpha < 1.0;
    return result;
  }

  static MaxStepResult contact(double alpha)
  {
    MaxStepResult result;
    result.alpha = alpha;
    result.contactAlpha = alpha;
    result.contactClamped = alpha < 1.0;
    return result;
  }
};

inline MaxStepResult mergeMaxStepResults(const MaxStepResult &lhs, const MaxStepResult &rhs)
{
  MaxStepResult result;
  result.alpha = std::min(lhs.alpha, rhs.alpha);
  result.materialAlpha = std::min(lhs.materialAlpha, rhs.materialAlpha);
  result.contactAlpha = std::min(lhs.contactAlpha, rhs.contactAlpha);
  result.materialClamped = lhs.materialClamped || rhs.materialClamped;
  result.contactClamped = lhs.contactClamped || rhs.contactClamped;
  return result;
}

struct SolveDiagnostics
{
  std::int64_t materialClampCount = 0;
  std::int64_t contactClampCount = 0;
  double minFeasibleAlpha = 1.0;
  double minMaterialFeasibleAlpha = 1.0;
  double minContactFeasibleAlpha = 1.0;
  double minLineSearchAlpha = 1.0;
  double minEffectiveAlpha = 1.0;
  double currentMaterialAlpha = 1.0;
  double currentContactAlpha = 1.0;

  void reset()
  {
    *this = SolveDiagnostics{};
  }

  void recordMaxStep(const MaxStepResult &result)
  {
    currentMaterialAlpha = result.materialAlpha;
    currentContactAlpha = result.contactAlpha;
    minFeasibleAlpha = std::min(minFeasibleAlpha, result.alpha);
    minMaterialFeasibleAlpha = std::min(minMaterialFeasibleAlpha, result.materialAlpha);
    minContactFeasibleAlpha = std::min(minContactFeasibleAlpha, result.contactAlpha);

    if (result.materialClamped)
      materialClampCount += 1;
    if (result.contactClamped)
      contactClampCount += 1;
  }

  void recordLineSearch(double feasibleAlpha, double lineSearchAlpha, double effectiveAlpha)
  {
    minFeasibleAlpha = std::min(minFeasibleAlpha, feasibleAlpha);
    minLineSearchAlpha = std::min(minLineSearchAlpha, lineSearchAlpha);
    minEffectiveAlpha = std::min(minEffectiveAlpha, effectiveAlpha);
  }
};

}  // namespace pgo::NonlinearOptimization

namespace pgo
{
using NonlinearOptimization::MaxStepResult;
using NonlinearOptimization::SolveDiagnostics;
}  // namespace pgo
