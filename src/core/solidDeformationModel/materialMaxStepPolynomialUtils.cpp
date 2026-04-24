

#include "materialMaxStepPolynomialUtils.h"

#include <algorithm>
#include <cmath>

namespace pgo::SolidDeformationModel
{
namespace ES = pgo::EigenSupport;

namespace
{
constexpr double kBoundaryEvalTol = 1e-12;
}

double detFromColumns(const ES::V3d &c0, const ES::V3d &c1, const ES::V3d &c2)
{
  ES::M3d M;
  M.col(0) = c0;
  M.col(1) = c1;
  M.col(2) = c2;
  return M.determinant();
}

BasicAlgorithms::CubicPolynomial buildDeterminantCubicFromAffineMatrixPath(const double AData[9], const double BData[9], double eps)
{
  const Eigen::Map<const ES::M3d> A(AData);
  const Eigen::Map<const ES::M3d> B(BData);

  const ES::V3d a0 = A.col(0);
  const ES::V3d a1 = A.col(1);
  const ES::V3d a2 = A.col(2);
  const ES::V3d b0 = B.col(0);
  const ES::V3d b1 = B.col(1);
  const ES::V3d b2 = B.col(2);

  BasicAlgorithms::CubicPolynomial poly;
  poly.c0 = detFromColumns(a0, a1, a2) - eps;
  poly.c1 = detFromColumns(b0, a1, a2) + detFromColumns(a0, b1, a2) + detFromColumns(a0, a1, b2);
  poly.c2 = detFromColumns(b0, b1, a2) + detFromColumns(b0, a1, b2) + detFromColumns(a0, b1, b2);
  poly.c3 = detFromColumns(b0, b1, b2);
  return poly;
}

double applyMaterialMaxStepSafetyClamp(double alpha)
{
  if (alpha >= 1.0) {
    return 1.0;
  }

  return std::max(kMaterialMaxStepMinClamp, std::min(1.0, alpha * kMaterialMaxStepInteriorSafety));
}

ConservativeFeasibleAlphaResult findConservativeFeasibleAlpha(const BasicAlgorithms::CubicPolynomial &poly, double eps)
{
  ConservativeFeasibleAlphaResult result;
  result.phi0 = poly.eval(0.0) + eps;
  if (result.phi0 <= eps) {
    result.alpha = kMaterialMaxStepMinClamp;
    result.illegalInitialState = true;
    return result;
  }

  std::vector<double> cuts;
  cuts.reserve(4);
  cuts.push_back(0.0);
  std::vector<double> criticalPoints = BasicAlgorithms::findCriticalPointsInUnitInterval(poly);
  cuts.insert(cuts.end(), criticalPoints.begin(), criticalPoints.end());
  cuts.push_back(1.0);

  std::sort(cuts.begin(), cuts.end());
  cuts.erase(std::unique(cuts.begin(), cuts.end(), [](double lhs, double rhs) {
    return std::abs(lhs - rhs) <= kBoundaryEvalTol;
  }), cuts.end());

  for (size_t i = 1; i < cuts.size(); i++) {
    const double left = cuts[i - 1];
    const double right = cuts[i];
    const double gr = poly.eval(right);

    if (gr > kBoundaryEvalTol) {
      continue;
    }

    if (std::abs(gr) <= kBoundaryEvalTol) {
      result.alpha = applyMaterialMaxStepSafetyClamp(right);
      return result;
    }

    result.alpha = applyMaterialMaxStepSafetyClamp(BasicAlgorithms::findFirstBoundaryRootByBisection(poly, left, right));
    return result;
  }

  result.alpha = 1.0;
  return result;
}
}
