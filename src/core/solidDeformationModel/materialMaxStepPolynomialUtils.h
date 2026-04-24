

#pragma once

#include "EigenDef.h"
#include "polynomialRootUtils.h"

#include <limits>

namespace pgo::SolidDeformationModel
{
inline constexpr double kMaterialMaxStepInteriorSafety = 0.99;
inline constexpr double kMaterialMaxStepMinClamp = 1e-12;
inline constexpr double kTetRelativeDetEps = 1e-8;
inline constexpr double kCubicRelativeDetEps = 1e-8;

struct ConservativeFeasibleAlphaResult
{
  double alpha = 1.0;
  bool illegalInitialState = false;
  double phi0 = std::numeric_limits<double>::quiet_NaN();
};

double detFromColumns(const EigenSupport::V3d &c0, const EigenSupport::V3d &c1, const EigenSupport::V3d &c2);
BasicAlgorithms::CubicPolynomial buildDeterminantCubicFromAffineMatrixPath(const double A[9], const double B[9], double eps);
double applyMaterialMaxStepSafetyClamp(double alpha);
ConservativeFeasibleAlphaResult findConservativeFeasibleAlpha(const BasicAlgorithms::CubicPolynomial &poly, double eps);
}
