/*
author: Bohan Wang
copyright to USC
*/

#pragma once

#include "EigenSupport.h"

#include <memory>

namespace pgo
{
namespace ES = EigenSupport;

namespace NonlinearOptimization
{
class NaturalCubicSpline2DAsNonlinearConstraints;

class NaturalCubicSplineFitting
{
public:
  NaturalCubicSplineFitting(const EigenSupport::VXd &xVals, const EigenSupport::VXd &yVals, int numNodes, double leftDeriv, double rightDeriv, double minXDist = 0.1);
  void buildSpline(const double* cpX, const double* cpY);
  void buildSpline(const double *cpX, const double *cpY, EigenSupport::VXd &param) const;
  int fit(const char *solverConfigFilename);

  double func(double x) const;

  int getNumControlPoints() const { return numControlPoints; }
  EigenSupport::V2d getControlPoint(int i) const { return EigenSupport::V2d(controlPointsX[i], controlPointsY[i]); }

  int getNumFunctions() const { return numControlPoints - 1; }
  EigenSupport::V4d getFunctionCoeff(int i) const { return parameters.segment<4>(i * 4); }

protected:
  EigenSupport::VXd xVals;
  EigenSupport::VXd yVals;

  int numControlPoints;
  double leftDerivative, rightDerivative;
  double minXDist;

  int nAll;

  EigenSupport::VXd controlPointsX, controlPointsY;
  EigenSupport::VXd parameters;

  std::shared_ptr<NaturalCubicSpline2DAsNonlinearConstraints> splineC;
};
}  // namespace NonlinearOptimization
}  // namespace pgo
