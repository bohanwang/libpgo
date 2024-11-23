/*
author: Bohan Wang
copyright to USC,MIT
*/

#include "elasticModel1DCubicSpline.h"

#include "naturalCubicSplineDerivatives.h"

#include "pgoLogging.h"

#include <memory>

namespace pgo::SolidDeformationModel
{
class ElasticModel1DCubicSplineImpl
{
public:
  std::shared_ptr<NonlinearOptimization::NaturalCubicSpline2DWithParameterDerivatives> spline;
  ES::VXd xValues;
  double c;
};
}  // namespace pgo::SolidDeformationModel

using namespace pgo::SolidDeformationModel;

ElasticModel1DCubicSpline::ElasticModel1DCubicSpline(double coeff, int nPoints, double xLeft, double xRight)
{
  impl = new ElasticModel1DCubicSplineImpl;
  impl->xValues.resize(nPoints);

  int nGaps = nPoints - 1;
  double detlaX = (xRight - xLeft) / nGaps;

  for (int i = 0; i < nPoints; i++)
    impl->xValues[i] = xLeft + i * detlaX;

  impl->spline = std::make_shared<NonlinearOptimization::NaturalCubicSpline2DWithParameterDerivatives>(nPoints, impl->xValues.data());
  impl->c = coeff;
}

ElasticModel1DCubicSpline::~ElasticModel1DCubicSpline()
{
  delete impl;
}

double ElasticModel1DCubicSpline::compute_psi(const double *param, const double F[9],
  const double[], const double[], const double S[]) const
{
  double x = F[0] - F[1];
  double psi = impl->spline->y(x, param) * impl->c;
  //if (psi < -1e-6) {
  //  ES::VXd paramV = ES::Mp<const ES::VXd>(param, impl->xValues.size());
  //  LGI << x << ",--" << paramV.transpose() << "--," << psi;
  //}

  return psi;
}

void ElasticModel1DCubicSpline::compute_P(const double *param, const double F[9],
  const double[], const double[], const double[], double P[]) const
{
  double x = F[0] - F[1];
  P[0] = impl->spline->dy_dx(x, param) * impl->c;
}

void ElasticModel1DCubicSpline::compute_dPdF(const double *param, const double F[9],
  const double[], const double[], const double[], double dPdF[]) const
{
  double x = F[0] - F[1];
  dPdF[0] = impl->spline->d2y_dx2(x, param) * impl->c;
}

int ElasticModel1DCubicSpline::getNumParameters() const
{
  return (int)impl->xValues.size();
}

double ElasticModel1DCubicSpline::compute_dpsi_dparam(const double *param, int i, const double F[],
  const double[], const double[], const double[]) const
{
  double x = F[0] - F[1];
  return impl->spline->dy_dparam(x, param, i) * impl->c;
}
// compute the 2nd order derivative with respect to the (i-th, j-th) parameter
double
ElasticModel1DCubicSpline::compute_d2psi_dparam2(const double *param, int i, int j,
  const double F[], const double[], const double[], const double[]) const
{
  double x = F[0] - F[1];
  return impl->spline->d2y_dparam2(x, param, i, j) * impl->c;
}
// compute the 2nd order derivative with respect to the i-th parameter and F
void ElasticModel1DCubicSpline::compute_dP_dparam(const double *param, int i, const double F[],
  const double[], const double[], const double[], double *ret) const
{
  double x = F[0] - F[1];
  ret[0] = impl->spline->d2y_dparam_dx(x, param, i) * impl->c;
}
