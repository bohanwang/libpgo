/*
author: Bohan Wang
copyright to USC,MIT
*/

#include "elasticModel1DQuadratic.h"

using namespace pgo::SolidDeformationModel;

double ElasticModel1DQuadratic::compute_psi(const double *param, const double F[9],
  const double[9], const double[9], const double[3]) const
{
  return param[0] * coeff * (F[0] - F[1]) * (F[0] - F[1]) * 0.5;
}

void ElasticModel1DQuadratic::compute_P(const double *param, const double F[9],
  const double[9], const double[9], const double[3], double P[9]) const
{
  P[0] = param[0] * coeff * (F[0] - F[1]);
}

void ElasticModel1DQuadratic::compute_dPdF(const double *param, const double F[9],
  const double[], const double[], const double[], double dPdF[81]) const
{
  dPdF[0] = param[0] * coeff;
}

double ElasticModel1DQuadratic::compute_dpsi_dparam(const double *, int, const double F[9],
  const double[9], const double[9], const double[3]) const
{
  return coeff * (F[0] - F[1]) * (F[0] - F[1]) * 0.5;
}
// compute the 2nd order derivative with respect to the (i-th, j-th) parameter
double ElasticModel1DQuadratic::compute_d2psi_dparam2(const double *, int, int,
  const double[], const double[], const double[], const double[]) const
{
  return 0;
}
// compute the 2nd order derivative with respect to the i-th parameter and F
void ElasticModel1DQuadratic::compute_dP_dparam(const double *, int, const double F[9],
  const double[], const double[], const double[], double *ret) const
{
  ret[0] = coeff * (F[0] - F[1]);
}
