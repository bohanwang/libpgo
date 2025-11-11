/*
author: Bohan Wang
copyright to USC,MIT,NUS
*/

#pragma once

#include "elasticModel.h"

namespace pgo
{
namespace SolidDeformationModel
{
class ElasticModel1DQuadratic : public ElasticModel
{
public:
  ElasticModel1DQuadratic(double c):
    coeff(c)
  {
    has3rdOrderDerivative = false;
  }

  virtual ~ElasticModel1DQuadratic() {}
  // F will be always two numbers, cur value, rest value
  virtual double compute_psi(const double *param, const double F[9],
    const double U[9], const double V[9], const double S[3]) const override;
  virtual void compute_P(const double *param, const double F[9],
    const double U[9], const double V[9], const double S[3], double P[9]) const override;
  virtual void compute_dPdF(const double *param, const double F[9],
    const double U[9], const double V[9], const double S[3], double dPdF[81]) const override;

  virtual int getNumParameters() const { return 1; };
  // compute the 1st order derivative with respect to the i-th parameter
  virtual double compute_dpsi_dparam(const double *param, int i, const double F[9],
    const double U[9], const double V[9], const double S[3]) const override;
  // compute the 2nd order derivative with respect to the (i-th, j-th) parameter
  virtual double compute_d2psi_dparam2(const double *param, int i, int j,
    const double F[9], const double U[9], const double V[9], const double S[3]) const override;
  // compute the 2nd order derivative with respect to the i-th parameter and F
  virtual void compute_dP_dparam(const double *param, int i, const double F[9],
    const double U[9], const double V[9], const double S[3], double *ret) const override;
protected:
  double coeff;
};

}  // namespace SolidDeformationModel
}  // namespace pgo