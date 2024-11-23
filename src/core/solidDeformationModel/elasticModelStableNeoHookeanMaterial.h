/*
author: Bohan Wang
copyright to USC,MIT
*/

#pragma once

#include "elasticModel.h"

namespace pgo
{
namespace SolidDeformationModel
{
class ElasticModelStableNeoHookeanMaterial : public ElasticModel
{
public:
  ElasticModelStableNeoHookeanMaterial(double mu, double lambda);
  virtual ~ElasticModelStableNeoHookeanMaterial();

  void enforceSPD(bool enforce) { enforceSPD_ = enforce ? 1 : 0; }

  virtual double compute_psi(const double *param, const double F[9],
    const double U[9], const double V[9], const double S[3]) const override;
  virtual void compute_P(const double *param, const double F[9],
    const double U[9], const double V[9], const double S[3], double P[9]) const override;
  virtual void compute_dPdF(const double *param, const double F[9],
    const double U[9], const double V[9], const double S[3], double dPdFOut[81]) const override;

  void setMaterial(double mu_, double lambda_);

protected:
  double _mu, _lambda, _ratio;
  int enforceSPD_ = 0;
};
}  // namespace SolidDeformationModel
}  // namespace pgo