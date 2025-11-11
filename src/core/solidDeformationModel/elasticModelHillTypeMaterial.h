/*
author: Bohan Wang
copyright to USC,MIT,NUS
*/

#pragma once

#include "elasticModel3DDeformationGradient.h"

namespace pgo
{
namespace SolidDeformationModel
{
class ElasticModelHillTypeMaterial : public ElasticModel3DDeformationGradient
{
public:
  ElasticModelHillTypeMaterial(double shapeParam, double maximalContractionForce, double optimalLengthRatio, const double fiberDirection[3]);
  virtual ~ElasticModelHillTypeMaterial() {}

  void enforceSPD(bool enforce) { enforceSPD_ = enforce ? 1 : 0; }

  virtual double compute_psi(const double *param, const double F[9],
    const double U[9], const double V[9], const double S[3]) const override;
  virtual void compute_P(const double *param, const double F[9],
    const double U[9], const double V[9], const double S[3], double P[9]) const override;
  virtual void compute_dPdF(const double *param, const double F[9],
    const double U[9], const double V[9], const double S[3], double dPdFOut[81]) const override;

  virtual int getNumParameters() const override { return 1; }
  virtual double compute_dpsi_dparam(const double *param, int i, const double F[9],
    const double U[9], const double V[9], const double S[3]) const override;
  virtual double compute_d2psi_dparam2(const double *param, int i, int j, const double F[9],
    const double U[9], const double V[9], const double S[3]) const override;
  virtual void compute_dP_dparam(const double *param, int i, const double F[9],
    const double U[9], const double V[9], const double S[3], double *ret) const override;

protected:
  double compute_length(const double F[9], double Fd[] = nullptr) const;
  void compute_dldF(const double F[9], const double Fd[], double dldF[9]) const;
  void compute_d2ldF2(const double F[9], const double Fd[], double d2ldF2[81]) const;

  double gamma;
  double maxf;
  double lo;
  double fiberDirection[3];

  double sqrt_gamma;
  double erf_sqrt_gamma;
  double sqrt_pi;
  double dFddT_dF[81];

  int enforceSPD_ = 1;
};

}  // namespace SolidDeformationModel
}  // namespace pgo