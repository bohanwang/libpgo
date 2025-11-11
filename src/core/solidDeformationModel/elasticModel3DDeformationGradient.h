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
class ElasticModel3DDeformationGradient : public ElasticModel
{
public:
  ElasticModel3DDeformationGradient() {}
  virtual ~ElasticModel3DDeformationGradient() {}

  virtual double compute_psi(const double *param, const double F[9],
    const double U[9], const double V[9], const double S[3]) const = 0;
  virtual void compute_P(const double *param, const double F[9],
    const double U[9], const double V[9], const double S[3], double P[9]) const = 0;
  virtual void compute_dPdF(const double *param, const double F[9],
    const double U[9], const double V[9], const double S[3], double dPdF[81]) const = 0;

  int getNumParameters() const override { return 0; }

  // compute the 1st order derivative with respect to the i-th parameter
  virtual double compute_dpsi_dparam(const double *param, int i, const double F[9],
    const double U[9], const double V[9], const double S[3]) const;
  // compute the 2nd order derivative with respect to the (i-th, j-th) parameter
  virtual double compute_d2psi_dparam2(const double *param, int i, int j,
    const double F[9], const double U[9], const double V[9], const double S[3]) const;
  // compute the 2nd order derivative with respect to the i-th parameter and F
  virtual void compute_dP_dparam(const double *param, int i, const double F[9],
    const double U[9], const double V[9], const double S[3], double *ret) const;

  // tensor[k * 81 + j * 9 + i] should containt d3psi / (dFi dFj dFk)
  virtual void compute_d2PdF2(const double *param, const double F[9],
    const double U[9], const double V[9], const double S[3], double d3psi_dFiFjFk[729]) const;
  virtual void compute_d2Pdparam2(const double *param, int i, int j, const double F[9],
    const double U[9], const double V[9], const double S[3], double d2p_dparam2[9]) const;
  virtual void compute_d2PdFdparam(const double *param, int i, const double F[9],
    const double U[9], const double V[9], const double S[3], double d2P_dFdparam[81]) const;

  bool Has3rdOrderDerivative() const { return has3rdOrderDerivative; }

protected:
  bool has3rdOrderDerivative = false;
};

inline double ElasticModel3DDeformationGradient::compute_dpsi_dparam(const double *, int,
  const double[9], const double[9], const double[9], const double[3]) const
{
  return 0;
}

inline double ElasticModel3DDeformationGradient::compute_d2psi_dparam2(const double *, int, int,
  const double[9], const double[9], const double[9], const double[3]) const
{
  return 0;
}

inline void ElasticModel3DDeformationGradient::compute_dP_dparam(const double *, int,
  const double[9], const double[9], const double[9], const double[3], double *) const
{
}

inline void ElasticModel3DDeformationGradient::compute_d2PdF2(const double *, const double[9],
  const double[9], const double[9], const double[3], double[729]) const
{
}

inline void ElasticModel3DDeformationGradient::compute_d2Pdparam2(const double *, int, int, const double[9],
  const double[9], const double[9], const double[3], double[9]) const
{
}

inline void ElasticModel3DDeformationGradient::compute_d2PdFdparam(const double *, int, const double[9],
  const double[9], const double[9], const double[3], double[81]) const
{
}

}  // namespace SolidDeformationModel
}  // namespace pgo