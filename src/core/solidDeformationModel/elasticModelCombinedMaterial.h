/*
author: Bohan Wang
copyright to USC,MIT,NUS
*/

#pragma once

#include "elasticModel3DDeformationGradient.h"

#include <stdexcept>
#include <utility>

namespace pgo
{
namespace SolidDeformationModel
{
template<int count>
class ElasticModelCombinedMaterial : public ElasticModel3DDeformationGradient
{
public:
  template<typename... T>
  ElasticModelCombinedMaterial(const T... materials);
  virtual ~ElasticModelCombinedMaterial() {}

  virtual double compute_psi(const double *param, const double F[9],
    const double U[9], const double V[9], const double S[3]) const override;
  virtual void compute_P(const double *param, const double F[9],
    const double U[9], const double V[9], const double S[3], double P[9]) const override;
  virtual void compute_dPdF(const double *param, const double F[9],
    const double U[9], const double V[9], const double S[3], double dPdFOut[81]) const override;

  virtual int getNumParameters() const override { return numTotalParameters; }
  virtual double compute_dpsi_dparam(const double *param, int i, const double F[9],
    const double U[9], const double V[9], const double S[3]) const override;
  virtual double compute_d2psi_dparam2(const double *param, int i, int j, const double F[9],
    const double U[9], const double V[9], const double S[3]) const override;
  virtual void compute_dP_dparam(const double *param, int i, const double F[9],
    const double U[9], const double V[9], const double S[3], double *ret) const override;

  virtual void compute_d2PdF2(const double *param, const double F[9],
    const double U[9], const double V[9], const double S[3], double d3psi_dFiFjFk[729]) const override;
  virtual void compute_d2Pdparam2(const double *param, int i, int j, const double F[9],
    const double U[9], const double V[9], const double S[3], double d2p_dparam2[9]) const override;
  virtual void compute_d2PdFdparam(const double *param, int i, const double F[9],
    const double U[9], const double V[9], const double S[3], double d2P_dFdparam[81]) const override;

  const ElasticModel3DDeformationGradient *getMaterial(int id) const { return materials[id]; }

protected:
  const ElasticModel3DDeformationGradient *materials[count];
  int parameterOffsets[count + 1];
  int numTotalParameters;
};

template<int count>
template<typename... T>
inline ElasticModelCombinedMaterial<count>::ElasticModelCombinedMaterial(const T... mats)
{
  static_assert(count > 0);
  static_assert(
    std::is_same<
      std::integer_sequence<bool, true, std::is_convertible<std::decay_t<T>, ElasticModel3DDeformationGradient *>::value...>,
      std::integer_sequence<bool, std::is_convertible<std::decay_t<T>, const ElasticModel3DDeformationGradient *>::value..., true>>::value,
    "T should be const MuscleElasticModel*");

  // std::cout << __PRETTY_FUNCTION__ << std::endl;

  const ElasticModel3DDeformationGradient *mats2[] = { mats... };

  int offset = 0;
  for (int i = 0; i < count; i++) {
    materials[i] = mats2[i];
    parameterOffsets[i] = offset;
    offset += materials[i]->getNumParameters();
  }

  parameterOffsets[count] = offset;
  numTotalParameters = offset;

  has3rdOrderDerivative = true;
  for (int i = 0; i < count; i++) {
    if (materials[i]->Has3rdOrderDerivative() == false) {
      has3rdOrderDerivative = false;
      break;
    }
  }
}

template<int count>
inline double ElasticModelCombinedMaterial<count>::compute_psi(const double *param,
  const double F[9], const double U[9], const double V[9], const double S[3]) const
{
  double energy = 0;
  for (int i = 0; i < count; i++) {
    energy += materials[i]->compute_psi(param + parameterOffsets[i], F, U, V, S);
  }

  return energy;
}

template<int count>
inline void ElasticModelCombinedMaterial<count>::compute_P(const double *param,
  const double F[9], const double U[9], const double V[9], const double S[3], double P[9]) const
{
  for (int i = 0; i < 9; i++) {
    P[i] = 0;
  }

  for (int i = 0; i < count; i++) {
    double tempP[9];
    materials[i]->compute_P(param + parameterOffsets[i], F, U, V, S, tempP);

    for (int j = 0; j < 9; j++) {
      P[j] += tempP[j];
    }
  }
}

template<int count>
inline void ElasticModelCombinedMaterial<count>::compute_dPdF(const double *param,
  const double F[9], const double U[9], const double V[9], const double S[3], double dPdFOut[81]) const
{
  for (int i = 0; i < 81; i++) {
    dPdFOut[i] = 0;
  }

  for (int i = 0; i < count; i++) {
    double tempdPdF[81];
    materials[i]->compute_dPdF(param + parameterOffsets[i], F, U, V, S, tempdPdF);

    for (int j = 0; j < 81; j++) {
      dPdFOut[j] += tempdPdF[j];
    }
  }
}

template<int count>
inline double ElasticModelCombinedMaterial<count>::compute_dpsi_dparam(const double *param, int i,
  const double F[9], const double U[9], const double V[9], const double S[3]) const
{
  int mi = 0;
  for (; mi < count; mi++) {
    if (i >= parameterOffsets[mi] && i < parameterOffsets[mi + 1])
      break;
  }

  if (mi >= count)
    return 0.0;

  return materials[mi]->compute_dpsi_dparam(param + parameterOffsets[mi], i - parameterOffsets[mi], F, U, V, S);
}

template<int count>
inline double ElasticModelCombinedMaterial<count>::compute_d2psi_dparam2(const double *param, int i, int j,
  const double F[9], const double U[9], const double V[9], const double S[3]) const
{
  int mi = 0;
  for (; mi < count; mi++) {
    if (i >= parameterOffsets[mi] && i < parameterOffsets[mi + 1])
      break;
  }

  if (mi >= count)
    return 0.0;

  int mj = 0;
  for (; mj < count; mj++) {
    if (j >= parameterOffsets[mj] && j < parameterOffsets[mj + 1])
      break;
  }

  if (mj >= count)
    return 0.0;

  if (mi != mj)
    abort();

  return materials[mi]->compute_d2psi_dparam2(param + parameterOffsets[mi], i - parameterOffsets[mi], j - parameterOffsets[mi],
    F, U, V, S);
}

template<int count>
inline void ElasticModelCombinedMaterial<count>::compute_dP_dparam(const double *param, int i,
  const double F[9], const double U[9], const double V[9], const double S[3], double *ret) const
{
  int mi = 0;
  for (; mi < count; mi++) {
    if (i >= parameterOffsets[mi] && i < parameterOffsets[mi + 1])
      break;
  }

  if (mi >= count) {
    for (int i = 0; i < 9; i++)
      ret[i] = 0.0;

    return;
  }

  materials[mi]->compute_dP_dparam(param + parameterOffsets[mi], i - parameterOffsets[mi], F, U, V, S, ret);
}

template<int count>
inline void ElasticModelCombinedMaterial<count>::compute_d2PdF2(const double *param, const double F[9],
  const double U[9], const double V[9], const double S[3], double d3psi_dFiFjFk[729]) const
{
  for (int i = 0; i < 729; i++) {
    d3psi_dFiFjFk[i] = 0;
  }

  for (int i = 0; i < count; i++) {
    double temp[729];
    materials[i]->compute_d2PdF2(param + parameterOffsets[i], F, U, V, S, temp);

    for (int j = 0; j < 729; j++) {
      d3psi_dFiFjFk[j] += temp[j];
    }
  }
}

template<int count>
inline void ElasticModelCombinedMaterial<count>::compute_d2Pdparam2(const double *param, int i, int j, const double F[9],
  const double U[9], const double V[9], const double S[3], double d2p_dparam2[9]) const
{
  int mi = 0;
  for (; mi < count; mi++) {
    if (i >= parameterOffsets[mi] && i < parameterOffsets[mi + 1])
      break;
  }

  for (int i = 0; i < 9; i++)
    d2p_dparam2[i] = 0;

  if (mi >= count)
    return;

  int mj = 0;
  for (; mj < count; mj++) {
    if (j >= parameterOffsets[mj] && j < parameterOffsets[mj + 1])
      break;
  }

  if (mj >= count)
    return;

  if (mi != mj)
    return;

  materials[mi]->compute_d2Pdparam2(param + parameterOffsets[mi], i - parameterOffsets[mi], j - parameterOffsets[mi],
    F, U, V, S, d2p_dparam2);
}

template<int count>
inline void ElasticModelCombinedMaterial<count>::compute_d2PdFdparam(const double *param, int i, const double F[9],
  const double U[9], const double V[9], const double S[3], double d2P_dFdparam[81]) const
{
  int mi = 0;
  for (; mi < count; mi++) {
    if (i >= parameterOffsets[mi] && i < parameterOffsets[mi + 1])
      break;
  }

  if (mi >= count) {
    for (int i = 0; i < 81; i++)
      d2P_dFdparam[i] = 0.0;

    return;
  }

  materials[mi]->compute_d2PdFdparam(param + parameterOffsets[mi], i - parameterOffsets[mi], F, U, V, S, d2P_dFdparam);
}

template<>
class ElasticModelCombinedMaterial<-1> : public ElasticModel3DDeformationGradient
{
public:
  ElasticModelCombinedMaterial(int numMaterials, const ElasticModel3DDeformationGradient *const *mats);
  virtual ~ElasticModelCombinedMaterial() { delete[] materials; }

  virtual double compute_psi(const double *param, const double F[9],
    const double U[9], const double V[9], const double S[3]) const override;
  virtual void compute_P(const double *param, const double F[9],
    const double U[9], const double V[9], const double S[3], double P[9]) const override;
  virtual void compute_dPdF(const double *param, const double F[9],
    const double U[9], const double V[9], const double S[3], double dPdFOut[81]) const override;

  virtual int getNumParameters() const override { return numTotalParameters; }
  virtual double compute_dpsi_dparam(const double *param, int i, const double F[9],
    const double U[9], const double V[9], const double S[3]) const override;
  virtual double compute_d2psi_dparam2(const double *param, int i, int j, const double F[9],
    const double U[9], const double V[9], const double S[3]) const override;
  virtual void compute_dP_dparam(const double *param, int i, const double F[9],
    const double U[9], const double V[9], const double S[3], double *ret) const override;

  virtual void compute_d2PdF2(const double *param, const double F[9],
    const double U[9], const double V[9], const double S[3], double d3psi_dFiFjFk[729]) const override;
  virtual void compute_d2Pdparam2(const double *param, int i, int j, const double F[9],
    const double U[9], const double V[9], const double S[3], double d2p_dparam2[9]) const override;
  virtual void compute_d2PdFdparam(const double *param, int i, const double F[9],
    const double U[9], const double V[9], const double S[3], double d2P_dFdparam[81]) const override;

protected:
  const ElasticModel3DDeformationGradient **materials;
  int count;

  int *parameterOffsets;
  int numTotalParameters;
};

inline ElasticModelCombinedMaterial<-1>::ElasticModelCombinedMaterial(int numMaterials, const ElasticModel3DDeformationGradient *const *mats):
  count(numMaterials)
{
  materials = new const ElasticModel3DDeformationGradient *[numMaterials];
  parameterOffsets = new int[numMaterials + 1];

  int offset = 0;
  for (int i = 0; i < count; i++) {
    materials[i] = mats[i];
    parameterOffsets[i] = offset;
    offset += materials[i]->getNumParameters();
  }

  parameterOffsets[count] = offset;
  numTotalParameters = offset;

  has3rdOrderDerivative = true;
  for (int i = 0; i < count; i++) {
    if (materials[i]->Has3rdOrderDerivative() == false) {
      has3rdOrderDerivative = false;
      break;
    }
  }
}

inline double ElasticModelCombinedMaterial<-1>::compute_psi(const double *param,
  const double F[9], const double U[9], const double V[9], const double S[3]) const
{
  double energy = 0;
  for (int i = 0; i < count; i++) {
    energy += materials[i]->compute_psi(param + parameterOffsets[i], F, U, V, S);
  }

  return energy;
}

inline void ElasticModelCombinedMaterial<-1>::compute_P(const double *param,
  const double F[9], const double U[9], const double V[9], const double S[3], double P[9]) const
{
  for (int i = 0; i < 9; i++) {
    P[i] = 0;
  }

  for (int i = 0; i < count; i++) {
    double tempP[9];
    materials[i]->compute_P(param + parameterOffsets[i], F, U, V, S, tempP);

    for (int j = 0; j < 9; j++) {
      P[j] += tempP[j];
    }
  }
}

inline void ElasticModelCombinedMaterial<-1>::compute_dPdF(const double *param,
  const double F[9], const double U[9], const double V[9], const double S[3], double dPdFOut[81]) const
{
  for (int i = 0; i < 81; i++) {
    dPdFOut[i] = 0;
  }

  for (int i = 0; i < count; i++) {
    double tempdPdF[81];
    materials[i]->compute_dPdF(param + parameterOffsets[i], F, U, V, S, tempdPdF);

    for (int j = 0; j < 81; j++) {
      dPdFOut[j] += tempdPdF[j];
    }
  }
}

inline double ElasticModelCombinedMaterial<-1>::compute_dpsi_dparam(const double *param, int i,
  const double F[9], const double U[9], const double V[9], const double S[3]) const
{
  int mi = 0;
  for (; mi < count; mi++) {
    if (i >= parameterOffsets[mi] && i < parameterOffsets[mi + 1])
      break;
  }

  if (mi >= count)
    return 0.0;

  return materials[mi]->compute_dpsi_dparam(param + parameterOffsets[mi], i - parameterOffsets[mi], F, U, V, S);
}

inline double ElasticModelCombinedMaterial<-1>::compute_d2psi_dparam2(const double *param, int i, int j,
  const double F[9], const double U[9], const double V[9], const double S[3]) const
{
  int mi = 0;
  for (; mi < count; mi++) {
    if (i >= parameterOffsets[mi] && i < parameterOffsets[mi + 1])
      break;
  }

  if (mi >= count)
    return 0.0;

  int mj = 0;
  for (; mj < count; mj++) {
    if (j >= parameterOffsets[mj] && j < parameterOffsets[mj + 1])
      break;
  }

  if (mj >= count)
    return 0.0;

  if (mi != mj)
    throw std::runtime_error("mi should be equal to mj");

  return materials[mi]->compute_d2psi_dparam2(param + parameterOffsets[mi], i - parameterOffsets[mi], j - parameterOffsets[mi],
    F, U, V, S);
}

inline void ElasticModelCombinedMaterial<-1>::compute_dP_dparam(const double *param, int i,
  const double F[9], const double U[9], const double V[9], const double S[3], double *ret) const
{
  int mi = 0;
  for (; mi < count; mi++) {
    if (i >= parameterOffsets[mi] && i < parameterOffsets[mi + 1])
      break;
  }

  if (mi >= count) {
    for (int i = 0; i < 9; i++)
      ret[i] = 0.0;

    return;
  }

  materials[mi]->compute_dP_dparam(param + parameterOffsets[mi], i - parameterOffsets[mi], F, U, V, S, ret);
}

inline void ElasticModelCombinedMaterial<-1>::compute_d2PdF2(const double *param, const double F[9],
  const double U[9], const double V[9], const double S[3], double d3psi_dFiFjFk[729]) const
{
  for (int i = 0; i < 729; i++) {
    d3psi_dFiFjFk[i] = 0;
  }

  for (int i = 0; i < count; i++) {
    double temp[729];
    materials[i]->compute_d2PdF2(param + parameterOffsets[i], F, U, V, S, temp);

    for (int j = 0; j < 729; j++) {
      d3psi_dFiFjFk[j] += temp[j];
    }
  }
}

inline void ElasticModelCombinedMaterial<-1>::compute_d2Pdparam2(const double *param, int i, int j, const double F[9],
  const double U[9], const double V[9], const double S[3], double d2p_dparam2[9]) const
{
  int mi = 0;
  for (; mi < count; mi++) {
    if (i >= parameterOffsets[mi] && i < parameterOffsets[mi + 1])
      break;
  }

  for (int i = 0; i < 9; i++)
    d2p_dparam2[i] = 0;

  if (mi >= count)
    return;

  int mj = 0;
  for (; mj < count; mj++) {
    if (j >= parameterOffsets[mj] && j < parameterOffsets[mj + 1])
      break;
  }

  if (mj >= count)
    return;

  if (mi != mj)
    return;

  materials[mi]->compute_d2Pdparam2(param + parameterOffsets[mi], i - parameterOffsets[mi], j - parameterOffsets[mi],
    F, U, V, S, d2p_dparam2);
}

inline void ElasticModelCombinedMaterial<-1>::compute_d2PdFdparam(const double *param, int i, const double F[9],
  const double U[9], const double V[9], const double S[3], double d2P_dFdparam[81]) const
{
  int mi = 0;
  for (; mi < count; mi++) {
    if (i >= parameterOffsets[mi] && i < parameterOffsets[mi + 1])
      break;
  }

  if (mi >= count) {
    for (int i = 0; i < 81; i++)
      d2P_dFdparam[i] = 0.0;

    return;
  }

  materials[mi]->compute_d2PdFdparam(param + parameterOffsets[mi], i - parameterOffsets[mi], F, U, V, S, d2P_dFdparam);
}

}  // namespace SolidDeformationModel
}  // namespace pgo