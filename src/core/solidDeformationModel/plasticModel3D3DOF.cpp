/*
author: Bohan Wang
copyright to USC,MIT
*/

#include "plasticModel3D3DOF.h"

#include "EigenSupport.h"

#include <cstring>

namespace ES = pgo::EigenSupport; 

using Map3 = Eigen::Map<ES::M3d>;
using MapC3 = Eigen::Map<const ES::M3d>;

using namespace pgo::SolidDeformationModel;

PlasticModel3D3DOF::PlasticModel3D3DOF(const double R_[9]):
  PlasticModel(3)
{
  std::memcpy(R, R_, sizeof(double) * 9);
  (Map3(RT)) = Map3(R).transpose();
}

void PlasticModel3D3DOF::setR(const double R_[9])
{
  std::memcpy(R, R_, sizeof(double) * 9);
  (Map3(RT)) = Map3(R).transpose();
}

void PlasticModel3D3DOF::computeA(const double *param, double A[9]) const
{
  ES::M3d S = ES::V3d(param[0], param[1], param[2]).asDiagonal();
  (Map3(A)) = (MapC3(RT)) * S * (MapC3(R));
}

void PlasticModel3D3DOF::computeAInv(const double *param, double AInv[9]) const
{
  ES::M3d SInv = ES::M3d::Identity();
  for (int i = 0; i < 3; i++) {
    if (param[i] < zeroThreshold)
      SInv(i, i) = 1.0 / zeroThreshold;
    else {
      SInv(i, i) = 1.0 / param[i];
    }
  }

  (Map3(AInv)) = (MapC3(RT)) * SInv * (MapC3(R));
}

double PlasticModel3D3DOF::compute_detA(const double *param) const
{
  return param[0] * param[1] * param[2];
}

void PlasticModel3D3DOF::compute_ddetA_da(const double *param, double *ddetA_da, int) const
{
  ddetA_da[0] = param[1] * param[2];
  ddetA_da[1] = param[0] * param[2];
  ddetA_da[2] = param[0] * param[1];
}

void PlasticModel3D3DOF::compute_d2detA_da2(const double *param, double *d2detA_da2, int leadingDim) const
{
#define ELIDX(i, j) ((j)*leadingDim + (i))

  d2detA_da2[ELIDX(0, 0)] = 0;
  d2detA_da2[ELIDX(0, 1)] = param[2];
  d2detA_da2[ELIDX(0, 2)] = param[1];

  d2detA_da2[ELIDX(1, 0)] = param[2];
  d2detA_da2[ELIDX(1, 1)] = 0;
  d2detA_da2[ELIDX(1, 2)] = param[0];

  d2detA_da2[ELIDX(2, 0)] = param[1];
  d2detA_da2[ELIDX(2, 1)] = param[0];
  d2detA_da2[ELIDX(2, 2)] = 0;

#undef ELIDX
}

void PlasticModel3D3DOF::compute_dAInv_da(const double *param, int pi, double ret[9]) const
{
  ES::M3d dSInv_dai = ES::M3d::Zero();

  double val = param[pi];
  if (val < zeroThreshold)
    val = zeroThreshold;

  dSInv_dai(pi, pi) = -1.0 / (val * val);
  (Map3(ret)) = (MapC3(RT)) * dSInv_dai * (MapC3(R));
}

void PlasticModel3D3DOF::compute_d2AInv_da2(const double *param, int pi, int pj, double ret[9]) const
{
  ES::M3d d2SInv_dai_daj = ES::M3d::Zero();

  if (pi == pj) {
    double val = param[pi];
    if (val < zeroThreshold)
      val = zeroThreshold;

    d2SInv_dai_daj(pi, pj) = 2.0 / (val * val * val);
  }

  (Map3(ret)) = (MapC3(RT)) * d2SInv_dai_daj * (MapC3(R));
}

void PlasticModel3D3DOF::projectParam(double *param, double zeroThreshold) const
{
  if (param[0] < zeroThreshold)
    param[0] = zeroThreshold;

  if (param[1] < zeroThreshold)
    param[1] = zeroThreshold;

  if (param[2] < zeroThreshold)
    param[2] = zeroThreshold;
}

void PlasticModel3D3DOF::compute_dparamfull_dparamsub(const double *, const double *basis, int numHandles, double *dpf_dps) const
{
  Eigen::Map<ES::MXd> B(dpf_dps, 3, 3 * numHandles);
  for (int i = 0; i < numHandles; i++) {
    B.block<3, 3>(0, i * 3) = ES::M3d::Identity() * basis[i];
  }
}

void PlasticModel3D3DOF::compute_d2paramfull_dparamsub2(const double *param, const double *basis, int numHandles, int pi, double *d2a_dz2) const
{
  memset(d2a_dz2, 0, (3 * numHandles) * (3 * numHandles) * sizeof(double));
}