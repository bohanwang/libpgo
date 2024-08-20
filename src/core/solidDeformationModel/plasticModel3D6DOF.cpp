/*
author: Bohan Wang
copyright to USC,MIT
*/

#include "plasticModel3D6DOF.h"

#include "EigenSupport.h"
#include "determinantDerivatives.h"

#include <cstring>

namespace ES = pgo::EigenSupport;
using Map3 = Eigen::Map<ES::M3d>;
using MapC3 = Eigen::Map<const ES::M3d>;

using namespace pgo::SolidDeformationModel;
using namespace pgo::NonlinearOptimization;

PlasticModel3D6DOF::PlasticModel3D6DOF(const double R_[9]):
  PlasticModel(6)
{
  if (R_) {
    std::memcpy(R, R_, sizeof(double) * 9);
    (Map3(RT)) = Map3(R).transpose();
  }
  else {
    (Map3(R)) = ES::M3d::Identity();
    (Map3(RT)) = ES::M3d::Identity();
  }
}

inline ES::M3d getS(const double *param)
{
  ES::M3d S;
  S << param[0], param[1], param[2],
    param[1], param[3], param[4],
    param[2], param[4], param[5];

  return S;
}

void PlasticModel3D6DOF::computeA(const double *param, double A[9]) const
{
  ES::M3d S = getS(param);
  (Map3(A)) = (MapC3(RT)) * S * (MapC3(R));
}

void PlasticModel3D6DOF::computeAInv(const double *param, double AInv[9]) const
{
  ES::M3d S = getS(param);
  ES::M3d SInv = S.fullPivLu().inverse();
  (Map3(AInv)) = (MapC3(RT)) * SInv * (MapC3(R));
}

double PlasticModel3D6DOF::compute_detA(const double *param) const
{
  ES::M3d S = getS(param);
  return S.determinant();
}

void PlasticModel3D6DOF::compute_ddetA_da(const double *param, double *ddetA_da, int) const
{
  Determinant::Dim3::ddetA_dA_sym(param, ddetA_da);
}

void PlasticModel3D6DOF::compute_d2detA_da2(const double *param, double *d2detA_da2, int leadingDim) const
{
  double deriv[36];
  Determinant::Dim3::d2detA_dA2_sym(param, deriv);

  (Eigen::Map<ES::MXd>(d2detA_da2, leadingDim, leadingDim)).block(0, 0, 6, 6) = Eigen::Map<Eigen::Matrix<double, 6, 6>>(deriv);
}

void PlasticModel3D6DOF::compute_dAInv_da(const double *param, int pi, double ret[9]) const
{
  ES::M3d S = getS(param);

  const double ei[6][9] = {
    { 1, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 1, 0, 1, 0, 0, 0, 0, 0 },
    { 0, 0, 1, 0, 0, 0, 1, 0, 0 },
    { 0, 0, 0, 0, 1, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0, 1, 0, 1, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 1 }
  };

  ES::M3d SInv = S.fullPivLu().inverse();
  ES::M3d Z = -SInv * MapC3(ei[pi]) * SInv;
  (Map3(ret)) = (MapC3(RT)) * Z * (MapC3(R));
}

void PlasticModel3D6DOF::compute_d2AInv_da2(const double *param, int pi, int pj, double ret[9]) const
{
  ES::M3d S = getS(param);

  const double ei[6][9] = {
    { 1, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 1, 0, 1, 0, 0, 0, 0, 0 },
    { 0, 0, 1, 0, 0, 0, 1, 0, 0 },
    { 0, 0, 0, 0, 1, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0, 1, 0, 1, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 1 }
  };

  ES::M3d SInv = S.fullPivLu().inverse();

  ES::M3d Z = SInv * MapC3(ei[pi]) * SInv * MapC3(ei[pj]) * SInv +
    SInv * MapC3(ei[pj]) * SInv * MapC3(ei[pi]) * SInv;

  (Map3(ret)) = (MapC3(RT)) * Z * (MapC3(R));
}

void PlasticModel3D6DOF::projectParam(double *param, double zeroThreshold) const
{
  ES::M3d Fp;
  computeA(param, Fp.data());

  Eigen::SelfAdjointEigenSolver<ES::M3d> eig(Fp);
  ES::M3d R = eig.eigenvectors();
  ES::V3d S = eig.eigenvalues();
  ES::V3d Struncate = S.cwiseMax(ES::V3d::Constant(zeroThreshold));

  ES::M3d FpPrime = R * Struncate.asDiagonal() * R.transpose();

  param[0] = FpPrime(0, 0);
  param[1] = FpPrime(0, 1);
  param[2] = FpPrime(0, 2);

  param[3] = FpPrime(1, 1);
  param[4] = FpPrime(1, 2);

  param[5] = FpPrime(2, 2);
}

void PlasticModel3D6DOF::compute_dparamfull_dparamsub(const double *, const double *basis, int numHandles, double *dpf_dps) const
{
  Eigen::Map<ES::MXd> B(dpf_dps, 6, 6 * numHandles);
  for (int i = 0; i < numHandles; i++) {
    B.block<6, 6>(0, i * 6) = Eigen::Matrix<double, 6, 6>::Identity() * basis[i];
  }
}

void PlasticModel3D6DOF::compute_d2paramfull_dparamsub2(const double *param, const double *basis, int numHandles, int pi, double *d2a_dz2) const
{
  memset(d2a_dz2, 0, (6 * numHandles) * (6 * numHandles) * sizeof(double));
}

void PlasticModel3D6DOF::computeR(const double *param, double ROut[9]) const
{
  ES::M3d Fp;
  computeA(param, Fp.data());

  Eigen::SelfAdjointEigenSolver<ES::M3d> eig(Fp);
  ES::M3d R = eig.eigenvectors();

  (Eigen::Map<ES::M3d>(ROut)) = R;
}