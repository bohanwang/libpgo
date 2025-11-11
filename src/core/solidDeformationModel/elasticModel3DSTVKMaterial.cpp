#include "elasticModel3DSTVKMaterial.h"

#include "polarDecompositionDerivatives.h"

#include <cmath>
#include <cstring>
#include <iostream>

namespace ES = pgo::EigenSupport;

using namespace pgo::SolidDeformationModel;

inline void crossProductMatrix(const ES::V3d &v, ES::M3d &K)
{
  K << 0, -v(2), v(1),
    v(2), 0, -v(0),
    -v(1), v(0), 0;
}

ElasticModel3DSTVKMaterial::ElasticModel3DSTVKMaterial(double _mu, double _lambda):
  mu(_mu), lambda(_lambda)
{
  C[0] << 0, -1, 0,
    1, 0, 0,
    0, 0, 0;
  C[1] << 0, 0, 0,
    0, 0, 1,
    0, -1, 0;
  C[2] << 0, 0, 1,
    0, 0, 0,
    -1, 0, 0;
}

ElasticModel3DSTVKMaterial::~ElasticModel3DSTVKMaterial()
{
}

ES::V3d ElasticModel3DSTVKMaterial::computeLowerInvariance(const ES::M3d &F, ES::M3d &R,EigenSupport::V3d *S_, EigenSupport::M3d *U_, EigenSupport::M3d *V_) const
{
  Eigen::JacobiSVD<ES::M3d, Eigen::NoQRPreconditioner> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);

  ES::M3d U = svd.matrixU();
  ES::M3d V = svd.matrixV();
  ES::V3d S = svd.singularValues();

  if (U.determinant() < 0.0) {
    U.col(2) *= -1.0;
    S(2) *= -1.0;
  }
  if (V.determinant() < 0.0) {
    V.col(2) *= -1.0;
    S(2) *= -1.0;
  }

  // U VT V S VT
  R = U * V.transpose();
  ES::M3d SS = V * S.asDiagonal() * V.transpose();

  ES::V3d i(SS.trace(), (SS * SS).trace(), SS.determinant());

  if (U_ != nullptr) {
    *U_ = U;
  }
  
  if (V_ != nullptr) {
    *V_ = V;
  }

  if (S_ != nullptr) {
    *S_ = S;
  }

  return i;
}

double ElasticModel3DSTVKMaterial::compute_psi(const double *param, const double _F[9], const double _U[], const double _V[], const double _S[]) const
{
  ES::M3d F = Eigen::Map<const ES::M3d>(_F), R;
  ES::V3d i = computeLowerInvariance(F, R);

  double I1 = i[0];
  double I2 = i[1];
  double I3 = i[2];

  double energy = lambda / 8.0 * (I2 - 3) * (I2 - 3) + mu / 8.0 * (8 * I1 * I3 + I2 * I2 + 2 * I1 * I1 * I2 - 4 * I2 - I1 * I1 * I1 * I1 + 6);

  return energy;
}

void ElasticModel3DSTVKMaterial::compute_P(const double *param, const double _F[9], const double _U[], const double _V[], const double _S[], double P[9]) const
{
  ES::M3d F = Eigen::Map<const ES::M3d>(_F), R;
  ES::V3d i = computeLowerInvariance(F, R);

  double I1 = i[0];
  double I2 = i[1];
  double I3 = i[2];

  double dpsi_dI1 = mu * I3 + 0.5 * mu * I1 * I2 - 0.5 * mu * I1 * I1 * I1;
  double dpsi_dI2 = 0.25 * lambda * (I2 - 3) + 0.25 * mu * I2 + 0.25 * mu * I1 * I1 - 0.5 * mu;
  double dpsi_dI3 = mu * I1;

  ES::V9d dI1_dF, dI2_dF, dI3_dF;
  for (int i = 0; i < 3; i++) {
    dI1_dF.segment<3>(3 * i) = R.col(i);
    dI2_dF.segment<3>(3 * i) = 2 * F.col(i);
  }
  dI3_dF.segment<3>(0) = F.col(1).cross(F.col(2));
  dI3_dF.segment<3>(3) = F.col(2).cross(F.col(0));
  dI3_dF.segment<3>(6) = F.col(0).cross(F.col(1));

  ES::V9d g = dpsi_dI1 * dI1_dF + dpsi_dI2 * dI2_dF + dpsi_dI3 * dI3_dF;

  (Eigen::Map<ES::V9d>(P)) = g;
}

void ElasticModel3DSTVKMaterial::compute_dPdF(const double *param, const double _F[9], const double _U[], const double _V[], const double _S[], double dPdF[81]) const
{
  ES::M3d F = Eigen::Map<const ES::M3d>(_F), R, U, V;
  ES::V3d s;
  ES::V3d i = computeLowerInvariance(F, R, &s, &U, &V);

  double I1 = i[0];
  double I2 = i[1];
  double I3 = i[2];

  double dpsi_dI1 = mu * I3 + 0.5 * mu * I1 * I2 - 0.5 * mu * I1 * I1 * I1;
  double dpsi_dI2 = 0.25 * lambda * (I2 - 3) + 0.25 * mu * I2 + 0.25 * mu * I1 * I1 - 0.5 * mu;
  double dpsi_dI3 = mu * I1;

  ES::V9d dI1_dF, dI2_dF, dI3_dF;
  for (int i = 0; i < 3; i++) {
    dI1_dF.segment<3>(3 * i) = R.col(i);
    dI2_dF.segment<3>(3 * i) = 2 * F.col(i);
  }
  dI3_dF.segment<3>(0) = F.col(1).cross(F.col(2));
  dI3_dF.segment<3>(3) = F.col(2).cross(F.col(0));
  dI3_dF.segment<3>(6) = F.col(0).cross(F.col(1));

  double eigv[3];
  eigv[0] = 2.0 / (s(0) + s(1));
  eigv[1] = 2.0 / (s(1) + s(2));
  eigv[2] = 2.0 / (s(2) + s(0));

  ES::M3d Q[3];
  Q[0] = 1.0 / sqrt(2) * U * C[0] * V.transpose();
  Q[1] = 1.0 / sqrt(2) * U * C[1] * V.transpose();
  Q[2] = 1.0 / sqrt(2) * U * C[2] * V.transpose();

  ES::V9d q[3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      q[i].segment<3>(3 * j) = Q[i].col(j);
    }
  }

  ES::M9d d2I1_dF2;
  d2I1_dF2.setZero();

  for (int i = 0; i < 3; i++) {
    d2I1_dF2 += eigv[i] * q[i] * q[i].transpose();
  }

  ES::M9d d2I2_dF2 = 2.0 * ES::M9d::Identity();

  ES::M9d d2I3_dF2;
  d2I3_dF2.setZero();

  ES::M3d f0_hat, f1_hat, f2_hat;
  crossProductMatrix(F.col(0), f0_hat);
  crossProductMatrix(F.col(1), f1_hat);
  crossProductMatrix(F.col(2), f2_hat);
  d2I3_dF2.block<3, 3>(0, 3) = -f2_hat;
  d2I3_dF2.block<3, 3>(0, 6) = f1_hat;
  d2I3_dF2.block<3, 3>(3, 0) = f2_hat;
  d2I3_dF2.block<3, 3>(3, 6) = -f0_hat;
  d2I3_dF2.block<3, 3>(6, 0) = -f1_hat;
  d2I3_dF2.block<3, 3>(6, 3) = f0_hat;

  ES::M9d H;
  H.setZero();
  H += dpsi_dI1 * d2I1_dF2 + dpsi_dI2 * d2I2_dF2 + dpsi_dI3 * d2I3_dF2;

  double d2psi_dI1dI1 = 0.5 * mu * I2 - 1.5 * mu * I1 * I1;
  double d2psi_dI1dI2 = 0.5 * mu * I1;
  double d2psi_dI1dI3 = mu;
  double d2psi_dI2dI2 = 0.25 * lambda + 0.25 * mu;

  H += d2psi_dI1dI1 * dI1_dF * dI1_dF.transpose();
  H += d2psi_dI1dI2 * (dI1_dF * dI2_dF.transpose() + dI2_dF * dI1_dF.transpose());
  H += d2psi_dI1dI3 * (dI1_dF * dI3_dF.transpose() + dI3_dF * dI1_dF.transpose());
  H += d2psi_dI2dI2 * dI2_dF * dI2_dF.transpose();

  (Eigen::Map<ES::M9d>(dPdF)) = H;
}
