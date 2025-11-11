#include "elasticModel2DFundamentalFormsSTVK.h"

#include "EigenDef.h"

using namespace pgo;
using namespace SolidDeformationModel;

namespace ES = pgo::EigenSupport;

static inline ES::M4d KroneckerProduct2(const ES::M2d &A, const ES::M2d &B)
{
  ES::M4d ret;

  for (int col = 0; col < 2; col++) {
    for (int row = 0; row < 2; row++) {
      ret.block<2, 2>(row * 2, col * 2) = A(row, col) * B;
    }
  }

  return ret;
}

static inline ES::V4d vecCM(const ES::M2d &X)
{  // column-major vec
  ES::V4d v;
  v << X(0, 0), X(1, 0), X(0, 1), X(1, 1);
  return v;
}

double ElasticModel2DFundamentalFormsSTVK::compute_psi_a(const double *param, const double a_[4], const double abar_[4]) const
{
  double Es = param[0];
  double nu_s = param[1];
  double h = param[4];

  double lameAlpha = Es * nu_s / ((1 + nu_s) * (1 - 2 * nu_s));
  double lameBeta = Es / (2 * (1 + nu_s));

  ES::M2d a = ES::Mp<const ES::M2d>(a_);
  ES::M2d abar = ES::Mp<const ES::M2d>(abar_);
  ES::M2d a_bar_inv = abar.fullPivHouseholderQr().inverse();
  ES::M2d M = a_bar_inv * a - ES::M2d::Identity();

  double trace_M = M.trace();
  double Wc = 0.5 * lameAlpha * (trace_M * trace_M) + lameBeta * (M * M).trace();

  return Wc * h;
}

void ElasticModel2DFundamentalFormsSTVK::compute_dpsi_da(const double *param, const double a_[4], const double abar_[4], double da_[4]) const
{
  double Es = param[0];
  double nu_s = param[1];
  double h = param[4];

  double lameAlpha = Es * nu_s / ((1 + nu_s) * (1 - 2 * nu_s));
  double lameBeta = Es / (2 * (1 + nu_s));

  ES::M2d a = ES::Mp<const ES::M2d>(a_);
  ES::M2d abar = ES::Mp<const ES::M2d>(abar_);
  ES::M2d a_bar_inv = abar.fullPivHouseholderQr().inverse();
  ES::M2d M = a_bar_inv * a - ES::M2d::Identity();

  double trace_M = M.trace();
  ES::M2d I = ES::M2d::Identity();
  ES::M2d S = (lameAlpha * trace_M * I + 2.0 * lameBeta * M) * a_bar_inv;

  (ES::Mp<ES::M2d>(da_)) = S * h;
}

void ElasticModel2DFundamentalFormsSTVK::compute_d2psi_da2(const double *param, const double a_[4], const double abar_[4], double da2_[16]) const
{
  double Es = param[0];
  double nu_s = param[1];
  double h = param[4];

  double lameAlpha = Es * nu_s / ((1 + nu_s) * (1 - 2 * nu_s));
  double lameBeta = Es / (2 * (1 + nu_s));

  ES::M2d a = ES::Mp<const ES::M2d>(a_);
  ES::M2d abar = ES::Mp<const ES::M2d>(abar_);
  ES::M2d a_bar_inv = abar.fullPivHouseholderQr().inverse();
  ES::M2d M = a_bar_inv * a - ES::M2d::Identity();

  // E = 0.5 alpha (trace(M))^2 + beta trace(M^2)
  // dE/dM = alpha trace(M) I + 2 beta M
  // d2E/dM2 = alpha I x I + 2 beta I

  // dE/da = dE/dM : dM/da
  // d2E/da2 = (dM/da)^T : d2E/dM2 : dM/da + dE/dM : d2M/da2

  ES::M2d I;
  I.setIdentity();

  ES::M4d d2EdM2;
  d2EdM2 << (lameAlpha + 2 * lameBeta), 0, 0, lameAlpha,
    0, 0, 2 * lameBeta, 0,
    0, 2 * lameBeta, 0, 0,
    lameAlpha, 0, 0, lameAlpha + 2 * lameBeta;

  ES::M4d dMda = KroneckerProduct2(I, a_bar_inv);
  ES::M4d d2Eda2 = dMda.transpose() * d2EdM2 * dMda;

  (ES::Mp<ES::M4d>(da2_)) = d2Eda2 * h;
}

double ElasticModel2DFundamentalFormsSTVK::compute_psi_b(const double *param, const double b_[4], const double abar_[4], const double bbar_[4]) const
{
  double Es = param[2];
  double nu_s = param[3];
  double h = param[4];

  double lameAlpha = Es * nu_s / ((1 + nu_s) * (1 - 2 * nu_s));
  double lameBeta = Es / (2 * (1 + nu_s));

  ES::M2d b = ES::Mp<const ES::M2d>(b_);
  ES::M2d abar = ES::Mp<const ES::M2d>(abar_);
  ES::M2d bbar = ES::Mp<const ES::M2d>(bbar_);
  ES::M2d a_bar_inv = abar.fullPivHouseholderQr().inverse();
  ES::M2d M = a_bar_inv * (b - bbar);

  double trace_M = M.trace();
  double Wc = 0.5 * lameAlpha * (trace_M * trace_M) + lameBeta * (M * M).trace();

  return Wc * h * h * h / 12;
}

void ElasticModel2DFundamentalFormsSTVK::compute_dpsi_db(const double *param, const double b_[4], const double abar_[4], const double bbar_[4], double db_[4]) const
{
  double Es = param[2];
  double nu_s = param[3];
  double h = param[4];

  double lameAlpha = Es * nu_s / ((1 + nu_s) * (1 - 2 * nu_s));
  double lameBeta = Es / (2 * (1 + nu_s));

  ES::M2d b = ES::Mp<const ES::M2d>(b_);
  ES::M2d abar = ES::Mp<const ES::M2d>(abar_);
  ES::M2d bbar = ES::Mp<const ES::M2d>(bbar_);
  ES::M2d a_bar_inv = abar.fullPivHouseholderQr().inverse();
  ES::M2d M = a_bar_inv * (b - bbar);

  double trace_M = M.trace();
  ES::M2d S = lameAlpha * trace_M * a_bar_inv + 2.0 * lameBeta * M * a_bar_inv;
  (ES::Mp<ES::M2d>(db_)) = S * h * h * h / 12;
}

void ElasticModel2DFundamentalFormsSTVK::compute_d2psi_db2(const double *param, const double b_[4], const double abar_[4], const double bbar_[4], double db2_[16]) const
{
  double Es = param[2];
  double nu_s = param[3];
  double h = param[4];

  double lameAlpha = Es * nu_s / ((1 + nu_s) * (1 - 2 * nu_s));
  double lameBeta = Es / (2 * (1 + nu_s));

  ES::M2d b = ES::Mp<const ES::M2d>(b_);
  ES::M2d abar = ES::Mp<const ES::M2d>(abar_);
  ES::M2d bbar = ES::Mp<const ES::M2d>(bbar_);
  ES::M2d a_bar_inv = abar.fullPivHouseholderQr().inverse();
  ES::M2d M = a_bar_inv * (b - bbar);

  ES::M2d I;
  I.setIdentity();

  ES::M4d d2EdM2;
  d2EdM2 << (lameAlpha + 2 * lameBeta), 0, 0, lameAlpha,
    0, 0, 2 * lameBeta, 0,
    0, 2 * lameBeta, 0, 0,
    lameAlpha, 0, 0, lameAlpha + 2 * lameBeta;

  ES::M4d dMdb = KroneckerProduct2(I, a_bar_inv);
  ES::M4d d2Edb2 = dMdb.transpose() * d2EdM2 * dMdb;

  (ES::Mp<ES::M4d>(db2_)) = d2Edb2 * h * h * h / 12;
}

void ElasticModel2DFundamentalFormsSTVK::compute_d2psi_dadabar(const double *param, const double a_[4], const double abar_[4], double dadabar_[16]) const
{
  double Es = param[0];
  double nu_s = param[1];
  double h = param[4];

  double lameAlpha = Es * nu_s / ((1 + nu_s) * (1 - 2 * nu_s));
  double lameBeta = Es / (2 * (1 + nu_s));

  // vec (ABX) = (X^T ⊗ A) vec(B)

  // d(abar^{-1]}) abar + abar^{-1} d(abar) = 0
  // d(abar^{-1)}) = - abar^{-1} d(abar) abar^{-1}

  // dM/dabar = d(abar^{-1} a - I) /dabar = d(z a - I) / dz * dz/dabar
  //          = (a^T ⊗ I) * d(abar^{-1})/dabar
  //          = - (a^T ⊗ I) * (abar^{-T} ⊗ abar^{-1})
  //

  // ES::M2d S = (lameAlpha * trace_M * I + 2.0 * lameBeta * M) * a_bar_inv;
  // dS/dabar = (c1 d trace(M)/dM * dM/dabar * I + c2 dM/dabar) * a_bar_inv + S * d(a_bar_inv)/dabar
  // dS/dabar = (c1 * I : dM/dabar + c2 dM/dabar) * a_bar_inv + S * d(a_bar_inv)/dabar

  ES::M2d a = ES::Mp<const ES::M2d>(a_);
  ES::M2d abar = ES::Mp<const ES::M2d>(abar_);
  ES::M2d a_bar_inv = abar.fullPivHouseholderQr().inverse();
  ES::M2d M = a_bar_inv * a - ES::M2d::Identity();
  ES::M2d I = ES::M2d::Identity();

  ES::M4d dadabar;
  dadabar.setZero();
  for (int i = 0; i < 4; i++) {
    ES::M2d dabar_i;
    dabar_i.setZero();
    dabar_i.data()[i] = 1.0;

    ES::M2d dabar_inv_i = -a_bar_inv * dabar_i * a_bar_inv;
    ES::M2d dM_dabar_i = dabar_inv_i * a;

    ES::M2d dS_dabar_i = (lameAlpha * dM_dabar_i.trace() * I + 2.0 * lameBeta * dM_dabar_i) * a_bar_inv +
      (lameAlpha * M.trace() * I + 2.0 * lameBeta * M) * dabar_inv_i;
    dadabar.col(i) = vecCM(dS_dabar_i);
  }

  (ES::Mp<ES::M4d>(dadabar_)) = dadabar * h;
}

void ElasticModel2DFundamentalFormsSTVK::compute_d2psi_db_dabar(const double *param, const double b_[4], const double abar_[4], const double bbar_[4], double dbdabar_[16]) const
{
  double Es = param[2];
  double nu_s = param[3];
  double h = param[4];

  double lameAlpha = Es * nu_s / ((1 + nu_s) * (1 - 2 * nu_s));
  double lameBeta = Es / (2 * (1 + nu_s));

  ES::M2d b = ES::Mp<const ES::M2d>(b_);
  ES::M2d abar = ES::Mp<const ES::M2d>(abar_);
  ES::M2d bbar = ES::Mp<const ES::M2d>(bbar_);
  ES::M2d a_bar_inv = abar.fullPivHouseholderQr().inverse();
  ES::M2d a_bar_inv_T = a_bar_inv.transpose();
  ES::M2d M = a_bar_inv * (b - bbar);
  ES::M2d I = ES::M2d::Identity();

  ES::M4d dbdabar;
  dbdabar.setZero();

  for (int i = 0; i < 4; i++) {
    ES::M2d dabar_i;
    dabar_i.setZero();
    dabar_i.data()[i] = 1.0;

    ES::M2d dabar_inv_i = -a_bar_inv * dabar_i * a_bar_inv;
    ES::M2d dM_dabar_i = dabar_inv_i * (b - bbar);

    ES::M2d dS_dabar_i = (lameAlpha * dM_dabar_i.trace() * I + 2.0 * lameBeta * dM_dabar_i) * a_bar_inv +
      (lameAlpha * M.trace() * I + 2.0 * lameBeta * M) * dabar_inv_i;
    dbdabar.col(i) = vecCM(dS_dabar_i);
  }

  (ES::Mp<ES::M4d>(dbdabar_)) = dbdabar * h * h * h / 12;
}

void ElasticModel2DFundamentalFormsSTVK::compute_d2psi_db_dbbar(const double *param, const double b_[4], const double abar_[4], const double bbar_[4], double dbdbbar_[16]) const
{
  double Es = param[2];
  double nu_s = param[3];
  double h = param[4];

  double lameAlpha = Es * nu_s / ((1 + nu_s) * (1 - 2 * nu_s));
  double lameBeta = Es / (2 * (1 + nu_s));

  ES::M2d b = ES::Mp<const ES::M2d>(b_);
  ES::M2d abar = ES::Mp<const ES::M2d>(abar_);
  ES::M2d bbar = ES::Mp<const ES::M2d>(bbar_);
  ES::M2d a_bar_inv = abar.fullPivHouseholderQr().inverse();
  ES::M2d a_bar_inv_T = a_bar_inv.transpose();
  ES::M2d M = a_bar_inv * (b - bbar);
  ES::M2d I = ES::M2d::Identity();

  ES::M4d dbdbbar;
  dbdbbar.setZero();
  for (int i = 0; i < 4; i++) {
    ES::M2d dbbar_i;
    dbbar_i.setZero();
    dbbar_i.data()[i] = 1.0;

    ES::M2d dM_dbbar_i = -a_bar_inv * dbbar_i;
    ES::M2d dS_dbbar_i = (lameAlpha * dM_dbbar_i.trace() * I + 2.0 * lameBeta * dM_dbbar_i) * a_bar_inv;

    dbdbbar.col(i) = vecCM(dS_dbbar_i);
  }

  (ES::Mp<ES::M4d>(dbdbbar_)) = dbdbbar * h * h * h / 12;
  // this route is not executed yet
}

void ElasticModel2DFundamentalFormsSTVK::compute_d2psi_da_dparam(const double *param, const double a_[4], const double abar_[4], double d2psi_dadparam[/*4 x numParams*/]) const
{
  double Es = param[0];
  double nu_s = param[1];
  double h = param[4];

  double lameAlpha = Es * nu_s / ((1 + nu_s) * (1 - 2 * nu_s));
  double lameBeta = Es / (2 * (1 + nu_s));

  ES::M2d a = ES::Mp<const ES::M2d>(a_);
  ES::M2d abar = ES::Mp<const ES::M2d>(abar_);
  ES::M2d a_bar_inv = abar.fullPivHouseholderQr().inverse();
  ES::M2d M = a_bar_inv * a - ES::M2d::Identity();

  double trace_M = M.trace();
  ES::M2d I = ES::M2d::Identity();
  ES::M2d S = (lameAlpha * trace_M * I + 2.0 * lameBeta * M) * a_bar_inv;

  Eigen::Map<ES::MXd> d2psi_dadparam_map(d2psi_dadparam, 4, 5);
  d2psi_dadparam_map.setZero();

  ES::M2d dS_dalpha = (trace_M * I) * a_bar_inv;
  ES::M2d dS_dbeta = (2.0 * M) * a_bar_inv;

  ES::M2d dS_dE = dS_dalpha * (nu_s / ((1 + nu_s) * (1 - 2 * nu_s))) + dS_dbeta * (1.0 / (2 * (1 + nu_s)));
  ES::M2d dS_dnu = dS_dalpha * (Es * (1 + 2 * nu_s * nu_s) / ((1 + nu_s) * (1 + nu_s) * (1 - 2 * nu_s) * (1 - 2 * nu_s)));
  dS_dnu += dS_dbeta * (-Es / (2 * (1 + nu_s) * (1 + nu_s)));

  d2psi_dadparam_map.col(0) = vecCM(dS_dE) * h;
  d2psi_dadparam_map.col(1) = vecCM(dS_dnu) * h;
  d2psi_dadparam_map.col(2).setZero();
  d2psi_dadparam_map.col(3).setZero();
  d2psi_dadparam_map.col(4) = vecCM(S);
}

void ElasticModel2DFundamentalFormsSTVK::compute_d2psi_db_dparam(const double *param, const double b_[4], const double abar_[4], const double bbar_[4], double d2psi_dbdparam[/*4 x numParams*/]) const
{
  double Es = param[2];
  double nu_s = param[3];
  double h = param[4];

  double lameAlpha = Es * nu_s / ((1 + nu_s) * (1 - 2 * nu_s));
  double lameBeta = Es / (2 * (1 + nu_s));

  ES::M2d b = ES::Mp<const ES::M2d>(b_);
  ES::M2d abar = ES::Mp<const ES::M2d>(abar_);
  ES::M2d bbar = ES::Mp<const ES::M2d>(bbar_);
  ES::M2d a_bar_inv = abar.fullPivHouseholderQr().inverse();
  ES::M2d M = a_bar_inv * (b - bbar);

  double trace_M = M.trace();
  ES::M2d I = ES::M2d::Identity();
  ES::M2d S = (lameAlpha * trace_M * I + 2.0 * lameBeta * M) * a_bar_inv;

  Eigen::Map<ES::MXd> d2psi_dbdparam_map(d2psi_dbdparam, 4, 5);
  d2psi_dbdparam_map.setZero();

  ES::M2d dS_dalpha = (trace_M * I) * a_bar_inv;
  ES::M2d dS_dbeta = (2.0 * M) * a_bar_inv;

  ES::M2d dS_dE = dS_dalpha * (nu_s / ((1 + nu_s) * (1 - 2 * nu_s))) + dS_dbeta * (1.0 / (2 * (1 + nu_s)));
  ES::M2d dS_dnu = dS_dalpha * (Es * (1 + 2 * nu_s * nu_s) / ((1 + nu_s) * (1 + nu_s) * (1 - 2 * nu_s) * (1 - 2 * nu_s)));
  dS_dnu += dS_dbeta * (-Es / (2 * (1 + nu_s) * (1 + nu_s)));

  d2psi_dbdparam_map.col(0).setZero();
  d2psi_dbdparam_map.col(1).setZero();
  d2psi_dbdparam_map.col(2) = vecCM(dS_dE) * h * h * h / 12;
  d2psi_dbdparam_map.col(3) = vecCM(dS_dnu) * h * h * h / 12;
  d2psi_dbdparam_map.col(4) = vecCM(S) * 3 * h * h / 12;
}