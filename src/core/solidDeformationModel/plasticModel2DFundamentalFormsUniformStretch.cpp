#include "plasticModel2DFundamentalFormsUniformStretch.h"

#include "EigenSupport.h"

pgo::SolidDeformationModel::PlasticModel2DFundamentalFormsUniformStretch::PlasticModel2DFundamentalFormsUniformStretch()
{
}

void pgo::SolidDeformationModel::PlasticModel2DFundamentalFormsUniformStretch::compute_abar(const double *params, double *a) const
{
  double s = params[0];
  EigenSupport::M2d S;
  S << s, 0, 0, s;

  (EigenSupport::Mp<EigenSupport::M2d>(a)) = S * abar * S;
}

void pgo::SolidDeformationModel::PlasticModel2DFundamentalFormsUniformStretch::compute_bbar(const double *params, double *b) const
{
  (EigenSupport::Mp<EigenSupport::M2d>(b)) = bbar;
}

void pgo::SolidDeformationModel::PlasticModel2DFundamentalFormsUniformStretch::compute_tbar(const double *params, double *t) const
{
  double s = params[0];
  EigenSupport::M2d S;
  S << s, 0, 0, s;

  (EigenSupport::Mp<EigenSupport::M3d>(t)) = tbar;
  (EigenSupport::Mp<EigenSupport::M3d>(t)).leftCols<2>() *= S;
}

double pgo::SolidDeformationModel::PlasticModel2DFundamentalFormsUniformStretch::computeArea(const double *params) const
{
  return areaRest * params[0] * params[0];
}

void pgo::SolidDeformationModel::PlasticModel2DFundamentalFormsUniformStretch::compute_dtbar_inv_dparam(const double *params, int j, double *dtbar_da) const
{
  EigenSupport::M3d S;
  S << 1, 0, 0,
    0, 1, 0,
    0, 0, 0;

  (EigenSupport::Mp<EigenSupport::M3d>(dtbar_da)) = tbar * S;
}

void pgo::SolidDeformationModel::PlasticModel2DFundamentalFormsUniformStretch::compute_dqbar_dparam(const double *params, int j, double *dqbar_da) const
{
  EigenSupport::M3d zero = EigenSupport::M3d::Zero();
  (EigenSupport::Mp<EigenSupport::M3d>(dqbar_da)) = zero;
}

void pgo::SolidDeformationModel::PlasticModel2DFundamentalFormsUniformStretch::compute_dK_dparam(const double *params, double *dK_da) const
{
  // K = det(L) = det(I_inv * II) = det(I_inv) * det(II) = det(S^-1 * Ibar^-1 * S^-1) * det(IIbar) = (1/s^4) * det(Ibar^-1) * det(IIbar)
  // H = trace(L) / 2

  double detAbar = abar.determinant();
  double detBbar = bbar.determinant();
  double s = params[0];

  double dK_ds = -4.0 * (1.0 / (s * s * s * s * s)) / detAbar * detBbar;
  *(dK_da) = dK_ds;
}

void pgo::SolidDeformationModel::PlasticModel2DFundamentalFormsUniformStretch::compute_dH_dparam(const double *params, double *dH_da) const
{
  // K = det(L) = det(I_inv * II) = det(I_inv) * det(II) = det(S^-1 * Ibar^-1 * S^-1) * det(IIbar) = (1/s^4) * det(Ibar^-1) * det(IIbar)
  // H = trace(L) / 2

  double s = params[0];
  EigenSupport::M2d I_inv = abar.inverse();
  EigenSupport::M2d II = bbar;
  EigenSupport::M2d L = I_inv * II;
  double dH_ds = -2.0 / (s * s * s) * L.trace() * 0.5;
  *(dH_da) = dH_ds;
}

void pgo::SolidDeformationModel::PlasticModel2DFundamentalFormsUniformStretch::compute_darea_dparam(const double *params, double *darea_da) const
{
  double s = params[0];
  double dA_ds = 2.0 * areaRest * s;
  *(darea_da) = dA_ds;
}

void pgo::SolidDeformationModel::PlasticModel2DFundamentalFormsUniformStretch::compute_dabar_dparam(const double *params, double *dabar_dparam) const
{
  double s = params[0];
  EigenSupport::M2d S;
  S << s, 0, 0, s;

  EigenSupport::M2d dS_ds;
  dS_ds << 1, 0, 0, 1;

  EigenSupport::M2d dAbar_ds = dS_ds * abar * S + S * abar * dS_ds;

  (EigenSupport::Mp<EigenSupport::M2d>(dabar_dparam)) = dAbar_ds;
}

void pgo::SolidDeformationModel::PlasticModel2DFundamentalFormsUniformStretch::compute_dbbar_dparam(const double *params, double *dbbar_dparam) const
{
  EigenSupport::M2d zero = EigenSupport::M2d::Zero();
  (EigenSupport::Mp<EigenSupport::M2d>(dbbar_dparam)) = zero;
}