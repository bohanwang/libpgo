#pragma once

#include "plasticModel.h"
#include "EigenDef.h"

namespace pgo
{
namespace SolidDeformationModel
{
class PlasticModel2DFundamentalForms : public PlasticModel
{
public:
  PlasticModel2DFundamentalForms() {}
  virtual ~PlasticModel2DFundamentalForms() {}

  void set_abar(const EigenSupport::M2d &abar_) { abar = abar_; }
  void set_bbar(const EigenSupport::M2d &bbar_) { bbar = bbar_; }
  void set_tbar(const EigenSupport::M3d &tbar_) { tbar = tbar_; }
  void set_qbar(const EigenSupport::M3d &qbar_) { qbar = qbar_; }
  void setArea(double a) { areaRest = a; }

  virtual int getNumParameters() const { return 0; }

  virtual void compute_abar(const double *params, double *a) const { (Eigen::Map<EigenSupport::M2d>(a)) = abar; }
  virtual void compute_bbar(const double *params, double *b) const { (Eigen::Map<EigenSupport::M2d>(b)) = bbar; }

  virtual void compute_tbar(const double *params, double *t) const { (Eigen::Map<EigenSupport::M3d>(t)) = tbar; }
  virtual void compute_qbar(const double *params, double *q) const { (Eigen::Map<EigenSupport::M3d>(q)) = qbar; }

  virtual double computeArea(const double *params) const { return areaRest; }

  virtual void compute_dtbar_inv_dparam(const double *params, int j, double *dtbar_da) const {}
  virtual void compute_dqbar_dparam(const double *params, int j, double *dqbar_da) const {}
  virtual void compute_dK_dparam(const double *params, double *dK_da) const {}
  virtual void compute_dH_dparam(const double *params, double *dH_da) const {}
  virtual void compute_darea_dparam(const double *params, double *darea_da) const {}

  virtual void compute_dabar_dparam(const double *params, double *dabar_dparam) const {}
  virtual void compute_dbbar_dparam(const double *params, double *dbbar_dparam) const {}

protected:
  EigenSupport::M2d abar = EigenSupport::M2d::Identity(), bbar = EigenSupport::M2d::Identity();
  EigenSupport::M3d tbar = EigenSupport::M3d::Identity(), qbar = EigenSupport::M3d::Identity();
  double areaRest = 1.0;
};
}  // namespace SolidDeformationModel
}  // namespace pgo