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
class ElasticModel2DFundamentalForms : public ElasticModel
{
public:
  ElasticModel2DFundamentalForms() {}
  virtual ~ElasticModel2DFundamentalForms() {}

  virtual double compute_psi_a(const double *param, const double a[4], const double abar[4]) const = 0;
  virtual double compute_psi_b(const double *param, const double b[4], const double abar[4], const double bbar[4]) const = 0;

  virtual void compute_dpsi_da(const double *param, const double a[4], const double abar[4], double da[4]) const = 0;
  virtual void compute_dpsi_db(const double *param, const double b[4], const double abar[4], const double bbar[4], double db[4]) const = 0;

  virtual void compute_d2psi_da2(const double *param, const double a[4], const double abar[4], double da2[16]) const = 0;
  virtual void compute_d2psi_db2(const double *param, const double b[4], const double abar[4], const double bbar[4], double db2[16]) const = 0;

  virtual void compute_d2psi_dadabar(const double *param, const double a[4], const double abar[4], double dadabar[16]) const {}
  virtual void compute_d2psi_db_dabar(const double *param, const double b[4], const double abar[4], const double bbar[4], double dbdabar[16]) const {}
  virtual void compute_d2psi_db_dbbar(const double *param, const double b[4], const double abar[4], const double bbar[4], double dbdbbar[16]) const {}

  virtual void compute_d2psi_da_dparam(const double *param, const double a[4], const double abar[4], double d2psi_dadparam[/*4 x numParams*/]) const {}
  virtual void compute_d2psi_db_dparam(const double *param, const double b[4], const double abar[4], const double bbar[4], double d2psi_dbdparam[/*4 x numParams*/]) const {}

  int getNumParameters() const override { return 0; };
};

}  // namespace SolidDeformationModel
}  // namespace pgo