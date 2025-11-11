#pragma once

#include "elasticModel3DDeformationGradient.h"
#include "EigenDef.h"

namespace pgo
{
namespace SolidDeformationModel
{
class ElasticModel3DMooneyRivlin : public ElasticModel3DDeformationGradient
{
public:
  ElasticModel3DMooneyRivlin(int N, const double *Cpq, int M, const double *D);
  virtual ~ElasticModel3DMooneyRivlin() {}

  int getNumParameters() const override { return 0; }
  void updateParameters(int N, const double *Cpq, int M, const double *D);

  virtual double compute_psi(const double *param, const double _F[9], const double _U[], const double _V[], const double _S[]) const override;
  virtual void compute_P(const double *param, const double _F[9], const double _U[], const double _V[], const double _S[], double P[9]) const override;
  virtual void compute_dPdF(const double *param, const double _F[9], const double _U[], const double _V[], const double _S[], double dPdF[81]) const override;  

protected:
  int N, M;

  EigenSupport::MXd Cpq;
  EigenSupport::VXd D;
};

}  // namespace SolidDeformationModel
}  // namespace pgo