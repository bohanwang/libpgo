#pragma once

#include "elasticModel3DDeformationGradient.h"
#include "EigenDef.h"

namespace pgo
{
namespace SolidDeformationModel
{
class ElasticModel3DSTVKMaterial : public ElasticModel3DDeformationGradient
{
public:
  ElasticModel3DSTVKMaterial(double _mu, double _lambda);
  virtual ~ElasticModel3DSTVKMaterial();

  int getNumParameters() const override { return 0; }
  void setMaterialParameters(double _mu, double _lambda)   {     mu = _mu;     lambda = _lambda;   }

  virtual double compute_psi(const double *param, const double _F[9], const double _U[], const double _V[], const double _S[]) const override;
  virtual void compute_P(const double *param, const double _F[9], const double _U[], const double _V[], const double _S[], double P[9]) const override;         // d||P||/dP  d||Q||/dQ
  virtual void compute_dPdF(const double *param, const double _F[9], const double _U[], const double _V[], const double _S[], double dPdF[81]) const override;  // d2||P||/dP2  d2||Q||/dQ2

protected:
  EigenSupport::V3d computeLowerInvariance(const EigenSupport::M3d &F, EigenSupport::M3d &R, 
    EigenSupport::V3d *S = nullptr, EigenSupport::M3d *U = nullptr, EigenSupport::M3d *V = nullptr) const;

  double mu, lambda;
  EigenSupport::M3d C[3];
};

}  // namespace SolidDeformationModel
}  // namespace pgo