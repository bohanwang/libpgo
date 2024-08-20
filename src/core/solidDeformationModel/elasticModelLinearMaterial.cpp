/*
author: Bohan Wang
copyright to USC,MIT
*/

#include "elasticModelLinearMaterial.h"

#include "EigenSupport.h"

using namespace pgo::SolidDeformationModel;
namespace ES = pgo::EigenSupport;

double ElasticModelLinearMaterial::compute_psi(const double *, const double F[9], const double[9], const double[9], const double[3]) const
{
  ES::Mp<const ES::M3d> FMap(F);
  ES::M3d strain = (FMap.transpose() + FMap) * 0.5 - ES::M3d::Identity();

  // E = mu eps : eps + lambda/2 trace^2(eps)
  double t = strain.trace();
  return strain.squaredNorm() * mu + t * t * lambda * 0.5;
}

void ElasticModelLinearMaterial::compute_P(const double *, const double F[9], const double[9], const double[9], const double[3], double P[9]) const
{
  // P = 2 mu eps + lambda trace(eps) I
  ES::Mp<const ES::M3d> FMap(F);
  ES::M3d strain = (FMap.transpose() + FMap) * 0.5 - ES::M3d::Identity();

  (ES::Mp<ES::M3d>(P)) = strain * 2 * mu + ES::M3d::Identity() * lambda * strain.trace();
}

void ElasticModelLinearMaterial::compute_dPdF(const double *, const double F[9], const double[9], const double[9], const double[3], double dPdFOut[81]) const
{
  // P = 2 mu eps + lambda trace(eps) I
  // P = mu (F + F^T) - 2mu I + lambda tr(F - I) I

  ES::Mp<ES::M9d> dPdF(dPdFOut);
  dPdF = ES::M9d::Identity() * mu;

  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      dPdF(i + 3 * j, j + 3 * i) += mu;
      dPdF(4 * i, 4 * j) += lambda;
    }
  }
}
