#pragma once

#include "deformationModel.h"

namespace pgo
{
namespace SolidDeformationModel
{
class CubicMeshDeformationModelInternal;

class CubicMeshDeformationModel : public DeformationModel
{
public:
  CubicMeshDeformationModel(const double restPositions[24], ElasticModel *elasticModel, PlasticModel *plasticModel);
  virtual ~CubicMeshDeformationModel();

  using DeformationModel::CacheData;

  virtual CacheData *allocateCacheData() const override;
  virtual void freeCacheData(CacheData *data) const override;
  virtual void prepareData(const double *x, const double *param, const double *materialParam, CacheData *cacheData) const override;

  void enableSPD(int enable) override;

  void vonMisesStress(const CacheData *cacheDataBase, int &nPt, double *stresses) const override;
  void maxStrain(const CacheData *cacheDataBase, int &nPt, double *stresses) const override;

  virtual double computeEnergy(const CacheData *cacheData) const override;
  virtual void compute_dE_dx(const CacheData *cacheData, double *grad) const override;
  virtual void compute_d2E_dx2(const CacheData *cacheData, double *hess) const override;

  void compute_dE_da(const CacheData *cacheData, double *grad) const;
  void compute_d2E_da2(const CacheData *cacheData, double *hess) const;
  virtual void compute_d2E_dxda(const CacheData *cacheData, double *hess) const override;

  void compute_dE_db(const CacheData *cacheData, double *grad) const;
  void compute_d2E_db2(const CacheData *cacheData, double *hess) const;
  virtual void compute_d2E_dxdb(const CacheData *cacheData, double *hess) const override;

  void compute_d2E_dadb(const CacheData *cacheData, double *hess) const;

  virtual int getNumVertices() const override { return 8; }
  virtual int getNumDOFs() const override { return 24; }
  virtual int getNumMaterialLocations() const override { return 8; }
  virtual LocalMaxStepResult computeLocalMaxStepSize(const double *x_local, const double *dx_local) const override;

  void computeF(const double *x, int materialLocationID, double F[9]) const;
  void computeFe(const CacheData *cacheData, int materialLocationID, double F[9]) const;
  void computeP(const CacheData *cacheData, int materialLocationID, double P[9]) const;
  void computedPdF(const CacheData *cacheData, int materialLocationID, double dPdF[81]) const;
  void computedFdx(const CacheData *cacheData, int materialLocationID, double *dFdx) const;
  void computeForceFromP(const CacheData *cacheData, int materialLocationID, const double P[9], double f[24]) const;

  double getWeightDetJ(int materialLocationID) const;

protected:
  CubicMeshDeformationModelInternal *ind;
  int numPlasticParams_ = 0;
  int numElasticParams_ = 0;
};
}  // namespace SolidDeformationModel
}  // namespace pgo
