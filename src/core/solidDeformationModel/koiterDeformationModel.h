#pragma once
/*
author: Bohan Wang
copyright to USC,MIT,NUS
*/

#include "deformationModel.h"

namespace pgo
{
namespace SolidDeformationModel
{
class KoiterDeformationModelInternal;
class KoiterDeformationModel : public SolidDeformationModel::DeformationModel
{
public:
  KoiterDeformationModel(const double X0[3], const double X1[3], const double X2[3],
    const double X3[3], const double X4[3], const double X5[3],
    ElasticModel *elasticModel, PlasticModel *plasticModel, int enforceSPD);

  virtual ~KoiterDeformationModel();

  CacheData *allocateCacheData() const override;
  void freeCacheData(CacheData *data) const override;
  void prepareData(const double *x, const double *param, const double *materialParam, CacheData *cacheData) const override;

  double computeEnergy(const CacheData *cacheData) const override;
  void compute_dE_dx(const CacheData *cacheData, double *grad) const override;
  void compute_d2E_dx2(const CacheData *cacheData, double *hess) const override;

  virtual void compute_dE_da(const CacheData *cacheData, double *grad) const {}
  virtual void compute_d2E_da2(const CacheData *cacheData, double *hess) const {}
  virtual void compute_d2E_dxda(const CacheData *cacheData, double *hess) const override;

  virtual void compute_dE_db(const CacheData *cacheData, double *grad) const {}
  virtual void compute_d2E_db2(const CacheData *cacheData, double *hess) const {}
  virtual void compute_d2E_dxdb(const CacheData *cacheData, double *hess) const override;

  virtual void compute_d2E_dadb(const CacheData *cacheData, double *hess) const {}

  virtual void compute_d3E_dx3(const CacheData *cacheData, double *tensor) const {}
  virtual void compute_d3E_dxdadx(const CacheData *cacheData, double *tensor) const {}
  virtual void compute_d3E_dxdada(const CacheData *cacheData, double *tensor) const {}

  virtual int getNumVertices() const override { return 6; }
  virtual int getNumDOFs() const override { return 18; }

protected:
  KoiterDeformationModelInternal *ind;
};
}  // namespace SolidDeformationModel
}  // namespace pgo
