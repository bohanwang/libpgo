/*
author: Bohan Wang
copyright to USC,MIT
*/

#pragma once

namespace pgo
{
namespace SolidDeformationModel
{
class ElasticModel;
class PlasticModel;

class DeformationModelCacheData
{
public:
  DeformationModelCacheData() {}
  virtual ~DeformationModelCacheData() {}
};

class DeformationModel
{
public:
  DeformationModel(const ElasticModel *elasticModel, const PlasticModel *plasticModel):
    em(elasticModel), pm(plasticModel) {}
  virtual ~DeformationModel() {}

  typedef DeformationModelCacheData CacheData;

  virtual CacheData *allocateCacheData() const = 0;
  virtual void freeCacheData(CacheData *data) const = 0;
  virtual void prepareData(const double *x, const double *param, const double *materialParam, CacheData *cacheData) const = 0;

  virtual void vonMisesStress(const CacheData *cacheDataBase, int &nPt, double *stresses) const {};
  virtual void maxStrain(const CacheData *cacheDataBase, int &nPt, double *stresses) const {};

  virtual double computeEnergy(const CacheData *cacheData) const = 0;
  virtual void compute_dE_dx(const CacheData *cacheData, double *grad) const = 0;
  virtual void compute_d2E_dx2(const CacheData *cacheData, double *hess) const = 0;

  virtual void compute_dE_da(const CacheData *cacheData, double *grad) const = 0;
  virtual void compute_d2E_da2(const CacheData *cacheData, double *hess) const = 0;
  virtual void compute_d2E_dxda(const CacheData *cacheData, double *hess) const = 0;

  virtual void compute_dE_db(const CacheData *cacheData, double *grad) const = 0;
  virtual void compute_d2E_db2(const CacheData *cacheData, double *hess) const = 0;
  virtual void compute_d2E_dxdb(const CacheData *cacheData, double *hess) const = 0;

  virtual void compute_d2E_dadb(const CacheData *cacheData, double *hess) const = 0;

  virtual void compute_d3E_dx3(const CacheData *cacheData, double *tensor) const = 0;
  virtual void compute_d3E_dxdadx(const CacheData *cacheData, double *tensor) const = 0;
  virtual void compute_d3E_dxdada(const CacheData *cacheData, double *tensor) const = 0;

  inline static double d3E_dx3_ijk(const double *tensor, int i, int j, int k, int dim) { return tensor[k * dim * dim + j * dim + i]; }
  inline static double &d3E_dx3_ijk(double *tensor, int i, int j, int k, int dim) { return tensor[k * dim * dim + j * dim + i]; }

  const PlasticModel *getPlasticModel() const { return pm; }
  const ElasticModel *getElasticModel() const { return em; }

  virtual int getNumVertices() const = 0;
  virtual int getNumDOFs() const = 0;

  // advanced routines
  virtual int getNumMaterialLocations() const { return 1; }
  virtual void computeF(const double *x, int materialLocationIDs, double F[9]) const {}
  virtual void computeP(const CacheData *cacheDataBase, int materialLocationIDs, double P[9]) const {}
  virtual void computedPdF(const CacheData *cacheDataBase, int materialLocationIDs, double dPdF[81]) const {}
  virtual void computedFdx(const CacheData *cacheDataBase, int materialLocationIDs, double *dFdx) const {}

protected:
  int numMaterialLocations = 1;

  const ElasticModel *em;
  const PlasticModel *pm;
};
}  // namespace SolidDeformationModel
}  // namespace pgo
