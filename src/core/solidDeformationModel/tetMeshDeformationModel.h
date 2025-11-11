/*
author: Bohan Wang
copyright to USC,MIT,NUS
*/

#pragma once

#include "deformationModel.h"

namespace pgo
{
namespace SolidDeformationModel
{
class TetMeshDeformationModelInternal;

class TetMeshDeformationModel : public DeformationModel
{
public:
  TetMeshDeformationModel(const double X0[3], const double X1[3], const double X2[3], const double X3[3],
    ElasticModel *elasticModel, PlasticModel *plasticModel);
    
  virtual ~TetMeshDeformationModel();

  using DeformationModel::CacheData;

  virtual CacheData *allocateCacheData() const override;
  virtual void freeCacheData(CacheData *data) const override;
  virtual void prepareData(const double *x, const double *param, const double *materialParam, CacheData *cacheData) const override;

  void vonMisesStress(const CacheData *cacheDataBase, int &nPt, double *stresses) const override;
  virtual void maxStrain(const CacheData *cacheDataBase, int &nPt, double *stresses) const override;

  virtual double computeEnergy(const CacheData *cacheData) const override;
  virtual void compute_dE_dx(const CacheData *cacheData, double *grad) const override;
  virtual void compute_d2E_dx2(const CacheData *cacheData, double *hess) const override;

  virtual void compute_dE_da(const CacheData *cacheData, double *grad) const override;
  virtual void compute_d2E_da2(const CacheData *cacheData, double *hess) const override;
  virtual void compute_d2E_dxda(const CacheData *cacheData, double *hess) const override;

  virtual void compute_dE_db(const CacheData *cacheData, double *grad) const override;
  virtual void compute_d2E_db2(const CacheData *cacheData, double *hess) const override;
  virtual void compute_d2E_dxdb(const CacheData *cacheData, double *hess) const override;

  virtual void compute_d2E_dadb(const CacheData *cacheData, double *hess) const override;

  // virtual void compute_d3E_dx3(const CacheData *cacheData, double *tensor) const override;
  // virtual void compute_d3E_dxdadx(const CacheData *cacheData, double *tensor) const override;
  // virtual void compute_d3E_dxdada(const CacheData *cacheData, double *tensor) const override;

  virtual int getNumVertices() const override { return 4; }
  virtual int getNumDOFs() const override { return 12; }

  // inline static double d3E_dx3_ijk(const double *tensor, int i, int j, int k) { return DeformationModel::d3E_dx3_ijk(tensor, i, j, k, 12); }
  // inline static double &d3E_dx3_ijk(double *tensor, int i, int j, int k) { return DeformationModel::d3E_dx3_ijk(tensor, i, j, k, 12); }

  void computeF(const double *x, int materialLocationIDs, double F[9]) const;
  void computeP(const CacheData *cacheDataBase, int materialLocationIDs, double P[9]) const;
  void computedPdF(const CacheData *cacheDataBase, int materialLocationIDs, double dPdF[81]) const;
  void computedFdx(const CacheData *cacheDataBase, int materialLocationIDs, double *dFdx) const;
  void computeForceFromP(const CacheData *cacheDataBase, const double P[9], double f[12]) const;

  static void computeDm(const double X[12], double Dm[9]);
  static void computeDs(const double x[12], double Ds[9]);
  static double computeVolume(const double X0[3], const double X1[3], const double X2[3], const double X3[3]);
  static void compute_dF_dx(const double DmInv[9], double dFdx[9 * 12]);

protected:
  TetMeshDeformationModelInternal *ind;
};
}  // namespace SolidDeformationModel
}  // namespace pgo
