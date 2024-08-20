/*
author: Bohan Wang
copyright to USC,MIT
*/

#pragma once

namespace pgo
{
namespace SolidDeformationModel
{
class HessianMatrixHandle;
class DeformationModelManager;
class DeformationModelAssemblerImpl;

enum class DeformationModelAssemblingType
{
  AT_ELASTIC,
  AT_ELASTIC_MULTI_PARAM,
  AT_ELASTIC_FULL_PARAM,
  AT_ELASTIC_MULTI_PLASTIC,
  AT_ELASTIC_FULL_PLASTIC,
};

class DeformationModelAssemblerCacheData;

class DeformationModelAssembler
{
public:
  DeformationModelAssembler(const DeformationModelManager *dm, DeformationModelAssemblingType at, const double *elementFlags = nullptr);
  virtual ~DeformationModelAssembler();

  void setFixedParameters(const double *param, int size);
  int getFixedParameters(double *param) const;

  virtual DeformationModelAssemblerCacheData *allocateCache() const = 0;
  virtual void freeCache(DeformationModelAssemblerCacheData *data) const = 0;

  virtual double computeEnergy(const double *x_param, DeformationModelAssemblerCacheData *cache) const = 0;
  virtual void computeGradient(const double *x_param, double *grad, DeformationModelAssemblerCacheData *cache) const = 0;
  virtual void computeHessian(const double *x_param, HessianMatrixHandle *hess, DeformationModelAssemblerCacheData *cache) const = 0;
  virtual void computeHessian(const double *x_param, double *hess, DeformationModelAssemblerCacheData *cache) const;
  virtual void compute3rdOrderTensor(const double *x_param, const double *lambda, HessianMatrixHandle *hess, DeformationModelAssemblerCacheData *cache) const;
  virtual int getNumDOFs() const = 0;

  const DeformationModelManager *getDeformationModelManager() const { return deformationModelManager; }
  const HessianMatrixHandle *getHessianTemplate() const;

protected:
  virtual void getFixedPlasticParameters(int ele, double *param) const = 0;
  virtual void getFixedElasticParameters(int ele, double *param) const = 0;
  virtual const double *getFixedDeformation() const = 0;

  const DeformationModelManager *deformationModelManager;
  DeformationModelAssemblingType assemblingType;
  int n3, nele, nvtx, neleVtx;

  DeformationModelAssemblerImpl *impl;

  const int enableSanityCheck = 1;
};
}  // namespace SolidDeformationModel
}  // namespace pgo