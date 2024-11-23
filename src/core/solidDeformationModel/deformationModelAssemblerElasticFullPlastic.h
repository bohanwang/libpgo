/*
author: Bohan Wang
copyright to USC,MIT
*/

#pragma once

#include "deformationModelAssembler.h"

namespace pgo
{
namespace SolidDeformationModel
{
class DeformationModelAssemblerElasticFullPlasticImpl;

class DeformationModelAssemblerElasticFullPlastic : public DeformationModelAssembler
{
public:
  DeformationModelAssemblerElasticFullPlastic(const DeformationModelManager *dm,
    const double *elementWeights, int numHandles, const double *elementFlags = nullptr);

  virtual DeformationModelAssemblerCacheData *allocateCache() const override;
  virtual void freeCache(DeformationModelAssemblerCacheData *data) const override;

  virtual double computeEnergy(const double *x_param, DeformationModelAssemblerCacheData *cache) const override;
  virtual void computeGradient(const double *x_param, double *grad, DeformationModelAssemblerCacheData *cache) const override;
  virtual void computeHessian(const double *x_param, HessianMatrixHandle *hess, DeformationModelAssemblerCacheData *cache) const override;

  virtual int getNumDOFs() const override { return n3 + nele * numPlasticParams; }

  void enableFullspaceScale(int enbale) { enableFullspaceScaling = enbale; }

protected:
  virtual void getFixedPlasticParameters(int ele, double *param) const override {}
  virtual void getFixedElasticParameters(int ele, double *param) const override;
  virtual const double *getFixedDeformation() const override { return nullptr; }

  void getPlasticParameters(const double *x_param, int ele, double *params) const;
  void getElasticParameters(const double *x_param, int ele, double *params) const;

  int numHandles;
  int numPlasticParams = 3;
  int numElasticParams = 0;
  int enableFullspaceScaling = 0;

  DeformationModelAssemblerElasticFullPlasticImpl *epimpl;
};
}  // namespace SolidDeformationModel
}  // namespace pgo
