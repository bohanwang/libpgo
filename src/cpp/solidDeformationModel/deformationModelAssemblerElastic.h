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
class DeformationModelAssemblerElasticImpl;

class DeformationModelAssemblerElastic : public DeformationModelAssembler
{
public:
  DeformationModelAssemblerElastic(const DeformationModelManager *dm,
    const double *elementWeights, int numHandles, const double *elementFlags = nullptr);

  virtual DeformationModelAssemblerCacheData *allocateCache() const override;
  virtual void freeCache(DeformationModelAssemblerCacheData *data) const override;

  virtual double computeEnergy(const double *x_param, DeformationModelAssemblerCacheData *cache) const override;
  virtual void computeGradient(const double *x_param, double *grad, DeformationModelAssemblerCacheData *cache) const override;
  virtual void computeHessian(const double *x_param, HessianMatrixHandle *hess, DeformationModelAssemblerCacheData *cache) const override;
  virtual int getNumDOFs() const override { return n3; }

  void enableFullspaceScale(int enbale) { enableFullspaceScaling = enbale; }

  void computeVonMisesStresses(const double *x_param, DeformationModelAssemblerCacheData *cache, double *elementStresses) const;
  void computeMaxStrains(const double *x_param, DeformationModelAssemblerCacheData *cache, double *elementStrain) const;

protected:
  virtual void getFixedPlasticParameters(int ele, double *param) const override;
  virtual void getFixedElasticParameters(int ele, double *param) const override;
  virtual const double *getFixedDeformation() const override { return nullptr; }

  int numHandles;
  int numPlasticParams = 3;
  int numElasticParams = 0;
  int enableFullspaceScaling = 0;

  DeformationModelAssemblerElasticImpl *eimpl;
};

}  // namespace SolidDeformationModel
}  // namespace pgo
