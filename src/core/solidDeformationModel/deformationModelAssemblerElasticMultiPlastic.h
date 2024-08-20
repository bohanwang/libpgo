/*
author: Bohan Wang
copyright to USC,MIT
*/

#pragma once

#include "deformationModelAssembler.h"
#include "simulationMesh.h"

namespace pgo
{
namespace SolidDeformationModel
{
class DeformationModelAssemblerElasticMultiPlasticImpl;

class DeformationModelAssemblerElasticMultiPlastic : public DeformationModelAssembler
{
public:
  DeformationModelAssemblerElasticMultiPlastic(const DeformationModelManager *dm,
    const double *elementWeights, int numHandles, const double *elementFlags = nullptr);

  virtual DeformationModelAssemblerCacheData *allocateCache() const override;
  virtual void freeCache(DeformationModelAssemblerCacheData *data) const override;

  virtual double computeEnergy(const double *x_param, DeformationModelAssemblerCacheData *cache) const override;
  virtual void computeGradient(const double *x_param, double *grad, DeformationModelAssemblerCacheData *cache) const override;
  virtual void computeHessian(const double *x_param, HessianMatrixHandle *hess, DeformationModelAssemblerCacheData *cache) const override;

  virtual int getNumDOFs() const override { return n3 + numHandles * numPlasticParams; }

  static void createWeightMatrix(const SimulationMesh *tetMesh, const int *elementIDs, int numElementIDs, double *W);

protected:
  virtual void getFixedPlasticParameters(int ele, double *param) const override {}
  virtual void getFixedElasticParameters(int ele, double *param) const override;
  virtual const double *getFixedDeformation() const override { return nullptr; }

  void getPlasticParameters(const double *x_param, int ele, double *params) const;

  int numHandles;
  int numPlasticParams = 3;
  int numElasticParams = 0;
  int enableFullspaceScaling = 0;

  DeformationModelAssemblerElasticMultiPlasticImpl *epimpl;
};
}  // namespace SolidDeformationModel
}  // namespace pgo
