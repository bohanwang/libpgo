/*
author: Bohan Wang
copyright to USC,MIT,NUS
*/

#pragma once

#include "deformationModelManager.h"
#include "EigenDef.h"

#include <memory>

namespace pgo
{
namespace SolidDeformationModel
{
class DeformationModelAssemblerCacheData;

class DeformationModelAssembler
{
public:
  struct MaterialMaxStepObservation
  {
    double alpha = 1.0;
    bool hasIllegalInitialState = false;
    int limitingElementId = -1;
    int limitingLocationId = -1;
  };

  DeformationModelAssembler(std::shared_ptr<const DeformationModelManager> dm, const double *elementFlags = nullptr);
  virtual ~DeformationModelAssembler();

  double computeEnergy(const double *x, const double *plasticParams, const double *elasticParams) const;
  MaterialMaxStepObservation computeMaxStepObservation(const double *x, const double *dx) const;
  double computeMaxStepSize(const double *x, const double *dx) const;
  void computeGradient(const double *x, const double *plasticParams, const double *elasticParams, double *grad) const;
  void computeHessian(const double *x, const double *plasticParams, const double *elasticParams, EigenSupport::SpMatD &hess) const;

  void compute_df_da(const double *x, const double *plasticParams, const double *elasticParams, EigenSupport::SpMatD &hess) const;
  void compute_df_db(const double *x, const double *plasticParams, const double *elasticParams, EigenSupport::SpMatD &hess) const;

  void computeVonMisesStresses(const double *x, const double *plasticParams, const double *elasticParams, double *elementStresses) const;
  void computeMaxStrains(const double *x, const double *plasticParams, const double *elasticParams, double *elementStrain) const;

  int getNumDOFs() const { return n3; }

  std::shared_ptr<const DeformationModelManager> getDeformationModelManager() const { return deformationModelManager; }
  const EigenSupport::SpMatD &getHessianTemplate() const { return KTemplate; }
  const EigenSupport::SpMatD &get_dfda_Template() const { return dfdaTemplate; }
  const EigenSupport::SpMatD &get_dfdb_Template() const { return dfdbTemplate; }

  // Deterministic assembly (opt-in; default off). The default paths accumulate
  // element contributions with atomic fetch_add and per-thread partial sums, so
  // the summation order depends on thread scheduling and results differ in the
  // last ULPs between otherwise identical runs. When enabled, element
  // contributions are computed in parallel into per-element buffers and then
  // summed per destination in a fixed order, which is bit-reproducible while
  // staying fully parallel.
  void setDeterministic(bool on) const;
  bool isDeterministic() const { return deterministic; }

protected:
  void getPlasticParameters(int ele, const double *paramsAll, double *param) const;
  void getElasticParameters(int ele, const double *paramsAll, double *param) const;

  std::shared_ptr<const DeformationModelManager> deformationModelManager;
  DeformationModelAssemblerCacheData *data;

  int n3, nele, nvtx, neleVtx, localDOFs;
  int numElasticParams = 0;
  int numPlasticParams = 0;

  typedef Eigen::Matrix<std::ptrdiff_t, Eigen::Dynamic, Eigen::Dynamic> DynamicIndexMatrix;

  EigenSupport::VXd restPositions;
  EigenSupport::SpMatD KTemplate, dfdaTemplate, dfdbTemplate;
  std::vector<DynamicIndexMatrix> elementKInverseIndices, element_dfda_InverseIndices, element_dfdb_InverseIndices;

  std::vector<double> elementFlags;
  std::vector<const DeformationModel *> femModels;

  const int enableSanityCheck = 1;

  // deterministic-assembly state (mutable: the compute* methods are const)
  void buildDeterministicMaps() const;
  mutable bool deterministic = false;
  mutable bool detMapsBuilt = false;
  mutable std::vector<double> elemEnergyBuf, elemGradBuf, elemHessBuf;
  // destination -> contribution lists (CSR): gradient by global DOF, Hessian by
  // position in the sparse value array
  mutable std::vector<std::ptrdiff_t> gradScatterStart, gradScatterSrc;
  mutable std::vector<std::ptrdiff_t> hessScatterStart, hessScatterSrc;
};
}  // namespace SolidDeformationModel
}  // namespace pgo
