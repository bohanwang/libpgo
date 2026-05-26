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

  DeformationModelAssembler(std::unique_ptr<const DeformationModelManager> dm, const double *elementFlags = nullptr);
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

  const DeformationModelManager &getDeformationModelManager() const { return *deformationModelManager; }
  const EigenSupport::SpMatD &getHessianTemplate() const { return KTemplate; }
  const EigenSupport::SpMatD &get_dfda_Template() const { return dfdaTemplate; }
  const EigenSupport::SpMatD &get_dfdb_Template() const { return dfdbTemplate; }

protected:
  void getPlasticParameters(int ele, const double *paramsAll, double *param) const;
  void getElasticParameters(int ele, const double *paramsAll, double *param) const;

  std::unique_ptr<const DeformationModelManager> deformationModelManager;
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
};
}  // namespace SolidDeformationModel
}  // namespace pgo
