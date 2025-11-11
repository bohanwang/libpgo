/*
author: Bohan Wang
copyright to USC,MIT,NUS
*/

#pragma once

#include "potentialEnergy.h"

#include <vector>

namespace pgo
{
namespace SolidDeformationModel
{
class DeformationModelAssembler;

class DeformationModelEnergy : public NonlinearOptimization::PotentialEnergy
{
public:
  DeformationModelEnergy(std::shared_ptr<DeformationModelAssembler> fma, const EigenSupport::VXd *restPosition = nullptr, int offset = 0);
  virtual ~DeformationModelEnergy();

  virtual double func(EigenSupport::ConstRefVecXd x) const override;
  virtual void gradient(EigenSupport::ConstRefVecXd x, EigenSupport::RefVecXd grad) const override;
  virtual void hessian(EigenSupport::ConstRefVecXd x, EigenSupport::SpMatD &hess) const override;
  virtual void createHessian(EigenSupport::SpMatD &hess) const override;
  virtual void getDOFs(std::vector<int> &dofs) const override { dofs = this->allDOFs; }
  virtual int getNumDOFs() const override { return (int)allDOFs.size(); }

  void setElasticParams(const EigenSupport::ConstRefVecXd elasticParams) { this->elasticParams = elasticParams; }
  void setPlasticParams(const EigenSupport::ConstRefVecXd plasticParams) { this->plasticParams = plasticParams; }

protected:
  std::shared_ptr<DeformationModelAssembler> forceModelAssembler;

  std::vector<int> allDOFs;
  EigenSupport::VXd restPosition;
  EigenSupport::VXd elasticParams;
  EigenSupport::VXd plasticParams;
};
}  // namespace SolidDeformationModel
}  // namespace pgo
