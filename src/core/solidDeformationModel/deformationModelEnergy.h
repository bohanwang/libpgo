/*
author: Bohan Wang
copyright to USC,MIT
*/

#pragma once

#include "potentialEnergy.h"

#include <vector>

namespace pgo
{
namespace SolidDeformationModel
{
class DeformationModelAssembler;
class DeformationModelAssemblerCacheData;

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

  void computeVonMisesStresses(EigenSupport::ConstRefVecXd x, EigenSupport::RefVecXd elementStresses) const;
  void computeMaxStrains(EigenSupport::ConstRefVecXd x, EigenSupport::RefVecXd elementStrain) const;

protected:
  std::shared_ptr<DeformationModelAssembler> forceModelAssembler;
  DeformationModelAssemblerCacheData *data;

  std::vector<int> allDOFs;
  EigenSupport::VXd restPosition;
};
}  // namespace SolidDeformationModel
}  // namespace pgo
