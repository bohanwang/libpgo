/*
author: Bohan Wang
copyright to USC,MIT,NUS
*/

#pragma once

#include "potentialEnergy.h"

#include <cstdint>
#include <memory>
#include <vector>

namespace pgo
{
namespace SolidDeformationModel
{
class DeformationModelAssembler;

class DeformationModelEnergy : public NonlinearOptimization::PotentialEnergy
{
public:
  DeformationModelEnergy(std::unique_ptr<DeformationModelAssembler> fma, const EigenSupport::VXd *restPosition = nullptr, int offset = 0);
  virtual ~DeformationModelEnergy();

  // Borrow the owned assembler (e.g. for stress queries).
  const DeformationModelAssembler &assembler() const { return *forceModelAssembler; }

  virtual double func(EigenSupport::ConstRefVecXd x) const override;
  virtual void gradient(EigenSupport::ConstRefVecXd x, EigenSupport::RefVecXd grad) const override;
  virtual void hessian(EigenSupport::ConstRefVecXd x, EigenSupport::SpMatD &hess) const override;
  virtual void createHessian(EigenSupport::SpMatD &hess) const override;
  virtual void getDOFs(std::vector<int> &dofs) const override { dofs = this->allDOFs; }
  virtual int getNumDOFs() const override { return (int)allDOFs.size(); }

  virtual NonlinearOptimization::MaxStepResult computeMaxStepLimit(EigenSupport::ConstRefVecXd x, EigenSupport::ConstRefVecXd dx) const override;

  void setElasticParams(const EigenSupport::ConstRefVecXd elasticParams) { this->elasticParams = elasticParams; }
  void setPlasticParams(const EigenSupport::ConstRefVecXd plasticParams) { this->plasticParams = plasticParams; }
  void setEnableMaterialMaxStep(bool enable) { enableMaterialMaxStep_ = enable; }
  bool isMaterialMaxStepEnabled() const { return enableMaterialMaxStep_; }
protected:
  std::unique_ptr<DeformationModelAssembler> forceModelAssembler;

  std::vector<int> allDOFs;
  EigenSupport::VXd restPosition;
  EigenSupport::VXd elasticParams;
  EigenSupport::VXd plasticParams;
  bool enableMaterialMaxStep_ = true;
};
}  // namespace SolidDeformationModel
}  // namespace pgo
