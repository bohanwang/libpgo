/*
author: Bohan Wang
copyright to USC, MIT
*/

#pragma once

#include "EigenSupport.h"
#include "potentialEnergy.h"

#include <vector>

namespace pgo
{
namespace ConstraintPotentialEnergies
{
class PotentialEnergyAligningMeshConnectivity : public NonlinearOptimization::PotentialEnergy
{
public:
  PotentialEnergyAligningMeshConnectivity(const EigenSupport::SpMatD &hessianBase_):
    hessianBase(hessianBase_)
  {
    allDOFs.resize(hessianBase.rows());
    for (int i = 0; i < (int)allDOFs.size(); i++)
      allDOFs[i] = i;
  }

  void setDOFs(const std::vector<int> &dofs);

  virtual ~PotentialEnergyAligningMeshConnectivity() {}

  virtual void createHessian(EigenSupport::SpMatD &h) const final { h = hessianBase; }
  virtual void getDOFs(std::vector<int> &adofs) const final { adofs = allDOFs; }
  virtual int getNumDOFs() const final { return (int)allDOFs.size(); }

protected:
  const EigenSupport::SpMatD &hessianBase;
  std::vector<int> allDOFs;
};

}  // namespace ConstraintPotentialEnergies
}  // namespace pgo