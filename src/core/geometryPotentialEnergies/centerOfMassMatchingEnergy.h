/*
author: Minghao Guo
*/
#pragma once

#include "potentialEnergy.h"
#include "tetMeshGeo.h"

#include <array>

namespace pgo
{
namespace PredefinedPotentialEnergies
{
class CenterOfMassMatchingEnergy : public NonlinearOptimization::PotentialEnergy
{
public:
  CenterOfMassMatchingEnergy(const Mesh::TetMeshGeo &tetMesh, const EigenSupport::V3d &tgtCoM, const EigenSupport::M3d &projMat);
  virtual ~CenterOfMassMatchingEnergy();

  virtual double func(EigenSupport::ConstRefVecXd x) const override;
  virtual void gradient(EigenSupport::ConstRefVecXd x, EigenSupport::RefVecXd grad) const override;
  virtual void hessian(EigenSupport::ConstRefVecXd x, EigenSupport::SpMatD &hess) const override {}

  void compute_com_and_energy_and_grad(EigenSupport::ConstRefVecXd x, EigenSupport::V3d &com, double &energy, EigenSupport::RefVecXd grad) const;

  virtual void createHessian(EigenSupport::SpMatD &hess) const { hess = hessTemplate; }
  virtual void getDOFs(std::vector<int> &dofs) const { dofs = allDOFs; }
  virtual int getNumDOFs() const { return (int)allDOFs.size(); }
  virtual int hasHessian() const { return 0; }

  const EigenSupport::V3d &getTgtCoM() const { return m_tgtCoM; }
  const EigenSupport::M3d &getProjMat() const { return m_projMat; }

protected:
  EigenSupport::MXi m_tet;
  EigenSupport::V3d m_tgtCoM;
  EigenSupport::M3d m_projMat;

  EigenSupport::SpMatD hessTemplate;
  std::vector<int> allDOFs;
};

}  // namespace PredefinedPotentialEnergies
}  // namespace pgo
