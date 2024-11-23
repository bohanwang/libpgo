/*
author: Bohan Wang
copyright to USC
*/

#pragma once

#include "potentialEnergy.h"

#include <cstdint>
#include <memory>

namespace pgo
{
namespace NonlinearOptimization
{
class PotentialEnergiesBuffer;

class PotentialEnergies : public PotentialEnergy
{
public:
  PotentialEnergies(int n);
  PotentialEnergies(int n, std::shared_ptr<PotentialEnergiesBuffer> buffer);

  virtual ~PotentialEnergies();

  void setDOFs(const std::vector<int> &dofs);

  void addPotentialEnergy(PotentialEnergy_p energy, double coeff = 1.0)
  {
    potentialEnergies.push_back(energy);
    energyCoeffs.push_back(coeff);
  }
  void init();

  virtual double func(EigenSupport::ConstRefVecXd x) const override;
  virtual void gradient(EigenSupport::ConstRefVecXd x, EigenSupport::RefVecXd grad) const override;
  virtual void hessian(EigenSupport::ConstRefVecXd x, EigenSupport::SpMatD &hess) const override;
  virtual void hessianVector(EigenSupport::ConstRefVecXd x, EigenSupport::ConstRefVecXd vec, EigenSupport::RefVecXd hessVec) const override;
  virtual void createHessian(EigenSupport::SpMatD &hess) const override { hess = hessianAll; }

  virtual void getDOFs(std::vector<int> &dofs) const override { dofs = allDOFs; }
  virtual int getNumDOFs() const override { return nAll; }

  const EigenSupport::SpMatD &getHessianTemplate() const { return hessianAll; }
  int getNNZHessain() const { return (int)hessianAll.nonZeros(); }
  void setEnergyCoeffs(int ith, double val) { energyCoeffs[ith] = val; }

  virtual int isQuadratic() const override { return isQuadraticEnergy; }
  virtual int hasHessianVector() const override { return hasHessianVectorProduct; }

  void printEnergy(EigenSupport::ConstRefVecXd x) const;

protected:
  void mapx(EigenSupport::ConstRefVecXd x, const std::vector<int> &dofs, EigenSupport::RefVecXd xlocal) const;

  std::vector<PotentialEnergy_p> potentialEnergies;
  std::vector<EigenSupport::SpMatI, Eigen::aligned_allocator<EigenSupport::SpMatI>> hessianMatrixMappings;
  std::vector<std::vector<int>> energyDOFs;
  EigenSupport::SpMatD hessianAll;

  std::vector<double> energyCoeffs;
  std::vector<int> allDOFs;
  int nAll;
  int isQuadraticEnergy = 0;
  int hasHessianVectorProduct = 0;

  std::shared_ptr<PotentialEnergiesBuffer> buffer;
};

typedef std::shared_ptr<PotentialEnergies> PotentialEnergies_p;
typedef std::shared_ptr<const PotentialEnergies> PotentialEnergies_const_p;
}  // namespace NonlinearOptimization
}  // namespace pgo