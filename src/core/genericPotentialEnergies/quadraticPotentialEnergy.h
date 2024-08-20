/*
author: Bohan Wang
copyright to USC, MIT
*/

#pragma once

#include "potentialEnergy.h"

#include <numeric>
#include <vector>

namespace pgo
{
namespace PredefinedPotentialEnergies
{
class QuadraticPotentialEnergyCache;

class QuadraticPotentialEnergy : public NonlinearOptimization::PotentialEnergy
{
public:
  // the energy is 1/2 xT A x
  QuadraticPotentialEnergy(const EigenSupport::SpMatD &A_);
  // the energy is 1/2 xT AT A x
  QuadraticPotentialEnergy(const EigenSupport::SpMatD &A_, int inParentheses);
  // the energy is 1/2 xT AT W A x
  QuadraticPotentialEnergy(const EigenSupport::SpMatD &A_, const double *W, int inParentheses);
  // the energy is 1/2 xT A x + bT x
  QuadraticPotentialEnergy(const EigenSupport::SpMatD &A_, const EigenSupport::VXd &b_);
  // the energy is 1/2 (A x + b)^2
  // = 1/2 x^T (A^T A) x + b^T A x + 1/2 b^T b
  QuadraticPotentialEnergy(const EigenSupport::SpMatD &A_, const EigenSupport::VXd &b_, int inParentheses);
  // the energy is 1/2 (A x + b)^T W (Ax + b)
  // = 1/2 x^T A^T W A x + 1/2 x^T A^T W b + 1/2 b^T W Ax + 1/2 b^T W b
  // = 1/2 x^T (A^T W A) x + b^T W Ax + 1/2 b^T W b
  QuadraticPotentialEnergy(const EigenSupport::SpMatD &A_, const EigenSupport::VXd &b_, const double *W, int inParentheses);

  void setDOFs(const std::vector<int> &dofs);

  virtual double func(EigenSupport::ConstRefVecXd x) const override;
  virtual void gradient(EigenSupport::ConstRefVecXd x, EigenSupport::RefVecXd grad) const override;
  virtual void hessian(EigenSupport::ConstRefVecXd x, EigenSupport::SpMatD &hess) const override;
  virtual void hessianVector(EigenSupport::ConstRefVecXd x, EigenSupport::ConstRefVecXd vec, EigenSupport::RefVecXd hVec) const override;

  virtual void createHessian(EigenSupport::SpMatD &hess) const override { hess = A; }
  virtual void getDOFs(std::vector<int> &dofs) const override { dofs = allDOFs; }
  virtual int getNumDOFs() const override { return (int)allDOFs.size(); }

  void gradientComponent(EigenSupport::SpMatD *A, EigenSupport::VXd *b) const;

  virtual int isQuadratic() const override { return 1; }
  virtual int hasHessianVector() const override { return 1; }

protected:
  const EigenSupport::SpMatD &A;
  const EigenSupport::VXd *b;
  int isInParentheses;

  std::vector<int> allDOFs;
  EigenSupport::SpMatD ATA;
  EigenSupport::VXd bTA;
  double c = 0.0;

  std::shared_ptr<QuadraticPotentialEnergyCache> cache;
};
}  // namespace NonlinearOptimization
}  // namespace pgo