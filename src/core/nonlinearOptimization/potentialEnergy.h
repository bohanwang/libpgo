/*
author: Bohan Wang
copyright to USC
*/

#pragma once

#include "EigenDef.h"

#include <memory>

namespace pgo
{
namespace NonlinearOptimization
{
class PotentialEnergy
{
public:
  PotentialEnergy() {}
  virtual ~PotentialEnergy() {}

  virtual double func(EigenSupport::ConstRefVecXd x) const = 0;
  virtual void gradient(EigenSupport::ConstRefVecXd x, EigenSupport::RefVecXd grad) const = 0;
  virtual void hessian(EigenSupport::ConstRefVecXd x, EigenSupport::SpMatD &hess) const = 0;
  virtual void hessianVector(EigenSupport::ConstRefVecXd x, EigenSupport::ConstRefVecXd vec, EigenSupport::RefVecXd hessVec) const;

  virtual void createHessian(EigenSupport::SpMatD &hess) const = 0;

  virtual void getDOFs(std::vector<int> &dofs) const = 0;
  virtual int getNumDOFs() const = 0;

  virtual int isQuadratic() const { return 0; }
  virtual int hasHessianVector() const { return 0; }
};

typedef std::shared_ptr<PotentialEnergy> PotentialEnergy_p;
typedef std::shared_ptr<const PotentialEnergy> PotentialEnergy_const_p;

inline void PotentialEnergy::hessianVector(EigenSupport::ConstRefVecXd, EigenSupport::ConstRefVecXd, EigenSupport::RefVecXd) const
{
}

}  // namespace NonlinearOptimization
}  // namespace pgo
