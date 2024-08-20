/*
author: Bohan Wang
copyright to USC, MIT
*/

#pragma once

#include "EigenSupport.h"

#include <memory>
#include <vector>

namespace pgo
{
namespace NonlinearOptimization
{
class ConstraintFunction
{
public:
  ConstraintFunction() {}
  virtual ~ConstraintFunction() {}

  virtual double func(EigenSupport::ConstRefVecXd x) const = 0;
  virtual void jacobian(EigenSupport::ConstRefVecXd x, EigenSupport::SpMatD &jac) const = 0;
  virtual void hessian(EigenSupport::ConstRefVecXd x, EigenSupport::SpMatD &hess) const = 0;

  virtual void createJacobian(EigenSupport::SpMatD &jac) const = 0;
  virtual void createHessian(EigenSupport::SpMatD &hess) const = 0;

  virtual const std::vector<int> &getDOFs() const { return dofs; }

protected:
  std::vector<int> dofs;
};

typedef std::shared_ptr<ConstraintFunction> ConstraintFunction_p;
typedef std::shared_ptr<const ConstraintFunction> ConstraintFunction_const_p;

}  // namespace NonlinearOptimization
}  // namespace pgo