/*
author: Bohan Wang
copyright to USC
*/

#pragma once

#include "constraintFunctions.h"

#include "EigenSupport.h"

#include <memory>
#include <vector>

namespace pgo
{
namespace ES = EigenSupport;

namespace NonlinearOptimization
{
class LinearConstraintFunctions : public ConstraintFunctions
{
public:
  LinearConstraintFunctions(const EigenSupport::SpMatD &C, EigenSupport::ConstRefVecXd d);

  void setd(EigenSupport::ConstRefMatXd &newd) { d.noalias() = newd; }

  virtual void func(EigenSupport::ConstRefVecXd x, EigenSupport::RefVecXd g) const override;
  virtual void jacobian(EigenSupport::ConstRefVecXd x, EigenSupport::SpMatD &jac) const override;
  virtual void hessian(EigenSupport::ConstRefVecXd x, EigenSupport::ConstRefVecXd lambda, EigenSupport::SpMatD &hess) const override;
  virtual void hessianVector(EigenSupport::ConstRefVecXd x, EigenSupport::ConstRefVecXd lambda, EigenSupport::ConstRefVecXd vec, EigenSupport::RefVecXd hessVec) const override;

  virtual void createJacobian(EigenSupport::SpMatD &jac) const override { jac = jacConst; }
  virtual void createHessian(EigenSupport::SpMatD &hess) const override { hess = lambdahZero; }
  virtual int getNumConstraints() const { return (int)jacConst.rows(); }

  virtual const EigenSupport::SpMatD &getlambdaHessianTemplate() const override { return lambdahZero; }
  virtual const EigenSupport::SpMatD &getJacobianTemplate() const override { return jacConst; }
  virtual int getNNZJacobian() const { return (int)jacConst.nonZeros(); }

  virtual bool isQuadratic() const override { return false; }
  virtual bool isLinear() const override { return true; }
  virtual bool hasHessianVector() const override { return true; }

protected:
  const EigenSupport::SpMatD &jacConst;
  EigenSupport::SpMatD lambdahZero;
  EigenSupport::VXd d;
};

typedef std::shared_ptr<LinearConstraintFunctions> LinearConstraintFunctions_p;
typedef std::shared_ptr<const LinearConstraintFunctions> LinearConstraintFunctions_const_p;

}  // namespace NonlinearOptimization
}  // namespace pgo