/*
author: Bohan Wang
copyright to USC, MIT
*/

#pragma once

#include "EigenDef.h"

#include <memory>
#include <vector>

namespace pgo
{
namespace NonlinearOptimization
{
class ConstraintFunctions
{
public:
  ConstraintFunctions(int nAll);
  virtual ~ConstraintFunctions();

  virtual void func(EigenSupport::ConstRefVecXd x, EigenSupport::RefVecXd g) const = 0;
  virtual void jacobian(EigenSupport::ConstRefVecXd x, EigenSupport::SpMatD &jac) const = 0;
  virtual void hessian(EigenSupport::ConstRefVecXd x, EigenSupport::ConstRefVecXd lambda, EigenSupport::SpMatD &hess) const = 0;
  virtual void hessianVector(EigenSupport::ConstRefVecXd x, EigenSupport::ConstRefVecXd lambda, EigenSupport::ConstRefVecXd vec, EigenSupport::RefVecXd hessVec) const;

  virtual void createJacobian(EigenSupport::SpMatD &jac) const { jac = jacobianTemplate; }
  virtual void createHessian(EigenSupport::SpMatD &hess) const { hess = lambdahTemplate; }
  virtual int getNumConstraints() const { return (int)jacobianTemplate.rows(); }

  virtual const EigenSupport::SpMatD &getlambdaHessianTemplate() const { return lambdahTemplate; }
  virtual const EigenSupport::SpMatD &getJacobianTemplate() const { return jacobianTemplate; }
  virtual int getNNZJacobian() const { return (int)jacobianTemplate.nonZeros(); }

  virtual bool isLinear() const { return false; }
  virtual bool isQuadratic() const { return false; }
  virtual bool hasHessianVector() const { return false; }

protected:
  void approxHessianVector(EigenSupport::ConstRefVecXd x, EigenSupport::ConstRefVecXd lambda, EigenSupport::ConstRefVecXd vec, double eps,
    EigenSupport::SpMatD &jacBuf0, EigenSupport::SpMatD &jacBuf1, EigenSupport::VXd &xleft, EigenSupport::VXd &xright, EigenSupport::RefVecXd hessVec) const;

  EigenSupport::SpMatD jacobianTemplate, lambdahTemplate;
  int nAll;
};

typedef std::shared_ptr<ConstraintFunctions> ConstraintFunctions_p;
typedef std::shared_ptr<const ConstraintFunctions> ConstraintFunctions_const_p;

inline void ConstraintFunctions::hessianVector(EigenSupport::ConstRefVecXd, EigenSupport::ConstRefVecXd, EigenSupport::ConstRefVecXd, EigenSupport::RefVecXd) const
{
}

}  // namespace NonlinearOptimization
}  // namespace pgo