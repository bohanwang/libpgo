/*
author: Bohan Wang
copyright to USC
*/

#include "linearConstraintFunctions.h"
#include "EigenSupport.h"

using namespace pgo;
using namespace pgo::NonlinearOptimization;

LinearConstraintFunctions::LinearConstraintFunctions(const ES::SpMatD &C, ES::ConstRefVecXd d_):
  ConstraintFunctions(static_cast<int>(C.cols())), jacConst(C), d(d_)
{
  lambdahZero.resize(C.cols(), C.cols());
}

void LinearConstraintFunctions::func(ES::ConstRefVecXd x, ES::RefVecXd g) const
{
  EigenSupport::mv(jacConst, x, g);
  g.noalias() += d;
}

void LinearConstraintFunctions::jacobian(ES::ConstRefVecXd, ES::SpMatD &jac) const
{
  memcpy(jac.valuePtr(), jacConst.valuePtr(), jacConst.nonZeros() * sizeof(double));
}

void LinearConstraintFunctions::hessian(ES::ConstRefVecXd, ES::ConstRefVecXd, ES::SpMatD &hess) const
{
  if (hess.valuePtr())
    memset(hess.valuePtr(), 0, sizeof(double) * hess.nonZeros());
}

void LinearConstraintFunctions::hessianVector(ES::ConstRefVecXd, ES::ConstRefVecXd, ES::ConstRefVecXd, ES::RefVecXd hessVec) const
{
  hessVec.setZero();
}
