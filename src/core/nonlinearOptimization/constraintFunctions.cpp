/*
author: Bohan Wang
copyright to USC, MIT
*/

#include "constraintFunctions.h"
#include "EigenSupport.h"

using namespace pgo;
using namespace pgo::NonlinearOptimization;
namespace ES = pgo::EigenSupport;

ConstraintFunctions::ConstraintFunctions(int na):
  nAll(na)
{
}

ConstraintFunctions::~ConstraintFunctions()
{
}

void ConstraintFunctions::approxHessianVector(ES::ConstRefVecXd x, ES::ConstRefVecXd lambda, ES::ConstRefVecXd vec, double eps,
  ES::SpMatD &jacBuf0, ES::SpMatD &jacBuf1, ES::VXd &xleft, ES::VXd &xright, ES::RefVecXd hessVec) const
{
  xleft.noalias() = x - eps * vec;
  xright.noalias() = x + eps * vec;

  jacobian(xleft, jacBuf0);
  ES::mv(jacBuf0, lambda, xleft, 1);

  jacobian(xright, jacBuf1);
  ES::mv(jacBuf1, lambda, xright, 1);

  hessVec.noalias() = (xright - xleft) * (0.5 / eps);
}