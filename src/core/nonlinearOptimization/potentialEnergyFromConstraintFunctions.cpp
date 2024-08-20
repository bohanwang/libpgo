/*
author: Bohan Wang
copyright to USC
*/

#include "potentialEnergyFromConstraintFunctions.h"
#include "EigenSupport.h"

#include <numeric>
#include <iostream>

using namespace pgo;
using namespace pgo::NonlinearOptimization;

PotentialEnergyConstraintFunctions::PotentialEnergyConstraintFunctions(int nAll, std::shared_ptr<const ConstraintFunctions> c):
  cnstt(c)
{
  g.resize(cnstt->getNumConstraints());
  cnstt->createJacobian(jac);
  cnstt->createHessian(lambdaHessian);

  allDOFs.resize(nAll);
  std::iota(allDOFs.begin(), allDOFs.end(), 0);

  if (lambdaHessian.nonZeros()) {
    ES::symbolicMm(jac, jac, JTJ, &mmData, 1);
    ES::mergeSparseMatrix(hessAll, JTJ, lambdaHessian);
    ES::small2Big(JTJ, hessAll, 0, 0, JTJMapping);
    ES::small2Big(lambdaHessian, hessAll, 0, 0, lambdaHessMapping);
  }
  else {
    ES::VXd zero(nAll);
    zero.setZero();

    cnstt->jacobian(zero, jac);
    ES::mm(jac, jac, JTJ, 1);
    hessAll = JTJ;
  }
}

PotentialEnergyConstraintFunctions::~PotentialEnergyConstraintFunctions()
{
  if (lambdaHessian.nonZeros()) {
    ES::destroySymbolicMmData(mmData);
  }
}

double PotentialEnergyConstraintFunctions::func(ES::ConstRefVecXd x) const
{
  g.setZero();
  cnstt->func(x, g);
  return g.dot(g) * 0.5;
}

void PotentialEnergyConstraintFunctions::gradient(ES::ConstRefVecXd x, ES::RefVecXd grad) const
{
  cnstt->jacobian(x, jac);
  cnstt->func(x, g);
  ES::mv(jac, g, grad, 1);
}

void PotentialEnergyConstraintFunctions::hessian(ES::ConstRefVecXd x, ES::SpMatD &hess) const
{
  memset(hess.valuePtr(), 0, hess.nonZeros() * sizeof(double));

  if (lambdaHessian.nonZeros()) {
    // g^2
    // dE = gT dgdx
    // d2E = dgdx^T dgdx + gT d2gdx
    cnstt->jacobian(x, jac);
    ES::mm(jac, jac, mmData, JTJ, 1);
    ES::addSmallToBig(1.0, JTJ, hess, 1.0, JTJMapping, 1);

    // std::cout << "Jac:\n";
    // std::cout << ES::MXd(jac) << std::endl;

    cnstt->func(x, g);
    cnstt->hessian(x, g, lambdaHessian);

    // std::cout << "g:\n"
    //           << g << std::endl;

    // std::cout << "x:\n"
    //           << x << std::endl;

    // std::cout << "H:\n";
    // std::cout << ES::MXd(lambdaHessian) << std::endl;

    ES::addSmallToBig(1.0, lambdaHessian, hess, 1.0, lambdaHessMapping, 1);
  }
  else {
    memcpy(hess.valuePtr(), JTJ.valuePtr(), sizeof(double) * JTJ.nonZeros());
  }
}