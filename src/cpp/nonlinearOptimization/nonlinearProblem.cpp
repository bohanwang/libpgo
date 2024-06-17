/*
author: Bohan Wang
copyright to USC
*/

#include "nonlinearProblem.h"
#include "EigenSupport.h"

#include "pgoLogging.h"

#include <cstring>
#include <iostream>

using namespace pgo;
using namespace pgo::NonlinearOptimization;
namespace ES = pgo::EigenSupport;

NonlinearProblem::NonlinearProblem(PotentialEnergy_const_p energy_, ConstraintFunctions_const_p constraints_):
  problem(energy_), constraints(constraints_)
{
  xinit.resize(problem->getNumDOFs());

  xlow.resize(problem->getNumDOFs());
  xhi.resize(problem->getNumDOFs());

  clow = Eigen::VectorXd::Zero(constraints->getNumConstraints());
  chi = Eigen::VectorXd::Zero(constraints->getNumConstraints());

  problem->createHessian(energyHessian);
  constraints->createJacobian(jac);
  constraints->createHessian(lambdaHessian);

  ES::mergeSparseMatrix(hessianAll, energyHessian, lambdaHessian);
  ES::small2Big(energyHessian, hessianAll, 0, 0, energyHessianMapping);
  ES::small2Big(lambdaHessian, hessianAll, 0, 0, lambdaHessianMapping);
}

NonlinearProblem::NonlinearProblem(PotentialEnergy_const_p energy_):
  problem(energy_)
{
  problem->createHessian(energyHessian);
  xinit.resize(problem->getNumDOFs());

  xlow.resize(problem->getNumDOFs());
  xhi.resize(problem->getNumDOFs());

  hessianAll = energyHessian;
  ES::small2Big(energyHessian, hessianAll, 0, 0, energyHessianMapping);
}

// NonlinearProblem::NonlinearProblem(PotentialEnergyDense_const_p energy_, ConstraintFunctionsDense_const_p constraints_):
//   problemDense(energy_), constraintsDense(constraints_)
// {
//   xinit.resize(problemDense->getNumDOFs());
//   xlow.resize(problemDense->getNumDOFs());
//   xhi.resize(problemDense->getNumDOFs());

//   clow = Eigen::VectorXd::Zero(constraintsDense->getNumConstraints());
//   chi = Eigen::VectorXd::Zero(constraintsDense->getNumConstraints());

//   energyHessianDense.resize(xinit.size(), xinit.size());
//   lambdaHessianDense.resize(xinit.size(), xinit.size());
//   jacDense.resize(constraintsDense->getNumConstraints(), xinit.size());

//   hessianAllDense.resize(xinit.size(), xinit.size());
// }

// NonlinearProblem::NonlinearProblem(PotentialEnergyDense_const_p energy_):
//   problemDense(energy_)
// {
//   xinit.resize(problemDense->getNumDOFs());

//   xlow.resize(problemDense->getNumDOFs());
//   xhi.resize(problemDense->getNumDOFs());

//   energyHessianDense.resize(xinit.size(), xinit.size());
//   hessianAllDense.resize(xinit.size(), xinit.size());
// }

NonlinearProblem::~NonlinearProblem()
{
}

int NonlinearProblem::getNumHessianNonZeros(int half) const
{
  if (half) {
    int inc = 0;
    for (Eigen::Index i = 0; i < hessianAll.outerSize(); i++) {
      for (Eigen::InnerIterator it(hessianAll, i); it; ++it) {
        if (it.row() > it.col())
          continue;
        inc++;
      }
    }

    return inc;
  }
  return (int)hessianAll.nonZeros();
}

ES::SpMatD &NonlinearProblem::hessian(EigenSupport::ConstRefVecXd x, EigenSupport::ConstRefVecXd lambda, double objScale)
{
  memset(energyHessian.valuePtr(), 0, sizeof(double) * energyHessian.nonZeros());
  memset(hessianAll.valuePtr(), 0, sizeof(double) * hessianAll.nonZeros());

  if (constraints && lambdaHessian.nonZeros())
    memset(lambdaHessian.valuePtr(), 0, sizeof(double) * lambdaHessian.nonZeros());

  if (objScale > 0) {
    if (constraints) {
      problem->hessian(x, energyHessian);
      ES::addSmallToBig(objScale, energyHessian, hessianAll, 0.0, energyHessianMapping, 1);

      if (lambdaHessian.nonZeros()) {
        constraints->hessian(x, lambda, lambdaHessian);
        ES::addSmallToBig(1.0, lambdaHessian, hessianAll, 1.0, lambdaHessianMapping, 1);
      }
    }
    else {
      problem->hessian(x, hessianAll);
      hessianAll *= objScale;
    }
  }
  else {
    if (constraints && lambdaHessian.nonZeros()) {
      constraints->hessian(x, lambda, lambdaHessian);
      ES::transferSmallToBig(lambdaHessian, hessianAll, lambdaHessianMapping, 1);
    }
  }

  return hessianAll;
}

// EigenSupport::MXd &NonlinearProblem::hessianDense(EigenSupport::ConstRefVecXd x, EigenSupport::ConstRefVecXd lambda, double objScale)
// {
//   energyHessianDense.setZero();
//   hessianAllDense.setZero();

//   if (constraintsDense) {
//     lambdaHessianDense.setZero();
//   }

//   if (objScale > 0) {
//     if (constraintsDense && constraintsDense->isQuadratic() == false && constraintsDense->isLinear() == false) {
//       problemDense->hessian(x, energyHessianDense);
//       constraintsDense->hessian(x, lambda, lambdaHessianDense);

//       hessianAllDense += energyHessianDense * objScale;
//       hessianAllDense += lambdaHessianDense;
//     }
//     else {
//       problemDense->hessian(x, hessianAllDense);
//       hessianAllDense *= objScale;
//     }
//   }
//   else {
//     if (constraintsDense && constraintsDense->isQuadratic() == false && constraintsDense->isLinear() == false) {
//       constraintsDense->hessian(x, lambda, hessianAllDense);
//     }
//   }

//   return hessianAllDense;
// }