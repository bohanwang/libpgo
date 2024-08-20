/*
author: Bohan Wang
copyright to USC
*/

#include "knitroProblem.h"

#include <knitro.h>

#include <cstring>
#include <iostream>

using namespace pgo;
using namespace pgo::NonlinearOptimization;

namespace ES = pgo::EigenSupport;

double KnitroProblem::inf()
{
  return KN_INFINITY;
}

KnitroProblem::KnitroProblem(PotentialEnergy_const_p energy_, ConstraintFunctions_const_p constraints_):
  NonlinearProblem(energy_, constraints_)
{
  for (int i = 0; i < problem->getNumDOFs(); i++) {
    xlow[i] = -inf();
    xhi[i] = inf();
  }

  for (int i = 0; i < constraints->getNumConstraints(); i++) {
    clow[i] = -inf();
    chi[i] = inf();
  }

  hVec.resize(hessianAll.rows());
  hVec1.resize(hessianAll.rows());

  xcur.resize(problem->getNumDOFs());
}

KnitroProblem::KnitroProblem(PotentialEnergy_const_p energy_):
  NonlinearProblem(energy_)
{
  for (int i = 0; i < problem->getNumDOFs(); i++) {
    xlow[i] = -inf();
    xhi[i] = inf();
  }

  hVec.resize(hessianAll.rows());

  xcur.resize(problem->getNumDOFs());
}


// KnitroProblem::KnitroProblem(PotentialEnergyDense_const_p energy_, ConstraintFunctionsDense_const_p constraints_):
//   NonlinearProblem(energy_, constraints_)
// {
//   for (int i = 0; i < problemDense->getNumDOFs(); i++) {
//     xlow[i] = -inf();
//     xhi[i] = inf();
//   }

//   for (int i = 0; i < constraintsDense->getNumConstraints(); i++) {
//     clow[i] = -inf();
//     chi[i] = inf();
//   }

//   hVec.resize(hessianAll.rows());
//   hVec1.resize(hessianAll.rows());

//   xcur.resize(problemDense->getNumDOFs());
// }

// KnitroProblem::KnitroProblem(PotentialEnergyDense_const_p energy_):
//   NonlinearProblem(energy_)
// {
//   for (int i = 0; i < problemDense->getNumDOFs(); i++) {
//     xlow[i] = -inf();
//     xhi[i] = inf();
//   }

//   hVec.resize(hessianAll.rows());

//   xcur.resize(problemDense->getNumDOFs());
// }

KnitroProblem::~KnitroProblem()
{
}

ES::VXd &KnitroProblem::hessianVector(EigenSupport::ConstRefVecXd x, EigenSupport::ConstRefVecXd lambda, EigenSupport::ConstRefVecXd vec, double objScale)
{
  memset(hVec.data(), 0, hVec.size() * sizeof(double));

  if (problem->hasHessianVector() && objScale > 0) {
    problem->hessianVector(x, vec, hVec);
    hVec *= objScale;
  }

  if (constraints && constraints->hasHessianVector()) {
    constraints->hessianVector(x, lambda, vec, hVec1);
    hVec += hVec1;
  }

  return hVec;
}

void KnitroProblem::setCurrentSolution(ES::ConstRefVecXd x, ES::ConstRefVecXd lambda)
{
  xcur.noalias() = x;
  lambdacur = lambda;
}

void KnitroProblem::iterationCallback()
{
  if (callbackFunc)
    callbackFunc(this);
}
