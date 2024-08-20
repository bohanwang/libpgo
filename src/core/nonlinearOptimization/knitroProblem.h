/*
author: Bohan Wang
copyright to USC
*/

#pragma once

#include "nonlinearProblem.h"

#include <memory>
#include <functional>

namespace pgo
{
namespace NonlinearOptimization
{

class KnitroProblem : public NonlinearProblem
{
public:
  KnitroProblem(PotentialEnergy_const_p energy, ConstraintFunctions_const_p constraints);
  KnitroProblem(PotentialEnergy_const_p energy);

  // KnitroProblem(PotentialEnergyDense_const_p energyDense, ConstraintFunctionsDense_const_p constraintsDense);
  // KnitroProblem(PotentialEnergyDense_const_p energy);

  virtual ~KnitroProblem();

  static double inf();
  EigenSupport::VXd &hessianVector(EigenSupport::ConstRefVecXd x, EigenSupport::ConstRefVecXd lambda, EigenSupport::ConstRefVecXd vec, double objScale);

  const EigenSupport::VXd &getCurrentx() const { return xcur; }
  const EigenSupport::VXd &getCurrentlambda() const { return lambdacur; }
  void setCurrentSolution(EigenSupport::ConstRefVecXd x, EigenSupport::ConstRefVecXd lambdacur);

  typedef std::function<void(KnitroProblem *const)> CallbackFunc;
  void setIterationCallback(CallbackFunc func) { callbackFunc = func; }
  void iterationCallback();

protected:
  EigenSupport::VXd hVec, hVec1;
  EigenSupport::VXd xcur, lambdacur;
  CallbackFunc callbackFunc;
};

}  // namespace NonlinearOptimization
}  // namespace pgo
