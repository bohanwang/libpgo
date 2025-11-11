/*
author: Bohan Wang
copyright to USC
*/

#pragma once

#include "potentialEnergy.h"
#include "constraintFunctions.h"

namespace pgo
{
namespace NonlinearOptimization
{
class NonlinearProblem
{
public:
  NonlinearProblem(PotentialEnergy_const_p energy, ConstraintFunctions_const_p constraints);
  NonlinearProblem(PotentialEnergy_const_p energy);

  // NonlinearProblem(PotentialEnergyDense_const_p energyDense, ConstraintFunctionsDense_const_p constraintsDense);
  // NonlinearProblem(PotentialEnergyDense_const_p energy);

  virtual ~NonlinearProblem();

  void setInit(EigenSupport::ConstRefVecXd x) { xinit = x; }
  void setRange(EigenSupport::ConstRefVecXd low, EigenSupport::ConstRefVecXd hi) { xlow = low, xhi = hi; }
  void setConstraintsRange(EigenSupport::ConstRefVecXd low, EigenSupport::ConstRefVecXd hi);

  const EigenSupport::VXd &getXLow() const { return xlow; }
  const EigenSupport::VXd &getXHi() const { return xhi; }

  const EigenSupport::VXd &getCLow() const { return clow; }
  const EigenSupport::VXd &getCHi() const { return chi; }

  const EigenSupport::VXd &getXInit() const { return xinit; }

  // bool isDense() const { return problemDense != nullptr; }
  int getn() const { return problem->getNumDOFs(); }
  int getm() const { return constraints ? constraints->getNumConstraints() : 0; }
  bool hasConstraints() const { return getm() > 0; }

  // sparse routine
  ConstraintFunctions_const_p getConstraintFunctions() const { return constraints; }
  PotentialEnergy_const_p getObjectiveFunction() const { return problem; }

  int getNumHessianNonZeros(int half) const;
  int getNumJacobianNonZeros() const { return constraints ? (int)jac.nonZeros() : 0; }
  EigenSupport::SpMatD &hessian(EigenSupport::ConstRefVecXd x, EigenSupport::ConstRefVecXd lambda, double objScale);
  EigenSupport::SpMatD &getJacBuffer() { return jac; }
  EigenSupport::SpMatD &getLambdaHessianBuffer() { return lambdaHessian; }
  EigenSupport::SpMatD &getEnergyHessianBuffer() { return energyHessian; }
  EigenSupport::SpMatD &getFinalHessianBuffer() { return hessianAll; }

  int isDense() const { return 0; }

  // dense routine
  // ConstraintFunctionsDense_const_p getDenseConstraintFunctions() const { return constraintsDense; }
  // PotentialEnergyDense_const_p getDenseObjectiveFunction() const { return problemDense; }
  // EigenSupport::MXd &hessianDense(EigenSupport::ConstRefVecXd x, EigenSupport::ConstRefVecXd lambda, double objScale);

protected:
  PotentialEnergy_const_p problem;
  ConstraintFunctions_const_p constraints;

  // PotentialEnergyDense_const_p problemDense;
  // ConstraintFunctionsDense_const_p constraintsDense;

  EigenSupport::VXd xlow, xhi, clow, chi;
  Eigen::VectorXd xinit;

  EigenSupport::SpMatD hessianAll, jac, lambdaHessian, energyHessian;
  EigenSupport::SpMatI energyHessianMapping;
  EigenSupport::SpMatI lambdaHessianMapping;

  EigenSupport::MXd hessianAllDense, lambdaHessianDense, energyHessianDense;
  EigenSupport::MXd jacDense;
};

inline void NonlinearProblem::setConstraintsRange(EigenSupport::ConstRefVecXd low, EigenSupport::ConstRefVecXd hi)
{
  if (constraints) {
    clow = low;
    chi = hi;
  }
}

}  // namespace NonlinearOptimization
}  // namespace pgo