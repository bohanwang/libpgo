/*
author: Bohan Wang
copyright to USC, MIT
*/

#pragma once

#include "constraintFunction.h"
#include "constraintFunctions.h"

#include "EigenSupport.h"

#include <memory>
#include <vector>

namespace pgo
{
namespace NonlinearOptimization
{
struct ConstraintFunctionsAssemblerBuffer;

class ConstraintFunctionsAssembler : public ConstraintFunctions
{
public:
  ConstraintFunctionsAssembler(int nAll);
  virtual ~ConstraintFunctionsAssembler();

  void addConstraint(ConstraintFunction_p constt);
  void addConstraint(std::shared_ptr<ConstraintFunctions> constts);
  void init();

  virtual void func(EigenSupport::ConstRefVecXd x, EigenSupport::RefVecXd g) const override;
  virtual void jacobian(EigenSupport::ConstRefVecXd x, EigenSupport::SpMatD &jac) const override;
  virtual void hessian(EigenSupport::ConstRefVecXd x, EigenSupport::ConstRefVecXd lambda, EigenSupport::SpMatD &hess) const override;
  virtual void hessianVector(EigenSupport::ConstRefVecXd x, EigenSupport::ConstRefVecXd lambda, EigenSupport::ConstRefVecXd vec, EigenSupport::RefVecXd hessVec) const override;

  virtual bool isLinear() const override { return false; }
  virtual bool isQuadratic() const override{ return false; }
  virtual bool hasHessianVector() const override { return hasHessianVectorRoutine; }

protected:
  std::vector<ConstraintFunction_p> constraints1D;
  std::vector<ConstraintFunctions_p> constraintsND;
  std::vector<EigenSupport::IDX> constraintNDOffsets;

  std::vector<EigenSupport::SpMatI> jacobian1DMappings, hessian1DMappings;
  std::vector<EigenSupport::SpMatI> jacobianNDMappings, hessianNDMappings;
  bool hasHessianVectorRoutine = false;

  std::shared_ptr<ConstraintFunctionsAssemblerBuffer> buf;
};

typedef std::shared_ptr<ConstraintFunctionsAssembler> ConstraintFunctionsAssembler_p;
typedef std::shared_ptr<const ConstraintFunctionsAssembler> ConstraintFunctionsAssembler_const_p;
}  // namespace NonlinearOptimization
}  // namespace pgo