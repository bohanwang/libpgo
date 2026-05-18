#pragma once

#include "EigenDef.h"
#include "ipc/core/surfaceIPCPairs.h"

#include <vector>

namespace pgo
{
namespace Contact
{
namespace CIPC
{

double computeSelfEnergy(
  EigenSupport::ConstRefVecXd dynPos,
  const SelfPairSet &pairs,
  int numVerts,
  double dhat,
  double kappa,
  double eps_ee);

void computeSelfGradient(
  EigenSupport::ConstRefVecXd dynPos,
  const SelfPairSet &pairs,
  int numVerts,
  double dhat,
  double kappa,
  double eps_ee,
  EigenSupport::RefVecXd grad);

void computeSelfHessian(
  EigenSupport::ConstRefVecXd dynPos,
  const SelfPairSet &pairs,
  int numVerts,
  double dhat,
  double kappa,
  double eps_ee,
  EigenSupport::SpMatD &hess);

void computeSelfAll(
  EigenSupport::ConstRefVecXd dynPos,
  const SelfPairSet &pairs,
  int numVerts,
  double dhat,
  double kappa,
  double eps_ee,
  double &energy,
  EigenSupport::VXd &grad,
  EigenSupport::SpMatD &hess);

}  // namespace CIPC
}  // namespace Contact
}  // namespace pgo
