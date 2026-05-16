#pragma once

#include "EigenDef.h"
#include "ipc/core/surfaceIPCPairs.h"
#include "ipc/external/obstacleSurface.h"

#include <vector>

namespace pgo
{
namespace Contact
{
namespace CIPC
{

double computeExternalEnergy(
  EigenSupport::ConstRefVecXd dynPos,
  const std::vector<ObstacleSurface> &obstacles,
  const ExternalPairSet &pairs,
  double dhat,
  double kappa,
  double eps_ee);

void computeExternalGradient(
  EigenSupport::ConstRefVecXd dynPos,
  const std::vector<ObstacleSurface> &obstacles,
  const ExternalPairSet &pairs,
  int numDynVerts,
  double dhat,
  double kappa,
  double eps_ee,
  EigenSupport::RefVecXd grad);

void computeExternalHessian(
  EigenSupport::ConstRefVecXd dynPos,
  const std::vector<ObstacleSurface> &obstacles,
  const ExternalPairSet &pairs,
  int numDynVerts,
  double dhat,
  double kappa,
  double eps_ee,
  EigenSupport::SpMatD &hess);

void computeExternalAll(
  EigenSupport::ConstRefVecXd dynPos,
  const std::vector<ObstacleSurface> &obstacles,
  const ExternalPairSet &pairs,
  int numDynVerts,
  double dhat,
  double kappa,
  double eps_ee,
  double &energy,
  EigenSupport::VXd &grad,
  EigenSupport::SpMatD &hess);

}  // namespace CIPC
}  // namespace Contact
}  // namespace pgo
