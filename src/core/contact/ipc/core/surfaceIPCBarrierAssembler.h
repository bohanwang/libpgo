/*
copyright to Bohan Wang
*/

#pragma once

#include "EigenDef.h"
#include "ipc/core/surfaceIPCPairs.h"
#include "ipc/external/obstacleSurface.h"

#include <memory>
#include <vector>

namespace pgo
{
namespace Contact
{
namespace CIPC
{

class SurfaceIPCBarrierAssembler
{
public:
  // --- Self pair assembly ---
  double computeSelfEnergy(
    EigenSupport::ConstRefVecXd dynPos,
    const SelfPairSet &pairs,
    int numVerts,
    double dhat,
    double kappa,
    double eps_ee) const;

  void computeSelfGradient(
    EigenSupport::ConstRefVecXd dynPos,
    const SelfPairSet &pairs,
    int numVerts,
    double dhat,
    double kappa,
    double eps_ee,
    EigenSupport::RefVecXd grad) const;

  void computeSelfHessian(
    EigenSupport::ConstRefVecXd dynPos,
    const SelfPairSet &pairs,
    int numVerts,
    double dhat,
    double kappa,
    double eps_ee,
    EigenSupport::SpMatD &hess) const;

  void computeSelfAll(
    EigenSupport::ConstRefVecXd dynPos,
    const SelfPairSet &pairs,
    int numVerts,
    double dhat,
    double kappa,
    double eps_ee,
    double &energy,
    EigenSupport::VXd &grad,
    EigenSupport::SpMatD &hess) const;

  // --- External pair assembly (dynamic-only block scatter) ---
  double computeExternalEnergy(
    EigenSupport::ConstRefVecXd dynPos,
    const std::vector<std::shared_ptr<ObstacleSurface>> &obstacles,
    const ExternalPairSet &pairs,
    double dhat,
    double kappa,
    double eps_ee) const;

  void computeExternalGradient(
    EigenSupport::ConstRefVecXd dynPos,
    const std::vector<std::shared_ptr<ObstacleSurface>> &obstacles,
    const ExternalPairSet &pairs,
    int numDynVerts,
    double dhat,
    double kappa,
    double eps_ee,
    EigenSupport::RefVecXd grad) const;

  void computeExternalHessian(
    EigenSupport::ConstRefVecXd dynPos,
    const std::vector<std::shared_ptr<ObstacleSurface>> &obstacles,
    const ExternalPairSet &pairs,
    int numDynVerts,
    double dhat,
    double kappa,
    double eps_ee,
    EigenSupport::SpMatD &hess) const;

  void computeExternalAll(
    EigenSupport::ConstRefVecXd dynPos,
    const std::vector<std::shared_ptr<ObstacleSurface>> &obstacles,
    const ExternalPairSet &pairs,
    int numDynVerts,
    double dhat,
    double kappa,
    double eps_ee,
    double &energy,
    EigenSupport::VXd &grad,
    EigenSupport::SpMatD &hess) const;
};

}  // namespace CIPC
}  // namespace Contact
}  // namespace pgo
