/*
copyright to Bohan Wang
*/

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

class SurfaceIPCBarrierAssembler
{
public:
  double computeEnergy(
    EigenSupport::ConstRefVecXd pos,
    const std::vector<PTPair> &ptPairs,
    const std::vector<EEPair> &eePairs,
    int numVerts,
    double dhat,
    double kappa,
    double eps_ee) const;

  void computeGradient(
    EigenSupport::ConstRefVecXd pos,
    const std::vector<PTPair> &ptPairs,
    const std::vector<EEPair> &eePairs,
    int numVerts,
    double dhat,
    double kappa,
    double eps_ee,
    EigenSupport::RefVecXd grad) const;

  void computeHessian(
    EigenSupport::ConstRefVecXd pos,
    const std::vector<PTPair> &ptPairs,
    const std::vector<EEPair> &eePairs,
    int numVerts,
    double dhat,
    double kappa,
    double eps_ee,
    EigenSupport::SpMatD &hess) const;

  void computeAll(
    EigenSupport::ConstRefVecXd x,
    const std::vector<PTPair> &ptPairs,
    const std::vector<EEPair> &eePairs,
    int numVerts,
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
