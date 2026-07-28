/*
copyright to Bohan Wang
*/

// Full implementation of Codimensional IPC collision handling.
// =============================================================================

#include "ipc/core/surfaceIPCCore.h"
#include "ipc/broadPhase/surfaceIPCSelfBroadPhase.h"
#include "ipc/core/surfaceIPCBarrierAssembler.h"
#include "ipc/core/surfaceIPCMaxStep.h"

#include <algorithm>
#include <stdexcept>

namespace pgo {
namespace Contact {
namespace CIPC {
SurfaceIPCCore::SurfaceIPCCore(const SurfaceIPCCore &other):
  dhat(other.dhat),
  kappa(other.kappa),
  eps_ee(other.eps_ee),
  slackness(other.slackness),
  topology_(other.topology_),
  ptPairs_(other.ptPairs_),
  eePairs_(other.eePairs_),
  hasPreparedState_(other.hasPreparedState_),
  preparedPositions_(other.preparedPositions_)
{
}

SurfaceIPCCore &SurfaceIPCCore::operator=(const SurfaceIPCCore &other)
{
  if (this == &other)
    return *this;

  dhat = other.dhat;
  kappa = other.kappa;
  eps_ee = other.eps_ee;
  slackness = other.slackness;
  topology_ = other.topology_;
  ptPairs_ = other.ptPairs_;
  eePairs_ = other.eePairs_;
  hasPreparedState_ = other.hasPreparedState_;
  preparedPositions_ = other.preparedPositions_;
  return *this;
}

void SurfaceIPCCore::setParameters(const Parameters &params)
{
  dhat = params.dhat;
  kappa = params.kappa;
  eps_ee = params.eps_ee;
  slackness = params.slackness;
  invalidatePreparedState();
}

SurfaceIPCCore::Parameters SurfaceIPCCore::getParameters() const
{
  Parameters params;
  params.dhat = dhat;
  params.kappa = kappa;
  params.eps_ee = eps_ee;
  params.slackness = slackness;
  return params;
}

// =========================================================================
//  CollisionIPC  —  mesh setup
// =========================================================================

void SurfaceIPCCore::setMesh(const MXd &V, const MXi &F)
{
  topology_.setMesh(V, F);
  invalidatePreparedState();
}

// =========================================================================
//  Spatial hash grid for O(n) broad-phase collision detection
// =========================================================================

// =========================================================================
//  Broad phase: find candidate PT and EE pairs using spatial hashing
//  Uses insert-then-query: insert one type, query with the other.
// =========================================================================
void SurfaceIPCCore::findCollisionPairs(const VXd &positions) const
{
  SurfaceIPCSelfBroadPhase().buildPairs(topology_, positions, dhat, ptPairs_, eePairs_);
}

void SurfaceIPCCore::invalidatePreparedState() const
{
  hasPreparedState_ = false;
  preparedPositions_.resize(0);
  ptPairs_.clear();
  eePairs_.clear();
}

bool SurfaceIPCCore::isPreparedFor(EigenSupport::ConstRefVecXd x_surf) const
{
  return hasPreparedState_ &&
    preparedPositions_.size() == x_surf.size() &&
    (preparedPositions_.array() == x_surf.array()).all();
}

void SurfaceIPCCore::prepareForSurfacePositions(EigenSupport::ConstRefVecXd x_surf) const
{
  preparedPositions_ = x_surf;
  findCollisionPairs(preparedPositions_);
  hasPreparedState_ = true;
}

void SurfaceIPCCore::requirePreparedState() const
{
  if (!hasPreparedState_)
    throw std::logic_error("SurfaceIPCCore prepared active pairs are missing. Call prepareForSurfacePositions() first.");
}

// =========================================================================
//  1)  Maximum step size  (CCD-based line search with spatial hashing)
// =========================================================================
double SurfaceIPCCore::computeMaxStepSize(EigenSupport::ConstRefVecXd x, EigenSupport::ConstRefVecXd dx) const
{
  const double alpha = SurfaceIPCMaxStep().compute(topology_, x, dx, dhat, slackness);
  return std::max(alpha, 1e-12);
}

// =========================================================================
//  2)  Energy
// =========================================================================
double SurfaceIPCCore::computeEnergy(EigenSupport::ConstRefVecXd pos) const
{
  prepareForSurfacePositions(pos);
  return computeEnergyWithPreparedPairs();
}

double SurfaceIPCCore::computeEnergyWithPreparedPairs() const
{
  requirePreparedState();
  return SurfaceIPCBarrierAssembler().computeEnergy(preparedPositions_, ptPairs_, eePairs_, topology_.numVerts, dhat, kappa, eps_ee);
}

// =========================================================================
//  2)  Gradient
// =========================================================================
void SurfaceIPCCore::computeGradient(EigenSupport::ConstRefVecXd pos, EigenSupport::RefVecXd grad) const
{
  prepareForSurfacePositions(pos);
  computeGradientWithPreparedPairs(grad);
}

void SurfaceIPCCore::computeGradientWithPreparedPairs(EigenSupport::RefVecXd grad) const
{
  requirePreparedState();
  SurfaceIPCBarrierAssembler().computeGradient(preparedPositions_, ptPairs_, eePairs_, topology_.numVerts, dhat, kappa, eps_ee, grad);
}

// =========================================================================
//  2)  Sparse Hessian
// =========================================================================
void SurfaceIPCCore::computeHessian(EigenSupport::ConstRefVecXd pos, SpMatD &hess) const
{
  prepareForSurfacePositions(pos);
  computeHessianWithPreparedPairs(hess);
}

void SurfaceIPCCore::computeHessianWithPreparedPairs(SpMatD &hess) const
{
  requirePreparedState();
  SurfaceIPCBarrierAssembler().computeHessian(preparedPositions_, ptPairs_, eePairs_, topology_.numVerts, dhat, kappa, eps_ee, hess);
}

// =========================================================================
//  Combined computation (single broad-phase pass)
// =========================================================================
void SurfaceIPCCore::computeAll(EigenSupport::ConstRefVecXd x,
  double &energy, VXd &grad, SpMatD &hess) const
{
  prepareForSurfacePositions(x);
  computeAllWithPreparedPairs(energy, grad, hess);
}

void SurfaceIPCCore::computeAllWithPreparedPairs(double &energy, VXd &grad, SpMatD &hess) const
{
  requirePreparedState();
  SurfaceIPCBarrierAssembler().computeAll(preparedPositions_, ptPairs_, eePairs_, topology_.numVerts, dhat, kappa, eps_ee, energy, grad, hess);
}

}  // namespace CIPC
}  // namespace Contact
}  // namespace pgo
