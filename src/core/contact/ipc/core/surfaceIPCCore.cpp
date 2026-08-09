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
#include <cmath>
#include <stdexcept>
#include <string>

namespace pgo
{
namespace Contact
{
namespace CIPC
{
SurfaceIPCCore::SurfaceIPCCore(const SurfaceIPCCore &other):
  dhat(other.dhat),
  kappa(other.kappa),
  eps_ee(other.eps_ee),
  slackness(other.slackness),
  projectHessianToPSD(other.projectHessianToPSD),
  hasMesh_(other.hasMesh_),
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
  projectHessianToPSD = other.projectHessianToPSD;
  hasMesh_ = other.hasMesh_;
  topology_ = other.topology_;
  ptPairs_ = other.ptPairs_;
  eePairs_ = other.eePairs_;
  hasPreparedState_ = other.hasPreparedState_;
  preparedPositions_ = other.preparedPositions_;
  return *this;
}

void SurfaceIPCCore::setParameters(const Parameters &params)
{
  if (!std::isfinite(params.dhat) || params.dhat <= 0.0)
    throw std::invalid_argument("SurfaceIPCCore::Parameters::dhat must be finite and positive.");
  if (!std::isfinite(params.kappa) || params.kappa <= 0.0)
    throw std::invalid_argument("SurfaceIPCCore::Parameters::kappa must be finite and positive.");
  if (!std::isfinite(params.eps_ee) || params.eps_ee < 0.0)
    throw std::invalid_argument("SurfaceIPCCore::Parameters::eps_ee must be finite and non-negative.");
  if (!std::isfinite(params.slackness) || params.slackness <= 0.0 || params.slackness > 1.0)
    throw std::invalid_argument("SurfaceIPCCore::Parameters::slackness must be finite and in (0, 1].");

  dhat = params.dhat;
  kappa = params.kappa;
  eps_ee = params.eps_ee;
  slackness = params.slackness;
  projectHessianToPSD = params.projectHessianToPSD;
  invalidatePreparedState();
}

SurfaceIPCCore::Parameters SurfaceIPCCore::getParameters() const
{
  Parameters params;
  params.dhat = dhat;
  params.kappa = kappa;
  params.eps_ee = eps_ee;
  params.slackness = slackness;
  params.projectHessianToPSD = projectHessianToPSD;
  return params;
}

// =========================================================================
//  CollisionIPC  —  mesh setup
// =========================================================================

void SurfaceIPCCore::setMesh(const MXd &V, const MXi &F)
{
  setMesh(V, F, std::vector<uint8_t>(static_cast<std::size_t>(V.rows()), uint8_t{ 1 }));
}

void SurfaceIPCCore::setMesh(const MXd &V, const MXi &F, const std::vector<uint8_t> &vertexIsDeformableMask)
{
  if (V.cols() != 3 || V.rows() <= 0)
    throw std::invalid_argument("SurfaceIPCCore mesh vertices must be a non-empty N x 3 matrix.");
  if (!V.allFinite())
    throw std::invalid_argument("SurfaceIPCCore mesh vertices must be finite.");
  if (F.cols() != 3)
    throw std::invalid_argument("SurfaceIPCCore mesh triangles must be an M x 3 index matrix.");
  if (F.size() > 0 && (F.minCoeff() < 0 || F.maxCoeff() >= V.rows()))
    throw std::invalid_argument("SurfaceIPCCore mesh triangles contain an out-of-range vertex index.");
  for (Eigen::Index fi = 0; fi < F.rows(); ++fi) {
    if (F(fi, 0) == F(fi, 1) || F(fi, 1) == F(fi, 2) || F(fi, 2) == F(fi, 0))
      throw std::invalid_argument("SurfaceIPCCore mesh triangles must reference three distinct vertices.");
  }

  topology_.setMesh(V, F, vertexIsDeformableMask);
  hasMesh_ = true;
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
  validateSurfaceState(x_surf, "surface positions");
  preparedPositions_ = x_surf;
  findCollisionPairs(preparedPositions_);
  hasPreparedState_ = true;
}

void SurfaceIPCCore::validateCollisionFreeSurfacePositions(
  EigenSupport::ConstRefVecXd x_surf) const
{
  validateSurfaceState(x_surf, "initial surface positions");
  invalidatePreparedState();
  findCollisionPairs(x_surf);

  for (const PTPair &pair : ptPairs_) {
    const double d2 = distance::computePTSqDist(
      vtx(x_surf, pair.p), vtx(x_surf, pair.t0), vtx(x_surf, pair.t1), vtx(x_surf, pair.t2));
    if (!std::isfinite(d2) || d2 <= 0.0) {
      throw std::invalid_argument(
        "IPC initial collision surface is not strictly collision-free: zero/non-finite point-triangle distance for point " +
        std::to_string(pair.p) + ".");
    }
  }

  for (const EEPair &pair : eePairs_) {
    const double d2 = distance::computeEESqDist(
      vtx(x_surf, pair.ea0), vtx(x_surf, pair.ea1), vtx(x_surf, pair.eb0), vtx(x_surf, pair.eb1));
    if (!std::isfinite(d2) || d2 <= 0.0) {
      throw std::invalid_argument(
        "IPC initial collision surface is not strictly collision-free: zero/non-finite edge-edge distance for edge (" +
        std::to_string(pair.ea0) + ", " + std::to_string(pair.ea1) + ").");
    }
  }

  invalidatePreparedState();
}

void SurfaceIPCCore::validateSurfaceState(EigenSupport::ConstRefVecXd x_surf, const char *argumentName) const
{
  if (!hasMesh_)
    throw std::logic_error("SurfaceIPCCore requires setMesh() before evaluation.");
  if (x_surf.size() != topology_.numSurfaceDOFs())
    throw std::invalid_argument(std::string("SurfaceIPCCore ") + argumentName + " size must equal 3 * numSurfaceVertices.");
  if (!x_surf.allFinite())
    throw std::invalid_argument(std::string("SurfaceIPCCore ") + argumentName + " must be finite.");
}

void SurfaceIPCCore::validateSurfaceGradient(EigenSupport::RefVecXd g_surf) const
{
  if (g_surf.size() != topology_.numSurfaceDOFs())
    throw std::invalid_argument("SurfaceIPCCore gradient size must equal 3 * numSurfaceVertices.");
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
  validateSurfaceState(x, "surface positions");
  validateSurfaceState(dx, "surface displacement");
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
  validateSurfaceGradient(grad);
  prepareForSurfacePositions(pos);
  computeGradientWithPreparedPairs(grad);
}

void SurfaceIPCCore::computeGradientWithPreparedPairs(EigenSupport::RefVecXd grad) const
{
  requirePreparedState();
  validateSurfaceGradient(grad);
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
  SurfaceIPCBarrierAssembler().computeHessian(preparedPositions_, ptPairs_, eePairs_, topology_.numVerts,
    dhat, kappa, eps_ee, projectHessianToPSD, hess);
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
  SurfaceIPCBarrierAssembler().computeAll(preparedPositions_, ptPairs_, eePairs_, topology_.numVerts,
    dhat, kappa, eps_ee, projectHessianToPSD, energy, grad, hess);
}

}  // namespace CIPC
}  // namespace Contact
}  // namespace pgo
