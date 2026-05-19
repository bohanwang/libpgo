/*
copyright to Bohan Wang
*/

// Full implementation of Codimensional IPC collision handling.
// =============================================================================

#include "ipc/core/surfaceIPCCore.h"
#include "scopedProfileSection.h"
#include "ipc/broadPhase/surfaceIPCBroadPhase.h"
#include "ipc/core/surfaceIPCSelfBarrierAssembler.h"
#include "ipc/core/surfaceIPCExternalBarrierAssembler.h"
#include "ipc/core/surfaceIPCMaxStep.h"
#include "ipc/profiling/surfaceIPCProfiling.h"

#include "pgoLogging.h"

#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <utility>

namespace pgo {
namespace Contact {
namespace CIPC {
using namespace pgo::EigenSupport;
static constexpr double kSmallContactAlphaWarnThreshold = 1e-2;

SurfaceIPCCore::SurfaceIPCCore(const SurfaceIPCCore &other):
  dhat(other.dhat),
  dhat_external(other.dhat_external),
  kappa(other.kappa),
  eps_ee(other.eps_ee),
  slackness(other.slackness),
  ccd_thickness(other.ccd_thickness),
  topology_(other.topology_),
  obstacles_(other.obstacles_),
  staticObstacles_(other.staticObstacles_)
{
}

SurfaceIPCCore &SurfaceIPCCore::operator=(const SurfaceIPCCore &other)
{
  if (this == &other)
    return *this;

  dhat = other.dhat;
  dhat_external = other.dhat_external;
  kappa = other.kappa;
  eps_ee = other.eps_ee;
  slackness = other.slackness;
  ccd_thickness = other.ccd_thickness;
  topology_ = other.topology_;
  obstacles_ = other.obstacles_;
  staticObstacles_ = other.staticObstacles_;
  return *this;
}

void SurfaceIPCCore::setParameters(const Parameters &params)
{
  dhat = params.dhat;
  dhat_external = params.dhat_external;
  kappa = params.kappa;
  eps_ee = params.eps_ee;
  slackness = params.slackness;
  ccd_thickness = params.ccd_thickness;
}

SurfaceIPCCore::Parameters SurfaceIPCCore::getParameters() const
{
  Parameters params;
  params.dhat = dhat;
  params.dhat_external = dhat_external;
  params.kappa = kappa;
  params.eps_ee = eps_ee;
  params.slackness = slackness;
  params.ccd_thickness = ccd_thickness;
  return params;
}

// =========================================================================
//  CollisionIPC  —  mesh setup
// =========================================================================

void SurfaceIPCCore::setMesh(const MXd &V, const MXi &F)
{
  topology_.setMesh(V, F);
}

// =========================================================================
//  Spatial hash grid for O(n) broad-phase collision detection
// =========================================================================

// =========================================================================
//  Broad phase: find candidate PT and EE pairs using spatial hashing
//  Uses insert-then-query: insert one type, query with the other.
// =========================================================================
SurfaceIPCActiveSet SurfaceIPCCore::buildActiveSet(EigenSupport::ConstRefVecXd x_surf) const
{
  Profiling::ScopedProfileSection scopedProfile(SurfaceIPCProfileSections::kBuildActiveSet);

  SurfaceIPCActiveSet activeSet;
  activeSet.positions = x_surf;
  buildSelfPairs(topology_, activeSet.positions, dhat, activeSet.selfPairs);

  if (!obstacles_.empty())
    buildExternalPairs(topology_, activeSet.positions, obstacles_, dhat_external, activeSet.externalPairs);

  if (auto logger = Logging::lgr(); logger && logger->should_log(spdlog::level::debug)) {
    const size_t selfTotal = activeSet.selfPairs.size();
    const size_t externalTotal = activeSet.externalPairs.size();
    SPDLOG_LOGGER_DEBUG(logger,
      "SurfaceIPCCore active pairs: selfPT={} selfEE={} selfTotal={} externalPT={} externalTP={} externalEE={} externalTotal={}",
      activeSet.selfPairs.ptPairs.size(), activeSet.selfPairs.eePairs.size(), selfTotal,
      activeSet.externalPairs.ptPairs.size(), activeSet.externalPairs.tpPairs.size(), activeSet.externalPairs.eePairs.size(), externalTotal);
  }

  return activeSet;
}

// =========================================================================
//  1)  Maximum step size  (CCD-based line search with spatial hashing)
// =========================================================================
NonlinearOptimization::MaxStepResult SurfaceIPCCore::computeMaxStepLimit(EigenSupport::ConstRefVecXd x, EigenSupport::ConstRefVecXd dx) const
{
  double alpha = computeSelfMaxStep(topology_, x, dx, dhat, slackness, ccd_thickness);
  alpha = std::min(alpha, computeExternalMaxStep(topology_, x, dx, obstacles_, dhat_external, slackness, ccd_thickness));

  const double clampedAlpha = std::max(alpha, 1e-12);

  if (clampedAlpha < 1.0) {
    if (clampedAlpha > 0.0 && clampedAlpha < kSmallContactAlphaWarnThreshold) {
      SPDLOG_LOGGER_WARN(Logging::lgr(),
        "IPC contact max step produced small contactFeasibleAlpha={} (slackness={}).",
        clampedAlpha, slackness);
    }

    if (auto logger = Logging::lgr(); logger && logger->should_log(spdlog::level::trace)) {
      SPDLOG_LOGGER_TRACE(logger,
        "IPC contact clamp: contactFeasibleAlpha={} slackness={}.",
        clampedAlpha, slackness);
    }
  }

  return NonlinearOptimization::MaxStepResult::contact(clampedAlpha);
}

// =========================================================================
//  2)  Energy
// =========================================================================
double SurfaceIPCCore::computeEnergy(EigenSupport::ConstRefVecXd pos) const
{
  Profiling::ScopedProfileSection scopedProfile(SurfaceIPCProfileSections::kEnergy);
  return computeEnergy(buildActiveSet(pos));
}

double SurfaceIPCCore::computeEnergy(const SurfaceIPCActiveSet &activeSet) const
{
  Profiling::ScopedProfileSection scopedProfile(SurfaceIPCProfileSections::kActiveSetEnergy);
  double e = computeSelfEnergy(activeSet.positions, activeSet.selfPairs, topology_.numVerts, dhat, kappa, eps_ee);
  if (!obstacles_.empty()) {
    e += computeExternalEnergy(
      activeSet.positions, obstacles_, activeSet.externalPairs, dhat_external, kappa, eps_ee);
  }
  return e;
}

// =========================================================================
//  2)  Gradient
// =========================================================================
void SurfaceIPCCore::computeGradient(EigenSupport::ConstRefVecXd pos, EigenSupport::RefVecXd grad) const
{
  Profiling::ScopedProfileSection scopedProfile(SurfaceIPCProfileSections::kGradient);
  computeGradient(buildActiveSet(pos), grad);
}

void SurfaceIPCCore::computeGradient(const SurfaceIPCActiveSet &activeSet, EigenSupport::RefVecXd grad) const
{
  Profiling::ScopedProfileSection scopedProfile(SurfaceIPCProfileSections::kActiveSetGradient);
  computeSelfGradient(activeSet.positions, activeSet.selfPairs, topology_.numVerts, dhat, kappa, eps_ee, grad);
  if (!obstacles_.empty()) {
    computeExternalGradient(
      activeSet.positions, obstacles_, activeSet.externalPairs, topology_.numVerts, dhat_external, kappa, eps_ee, grad);
  }
}

// =========================================================================
//  2)  Sparse Hessian
// =========================================================================
void SurfaceIPCCore::computeHessian(EigenSupport::ConstRefVecXd pos, SpMatD &hess) const
{
  Profiling::ScopedProfileSection scopedProfile(SurfaceIPCProfileSections::kHessian);
  computeHessian(buildActiveSet(pos), hess);
}

void SurfaceIPCCore::computeHessian(const SurfaceIPCActiveSet &activeSet, SpMatD &hess) const
{
  Profiling::ScopedProfileSection scopedProfile(SurfaceIPCProfileSections::kActiveSetHessian);
  computeSelfHessian(activeSet.positions, activeSet.selfPairs, topology_.numVerts, dhat, kappa, eps_ee, hess);
  if (!obstacles_.empty()) {
    computeExternalHessian(
      activeSet.positions, obstacles_, activeSet.externalPairs, topology_.numVerts, dhat_external, kappa, eps_ee, hess);
  }
  if (auto logger = Logging::lgr(); logger && logger->should_log(spdlog::level::debug))
    SPDLOG_LOGGER_DEBUG(logger, "# nonzeros in Hessian: {}", hess.nonZeros());
}

// =========================================================================
//  Combined computation (single broad-phase pass)
// =========================================================================
void SurfaceIPCCore::computeAll(EigenSupport::ConstRefVecXd x,
  double &energy, VXd &grad, SpMatD &hess) const
{
  Profiling::ScopedProfileSection scopedProfile(SurfaceIPCProfileSections::kCombined);
  computeAll(buildActiveSet(x), energy, grad, hess);
}

void SurfaceIPCCore::computeAll(const SurfaceIPCActiveSet &activeSet, double &energy, VXd &grad, SpMatD &hess) const
{
  Profiling::ScopedProfileSection scopedProfile(SurfaceIPCProfileSections::kActiveSetCombined);
  computeSelfAll(activeSet.positions, activeSet.selfPairs, topology_.numVerts, dhat, kappa, eps_ee, energy, grad, hess);
  if (!obstacles_.empty()) {
    double extEnergy = 0.0;
    computeExternalAll(
      activeSet.positions, obstacles_, activeSet.externalPairs, topology_.numVerts, dhat_external, kappa, eps_ee, extEnergy, grad, hess);
    energy += extEnergy;
  }
  if (auto logger = Logging::lgr(); logger && logger->should_log(spdlog::level::debug))
    SPDLOG_LOGGER_DEBUG(logger, "# nonzeros in Hessian: {}", hess.nonZeros());
}

// =========================================================================
//  Obstacle (external) registration
// =========================================================================
void SurfaceIPCCore::setObstacles(std::vector<ObstacleSurface> obstacles)
{
  obstacles_ = std::move(obstacles);
  staticObstacles_.assign(obstacles_.size(), false);
  for (std::size_t slot = 0; slot < obstacles_.size(); ++slot)
    obstacles_[slot].setObjectId(static_cast<int32_t>(slot));
}

void SurfaceIPCCore::markObstacleStatic(int32_t objectId)
{
  if (objectId < 0 || static_cast<std::size_t>(objectId) >= obstacles_.size())
    return;
  staticObstacles_[static_cast<std::size_t>(objectId)] = true;
  obstacles_[static_cast<std::size_t>(objectId)].update(0.0);
}

void SurfaceIPCCore::setObstacleTime(double t)
{
  for (std::size_t i = 0; i < obstacles_.size(); ++i) {
    if (!staticObstacles_[i])
      obstacles_[i].update(t);
  }
}

}  // namespace CIPC
}  // namespace Contact
}  // namespace pgo
