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
  topology_(other.topology_),
  preparedState_(other.preparedState_),
  obstacles_(other.obstacles_)
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
  topology_ = other.topology_;
  preparedState_ = other.preparedState_;
  obstacles_ = other.obstacles_;
  return *this;
}

void SurfaceIPCCore::setParameters(const Parameters &params)
{
  dhat = params.dhat;
  dhat_external = params.dhat_external;
  kappa = params.kappa;
  eps_ee = params.eps_ee;
  slackness = params.slackness;
  preparedState_.clear();
}

SurfaceIPCCore::Parameters SurfaceIPCCore::getParameters() const
{
  Parameters params;
  params.dhat = dhat;
  params.dhat_external = dhat_external;
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
  preparedState_.clear();
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
  buildSelfPairs(topology_, positions, dhat, preparedState_.selfPairs);

  preparedState_.externalPairs.clear();

  if (obstacles_.empty())
    return;

  buildExternalPairs(topology_, positions, obstacles_, dhat_external, preparedState_.externalPairs);
}

void SurfaceIPCCore::prepareForSurfacePositions(EigenSupport::ConstRefVecXd x_surf) const
{
  Profiling::ScopedProfileSection scopedProfile(SurfaceIPCProfileSections::kPrepareActivePairs);
  preparedState_.positions = x_surf;
  findCollisionPairs(preparedState_.positions);
  if (auto logger = Logging::lgr(); logger) {
    const size_t selfTotal = preparedState_.selfPairs.ptPairs.size() + preparedState_.selfPairs.eePairs.size();
    const size_t externalTotal = preparedState_.externalPairs.ptPairs.size() + preparedState_.externalPairs.tpPairs.size() + preparedState_.externalPairs.eePairs.size();
    SPDLOG_LOGGER_INFO(logger,
      "SurfaceIPCCore active pairs: selfPT={} selfEE={} selfTotal={} externalPT={} externalTP={} externalEE={} externalTotal={}",
      preparedState_.selfPairs.ptPairs.size(), preparedState_.selfPairs.eePairs.size(), selfTotal,
      preparedState_.externalPairs.ptPairs.size(), preparedState_.externalPairs.tpPairs.size(), preparedState_.externalPairs.eePairs.size(), externalTotal);
  }
  preparedState_.hasState = true;
}

// =========================================================================
//  1)  Maximum step size  (CCD-based line search with spatial hashing)
// =========================================================================
NonlinearOptimization::MaxStepResult SurfaceIPCCore::computeMaxStepLimit(EigenSupport::ConstRefVecXd x, EigenSupport::ConstRefVecXd dx) const
{
  double alpha = computeSelfMaxStep(topology_, x, dx, dhat, slackness);
  alpha = std::min(alpha, computeExternalMaxStep(topology_, x, dx, obstacles_, dhat_external, slackness));

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
  prepareForSurfacePositions(pos);
  return computeEnergyWithPreparedPairs();
}

double SurfaceIPCCore::computeEnergyWithPreparedPairs() const
{
  Profiling::ScopedProfileSection scopedProfile(SurfaceIPCProfileSections::kPreparedEnergy);
  if (!preparedState_.hasState)
    throw std::logic_error("SurfaceIPCCore prepared active pairs are missing. Call prepareForSurfacePositions() first.");
  double e = computeSelfEnergy(preparedState_.positions, preparedState_.selfPairs, topology_.numVerts, dhat, kappa, eps_ee);
  if (!obstacles_.empty()) {
    e += computeExternalEnergy(
      preparedState_.positions, obstacles_, preparedState_.externalPairs, dhat_external, kappa, eps_ee);
  }
  return e;
}

// =========================================================================
//  2)  Gradient
// =========================================================================
void SurfaceIPCCore::computeGradient(EigenSupport::ConstRefVecXd pos, EigenSupport::RefVecXd grad) const
{
  Profiling::ScopedProfileSection scopedProfile(SurfaceIPCProfileSections::kGradient);
  prepareForSurfacePositions(pos);
  computeGradientWithPreparedPairs(grad);
}

void SurfaceIPCCore::computeGradientWithPreparedPairs(EigenSupport::RefVecXd grad) const
{
  Profiling::ScopedProfileSection scopedProfile(SurfaceIPCProfileSections::kPreparedGradient);
  if (!preparedState_.hasState)
    throw std::logic_error("SurfaceIPCCore prepared active pairs are missing. Call prepareForSurfacePositions() first.");
  computeSelfGradient(preparedState_.positions, preparedState_.selfPairs, topology_.numVerts, dhat, kappa, eps_ee, grad);
  if (!obstacles_.empty()) {
    computeExternalGradient(
      preparedState_.positions, obstacles_, preparedState_.externalPairs, topology_.numVerts, dhat_external, kappa, eps_ee, grad);
  }
}

// =========================================================================
//  2)  Sparse Hessian
// =========================================================================
void SurfaceIPCCore::computeHessian(EigenSupport::ConstRefVecXd pos, SpMatD &hess) const
{
  Profiling::ScopedProfileSection scopedProfile(SurfaceIPCProfileSections::kHessian);
  prepareForSurfacePositions(pos);
  computeHessianWithPreparedPairs(hess);
}

void SurfaceIPCCore::computeHessianWithPreparedPairs(SpMatD &hess) const
{
  Profiling::ScopedProfileSection scopedProfile(SurfaceIPCProfileSections::kPreparedHessian);
  if (!preparedState_.hasState)
    throw std::logic_error("SurfaceIPCCore prepared active pairs are missing. Call prepareForSurfacePositions() first.");
  computeSelfHessian(preparedState_.positions, preparedState_.selfPairs, topology_.numVerts, dhat, kappa, eps_ee, hess);
  if (!obstacles_.empty()) {
    computeExternalHessian(
      preparedState_.positions, obstacles_, preparedState_.externalPairs, topology_.numVerts, dhat_external, kappa, eps_ee, hess);
  }
  if (auto logger = Logging::lgr(); logger)
    SPDLOG_LOGGER_INFO(logger, "# nonzeros in Hessian: {}", hess.nonZeros());
}

// =========================================================================
//  Combined computation (single broad-phase pass)
// =========================================================================
void SurfaceIPCCore::computeAll(EigenSupport::ConstRefVecXd x,
  double &energy, VXd &grad, SpMatD &hess) const
{
  Profiling::ScopedProfileSection scopedProfile(SurfaceIPCProfileSections::kCombined);
  prepareForSurfacePositions(x);
  computeAllWithPreparedPairs(energy, grad, hess);
}

void SurfaceIPCCore::computeAllWithPreparedPairs(double &energy, VXd &grad, SpMatD &hess) const
{
  if (!preparedState_.hasState)
    throw std::logic_error("SurfaceIPCCore prepared active pairs are missing. Call prepareForSurfacePositions() first.");
  computeSelfAll(preparedState_.positions, preparedState_.selfPairs, topology_.numVerts, dhat, kappa, eps_ee, energy, grad, hess);
  if (!obstacles_.empty()) {
    double extEnergy = 0.0;
    computeExternalAll(
      preparedState_.positions, obstacles_, preparedState_.externalPairs, topology_.numVerts, dhat_external, kappa, eps_ee, extEnergy, grad, hess);
    energy += extEnergy;
  }
  if (auto logger = Logging::lgr(); logger)
    SPDLOG_LOGGER_INFO(logger, "# nonzeros in Hessian: {}", hess.nonZeros());
}

// =========================================================================
//  Obstacle (external) registration
// =========================================================================
int32_t SurfaceIPCCore::addObstacleSurface(std::shared_ptr<ObstacleSurface> obs)
{
  if (!obs)
    throw std::invalid_argument("SurfaceIPCCore::addObstacleSurface: obstacle must not be null.");

  int32_t id = static_cast<int32_t>(obstacles_.size());
  obstacles_.push_back(std::move(obs));
  obstacles_.back()->setObjectId(id);
  preparedState_.clear();
  return id;
}

void SurfaceIPCCore::clearObstacleSurfaces()
{
  obstacles_.clear();
  preparedState_.clear();
}

void SurfaceIPCCore::updateObstacleStage(double tStart, double tEnd)
{
  for (auto &obs : obstacles_)
    obs->update(tStart, tEnd);
  preparedState_.clear();
}

}  // namespace CIPC
}  // namespace Contact
}  // namespace pgo
