/*
copyright to Bohan Wang
*/

// Full implementation of Codimensional IPC collision handling.
// =============================================================================

#include "ipc/core/surfaceIPCCore.h"
#include "scopedProfileSection.h"
#include "ipc/broadPhase/surfaceIPCSelfBroadPhase.h"
#include "ipc/core/surfaceIPCBarrierAssembler.h"
#include "ipc/core/surfaceIPCMaxStep.h"
#include "ipc/profiling/surfaceIPCProfiling.h"

#include "pgoLogging.h"

#include <atomic>
#include <algorithm>
#include <cstdint>

namespace pgo {
namespace Contact {
namespace CIPC {
static constexpr double kSmallContactAlphaWarnThreshold = 1e-2;

void updateMinAtomic(std::atomic<double> &target, double value)
{
  double current = target.load(std::memory_order_relaxed);
  while (value < current && !target.compare_exchange_weak(current, value, std::memory_order_relaxed)) {
  }
}

SurfaceIPCCore::SurfaceIPCCore(const SurfaceIPCCore &other):
  dhat(other.dhat),
  kappa(other.kappa),
  eps_ee(other.eps_ee),
  slackness(other.slackness),
  topology_(other.topology_),
  contactClampCount_(other.contactClampCount_.load(std::memory_order_relaxed)),
  minContactFeasibleAlphaThisSolve_(other.minContactFeasibleAlphaThisSolve_.load(std::memory_order_relaxed)),
  ptPairs_(other.ptPairs_),
  eePairs_(other.eePairs_)
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
  contactClampCount_.store(other.contactClampCount_.load(std::memory_order_relaxed), std::memory_order_relaxed);
  minContactFeasibleAlphaThisSolve_.store(other.minContactFeasibleAlphaThisSolve_.load(std::memory_order_relaxed), std::memory_order_relaxed);
  ptPairs_ = other.ptPairs_;
  eePairs_ = other.eePairs_;
  return *this;
}

void SurfaceIPCCore::setParameters(const Parameters &params)
{
  dhat = params.dhat;
  kappa = params.kappa;
  eps_ee = params.eps_ee;
  slackness = params.slackness;
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

void SurfaceIPCCore::resetContactMaxStepStats() const
{
  contactClampCount_.store(0, std::memory_order_relaxed);
  minContactFeasibleAlphaThisSolve_.store(1.0, std::memory_order_relaxed);
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
void SurfaceIPCCore::findCollisionPairs(const VXd &positions) const
{
  SurfaceIPCSelfBroadPhase().buildPairs(topology_, positions, dhat, ptPairs_, eePairs_);
}

// =========================================================================
//  1)  Maximum step size  (CCD-based line search with spatial hashing)
// =========================================================================
double SurfaceIPCCore::computeMaxStepSize(EigenSupport::ConstRefVecXd x, EigenSupport::ConstRefVecXd dx) const
{
  const double alpha = SurfaceIPCMaxStep().compute(topology_, x, dx, dhat, slackness);
  const double clampedAlpha = std::max(alpha, 1e-12);
  updateMinAtomic(minContactFeasibleAlphaThisSolve_, clampedAlpha);

  if (clampedAlpha < 1.0) {
    const std::int64_t clampCount = contactClampCount_.fetch_add(1, std::memory_order_relaxed) + 1;

    if (clampedAlpha > 0.0 && clampedAlpha < kSmallContactAlphaWarnThreshold) {
      SPDLOG_LOGGER_WARN(Logging::lgr(),
        "IPC contact max step produced small contactFeasibleAlpha={} (contactClampCount={}, slackness={}).",
        clampedAlpha, clampCount, slackness);
    }

    if (auto logger = Logging::lgr(); logger && logger->should_log(spdlog::level::trace)) {
      SPDLOG_LOGGER_TRACE(logger,
        "IPC contact clamp: contactFeasibleAlpha={} contactClampCount={} slackness={}.",
        clampedAlpha, clampCount, slackness);
    }
  }

  return clampedAlpha;
}

// =========================================================================
//  2)  Energy
// =========================================================================
double SurfaceIPCCore::computeEnergy(EigenSupport::ConstRefVecXd pos) const
{
  Profiling::ScopedProfileSection scopedProfile(SurfaceIPCProfileSections::kEnergy);
  findCollisionPairs(VXd(pos));
  return SurfaceIPCBarrierAssembler().computeEnergy(pos, ptPairs_, eePairs_, topology_.numVerts, dhat, kappa, eps_ee);
}

// =========================================================================
//  2)  Gradient
// =========================================================================
void SurfaceIPCCore::computeGradient(EigenSupport::ConstRefVecXd pos, EigenSupport::RefVecXd grad) const
{
  Profiling::ScopedProfileSection scopedProfile(SurfaceIPCProfileSections::kGradient);
  findCollisionPairs(VXd(pos));
  SurfaceIPCBarrierAssembler().computeGradient(pos, ptPairs_, eePairs_, topology_.numVerts, dhat, kappa, eps_ee, grad);
}

// =========================================================================
//  2)  Sparse Hessian
// =========================================================================
void SurfaceIPCCore::computeHessian(EigenSupport::ConstRefVecXd pos, SpMatD &hess) const
{
  Profiling::ScopedProfileSection scopedProfile(SurfaceIPCProfileSections::kHessian);
  findCollisionPairs(VXd(pos));
  SurfaceIPCBarrierAssembler().computeHessian(pos, ptPairs_, eePairs_, topology_.numVerts, dhat, kappa, eps_ee, hess);
}

// =========================================================================
//  Combined computation (single broad-phase pass)
// =========================================================================
void SurfaceIPCCore::computeAll(EigenSupport::ConstRefVecXd x,
  double &energy, VXd &grad, SpMatD &hess) const
{
  Profiling::ScopedProfileSection scopedProfile(SurfaceIPCProfileSections::kCombined);
  findCollisionPairs(x);
  SurfaceIPCBarrierAssembler().computeAll(x, ptPairs_, eePairs_, topology_.numVerts, dhat, kappa, eps_ee, energy, grad, hess);
}

}  // namespace CIPC
}  // namespace Contact
}  // namespace pgo
