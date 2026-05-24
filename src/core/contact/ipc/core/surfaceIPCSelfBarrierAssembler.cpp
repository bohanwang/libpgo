#include "surfaceIPCSelfBarrierAssembler.h"
#include "surfaceIPCBarrierKernels.h"

#include "scopedProfileSection.h"
#include "ipc/profiling/surfaceIPCProfiling.h"

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>

#include <atomic>
#include <cstdint>
#include <functional>
#include <stdexcept>
#include <vector>

namespace pgo {
namespace Contact {
namespace IPC {
using namespace pgo::EigenSupport;

static V3d vtx(ConstRefVecXd x, int i)
{
  return x.segment<3>(3 * i);
}

static void scatterSelfGrad(const V12d &local, const int idx[4], RefVecXd grad)
{
  double *gdata = grad.data();
  for (int i = 0; i < 4; ++i)
    if (idx[i] >= 0)
      for (int d = 0; d < 3; ++d)
        std::atomic_ref<double>(gdata[3 * idx[i] + d])
          .fetch_add(local[3 * i + d], std::memory_order_relaxed);
}

static void scatterSelfHessian(int pairIdx, const M12d &localH, const int idx[4],
  std::vector<TripletD> &triplets)
{
  TripletD *base = triplets.data() + 144 * pairIdx;
  int k = 0;
  for (int i = 0; i < 4; ++i) {
    int ri = (idx[i] >= 0) ? 3 * idx[i] : 0;
    bool vi = (idx[i] >= 0);
    for (int j = 0; j < 4; ++j) {
      int cj = (idx[j] >= 0) ? 3 * idx[j] : 0;
      bool valid = vi && (idx[j] >= 0);
      for (int di = 0; di < 3; ++di)
        for (int dj = 0; dj < 3; ++dj, ++k)
          base[k] = TripletD(ri + di, cj + dj,
            valid ? localH(3 * i + di, 3 * j + dj) : 0.0);
    }
  }
}

// =========================================================================
//  Self Energy
// =========================================================================

double computeSelfEnergy(
  EigenSupport::ConstRefVecXd dynPos,
  const SelfPairSet &pairs,
  int numVerts,
  double dhat,
  double kappa,
  double eps_ee)
{
  (void)numVerts;
  double dhat2 = dhat * dhat;

  // PT pairs
  double ptEnergy = tbb::parallel_reduce(
    tbb::blocked_range<int>(0, (int)pairs.ptPairs.size()), 0.0,
    [&](const tbb::blocked_range<int> &range, double localE) {
      for (int i = range.begin(); i < range.end(); ++i) {
        auto &pair = pairs.ptPairs[i];
        auto k = barrier_kernels::pointTriangle(
          vtx(dynPos, pair.p), vtx(dynPos, pair.t0), vtx(dynPos, pair.t1), vtx(dynPos, pair.t2),
          pair.weight, dhat2, kappa, false, false);
        if (k.active)
          localE += k.energy;
      }
      return localE;
    },
    std::plus<double>());

  // EE pairs
  double eeEnergy = tbb::parallel_reduce(
    tbb::blocked_range<int>(0, (int)pairs.eePairs.size()), 0.0,
    [&](const tbb::blocked_range<int> &range, double localE) {
      for (int i = range.begin(); i < range.end(); ++i) {
        auto &pair = pairs.eePairs[i];
        auto k = barrier_kernels::edgeEdge(
          vtx(dynPos, pair.ea0), vtx(dynPos, pair.ea1), vtx(dynPos, pair.eb0), vtx(dynPos, pair.eb1),
          pair.weight, dhat2, kappa, eps_ee, false, false);
        if (k.active)
          localE += k.energy;
      }
      return localE;
    },
    std::plus<double>());

  return ptEnergy + eeEnergy;
}

// =========================================================================
//  Self Gradient
// =========================================================================

void computeSelfGradient(
  EigenSupport::ConstRefVecXd dynPos,
  const SelfPairSet &pairs,
  int numVerts,
  double dhat,
  double kappa,
  double eps_ee,
  EigenSupport::RefVecXd grad)
{
  int n = 3 * numVerts;
  if (grad.size() != n)
    throw std::runtime_error("Gradient vector has wrong size");

  grad.setZero();

  double dhat2 = dhat * dhat;

  // PT pairs
  tbb::parallel_for(
    tbb::blocked_range<int>(0, (int)pairs.ptPairs.size()),
    [&](const tbb::blocked_range<int> &range) {
      for (int i = range.begin(); i < range.end(); ++i) {
        auto &pair = pairs.ptPairs[i];
        auto k = barrier_kernels::pointTriangle(
          vtx(dynPos, pair.p), vtx(dynPos, pair.t0), vtx(dynPos, pair.t1), vtx(dynPos, pair.t2),
          pair.weight, dhat2, kappa, true, false);
        if (!k.active)
          continue;
        int idx[4] = { pair.p, pair.t0, pair.t1, pair.t2 };
        scatterSelfGrad(k.gradient, idx, grad);
      }
    });

  // EE pairs
  tbb::parallel_for(
    tbb::blocked_range<int>(0, (int)pairs.eePairs.size()),
    [&](const tbb::blocked_range<int> &range) {
      for (int i = range.begin(); i < range.end(); ++i) {
        auto &pair = pairs.eePairs[i];
        auto k = barrier_kernels::edgeEdge(
          vtx(dynPos, pair.ea0), vtx(dynPos, pair.ea1), vtx(dynPos, pair.eb0), vtx(dynPos, pair.eb1),
          pair.weight, dhat2, kappa, eps_ee, true, false);
        if (!k.active)
          continue;
        int idx[4] = { pair.ea0, pair.ea1, pair.eb0, pair.eb1 };
        scatterSelfGrad(k.gradient, idx, grad);
      }
    });
}

// =========================================================================
//  Self Hessian
// =========================================================================

void computeSelfHessian(
  EigenSupport::ConstRefVecXd dynPos,
  const SelfPairSet &pairs,
  int numVerts,
  double dhat,
  double kappa,
  double eps_ee,
  SpMatD &hess)
{
  int n = 3 * numVerts;
  int nPT = (int)pairs.ptPairs.size();
  int nEE = (int)pairs.eePairs.size();
  int totalPairs = nPT + nEE;

  std::vector<TripletD> triplets(144 * totalPairs, TripletD(0, 0, 0.0));

  double dhat2 = dhat * dhat;

  // PT pairs
  tbb::parallel_for(
    tbb::blocked_range<int>(0, nPT),
    [&](const tbb::blocked_range<int> &range) {
      for (int i = range.begin(); i < range.end(); ++i) {
        auto &pair = pairs.ptPairs[i];
        auto k = barrier_kernels::pointTriangle(
          vtx(dynPos, pair.p), vtx(dynPos, pair.t0), vtx(dynPos, pair.t1), vtx(dynPos, pair.t2),
          pair.weight, dhat2, kappa, false, true);
        if (!k.active)
          continue;
        int idx[4] = { pair.p, pair.t0, pair.t1, pair.t2 };
        scatterSelfHessian(i, k.hessian, idx, triplets);
      }
    });

  // EE pairs
  tbb::parallel_for(
    tbb::blocked_range<int>(0, nEE),
    [&](const tbb::blocked_range<int> &range) {
      for (int i = range.begin(); i < range.end(); ++i) {
        auto &pair = pairs.eePairs[i];
        auto k = barrier_kernels::edgeEdge(
          vtx(dynPos, pair.ea0), vtx(dynPos, pair.ea1), vtx(dynPos, pair.eb0), vtx(dynPos, pair.eb1),
          pair.weight, dhat2, kappa, eps_ee, false, true);
        if (!k.active)
          continue;
        int idx[4] = { pair.ea0, pair.ea1, pair.eb0, pair.eb1 };
        scatterSelfHessian(nPT + i, k.hessian, idx, triplets);
      }
    });

  hess.resize(n, n);
  hess.setFromTriplets(triplets.begin(), triplets.end());
}

// =========================================================================
//  Self Combined (energy + gradient + hessian in single pass)
// =========================================================================

void computeSelfAll(
  EigenSupport::ConstRefVecXd dynPos,
  const SelfPairSet &pairs,
  int numVerts,
  double dhat,
  double kappa,
  double eps_ee,
  double &energy,
  VXd &grad,
  SpMatD &hess)
{
  int n = 3 * numVerts;
  int nPT = (int)pairs.ptPairs.size();
  int nEE = (int)pairs.eePairs.size();
  int totalPairs = nPT + nEE;

  Profiling::ScopedProfileSection scopedProfile(SurfaceIPCProfileSections::kActiveSetSelfCombined);
  Profiling::recordProfileCounter(SurfaceIPCProfileSections::kActiveSetSelfPTPairCount, static_cast<std::uint64_t>(nPT));
  Profiling::recordProfileCounter(SurfaceIPCProfileSections::kActiveSetSelfEEPairCount, static_cast<std::uint64_t>(nEE));

  energy = 0.0;
  grad.setZero(n);
  std::vector<TripletD> triplets(144 * totalPairs, TripletD(0, 0, 0.0));

  double dhat2 = dhat * dhat;

  // ---- PT pairs ----
  double ptEnergy = 0.0;
  {
    Profiling::ScopedProfileSection ptProfile(SurfaceIPCProfileSections::kActiveSetSelfPTCombined);
    ptEnergy = tbb::parallel_reduce(
      tbb::blocked_range<int>(0, nPT), 0.0,
      [&](const tbb::blocked_range<int> &range, double localE) {
        for (int i = range.begin(); i < range.end(); ++i) {
          auto &pair = pairs.ptPairs[i];
          auto k = barrier_kernels::pointTriangle(
            vtx(dynPos, pair.p), vtx(dynPos, pair.t0), vtx(dynPos, pair.t1), vtx(dynPos, pair.t2),
            pair.weight, dhat2, kappa, true, true);
          if (!k.active)
            continue;
          localE += k.energy;
          int idx[4] = { pair.p, pair.t0, pair.t1, pair.t2 };
          scatterSelfGrad(k.gradient, idx, grad);
          scatterSelfHessian(i, k.hessian, idx, triplets);
        }
        return localE;
      },
      std::plus<double>());
  }

  // ---- EE pairs ----
  double eeEnergy = 0.0;
  {
    Profiling::ScopedProfileSection eeProfile(SurfaceIPCProfileSections::kActiveSetSelfEECombined);
    eeEnergy = tbb::parallel_reduce(
      tbb::blocked_range<int>(0, nEE), 0.0,
      [&](const tbb::blocked_range<int> &range, double localE) {
        for (int i = range.begin(); i < range.end(); ++i) {
          auto &pair = pairs.eePairs[i];
          auto k = barrier_kernels::edgeEdge(
            vtx(dynPos, pair.ea0), vtx(dynPos, pair.ea1), vtx(dynPos, pair.eb0), vtx(dynPos, pair.eb1),
            pair.weight, dhat2, kappa, eps_ee, true, true);
          if (!k.active)
            continue;
          localE += k.energy;
          int idx[4] = { pair.ea0, pair.ea1, pair.eb0, pair.eb1 };
          scatterSelfGrad(k.gradient, idx, grad);
          scatterSelfHessian(nPT + i, k.hessian, idx, triplets);
        }
        return localE;
      },
      std::plus<double>());
  }

  energy = ptEnergy + eeEnergy;

  hess.resize(n, n);
  hess.setFromTriplets(triplets.begin(), triplets.end());
}

}  // namespace IPC
}  // namespace Contact
}  // namespace pgo
