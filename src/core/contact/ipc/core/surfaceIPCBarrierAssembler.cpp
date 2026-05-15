/*
copyright to Bohan Wang
*/

#include "surfaceIPCBarrierAssembler.h"

#include "../geometry/ipcBarrier.h"
#include "../geometry/ipcDistancePrimitives.h"
#include "../geometry/ipcHessianProjection.h"

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>

#include <atomic>
#include <functional>
#include <stdexcept>
#include <vector>

namespace pgo {
namespace Contact {
namespace CIPC {
using namespace pgo::EigenSupport;

static V3d vtx(ConstRefVecXd x, int i)
{
  return x.segment<3>(3 * i);
}

// =========================================================================
//  2)  Energy
// =========================================================================
double SurfaceIPCBarrierAssembler::computeEnergy(
  EigenSupport::ConstRefVecXd pos,
  const std::vector<PTPair> &ptPairs,
  const std::vector<EEPair> &eePairs,
  int numVerts,
  double dhat,
  double kappa,
  double eps_ee) const
{
  (void)numVerts;
  double ee_eps = eps_ee;
  double dhat2 = dhat * dhat;

  // PT pairs
  double ptEnergy = tbb::parallel_reduce(
    tbb::blocked_range<int>(0, (int)ptPairs.size()), 0.0,
    [&](const tbb::blocked_range<int> &range, double localE) {
      for (int i = range.begin(); i < range.end(); ++i) {
        auto &pair = ptPairs[i];
        V3d p = vtx(pos, pair.p);
        V3d t0 = vtx(pos, pair.t0);
        V3d t1 = vtx(pos, pair.t1);
        V3d t2 = vtx(pos, pair.t2);
        double d2 = distance::computePTSqDist(p, t0, t1, t2);
        if (d2 < dhat2 && d2 > 0.0)
          localE += pair.weight * kappa * barrier::b(d2, dhat2);
      }
      return localE;
    },
    std::plus<double>());

  // EE pairs
  double eeEnergy = tbb::parallel_reduce(
    tbb::blocked_range<int>(0, (int)eePairs.size()), 0.0,
    [&](const tbb::blocked_range<int> &range, double localE) {
      for (int i = range.begin(); i < range.end(); ++i) {
        auto &pair = eePairs[i];
        V3d ea0 = vtx(pos, pair.ea0);
        V3d ea1 = vtx(pos, pair.ea1);
        V3d eb0 = vtx(pos, pair.eb0);
        V3d eb1 = vtx(pos, pair.eb1);
        double d2 = distance::computeEESqDist(ea0, ea1, eb0, eb1);
        if (d2 < dhat2 && d2 > 0.0) {
          double m = 1.0;
          if (ee_eps > 0.0)
            m = distance::eeMollifier(ea0, ea1, eb0, eb1, ee_eps);
          localE += pair.weight * kappa * m * barrier::b(d2, dhat2);
        }
      }
      return localE;
    },
    std::plus<double>());

  return ptEnergy + eeEnergy;
}

// =========================================================================
//  2)  Gradient
// =========================================================================
void SurfaceIPCBarrierAssembler::computeGradient(
  EigenSupport::ConstRefVecXd pos,
  const std::vector<PTPair> &ptPairs,
  const std::vector<EEPair> &eePairs,
  int numVerts,
  double dhat,
  double kappa,
  double eps_ee,
  EigenSupport::RefVecXd grad) const
{
  int n = 3 * numVerts;
  if (grad.size() != n)
    throw std::runtime_error("Gradient vector has wrong size");

  grad.setZero();

  // Atomic scatter: accumulate a 12-vector into global gradient
  auto scatter = [&](const V12d &local, const int idx[4]) {
    double *gdata = grad.data();
    for (int i = 0; i < 4; ++i)
      if (idx[i] >= 0)
        for (int d = 0; d < 3; ++d)
          std::atomic_ref<double>(gdata[3 * idx[i] + d])
            .fetch_add(local[3 * i + d], std::memory_order_relaxed);
  };

  double dhat2 = dhat * dhat;

  // PT pairs
  tbb::parallel_for(
    tbb::blocked_range<int>(0, (int)ptPairs.size()),
    [&](const tbb::blocked_range<int> &range) {
      for (int i = range.begin(); i < range.end(); ++i) {
        auto &pair = ptPairs[i];
        V3d p = vtx(pos, pair.p);
        V3d t0 = vtx(pos, pair.t0);
        V3d t1 = vtx(pos, pair.t1);
        V3d t2 = vtx(pos, pair.t2);

        double d2 = distance::computePTSqDist(p, t0, t1, t2);
        if (d2 >= dhat2 || d2 <= 0.0)
          continue;

        V12d gd2 = distance::computePTSqDistGrad(p, t0, t1, t2);
        double coeff = pair.weight * kappa * barrier::dbds(d2, dhat2);
        V12d gE = coeff * gd2;

        int idx[4] = { pair.p, pair.t0, pair.t1, pair.t2 };
        scatter(gE, idx);
      }
    });

  // EE pairs
  double ee_eps = eps_ee;
  tbb::parallel_for(
    tbb::blocked_range<int>(0, (int)eePairs.size()),
    [&](const tbb::blocked_range<int> &range) {
      for (int i = range.begin(); i < range.end(); ++i) {
        auto &pair = eePairs[i];
        V3d ea0 = vtx(pos, pair.ea0);
        V3d ea1 = vtx(pos, pair.ea1);
        V3d eb0 = vtx(pos, pair.eb0);
        V3d eb1 = vtx(pos, pair.eb1);

        double d2 = distance::computeEESqDist(ea0, ea1, eb0, eb1);
        if (d2 >= dhat2 || d2 <= 0.0)
          continue;

        double wk = pair.weight * kappa;
        V12d gd2 = distance::computeEESqDistGrad(ea0, ea1, eb0, eb1);
        double b_val = barrier::b(d2, dhat2);
        double dbv = barrier::dbds(d2, dhat2);

        V12d gE;
        if (ee_eps > 0.0) {
          double m = distance::eeMollifier(ea0, ea1, eb0, eb1, ee_eps);
          V12d gm = distance::eeMollifierGrad(ea0, ea1, eb0, eb1, ee_eps);
          gE = wk * (gm * b_val + m * dbv * gd2);
        }
        else {
          gE = wk * dbv * gd2;
        }

        int idx[4] = { pair.ea0, pair.ea1, pair.eb0, pair.eb1 };
        scatter(gE, idx);
      }
    });

}

// =========================================================================
//  2)  Sparse Hessian
// =========================================================================
void SurfaceIPCBarrierAssembler::computeHessian(
  EigenSupport::ConstRefVecXd pos,
  const std::vector<PTPair> &ptPairs,
  const std::vector<EEPair> &eePairs,
  int numVerts,
  double dhat,
  double kappa,
  double eps_ee,
  SpMatD &hess) const
{
  int n = 3 * numVerts;
  int nPT = (int)ptPairs.size();
  int nEE = (int)eePairs.size();
  int totalPairs = nPT + nEE;

  // Preallocate exactly 144 triplet slots per pair (12x12 block).
  // Each pair writes to its own contiguous slice — no conflicts.
  std::vector<TripletD> triplets(144 * totalPairs, TripletD(0, 0, 0.0));

  // Write 144 triplets from a 12x12 local hessian into a preallocated slice.
  auto scatterH = [&](int pairIdx, const M12d &localH, const int idx[4]) {
    TripletD *base = triplets.data() + 144 * pairIdx;
    int k = 0;
    for (int i = 0; i < 4; ++i) {
      int ri = (idx[i] >= 0) ? 3 * idx[i] : 0;
      bool vi = (idx[i] >= 0);
      for (int j = 0; j < 4; ++j) {
        int cj = (idx[j] >= 0) ? 3 * idx[j] : 0;
        bool valid = vi && (idx[j] >= 0);

        for (int di = 0; di < 3; ++di) {
          for (int dj = 0; dj < 3; ++dj, ++k) {
            base[k] = TripletD(ri + di, cj + dj,
              valid ? localH(3 * i + di, 3 * j + dj) : 0.0);
          }
        }
      }
    }
  };

  double dhat2 = dhat * dhat;

  // PT pairs
  tbb::parallel_for(
    tbb::blocked_range<int>(0, nPT),
    [&](const tbb::blocked_range<int> &range) {
      for (int i = range.begin(); i < range.end(); ++i) {
        auto &pair = ptPairs[i];
        V3d p = vtx(pos, pair.p);
        V3d t0 = vtx(pos, pair.t0);
        V3d t1 = vtx(pos, pair.t1);
        V3d t2 = vtx(pos, pair.t2);

        double d2 = distance::computePTSqDist(p, t0, t1, t2);
        if (d2 >= dhat2 || d2 <= 0.0)
          continue;  // slot already zeroed

        V12d gd2 = distance::computePTSqDistGrad(p, t0, t1, t2);
        M12d Hd2 = distance::computePTSqDistHess(p, t0, t1, t2);
        double wk = pair.weight * kappa;
        double gp = barrier::dbds(d2, dhat2);
        double gpp = barrier::d2bds2(d2, dhat2);
        M12d localH = wk * (gpp * gd2 * gd2.transpose() + gp * Hd2);
        localH = projectToPSD(localH);

        int idx[4] = { pair.p, pair.t0, pair.t1, pair.t2 };
        scatterH(i, localH, idx);
      }
    });

  // EE pairs
  double ee_eps = eps_ee;
  tbb::parallel_for(
    tbb::blocked_range<int>(0, nEE),
    [&](const tbb::blocked_range<int> &range) {
      for (int i = range.begin(); i < range.end(); ++i) {
        auto &pair = eePairs[i];
        V3d ea0 = vtx(pos, pair.ea0);
        V3d ea1 = vtx(pos, pair.ea1);
        V3d eb0 = vtx(pos, pair.eb0);
        V3d eb1 = vtx(pos, pair.eb1);

        double d2 = distance::computeEESqDist(ea0, ea1, eb0, eb1);
        if (d2 >= dhat2 || d2 <= 0.0)
          continue;  // slot already zeroed

        double wk = pair.weight * kappa;
        V12d gd2 = distance::computeEESqDistGrad(ea0, ea1, eb0, eb1);
        M12d Hd2 = distance::computeEESqDistHess(ea0, ea1, eb0, eb1);
        double gp = barrier::dbds(d2, dhat2);
        double gpp = barrier::d2bds2(d2, dhat2);

        M12d localH;
        if (ee_eps > 0.0) {
          double m = distance::eeMollifier(ea0, ea1, eb0, eb1, ee_eps);
          V12d gm = distance::eeMollifierGrad(ea0, ea1, eb0, eb1, ee_eps);
          M12d Hm = distance::eeMollifierHess(ea0, ea1, eb0, eb1, ee_eps);
          double bv = barrier::b(d2, dhat2);
          double dbv = barrier::dbds(d2, dhat2);
          V12d gb = dbv * gd2;
          localH = wk * (bv * Hm + gm * gb.transpose() + gb * gm.transpose() + m * (gpp * gd2 * gd2.transpose() + gp * Hd2));
        }
        else {
          localH = wk * (gpp * gd2 * gd2.transpose() + gp * Hd2);
        }
        localH = projectToPSD(localH);

        int idx[4] = { pair.ea0, pair.ea1, pair.eb0, pair.eb1 };
        scatterH(nPT + i, localH, idx);
      }
    });

  hess.resize(n, n);
  hess.setFromTriplets(triplets.begin(), triplets.end());
}

// =========================================================================
//  Combined computation (single broad-phase pass)
// =========================================================================
void SurfaceIPCBarrierAssembler::computeAll(
  EigenSupport::ConstRefVecXd x,
  const std::vector<PTPair> &ptPairs,
  const std::vector<EEPair> &eePairs,
  int numVerts,
  double dhat,
  double kappa,
  double eps_ee,
  double &energy,
  VXd &grad,
  SpMatD &hess) const
{
  int n = 3 * numVerts;
  int nPT = (int)ptPairs.size();
  int nEE = (int)eePairs.size();
  int totalPairs = nPT + nEE;

  energy = 0.0;
  grad.setZero(n);
  std::vector<TripletD> triplets(144 * totalPairs, TripletD(0, 0, 0.0));

  // Atomic scatter for gradient
  auto scatterG = [&](const V12d &local, const int idx[4]) {
    double *gdata = grad.data();
    for (int i = 0; i < 4; ++i)
      if (idx[i] >= 0)
        for (int d = 0; d < 3; ++d)
          std::atomic_ref<double>(gdata[3 * idx[i] + d])
            .fetch_add(local[3 * i + d], std::memory_order_relaxed);
  };

  // Preallocated scatter for hessian
  auto scatterH = [&](int pairIdx, const M12d &localH, const int idx[4]) {
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
  };

  double ee_eps = eps_ee;
  double dhat2 = dhat * dhat;

  // ---- PT pairs ----
  double ptEnergy = tbb::parallel_reduce(
    tbb::blocked_range<int>(0, nPT), 0.0,
    [&](const tbb::blocked_range<int> &range, double localE) {
      for (int i = range.begin(); i < range.end(); ++i) {
        auto &pair = ptPairs[i];
        V3d p = vtx(x, pair.p);
        V3d t0 = vtx(x, pair.t0);
        V3d t1 = vtx(x, pair.t1);
        V3d t2 = vtx(x, pair.t2);

        double d2 = distance::computePTSqDist(p, t0, t1, t2);
        if (d2 >= dhat2 || d2 <= 0.0)
          continue;  // triplet slot already zeroed

        double wk = pair.weight * kappa;

        // energy
        localE += wk * barrier::b(d2, dhat2);

        // gradient
        V12d gd2 = distance::computePTSqDistGrad(p, t0, t1, t2);
        double gp = barrier::dbds(d2, dhat2);
        V12d gE = wk * gp * gd2;
        int idx[4] = { pair.p, pair.t0, pair.t1, pair.t2 };
        scatterG(gE, idx);

        // hessian
        M12d Hd2 = distance::computePTSqDistHess(p, t0, t1, t2);
        double gpp = barrier::d2bds2(d2, dhat2);
        M12d localH = wk * (gpp * gd2 * gd2.transpose() + gp * Hd2);
        localH = projectToPSD(localH);
        scatterH(i, localH, idx);
      }
      return localE;
    },
    std::plus<double>());

  // ---- EE pairs ----
  double eeEnergy = tbb::parallel_reduce(
    tbb::blocked_range<int>(0, nEE), 0.0,
    [&](const tbb::blocked_range<int> &range, double localE) {
      for (int i = range.begin(); i < range.end(); ++i) {
        auto &pair = eePairs[i];
        V3d ea0 = vtx(x, pair.ea0);
        V3d ea1 = vtx(x, pair.ea1);
        V3d eb0 = vtx(x, pair.eb0);
        V3d eb1 = vtx(x, pair.eb1);

        double d2 = distance::computeEESqDist(ea0, ea1, eb0, eb1);
        if (d2 >= dhat2 || d2 <= 0.0)
          continue;  // triplet slot already zeroed

        double wk = pair.weight * kappa;
        double m = 1.0;
        V12d gm;
        gm.setZero();
        M12d Hm;
        Hm.setZero();
        if (ee_eps > 0.0) {
          m = distance::eeMollifier(ea0, ea1, eb0, eb1, ee_eps);
          gm = distance::eeMollifierGrad(ea0, ea1, eb0, eb1, ee_eps);
          Hm = distance::eeMollifierHess(ea0, ea1, eb0, eb1, ee_eps);
        }

        double bv = barrier::b(d2, dhat2);

        // energy
        localE += wk * m * bv;

        // gradient
        V12d gd2 = distance::computeEESqDistGrad(ea0, ea1, eb0, eb1);
        double dbv = barrier::dbds(d2, dhat2);
        V12d gE = wk * (gm * bv + m * dbv * gd2);
        int idx[4] = { pair.ea0, pair.ea1, pair.eb0, pair.eb1 };
        scatterG(gE, idx);

        // hessian
        M12d Hd2 = distance::computeEESqDistHess(ea0, ea1, eb0, eb1);
        double gp = barrier::dbds(d2, dhat2);
        double gpp = barrier::d2bds2(d2, dhat2);
        V12d gb = dbv * gd2;
        M12d localH;
        if (ee_eps > 0.0) {
          localH = wk * (bv * Hm + gm * gb.transpose() + gb * gm.transpose() + m * (gpp * gd2 * gd2.transpose() + gp * Hd2));
        }
        else {
          localH = wk * (gpp * gd2 * gd2.transpose() + gp * Hd2);
        }
        localH = projectToPSD(localH);
        scatterH(nPT + i, localH, idx);
      }
      return localE;
    },
    std::plus<double>());

  energy = ptEnergy + eeEnergy;

  hess.resize(n, n);
  hess.setFromTriplets(triplets.begin(), triplets.end());
}

// =========================================================================
//  External pair barrier assembly (dynamic-only block scatter)
// =========================================================================

static const VXd &obsPositions(
  const std::vector<std::shared_ptr<ObstacleSurface>> &obstacles,
  int32_t obsId)
{
  return obstacles[static_cast<std::size_t>(obsId)]->currentPositions();
}

static V3d obsVtx(const VXd &pos, int i)
{
  return pos.segment<3>(3 * i);
}

static V3d dynVtx(ConstRefVecXd x, int i)
{
  return x.segment<3>(3 * i);
}

// Gradient scatter helpers — dynamic-only, static dispatch per pair type

static void scatterExternalPTGrad(const V12d &g_local, int dynVertex, RefVecXd g_surf)
{
  double *gd = g_surf.data();
  for (int d = 0; d < 3; ++d)
    std::atomic_ref<double>(gd[3 * dynVertex + d])
      .fetch_add(g_local[d], std::memory_order_relaxed);
}

static void scatterExternalTPGrad(const V12d &g_local, const std::array<int, 3> &dynTri, RefVecXd g_surf)
{
  double *gd = g_surf.data();
  for (int i = 0; i < 3; ++i)
    for (int d = 0; d < 3; ++d)
      std::atomic_ref<double>(gd[3 * dynTri[i] + d])
        .fetch_add(g_local[3 * (i + 1) + d], std::memory_order_relaxed);
}

static void scatterExternalEEGrad(const V12d &g_local, const std::array<int, 2> &dynEdge, RefVecXd g_surf)
{
  double *gd = g_surf.data();
  for (int i = 0; i < 2; ++i)
    for (int d = 0; d < 3; ++d)
      std::atomic_ref<double>(gd[3 * dynEdge[i] + d])
        .fetch_add(g_local[3 * i + d], std::memory_order_relaxed);
}

// Hessian scatter helpers — dynamic-only block, per-triplet preallocation

struct ExtHessianScatterState
{
  std::vector<TripletD> triplets;
  int pairCount = 0;
};

// PT: only 3x3 block at (0,0) → (dynVertex, dynVertex)
static void scatterExternalPTHessian(int tripletOffset, const M12d &localH, int dynVertex, ExtHessianScatterState &state)
{
  TripletD *base = state.triplets.data() + tripletOffset;
  int k = 0;
  int ri = 3 * dynVertex;
  int cj = 3 * dynVertex;
  for (int di = 0; di < 3; ++di)
    for (int dj = 0; dj < 3; ++dj, ++k)
      base[k] = TripletD(ri + di, cj + dj, localH(di, dj));
}

// TP: only 9x9 block at (1..3, 1..3) → (dynTri, dynTri)
static void scatterExternalTPHessian(int tripletOffset, const M12d &localH, const std::array<int, 3> &dynTri, ExtHessianScatterState &state)
{
  TripletD *base = state.triplets.data() + tripletOffset;
  int k = 0;
  for (int i = 0; i < 3; ++i) {
    int ri = 3 * dynTri[i];
    for (int j = 0; j < 3; ++j) {
      int cj = 3 * dynTri[j];
      for (int di = 0; di < 3; ++di)
        for (int dj = 0; dj < 3; ++dj, ++k)
          base[k] = TripletD(ri + di, cj + dj, localH(3 * (i + 1) + di, 3 * (j + 1) + dj));
    }
  }
}

// EE: only 6x6 block at (0..1, 0..1) → (dynEdge, dynEdge)
static void scatterExternalEEHessian(int tripletOffset, const M12d &localH, const std::array<int, 2> &dynEdge, ExtHessianScatterState &state)
{
  TripletD *base = state.triplets.data() + tripletOffset;
  int k = 0;
  for (int i = 0; i < 2; ++i) {
    int ri = 3 * dynEdge[i];
    for (int j = 0; j < 2; ++j) {
      int cj = 3 * dynEdge[j];
      for (int di = 0; di < 3; ++di)
        for (int dj = 0; dj < 3; ++dj, ++k)
          base[k] = TripletD(ri + di, cj + dj, localH(3 * i + di, 3 * j + dj));
    }
  }
}

double SurfaceIPCBarrierAssembler::computeExternalEnergy(
  ConstRefVecXd dynPos,
  const std::vector<std::shared_ptr<ObstacleSurface>> &obstacles,
  const std::vector<ExternalPTPair> &extPTPairs,
  const std::vector<ExternalTPPair> &extTPPairs,
  const std::vector<ExternalEEPair> &extEEPairs,
  double dhat,
  double kappa,
  double eps_ee) const
{
  (void)eps_ee;
  double dhat2 = dhat * dhat;

  double ptEnergy = tbb::parallel_reduce(
    tbb::blocked_range<int>(0, (int)extPTPairs.size()), 0.0,
    [&](const tbb::blocked_range<int> &range, double localE) {
      for (int i = range.begin(); i < range.end(); ++i) {
        auto &pair = extPTPairs[i];
        const VXd &obsP = obsPositions(obstacles, pair.obstacleObjectId);
        V3d p = dynVtx(dynPos, pair.dynVertex);
        V3d t0 = obsVtx(obsP, pair.obsTri[0]);
        V3d t1 = obsVtx(obsP, pair.obsTri[1]);
        V3d t2 = obsVtx(obsP, pair.obsTri[2]);
        double d2 = distance::computePTSqDist(p, t0, t1, t2);
        if (d2 < dhat2 && d2 > 0.0)
          localE += pair.weight * kappa * barrier::b(d2, dhat2);
      }
      return localE;
    },
    std::plus<double>());

  double tpEnergy = tbb::parallel_reduce(
    tbb::blocked_range<int>(0, (int)extTPPairs.size()), 0.0,
    [&](const tbb::blocked_range<int> &range, double localE) {
      for (int i = range.begin(); i < range.end(); ++i) {
        auto &pair = extTPPairs[i];
        const VXd &obsP = obsPositions(obstacles, pair.obstacleObjectId);
        V3d p = obsVtx(obsP, pair.obsVertex);
        V3d t0 = dynVtx(dynPos, pair.dynTri[0]);
        V3d t1 = dynVtx(dynPos, pair.dynTri[1]);
        V3d t2 = dynVtx(dynPos, pair.dynTri[2]);
        double d2 = distance::computePTSqDist(p, t0, t1, t2);
        if (d2 < dhat2 && d2 > 0.0)
          localE += pair.weight * kappa * barrier::b(d2, dhat2);
      }
      return localE;
    },
    std::plus<double>());

  double eeEnergy = tbb::parallel_reduce(
    tbb::blocked_range<int>(0, (int)extEEPairs.size()), 0.0,
    [&](const tbb::blocked_range<int> &range, double localE) {
      for (int i = range.begin(); i < range.end(); ++i) {
        auto &pair = extEEPairs[i];
        const VXd &obsP = obsPositions(obstacles, pair.obstacleObjectId);
        V3d ea0 = dynVtx(dynPos, pair.dynEdge[0]);
        V3d ea1 = dynVtx(dynPos, pair.dynEdge[1]);
        V3d eb0 = obsVtx(obsP, pair.obsEdge[0]);
        V3d eb1 = obsVtx(obsP, pair.obsEdge[1]);
        double d2 = distance::computeEESqDist(ea0, ea1, eb0, eb1);
        if (d2 < dhat2 && d2 > 0.0) {
          double m = 1.0;
          if (eps_ee > 0.0)
            m = distance::eeMollifier(ea0, ea1, eb0, eb1, eps_ee);
          localE += pair.weight * kappa * m * barrier::b(d2, dhat2);
        }
      }
      return localE;
    },
    std::plus<double>());

  return ptEnergy + tpEnergy + eeEnergy;
}

void SurfaceIPCBarrierAssembler::computeExternalGradient(
  ConstRefVecXd dynPos,
  const std::vector<std::shared_ptr<ObstacleSurface>> &obstacles,
  const std::vector<ExternalPTPair> &extPTPairs,
  const std::vector<ExternalTPPair> &extTPPairs,
  const std::vector<ExternalEEPair> &extEEPairs,
  int numDynVerts,
  double dhat,
  double kappa,
  double eps_ee,
  RefVecXd grad) const
{
  int n = 3 * numDynVerts;
  if (grad.size() != n)
    throw std::runtime_error("External gradient vector has wrong size");

  double dhat2 = dhat * dhat;

  // PT pairs
  tbb::parallel_for(
    tbb::blocked_range<int>(0, (int)extPTPairs.size()),
    [&](const tbb::blocked_range<int> &range) {
      for (int i = range.begin(); i < range.end(); ++i) {
        auto &pair = extPTPairs[i];
        const VXd &obsP = obsPositions(obstacles, pair.obstacleObjectId);
        V3d p = dynVtx(dynPos, pair.dynVertex);
        V3d t0 = obsVtx(obsP, pair.obsTri[0]);
        V3d t1 = obsVtx(obsP, pair.obsTri[1]);
        V3d t2 = obsVtx(obsP, pair.obsTri[2]);

        double d2 = distance::computePTSqDist(p, t0, t1, t2);
        if (d2 >= dhat2 || d2 <= 0.0)
          continue;

        V12d gd2 = distance::computePTSqDistGrad(p, t0, t1, t2);
        double coeff = pair.weight * kappa * barrier::dbds(d2, dhat2);
        V12d gE = coeff * gd2;
        scatterExternalPTGrad(gE, pair.dynVertex, grad);
      }
    });

  // TP pairs
  tbb::parallel_for(
    tbb::blocked_range<int>(0, (int)extTPPairs.size()),
    [&](const tbb::blocked_range<int> &range) {
      for (int i = range.begin(); i < range.end(); ++i) {
        auto &pair = extTPPairs[i];
        const VXd &obsP = obsPositions(obstacles, pair.obstacleObjectId);
        V3d p = obsVtx(obsP, pair.obsVertex);
        V3d t0 = dynVtx(dynPos, pair.dynTri[0]);
        V3d t1 = dynVtx(dynPos, pair.dynTri[1]);
        V3d t2 = dynVtx(dynPos, pair.dynTri[2]);

        double d2 = distance::computePTSqDist(p, t0, t1, t2);
        if (d2 >= dhat2 || d2 <= 0.0)
          continue;

        V12d gd2 = distance::computePTSqDistGrad(p, t0, t1, t2);
        double coeff = pair.weight * kappa * barrier::dbds(d2, dhat2);
        V12d gE = coeff * gd2;
        scatterExternalTPGrad(gE, pair.dynTri, grad);
      }
    });

  // EE pairs
  double ee_eps = eps_ee;
  tbb::parallel_for(
    tbb::blocked_range<int>(0, (int)extEEPairs.size()),
    [&](const tbb::blocked_range<int> &range) {
      for (int i = range.begin(); i < range.end(); ++i) {
        auto &pair = extEEPairs[i];
        const VXd &obsP = obsPositions(obstacles, pair.obstacleObjectId);
        V3d ea0 = dynVtx(dynPos, pair.dynEdge[0]);
        V3d ea1 = dynVtx(dynPos, pair.dynEdge[1]);
        V3d eb0 = obsVtx(obsP, pair.obsEdge[0]);
        V3d eb1 = obsVtx(obsP, pair.obsEdge[1]);

        double d2 = distance::computeEESqDist(ea0, ea1, eb0, eb1);
        if (d2 >= dhat2 || d2 <= 0.0)
          continue;

        double wk = pair.weight * kappa;
        V12d gd2 = distance::computeEESqDistGrad(ea0, ea1, eb0, eb1);
        double b_val = barrier::b(d2, dhat2);
        double dbv = barrier::dbds(d2, dhat2);

        V12d gE;
        if (ee_eps > 0.0) {
          double m = distance::eeMollifier(ea0, ea1, eb0, eb1, ee_eps);
          V12d gm = distance::eeMollifierGrad(ea0, ea1, eb0, eb1, ee_eps);
          gE = wk * (gm * b_val + m * dbv * gd2);
        }
        else {
          gE = wk * dbv * gd2;
        }

        scatterExternalEEGrad(gE, pair.dynEdge, grad);
      }
    });
}

void SurfaceIPCBarrierAssembler::computeExternalHessian(
  ConstRefVecXd dynPos,
  const std::vector<std::shared_ptr<ObstacleSurface>> &obstacles,
  const std::vector<ExternalPTPair> &extPTPairs,
  const std::vector<ExternalTPPair> &extTPPairs,
  const std::vector<ExternalEEPair> &extEEPairs,
  int numDynVerts,
  double dhat,
  double kappa,
  double eps_ee,
  SpMatD &hess) const
{
  int n = 3 * numDynVerts;
  int nPT = (int)extPTPairs.size();
  int nTP = (int)extTPPairs.size();
  int nEE = (int)extEEPairs.size();

  ExtHessianScatterState state;
  state.pairCount = nPT * 9 + nTP * 81 + nEE * 36;
  state.triplets.resize(state.pairCount, TripletD(0, 0, 0.0));

  double dhat2 = dhat * dhat;

  // PT pairs
  tbb::parallel_for(
    tbb::blocked_range<int>(0, nPT),
    [&](const tbb::blocked_range<int> &range) {
      for (int i = range.begin(); i < range.end(); ++i) {
        auto &pair = extPTPairs[i];
        const VXd &obsP = obsPositions(obstacles, pair.obstacleObjectId);
        V3d p = dynVtx(dynPos, pair.dynVertex);
        V3d t0 = obsVtx(obsP, pair.obsTri[0]);
        V3d t1 = obsVtx(obsP, pair.obsTri[1]);
        V3d t2 = obsVtx(obsP, pair.obsTri[2]);

        double d2 = distance::computePTSqDist(p, t0, t1, t2);
        if (d2 >= dhat2 || d2 <= 0.0)
          continue;

        V12d gd2 = distance::computePTSqDistGrad(p, t0, t1, t2);
        M12d Hd2 = distance::computePTSqDistHess(p, t0, t1, t2);
        double wk = pair.weight * kappa;
        double gp = barrier::dbds(d2, dhat2);
        double gpp = barrier::d2bds2(d2, dhat2);
        M12d localH = wk * (gpp * gd2 * gd2.transpose() + gp * Hd2);
        localH = projectToPSD(localH);
        scatterExternalPTHessian(9 * i, localH, pair.dynVertex, state);
      }
    });

  // TP pairs
  tbb::parallel_for(
    tbb::blocked_range<int>(0, nTP),
    [&](const tbb::blocked_range<int> &range) {
      for (int i = range.begin(); i < range.end(); ++i) {
        auto &pair = extTPPairs[i];
        const VXd &obsP = obsPositions(obstacles, pair.obstacleObjectId);
        V3d p = obsVtx(obsP, pair.obsVertex);
        V3d t0 = dynVtx(dynPos, pair.dynTri[0]);
        V3d t1 = dynVtx(dynPos, pair.dynTri[1]);
        V3d t2 = dynVtx(dynPos, pair.dynTri[2]);

        double d2 = distance::computePTSqDist(p, t0, t1, t2);
        if (d2 >= dhat2 || d2 <= 0.0)
          continue;

        V12d gd2 = distance::computePTSqDistGrad(p, t0, t1, t2);
        M12d Hd2 = distance::computePTSqDistHess(p, t0, t1, t2);
        double wk = pair.weight * kappa;
        double gp = barrier::dbds(d2, dhat2);
        double gpp = barrier::d2bds2(d2, dhat2);
        M12d localH = wk * (gpp * gd2 * gd2.transpose() + gp * Hd2);
        localH = projectToPSD(localH);
        scatterExternalTPHessian(9 * nPT + 81 * i, localH, pair.dynTri, state);
      }
    });

  // EE pairs
  double ee_eps = eps_ee;
  tbb::parallel_for(
    tbb::blocked_range<int>(0, nEE),
    [&](const tbb::blocked_range<int> &range) {
      for (int i = range.begin(); i < range.end(); ++i) {
        auto &pair = extEEPairs[i];
        const VXd &obsP = obsPositions(obstacles, pair.obstacleObjectId);
        V3d ea0 = dynVtx(dynPos, pair.dynEdge[0]);
        V3d ea1 = dynVtx(dynPos, pair.dynEdge[1]);
        V3d eb0 = obsVtx(obsP, pair.obsEdge[0]);
        V3d eb1 = obsVtx(obsP, pair.obsEdge[1]);

        double d2 = distance::computeEESqDist(ea0, ea1, eb0, eb1);
        if (d2 >= dhat2 || d2 <= 0.0)
          continue;

        double wk = pair.weight * kappa;
        V12d gd2 = distance::computeEESqDistGrad(ea0, ea1, eb0, eb1);
        M12d Hd2 = distance::computeEESqDistHess(ea0, ea1, eb0, eb1);
        double gp = barrier::dbds(d2, dhat2);
        double gpp = barrier::d2bds2(d2, dhat2);

        M12d localH;
        if (ee_eps > 0.0) {
          double m = distance::eeMollifier(ea0, ea1, eb0, eb1, ee_eps);
          V12d gm = distance::eeMollifierGrad(ea0, ea1, eb0, eb1, ee_eps);
          M12d Hm = distance::eeMollifierHess(ea0, ea1, eb0, eb1, ee_eps);
          double bv = barrier::b(d2, dhat2);
          double dbv = barrier::dbds(d2, dhat2);
          V12d gb = dbv * gd2;
          localH = wk * (bv * Hm + gm * gb.transpose() + gb * gm.transpose() + m * (gpp * gd2 * gd2.transpose() + gp * Hd2));
        }
        else {
          localH = wk * (gpp * gd2 * gd2.transpose() + gp * Hd2);
        }
        localH = projectToPSD(localH);
        scatterExternalEEHessian(9 * nPT + 81 * nTP + 36 * i, localH, pair.dynEdge, state);
      }
    });

  // Merge external Hessian into existing (additive)
  SpMatD hessExt(n, n);
  hessExt.setFromTriplets(state.triplets.begin(), state.triplets.end());
  if (hess.nonZeros() == 0) {
    hess = std::move(hessExt);
  }
  else {
    hess += hessExt;
  }
}

void SurfaceIPCBarrierAssembler::computeExternalAll(
  ConstRefVecXd dynPos,
  const std::vector<std::shared_ptr<ObstacleSurface>> &obstacles,
  const std::vector<ExternalPTPair> &extPTPairs,
  const std::vector<ExternalTPPair> &extTPPairs,
  const std::vector<ExternalEEPair> &extEEPairs,
  int numDynVerts,
  double dhat,
  double kappa,
  double eps_ee,
  double &energy,
  VXd &grad,
  SpMatD &hess) const
{
  int n = 3 * numDynVerts;
  int nPT = (int)extPTPairs.size();
  int nTP = (int)extTPPairs.size();
  int nEE = (int)extEEPairs.size();

  energy = 0.0;
  if (grad.size() != n)
    grad.setZero(n);

  ExtHessianScatterState hState;
  hState.pairCount = nPT * 9 + nTP * 81 + nEE * 36;
  hState.triplets.resize(hState.pairCount, TripletD(0, 0, 0.0));

  double dhat2 = dhat * dhat;

  // PT pairs
  double ptEnergy = tbb::parallel_reduce(
    tbb::blocked_range<int>(0, nPT), 0.0,
    [&](const tbb::blocked_range<int> &range, double localE) {
      for (int i = range.begin(); i < range.end(); ++i) {
        auto &pair = extPTPairs[i];
        const VXd &obsP = obsPositions(obstacles, pair.obstacleObjectId);
        V3d p = dynVtx(dynPos, pair.dynVertex);
        V3d t0 = obsVtx(obsP, pair.obsTri[0]);
        V3d t1 = obsVtx(obsP, pair.obsTri[1]);
        V3d t2 = obsVtx(obsP, pair.obsTri[2]);

        double d2 = distance::computePTSqDist(p, t0, t1, t2);
        if (d2 >= dhat2 || d2 <= 0.0)
          continue;

        double wk = pair.weight * kappa;
        localE += wk * barrier::b(d2, dhat2);

        V12d gd2 = distance::computePTSqDistGrad(p, t0, t1, t2);
        double gp = barrier::dbds(d2, dhat2);
        V12d gE = wk * gp * gd2;
        scatterExternalPTGrad(gE, pair.dynVertex, grad);

        M12d Hd2 = distance::computePTSqDistHess(p, t0, t1, t2);
        double gpp = barrier::d2bds2(d2, dhat2);
        M12d localH = wk * (gpp * gd2 * gd2.transpose() + gp * Hd2);
        localH = projectToPSD(localH);
        scatterExternalPTHessian(9 * i, localH, pair.dynVertex, hState);
      }
      return localE;
    },
    std::plus<double>());

  // TP pairs
  double tpEnergy = tbb::parallel_reduce(
    tbb::blocked_range<int>(0, nTP), 0.0,
    [&](const tbb::blocked_range<int> &range, double localE) {
      for (int i = range.begin(); i < range.end(); ++i) {
        auto &pair = extTPPairs[i];
        const VXd &obsP = obsPositions(obstacles, pair.obstacleObjectId);
        V3d p = obsVtx(obsP, pair.obsVertex);
        V3d t0 = dynVtx(dynPos, pair.dynTri[0]);
        V3d t1 = dynVtx(dynPos, pair.dynTri[1]);
        V3d t2 = dynVtx(dynPos, pair.dynTri[2]);

        double d2 = distance::computePTSqDist(p, t0, t1, t2);
        if (d2 >= dhat2 || d2 <= 0.0)
          continue;

        double wk = pair.weight * kappa;
        localE += wk * barrier::b(d2, dhat2);

        V12d gd2 = distance::computePTSqDistGrad(p, t0, t1, t2);
        double gp = barrier::dbds(d2, dhat2);
        V12d gE = wk * gp * gd2;
        scatterExternalTPGrad(gE, pair.dynTri, grad);

        M12d Hd2 = distance::computePTSqDistHess(p, t0, t1, t2);
        double gpp = barrier::d2bds2(d2, dhat2);
        M12d localH = wk * (gpp * gd2 * gd2.transpose() + gp * Hd2);
        localH = projectToPSD(localH);
        scatterExternalTPHessian(9 * nPT + 81 * i, localH, pair.dynTri, hState);
      }
      return localE;
    },
    std::plus<double>());

  // EE pairs
  double ee_eps = eps_ee;
  double eeEnergy = tbb::parallel_reduce(
    tbb::blocked_range<int>(0, nEE), 0.0,
    [&](const tbb::blocked_range<int> &range, double localE) {
      for (int i = range.begin(); i < range.end(); ++i) {
        auto &pair = extEEPairs[i];
        const VXd &obsP = obsPositions(obstacles, pair.obstacleObjectId);
        V3d ea0 = dynVtx(dynPos, pair.dynEdge[0]);
        V3d ea1 = dynVtx(dynPos, pair.dynEdge[1]);
        V3d eb0 = obsVtx(obsP, pair.obsEdge[0]);
        V3d eb1 = obsVtx(obsP, pair.obsEdge[1]);

        double d2 = distance::computeEESqDist(ea0, ea1, eb0, eb1);
        if (d2 >= dhat2 || d2 <= 0.0)
          continue;

        double wk = pair.weight * kappa;
        double m = 1.0;
        V12d gm;
        gm.setZero();
        M12d Hm;
        Hm.setZero();
        if (ee_eps > 0.0) {
          m = distance::eeMollifier(ea0, ea1, eb0, eb1, ee_eps);
          gm = distance::eeMollifierGrad(ea0, ea1, eb0, eb1, ee_eps);
          Hm = distance::eeMollifierHess(ea0, ea1, eb0, eb1, ee_eps);
        }

        double bv = barrier::b(d2, dhat2);
        localE += wk * m * bv;

        V12d gd2 = distance::computeEESqDistGrad(ea0, ea1, eb0, eb1);
        double dbv = barrier::dbds(d2, dhat2);
        V12d gE = wk * (gm * bv + m * dbv * gd2);
        scatterExternalEEGrad(gE, pair.dynEdge, grad);

        M12d Hd2 = distance::computeEESqDistHess(ea0, ea1, eb0, eb1);
        double gp = barrier::dbds(d2, dhat2);
        double gpp = barrier::d2bds2(d2, dhat2);
        V12d gb = dbv * gd2;
        M12d localH;
        if (ee_eps > 0.0) {
          localH = wk * (bv * Hm + gm * gb.transpose() + gb * gm.transpose() + m * (gpp * gd2 * gd2.transpose() + gp * Hd2));
        }
        else {
          localH = wk * (gpp * gd2 * gd2.transpose() + gp * Hd2);
        }
        localH = projectToPSD(localH);
        scatterExternalEEHessian(9 * nPT + 81 * nTP + 36 * i, localH, pair.dynEdge, hState);
      }
      return localE;
    },
    std::plus<double>());

  energy = ptEnergy + tpEnergy + eeEnergy;

  SpMatD hessExt(n, n);
  hessExt.setFromTriplets(hState.triplets.begin(), hState.triplets.end());
  if (hess.nonZeros() == 0) {
    hess = std::move(hessExt);
  }
  else {
    hess += hessExt;
  }
}

}  // namespace CIPC
}  // namespace Contact
}  // namespace pgo
