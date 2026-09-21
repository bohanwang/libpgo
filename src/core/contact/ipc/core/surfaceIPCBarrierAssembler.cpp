/*
copyright to Bohan Wang
*/

#include "surfaceIPCBarrierAssembler.h"

#include "../geometry/ipcBarrier.h"
#include "../geometry/ipcDistancePrimitives.h"
#include "../geometry/ipcHessianProjection.h"
#include "pgoLogging.h"

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>

#include <algorithm>
#include <atomic>
#include <functional>
#include <stdexcept>
#include <vector>

namespace pgo {
namespace Contact {
namespace CIPC {
using namespace pgo::EigenSupport;

namespace {

V3d vtx(ConstRefVecXd x, int i)
{
  return x.segment<3>(3 * i);
}

// Per-pair barrier kernels. Each optional output is computed only when its
// pointer is non-null. They return false for inactive pairs (outside the
// barrier support or at zero distance), leaving the outputs untouched.
bool evaluatePT(ConstRefVecXd pos, const PTPair &pair, double dhat2, double kappa,
  bool projectHessianToPSD, double *energy, V12d *grad, M12d *hess)
{
  const V3d p = vtx(pos, pair.p);
  const V3d t0 = vtx(pos, pair.t0);
  const V3d t1 = vtx(pos, pair.t1);
  const V3d t2 = vtx(pos, pair.t2);

  const double d2 = distance::computePTSqDist(p, t0, t1, t2);
  if (d2 >= dhat2 || d2 <= 0.0)
    return false;

  const double wk = pair.weight * kappa;
  if (energy)
    *energy = wk * barrier::b(d2, dhat2);

  if (grad || hess) {
    const V12d gd2 = distance::computePTSqDistGrad(p, t0, t1, t2);
    const double gp = barrier::dbds(d2, dhat2);
    if (grad)
      *grad = wk * gp * gd2;
    if (hess) {
      const M12d Hd2 = distance::computePTSqDistHess(p, t0, t1, t2);
      const double gpp = barrier::d2bds2(d2, dhat2);
      M12d localH = wk * (gpp * gd2 * gd2.transpose() + gp * Hd2);
      if (projectHessianToPSD)
        localH = projectToPSD(localH);
      *hess = localH;
    }
  }
  return true;
}

bool evaluateEE(ConstRefVecXd pos, const EEPair &pair, double dhat2, double kappa, double ee_eps,
  bool projectHessianToPSD, double *energy, V12d *grad, M12d *hess)
{
  const V3d ea0 = vtx(pos, pair.ea0);
  const V3d ea1 = vtx(pos, pair.ea1);
  const V3d eb0 = vtx(pos, pair.eb0);
  const V3d eb1 = vtx(pos, pair.eb1);

  const double d2 = distance::computeEESqDist(ea0, ea1, eb0, eb1);
  if (d2 >= dhat2 || d2 <= 0.0)
    return false;

  const double wk = pair.weight * kappa;
  const double bv = barrier::b(d2, dhat2);
  double m = 1.0;
  if (ee_eps > 0.0)
    m = distance::eeMollifier(ea0, ea1, eb0, eb1, ee_eps);

  if (energy)
    *energy = wk * m * bv;

  if (grad || hess) {
    const V12d gd2 = distance::computeEESqDistGrad(ea0, ea1, eb0, eb1);
    const double dbv = barrier::dbds(d2, dhat2);
    V12d gm;
    gm.setZero();
    if (ee_eps > 0.0)
      gm = distance::eeMollifierGrad(ea0, ea1, eb0, eb1, ee_eps);

    if (grad) {
      if (ee_eps > 0.0)
        *grad = wk * (gm * bv + m * dbv * gd2);
      else
        *grad = wk * dbv * gd2;
    }
    if (hess) {
      const M12d Hd2 = distance::computeEESqDistHess(ea0, ea1, eb0, eb1);
      const double gpp = barrier::d2bds2(d2, dhat2);
      M12d localH;
      if (ee_eps > 0.0) {
        const M12d Hm = distance::eeMollifierHess(ea0, ea1, eb0, eb1, ee_eps);
        const V12d gb = dbv * gd2;
        localH = wk * (bv * Hm + gm * gb.transpose() + gb * gm.transpose() + m * (gpp * gd2 * gd2.transpose() + dbv * Hd2));
      }
      else {
        localH = wk * (gpp * gd2 * gd2.transpose() + dbv * Hd2);
      }
      if (projectHessianToPSD)
        localH = projectToPSD(localH);
      *hess = localH;
    }
  }
  return true;
}

void ptIndices(const PTPair &pair, int idx[4])
{
  idx[0] = pair.p;
  idx[1] = pair.t0;
  idx[2] = pair.t1;
  idx[3] = pair.t2;
}

void eeIndices(const EEPair &pair, int idx[4])
{
  idx[0] = pair.ea0;
  idx[1] = pair.ea1;
  idx[2] = pair.eb0;
  idx[3] = pair.eb1;
}

// Non-deterministic scatter: accumulate a 12-vector into the global gradient
// with atomic adds (summation order follows thread scheduling).
void scatterAtomic(double *gdata, const V12d &local, const int idx[4])
{
  for (int i = 0; i < 4; ++i)
    for (int d = 0; d < 3; ++d)
      std::atomic_ref<double>(gdata[3 * idx[i] + d]).fetch_add(local[3 * i + d], std::memory_order_relaxed);
}

// Deterministic gather: each vertex sums its stored blocks in map order.
void gatherBlocks(const SurfaceIPCGradientGatherMap &gather, const std::vector<double> &blocks, RefVecXd grad)
{
  const int numVerts = gather.numVertices();
  tbb::parallel_for(tbb::blocked_range<int>(0, numVerts), [&](const tbb::blocked_range<int> &range) {
    for (int v = range.begin(); v < range.end(); ++v) {
      V3d s = V3d::Zero();
      for (std::ptrdiff_t k = gather.start[v]; k < gather.start[v + 1]; ++k)
        s += Eigen::Map<const V3d>(blocks.data() + 3 * gather.blocks[k]);
      grad.segment<3>(3 * v) = s;
    }
  });
}

void storeBlock(std::vector<double> &blocks, std::size_t pairIndex, const V12d &local)
{
  std::copy(local.data(), local.data() + 12, blocks.data() + 12 * pairIndex);
}

// Write 144 triplets from a 12x12 local hessian into the pair's preallocated
// slice. Inactive pairs keep their zero-valued (0, 0) placeholders.
void scatterHessian(std::vector<TripletD> &triplets, std::size_t pairIdx, const M12d &localH, const int idx[4])
{
  TripletD *base = triplets.data() + 144 * pairIdx;
  int k = 0;
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
      for (int di = 0; di < 3; ++di)
        for (int dj = 0; dj < 3; ++dj, ++k)
          base[k] = TripletD(3 * idx[i] + di, 3 * idx[j] + dj, localH(3 * i + di, 3 * j + dj));
}

void validateGatherMap(const SurfaceIPCGradientGatherMap *gather, std::size_t numPairs, int numVerts)
{
  if (!gather)
    return;
  if (gather->numVertices() != numVerts || gather->numBlocks() != 4 * numPairs)
    throw std::invalid_argument("SurfaceIPCGradientGatherMap was not built for the current pairs.");
}

}  // namespace

void SurfaceIPCGradientGatherMap::build(const std::vector<PTPair> &ptPairs, const std::vector<EEPair> &eePairs, int numVerts)
{
  start.assign(static_cast<std::size_t>(numVerts) + 1, 0);
  auto count = [&](const int idx[4]) {
    for (int s = 0; s < 4; ++s)
      start[idx[s] + 1]++;
  };
  int idx[4];
  for (const PTPair &pair : ptPairs) {
    ptIndices(pair, idx);
    count(idx);
  }
  for (const EEPair &pair : eePairs) {
    eeIndices(pair, idx);
    count(idx);
  }
  for (int v = 0; v < numVerts; ++v)
    start[v + 1] += start[v];

  blocks.assign(static_cast<std::size_t>(start[numVerts]), -1);
  std::vector<std::ptrdiff_t> cursor(start.begin(), start.end() - 1);
  std::ptrdiff_t pairIndex = 0;
  auto place = [&](const int idx[4]) {
    for (int s = 0; s < 4; ++s)
      blocks[cursor[idx[s]]++] = 4 * pairIndex + s;
    ++pairIndex;
  };
  for (const PTPair &pair : ptPairs) {
    ptIndices(pair, idx);
    place(idx);
  }
  for (const EEPair &pair : eePairs) {
    eeIndices(pair, idx);
    place(idx);
  }
}

// =========================================================================
//  Energy
// =========================================================================
double SurfaceIPCBarrierAssembler::computeEnergy(
  EigenSupport::ConstRefVecXd pos,
  const std::vector<PTPair> &ptPairs,
  const std::vector<EEPair> &eePairs,
  int numVerts,
  double dhat,
  double kappa,
  double eps_ee,
  bool deterministic) const
{
  (void)numVerts;
  const double dhat2 = dhat * dhat;
  const int nPT = static_cast<int>(ptPairs.size());
  const int nEE = static_cast<int>(eePairs.size());

  if (deterministic) {
    // per-pair energies, summed in pair order
    std::vector<double> pairEnergy(static_cast<std::size_t>(nPT + nEE), 0.0);
    tbb::parallel_for(tbb::blocked_range<int>(0, nPT), [&](const tbb::blocked_range<int> &range) {
      for (int i = range.begin(); i < range.end(); ++i)
        evaluatePT(pos, ptPairs[i], dhat2, kappa, false, &pairEnergy[i], nullptr, nullptr);
    });
    tbb::parallel_for(tbb::blocked_range<int>(0, nEE), [&](const tbb::blocked_range<int> &range) {
      for (int i = range.begin(); i < range.end(); ++i)
        evaluateEE(pos, eePairs[i], dhat2, kappa, eps_ee, false, &pairEnergy[nPT + i], nullptr, nullptr);
    });
    double energy = 0.0;
    for (double e : pairEnergy)
      energy += e;
    return energy;
  }

  const double ptEnergy = tbb::parallel_reduce(
    tbb::blocked_range<int>(0, nPT), 0.0,
    [&](const tbb::blocked_range<int> &range, double localE) {
      for (int i = range.begin(); i < range.end(); ++i) {
        double e = 0.0;
        if (evaluatePT(pos, ptPairs[i], dhat2, kappa, false, &e, nullptr, nullptr))
          localE += e;
      }
      return localE;
    },
    std::plus<double>());

  const double eeEnergy = tbb::parallel_reduce(
    tbb::blocked_range<int>(0, nEE), 0.0,
    [&](const tbb::blocked_range<int> &range, double localE) {
      for (int i = range.begin(); i < range.end(); ++i) {
        double e = 0.0;
        if (evaluateEE(pos, eePairs[i], dhat2, kappa, eps_ee, false, &e, nullptr, nullptr))
          localE += e;
      }
      return localE;
    },
    std::plus<double>());

  return ptEnergy + eeEnergy;
}

// =========================================================================
//  Gradient
// =========================================================================
void SurfaceIPCBarrierAssembler::computeGradient(
  EigenSupport::ConstRefVecXd pos,
  const std::vector<PTPair> &ptPairs,
  const std::vector<EEPair> &eePairs,
  int numVerts,
  double dhat,
  double kappa,
  double eps_ee,
  EigenSupport::RefVecXd grad,
  const SurfaceIPCGradientGatherMap *gather) const
{
  const int n = 3 * numVerts;
  if (grad.size() != n)
    throw std::runtime_error("Gradient vector has wrong size");
  const int nPT = static_cast<int>(ptPairs.size());
  const int nEE = static_cast<int>(eePairs.size());
  validateGatherMap(gather, static_cast<std::size_t>(nPT + nEE), numVerts);

  grad.setZero();
  const double dhat2 = dhat * dhat;

  if (gather) {
    std::vector<double> blocks(12 * static_cast<std::size_t>(nPT + nEE), 0.0);
    tbb::parallel_for(tbb::blocked_range<int>(0, nPT), [&](const tbb::blocked_range<int> &range) {
      for (int i = range.begin(); i < range.end(); ++i) {
        V12d gE;
        if (evaluatePT(pos, ptPairs[i], dhat2, kappa, false, nullptr, &gE, nullptr))
          storeBlock(blocks, i, gE);
      }
    });
    tbb::parallel_for(tbb::blocked_range<int>(0, nEE), [&](const tbb::blocked_range<int> &range) {
      for (int i = range.begin(); i < range.end(); ++i) {
        V12d gE;
        if (evaluateEE(pos, eePairs[i], dhat2, kappa, eps_ee, false, nullptr, &gE, nullptr))
          storeBlock(blocks, static_cast<std::size_t>(nPT) + i, gE);
      }
    });
    gatherBlocks(*gather, blocks, grad);
    return;
  }

  double *gdata = grad.data();
  tbb::parallel_for(tbb::blocked_range<int>(0, nPT), [&](const tbb::blocked_range<int> &range) {
    for (int i = range.begin(); i < range.end(); ++i) {
      V12d gE;
      if (!evaluatePT(pos, ptPairs[i], dhat2, kappa, false, nullptr, &gE, nullptr))
        continue;
      int idx[4];
      ptIndices(ptPairs[i], idx);
      scatterAtomic(gdata, gE, idx);
    }
  });
  tbb::parallel_for(tbb::blocked_range<int>(0, nEE), [&](const tbb::blocked_range<int> &range) {
    for (int i = range.begin(); i < range.end(); ++i) {
      V12d gE;
      if (!evaluateEE(pos, eePairs[i], dhat2, kappa, eps_ee, false, nullptr, &gE, nullptr))
        continue;
      int idx[4];
      eeIndices(eePairs[i], idx);
      scatterAtomic(gdata, gE, idx);
    }
  });
}

// =========================================================================
//  Sparse Hessian
// =========================================================================
void SurfaceIPCBarrierAssembler::computeHessian(
  EigenSupport::ConstRefVecXd pos,
  const std::vector<PTPair> &ptPairs,
  const std::vector<EEPair> &eePairs,
  int numVerts,
  double dhat,
  double kappa,
  double eps_ee,
  bool projectHessianToPSD,
  SpMatD &hess) const
{
  const int n = 3 * numVerts;
  const int nPT = static_cast<int>(ptPairs.size());
  const int nEE = static_cast<int>(eePairs.size());
  const double dhat2 = dhat * dhat;

  // Preallocate exactly 144 triplet slots per pair (12x12 block). Each pair
  // writes to its own contiguous slice, so the triplet order is fixed by the
  // pair order and setFromTriplets sums duplicates sequentially.
  std::vector<TripletD> triplets(144 * static_cast<std::size_t>(nPT + nEE), TripletD(0, 0, 0.0));

  tbb::parallel_for(tbb::blocked_range<int>(0, nPT), [&](const tbb::blocked_range<int> &range) {
    for (int i = range.begin(); i < range.end(); ++i) {
      M12d localH;
      if (!evaluatePT(pos, ptPairs[i], dhat2, kappa, projectHessianToPSD, nullptr, nullptr, &localH))
        continue;  // slot already zeroed
      int idx[4];
      ptIndices(ptPairs[i], idx);
      scatterHessian(triplets, i, localH, idx);
    }
  });
  tbb::parallel_for(tbb::blocked_range<int>(0, nEE), [&](const tbb::blocked_range<int> &range) {
    for (int i = range.begin(); i < range.end(); ++i) {
      M12d localH;
      if (!evaluateEE(pos, eePairs[i], dhat2, kappa, eps_ee, projectHessianToPSD, nullptr, nullptr, &localH))
        continue;  // slot already zeroed
      int idx[4];
      eeIndices(eePairs[i], idx);
      scatterHessian(triplets, static_cast<std::size_t>(nPT) + i, localH, idx);
    }
  });

  hess.resize(n, n);
  hess.setFromTriplets(triplets.begin(), triplets.end());
  if (auto logger = Logging::lgr(); logger)
    SPDLOG_LOGGER_INFO(logger, "# nonzeros in Hessian: {}", hess.nonZeros());
}

// =========================================================================
//  Combined computation
// =========================================================================
void SurfaceIPCBarrierAssembler::computeAll(
  EigenSupport::ConstRefVecXd x,
  const std::vector<PTPair> &ptPairs,
  const std::vector<EEPair> &eePairs,
  int numVerts,
  double dhat,
  double kappa,
  double eps_ee,
  bool projectHessianToPSD,
  double &energy,
  VXd &grad,
  SpMatD &hess,
  const SurfaceIPCGradientGatherMap *gather) const
{
  const int n = 3 * numVerts;
  const int nPT = static_cast<int>(ptPairs.size());
  const int nEE = static_cast<int>(eePairs.size());
  const std::size_t numPairs = static_cast<std::size_t>(nPT + nEE);
  validateGatherMap(gather, numPairs, numVerts);

  energy = 0.0;
  grad.setZero(n);
  std::vector<TripletD> triplets(144 * numPairs, TripletD(0, 0, 0.0));
  const double dhat2 = dhat * dhat;

  if (gather) {
    std::vector<double> pairEnergy(numPairs, 0.0);
    std::vector<double> blocks(12 * numPairs, 0.0);
    tbb::parallel_for(tbb::blocked_range<int>(0, nPT), [&](const tbb::blocked_range<int> &range) {
      for (int i = range.begin(); i < range.end(); ++i) {
        V12d gE;
        M12d localH;
        if (!evaluatePT(x, ptPairs[i], dhat2, kappa, projectHessianToPSD, &pairEnergy[i], &gE, &localH))
          continue;
        storeBlock(blocks, i, gE);
        int idx[4];
        ptIndices(ptPairs[i], idx);
        scatterHessian(triplets, i, localH, idx);
      }
    });
    tbb::parallel_for(tbb::blocked_range<int>(0, nEE), [&](const tbb::blocked_range<int> &range) {
      for (int i = range.begin(); i < range.end(); ++i) {
        V12d gE;
        M12d localH;
        const std::size_t pairIndex = static_cast<std::size_t>(nPT) + i;
        if (!evaluateEE(x, eePairs[i], dhat2, kappa, eps_ee, projectHessianToPSD, &pairEnergy[pairIndex], &gE, &localH))
          continue;
        storeBlock(blocks, pairIndex, gE);
        int idx[4];
        eeIndices(eePairs[i], idx);
        scatterHessian(triplets, pairIndex, localH, idx);
      }
    });
    for (double e : pairEnergy)
      energy += e;
    gatherBlocks(*gather, blocks, grad);
  }
  else {
    double *gdata = grad.data();
    const double ptEnergy = tbb::parallel_reduce(
      tbb::blocked_range<int>(0, nPT), 0.0,
      [&](const tbb::blocked_range<int> &range, double localE) {
        for (int i = range.begin(); i < range.end(); ++i) {
          double e = 0.0;
          V12d gE;
          M12d localH;
          if (!evaluatePT(x, ptPairs[i], dhat2, kappa, projectHessianToPSD, &e, &gE, &localH))
            continue;
          localE += e;
          int idx[4];
          ptIndices(ptPairs[i], idx);
          scatterAtomic(gdata, gE, idx);
          scatterHessian(triplets, i, localH, idx);
        }
        return localE;
      },
      std::plus<double>());
    const double eeEnergy = tbb::parallel_reduce(
      tbb::blocked_range<int>(0, nEE), 0.0,
      [&](const tbb::blocked_range<int> &range, double localE) {
        for (int i = range.begin(); i < range.end(); ++i) {
          double e = 0.0;
          V12d gE;
          M12d localH;
          if (!evaluateEE(x, eePairs[i], dhat2, kappa, eps_ee, projectHessianToPSD, &e, &gE, &localH))
            continue;
          localE += e;
          int idx[4];
          eeIndices(eePairs[i], idx);
          scatterAtomic(gdata, gE, idx);
          scatterHessian(triplets, static_cast<std::size_t>(nPT) + i, localH, idx);
        }
        return localE;
      },
      std::plus<double>());
    energy = ptEnergy + eeEnergy;
  }

  hess.resize(n, n);
  hess.setFromTriplets(triplets.begin(), triplets.end());
}

}  // namespace CIPC
}  // namespace Contact
}  // namespace pgo
