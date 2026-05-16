#include "surfaceIPCSelfBarrierAssembler.h"

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
  double ee_eps = eps_ee;
  double dhat2 = dhat * dhat;

  // PT pairs
  double ptEnergy = tbb::parallel_reduce(
    tbb::blocked_range<int>(0, (int)pairs.ptPairs.size()), 0.0,
    [&](const tbb::blocked_range<int> &range, double localE) {
      for (int i = range.begin(); i < range.end(); ++i) {
        auto &pair = pairs.ptPairs[i];
        V3d p = vtx(dynPos, pair.p);
        V3d t0 = vtx(dynPos, pair.t0);
        V3d t1 = vtx(dynPos, pair.t1);
        V3d t2 = vtx(dynPos, pair.t2);
        double d2 = distance::computePTSqDist(p, t0, t1, t2);
        if (d2 < dhat2 && d2 > 0.0)
          localE += pair.weight * kappa * barrier::b(d2, dhat2);
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
        V3d ea0 = vtx(dynPos, pair.ea0);
        V3d ea1 = vtx(dynPos, pair.ea1);
        V3d eb0 = vtx(dynPos, pair.eb0);
        V3d eb1 = vtx(dynPos, pair.eb1);
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
    tbb::blocked_range<int>(0, (int)pairs.ptPairs.size()),
    [&](const tbb::blocked_range<int> &range) {
      for (int i = range.begin(); i < range.end(); ++i) {
        auto &pair = pairs.ptPairs[i];
        V3d p = vtx(dynPos, pair.p);
        V3d t0 = vtx(dynPos, pair.t0);
        V3d t1 = vtx(dynPos, pair.t1);
        V3d t2 = vtx(dynPos, pair.t2);

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
    tbb::blocked_range<int>(0, (int)pairs.eePairs.size()),
    [&](const tbb::blocked_range<int> &range) {
      for (int i = range.begin(); i < range.end(); ++i) {
        auto &pair = pairs.eePairs[i];
        V3d ea0 = vtx(dynPos, pair.ea0);
        V3d ea1 = vtx(dynPos, pair.ea1);
        V3d eb0 = vtx(dynPos, pair.eb0);
        V3d eb1 = vtx(dynPos, pair.eb1);

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
        auto &pair = pairs.ptPairs[i];
        V3d p = vtx(dynPos, pair.p);
        V3d t0 = vtx(dynPos, pair.t0);
        V3d t1 = vtx(dynPos, pair.t1);
        V3d t2 = vtx(dynPos, pair.t2);

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
        auto &pair = pairs.eePairs[i];
        V3d ea0 = vtx(dynPos, pair.ea0);
        V3d ea1 = vtx(dynPos, pair.ea1);
        V3d eb0 = vtx(dynPos, pair.eb0);
        V3d eb1 = vtx(dynPos, pair.eb1);

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

        int idx[4] = { pair.ea0, pair.ea1, pair.eb0, pair.eb1 };
        scatterH(nPT + i, localH, idx);
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

  energy = 0.0;
  grad.setZero(n);
  std::vector<TripletD> triplets(144 * totalPairs, TripletD(0, 0, 0.0));

  auto scatterG = [&](const V12d &local, const int idx[4]) {
    double *gdata = grad.data();
    for (int i = 0; i < 4; ++i)
      if (idx[i] >= 0)
        for (int d = 0; d < 3; ++d)
          std::atomic_ref<double>(gdata[3 * idx[i] + d])
            .fetch_add(local[3 * i + d], std::memory_order_relaxed);
  };

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
        auto &pair = pairs.ptPairs[i];
        V3d p = vtx(dynPos, pair.p);
        V3d t0 = vtx(dynPos, pair.t0);
        V3d t1 = vtx(dynPos, pair.t1);
        V3d t2 = vtx(dynPos, pair.t2);

        double d2 = distance::computePTSqDist(p, t0, t1, t2);
        if (d2 >= dhat2 || d2 <= 0.0)
          continue;

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
        auto &pair = pairs.eePairs[i];
        V3d ea0 = vtx(dynPos, pair.ea0);
        V3d ea1 = vtx(dynPos, pair.ea1);
        V3d eb0 = vtx(dynPos, pair.eb0);
        V3d eb1 = vtx(dynPos, pair.eb1);

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

}  // namespace CIPC
}  // namespace Contact
}  // namespace pgo
