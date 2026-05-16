/*
copyright to Bohan Wang
*/

// Full implementation of Codimensional IPC collision handling.
// =============================================================================

#include "ipc/core/surfaceIPCCore.h"
#include "scopedProfileSection.h"
#include "ipc/broadPhase/surfaceIPCSelfBroadPhase.h"
#include "ipc/broadPhase/spatialHashGrid.h"
#include "ipc/core/surfaceIPCBarrierAssembler.h"
#include "ipc/core/surfaceIPCMaxStep.h"
#include "ipc/geometry/ipcCCD.h"
#include "ipc/geometry/ipcDistancePrimitives.h"
#include "ipc/profiling/surfaceIPCProfiling.h"

#include <tbb/blocked_range.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>

#include "pgoLogging.h"

#include <algorithm>
#include <cstdint>
#include <stdexcept>

namespace pgo {
namespace Contact {
namespace CIPC {
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
static V3d obsVtx(const VXd &pos, int i)
{
  return pos.segment<3>(3 * i);
}

void SurfaceIPCCore::findCollisionPairs(const VXd &positions) const
{
  SurfaceIPCSelfBroadPhase().buildPairs(topology_, positions, dhat, preparedState_.ptPairs, preparedState_.eePairs);

  preparedState_.externalPTPairs.clear();
  preparedState_.externalTPPairs.clear();
  preparedState_.externalEEPairs.clear();

  if (obstacles_.empty())
    return;

  const double inflate = dhat_external;
  const double dhat2 = dhat_external * dhat_external;

  auto getV = [&](int i) -> V3d {
    return positions.segment<3>(3 * i);
  };

  int nDynTri = (int)topology_.triangles.size();
  int nDynEdge = (int)topology_.edges.size();

  // Build dynamic-side AABBs
  std::vector<SpatialHashGrid::AABB> dynVertBox(topology_.numVerts);
  tbb::parallel_for(tbb::blocked_range<int>(0, topology_.numVerts),
    [&](const tbb::blocked_range<int> &r) {
      for (int vi = r.begin(); vi < r.end(); ++vi)
        dynVertBox[vi].init(getV(vi), inflate);
    });

  std::vector<SpatialHashGrid::AABB> dynTriBox(nDynTri);
  tbb::parallel_for(tbb::blocked_range<int>(0, nDynTri),
    [&](const tbb::blocked_range<int> &r) {
      for (int fi = r.begin(); fi < r.end(); ++fi) {
        auto &tri = topology_.triangles[fi];
        dynTriBox[fi].init(getV(tri[0]), inflate);
        dynTriBox[fi].expand(getV(tri[1]), inflate);
        dynTriBox[fi].expand(getV(tri[2]), inflate);
      }
    });

  std::vector<SpatialHashGrid::AABB> dynEdgeBox(nDynEdge);
  tbb::parallel_for(tbb::blocked_range<int>(0, nDynEdge),
    [&](const tbb::blocked_range<int> &r) {
      for (int ei = r.begin(); ei < r.end(); ++ei) {
        dynEdgeBox[ei].init(getV(topology_.edges[ei][0]), inflate);
        dynEdgeBox[ei].expand(getV(topology_.edges[ei][1]), inflate);
      }
    });

  for (const auto &obs : obstacles_) {
    const VXd &obsPos = obs->currentPositions();
    int nObsVert = (int)obsPos.size() / 3;
    int nObsTri = (int)obs->triangles().rows();
    int nObsEdge = (int)obs->uniqueEdges().rows();

    int32_t obsId = obs->objectId();

    // Precompute obstacle triangle areas and edge lengths from current positions
    std::vector<double> obsTriArea(nObsTri, 0.0);
    for (int fi = 0; fi < nObsTri; ++fi) {
      V3d v0 = obsVtx(obsPos, obs->triangles()(fi, 0));
      V3d v1 = obsVtx(obsPos, obs->triangles()(fi, 1));
      V3d v2 = obsVtx(obsPos, obs->triangles()(fi, 2));
      obsTriArea[fi] = 0.5 * (v1 - v0).cross(v2 - v0).norm();
    }
    std::vector<double> obsEdgeLen(nObsEdge, 0.0);
    for (int ei = 0; ei < nObsEdge; ++ei) {
      V3d e0 = obsVtx(obsPos, obs->uniqueEdges()(ei, 0));
      V3d e1 = obsVtx(obsPos, obs->uniqueEdges()(ei, 1));
      obsEdgeLen[ei] = (e1 - e0).norm();
    }

    // Build obstacle AABBs
    std::vector<SpatialHashGrid::AABB> obsVertBox(nObsVert);
    tbb::parallel_for(tbb::blocked_range<int>(0, nObsVert),
      [&](const tbb::blocked_range<int> &r) {
        for (int vi = r.begin(); vi < r.end(); ++vi)
          obsVertBox[vi].init(obsVtx(obsPos, vi), inflate);
      });

    std::vector<SpatialHashGrid::AABB> obsTriBox(nObsTri);
    tbb::parallel_for(tbb::blocked_range<int>(0, nObsTri),
      [&](const tbb::blocked_range<int> &r) {
        for (int fi = r.begin(); fi < r.end(); ++fi) {
          obsTriBox[fi].init(obsVtx(obsPos, obs->triangles()(fi, 0)), inflate);
          obsTriBox[fi].expand(obsVtx(obsPos, obs->triangles()(fi, 1)), inflate);
          obsTriBox[fi].expand(obsVtx(obsPos, obs->triangles()(fi, 2)), inflate);
        }
      });

    std::vector<SpatialHashGrid::AABB> obsEdgeBox(nObsEdge);
    tbb::parallel_for(tbb::blocked_range<int>(0, nObsEdge),
      [&](const tbb::blocked_range<int> &r) {
        for (int ei = r.begin(); ei < r.end(); ++ei) {
          obsEdgeBox[ei].init(obsVtx(obsPos, obs->uniqueEdges()(ei, 0)), inflate);
          obsEdgeBox[ei].expand(obsVtx(obsPos, obs->uniqueEdges()(ei, 1)), inflate);
        }
      });

    // Compute cell size for this obstacle
    double avgBoxDiag = 0.0;
    for (const auto &aabb : obsTriBox)
      avgBoxDiag += (aabb.hi - aabb.lo).norm();
    double cellSize = nObsTri > 0 ? std::max(avgBoxDiag / nObsTri, 1e-6) : std::max(1e-6, dhat_external);

    // ---- External PT: dyn vertex × obs triangle ----
    {
      SpatialHashGrid obsTriHash(nObsTri);
      obsTriHash.setCellSize(cellSize);
      for (int fi = 0; fi < nObsTri; ++fi)
        obsTriHash.insert(obsTriBox[fi], fi);

      tbb::enumerable_thread_specific<std::vector<int>> tls_visited(
        [nObsTri]() { return std::vector<int>(nObsTri, 0); });
      tbb::enumerable_thread_specific<std::vector<int>> tls_candidates;
      tbb::enumerable_thread_specific<std::vector<ExternalPTPair>> tls_pairs;

      tbb::parallel_for(tbb::blocked_range<int>(0, topology_.numVerts),
        [&](const tbb::blocked_range<int> &range) {
          auto &visited = tls_visited.local();
          auto &candidates = tls_candidates.local();
          auto &localPairs = tls_pairs.local();

          for (int vi = range.begin(); vi < range.end(); ++vi) {
            candidates.clear();
            obsTriHash.query(dynVertBox[vi], -1, visited, vi + 1, candidates);

            for (int fi : candidates) {
              if (!dynVertBox[vi].overlaps(obsTriBox[fi]))
                continue;

              V3d vp = getV(vi);
              V3d vt0 = obsVtx(obsPos, obs->triangles()(fi, 0));
              V3d vt1 = obsVtx(obsPos, obs->triangles()(fi, 1));
              V3d vt2 = obsVtx(obsPos, obs->triangles()(fi, 2));
              double d2 = distance::computePTSqDist(vp, vt0, vt1, vt2);
              if (d2 < dhat2 && d2 > 0.0) {
                double w = topology_.vertexArea[vi] * obsTriArea[fi];
                localPairs.push_back({ obsId, vi,
                  {{ obs->triangles()(fi, 0), obs->triangles()(fi, 1), obs->triangles()(fi, 2) }},
                  w });
              }
            }
          }
        });

      for (auto &lp : tls_pairs)
        preparedState_.externalPTPairs.insert(preparedState_.externalPTPairs.end(), lp.begin(), lp.end());
    }

    // ---- External TP: obs vertex × dyn triangle ----
    {
      SpatialHashGrid dynTriHash(nDynTri);
      dynTriHash.setCellSize(cellSize);
      for (int fi = 0; fi < nDynTri; ++fi)
        dynTriHash.insert(dynTriBox[fi], fi);

      tbb::enumerable_thread_specific<std::vector<int>> tls_visited(
        [nDynTri]() { return std::vector<int>(nDynTri, 0); });
      tbb::enumerable_thread_specific<std::vector<int>> tls_candidates;
      tbb::enumerable_thread_specific<std::vector<ExternalTPPair>> tls_pairs;

      tbb::parallel_for(tbb::blocked_range<int>(0, nObsVert),
        [&](const tbb::blocked_range<int> &range) {
          auto &visited = tls_visited.local();
          auto &candidates = tls_candidates.local();
          auto &localPairs = tls_pairs.local();

          for (int ovi = range.begin(); ovi < range.end(); ++ovi) {
            candidates.clear();
            dynTriHash.query(obsVertBox[ovi], -1, visited, ovi + 1, candidates);

            for (int fi : candidates) {
              if (!obsVertBox[ovi].overlaps(dynTriBox[fi]))
                continue;

              auto &tri = topology_.triangles[fi];
              V3d vp = obsVtx(obsPos, ovi);
              V3d vt0 = getV(tri[0]);
              V3d vt1 = getV(tri[1]);
              V3d vt2 = getV(tri[2]);
              double d2 = distance::computePTSqDist(vp, vt0, vt1, vt2);
              if (d2 < dhat2 && d2 > 0.0) {
                double w = 1.0 * topology_.triArea[fi];
                localPairs.push_back({ obsId,
                  {{ tri[0], tri[1], tri[2] }}, ovi, w });
              }
            }
          }
        });

      for (auto &lp : tls_pairs)
        preparedState_.externalTPPairs.insert(preparedState_.externalTPPairs.end(), lp.begin(), lp.end());
    }

    // ---- External EE: dyn edge × obs edge ----
    {
      SpatialHashGrid obsEdgeHash(nObsEdge);
      obsEdgeHash.setCellSize(cellSize);
      for (int ei = 0; ei < nObsEdge; ++ei)
        obsEdgeHash.insert(obsEdgeBox[ei], ei);

      tbb::enumerable_thread_specific<std::vector<int>> tls_visited(
        [nObsEdge]() { return std::vector<int>(nObsEdge, 0); });
      tbb::enumerable_thread_specific<std::vector<int>> tls_candidates;
      tbb::enumerable_thread_specific<std::vector<ExternalEEPair>> tls_pairs;

      tbb::parallel_for(tbb::blocked_range<int>(0, nDynEdge),
        [&](const tbb::blocked_range<int> &range) {
          auto &visited = tls_visited.local();
          auto &candidates = tls_candidates.local();
          auto &localPairs = tls_pairs.local();

          for (int ei = range.begin(); ei < range.end(); ++ei) {
            candidates.clear();
            obsEdgeHash.query(dynEdgeBox[ei], -1, visited, ei + 1, candidates);

            int a0 = topology_.edges[ei][0], a1 = topology_.edges[ei][1];
            for (int ej : candidates) {
              if (!dynEdgeBox[ei].overlaps(obsEdgeBox[ej]))
                continue;

              int b0 = obs->uniqueEdges()(ej, 0);
              int b1 = obs->uniqueEdges()(ej, 1);
              V3d va0 = getV(a0), va1 = getV(a1);
              V3d vb0 = obsVtx(obsPos, b0), vb1 = obsVtx(obsPos, b1);
              double d2 = distance::computeEESqDist(va0, va1, vb0, vb1);
              if (d2 < dhat2 && d2 > 0.0) {
                double w = topology_.edgeLength[ei] * obsEdgeLen[ej];
                localPairs.push_back({ obsId,
                  {{ a0, a1 }}, {{ b0, b1 }}, w });
              }
            }
          }
        });

      for (auto &lp : tls_pairs)
        preparedState_.externalEEPairs.insert(preparedState_.externalEEPairs.end(), lp.begin(), lp.end());
    }
  }
}

void SurfaceIPCCore::prepareForSurfacePositions(EigenSupport::ConstRefVecXd x_surf) const
{
  Profiling::ScopedProfileSection scopedProfile(SurfaceIPCProfileSections::kPrepareActivePairs);
  preparedState_.positions = x_surf;
  findCollisionPairs(preparedState_.positions);
  if (auto logger = Logging::lgr(); logger) {
    const size_t selfTotal = preparedState_.ptPairs.size() + preparedState_.eePairs.size();
    const size_t externalTotal = preparedState_.externalPTPairs.size() + preparedState_.externalTPPairs.size() + preparedState_.externalEEPairs.size();
    SPDLOG_LOGGER_INFO(logger,
      "SurfaceIPCCore active pairs: selfPT={} selfEE={} selfTotal={} externalPT={} externalTP={} externalEE={} externalTotal={}",
      preparedState_.ptPairs.size(), preparedState_.eePairs.size(), selfTotal,
      preparedState_.externalPTPairs.size(), preparedState_.externalTPPairs.size(), preparedState_.externalEEPairs.size(), externalTotal);
  }
  preparedState_.hasState = true;
}

// =========================================================================
//  1)  Maximum step size  (CCD-based line search with spatial hashing)
// =========================================================================
NonlinearOptimization::MaxStepResult SurfaceIPCCore::computeMaxStepLimit(EigenSupport::ConstRefVecXd x, EigenSupport::ConstRefVecXd dx) const
{
  double alpha = SurfaceIPCMaxStep().compute(topology_, x, dx, dhat, slackness);

  // External CCD — obstacle displacement = current - previous (fixed per stage)
  if (!obstacles_.empty()) {
    const VXd pos = VXd(x);

    auto getV = [&](int i) -> V3d {
      return pos.segment<3>(3 * i);
    };
    auto getdV = [&](int i) -> V3d {
      return dx.segment<3>(3 * i);
    };

    int nDynTri = (int)topology_.triangles.size();
    int nDynEdge = (int)topology_.edges.size();

    // Build dynamic swept AABBs
    std::vector<SpatialHashGrid::AABB> dynVertBox(topology_.numVerts);
    std::vector<SpatialHashGrid::AABB> dynTriBox(nDynTri);
    std::vector<SpatialHashGrid::AABB> dynEdgeBox(nDynEdge);

    tbb::parallel_for(tbb::blocked_range<int>(0, topology_.numVerts),
      [&](const tbb::blocked_range<int> &r) {
        for (int vi = r.begin(); vi < r.end(); ++vi) {
          V3d p0 = getV(vi), p1 = p0 + getdV(vi);
          dynVertBox[vi].init(p0, 0.0);
          dynVertBox[vi].expand(p1);
        }
      });

    tbb::parallel_for(tbb::blocked_range<int>(0, nDynTri),
      [&](const tbb::blocked_range<int> &r) {
        for (int fi = r.begin(); fi < r.end(); ++fi) {
          auto &tri = topology_.triangles[fi];
          V3d v0 = getV(tri[0]), v1 = getV(tri[1]), v2 = getV(tri[2]);
          V3d d0 = getdV(tri[0]), d1 = getdV(tri[1]), d2 = getdV(tri[2]);
          dynTriBox[fi].init(v0, 0.0);
          dynTriBox[fi].expand(v1);
          dynTriBox[fi].expand(v2);
          dynTriBox[fi].expand(v0 + d0);
          dynTriBox[fi].expand(v1 + d1);
          dynTriBox[fi].expand(v2 + d2);
        }
      });

    tbb::parallel_for(tbb::blocked_range<int>(0, nDynEdge),
      [&](const tbb::blocked_range<int> &r) {
        for (int ei = r.begin(); ei < r.end(); ++ei) {
          V3d a0 = getV(topology_.edges[ei][0]), a1 = getV(topology_.edges[ei][1]);
          V3d da0 = getdV(topology_.edges[ei][0]), da1 = getdV(topology_.edges[ei][1]);
          dynEdgeBox[ei].init(a0, 0.0);
          dynEdgeBox[ei].expand(a1);
          dynEdgeBox[ei].expand(a0 + da0);
          dynEdgeBox[ei].expand(a1 + da1);
        }
      });

    for (const auto &obs : obstacles_) {
      const VXd &obsCur = obs->currentPositions();
      const VXd &obsPrev = obs->previousPositions();
      int nObsVert = (int)obsCur.size() / 3;
      int nObsTri = (int)obs->triangles().rows();
      int nObsEdge = (int)obs->uniqueEdges().rows();

      auto obsV = [&](int i) -> V3d { return obsPrev.segment<3>(3 * i); };
      auto obsDisp = [&](int i) -> V3d { return obsCur.segment<3>(3 * i) - obsPrev.segment<3>(3 * i); };

      // Build obstacle swept AABBs
      std::vector<SpatialHashGrid::AABB> obsVertBox(nObsVert);
      std::vector<SpatialHashGrid::AABB> obsTriBox(nObsTri);
      std::vector<SpatialHashGrid::AABB> obsEdgeBox(nObsEdge);

      tbb::parallel_for(tbb::blocked_range<int>(0, nObsVert),
        [&](const tbb::blocked_range<int> &r) {
          for (int vi = r.begin(); vi < r.end(); ++vi) {
            V3d p0 = obsV(vi), p1 = p0 + obsDisp(vi);
            obsVertBox[vi].init(p0, 0.0);
            obsVertBox[vi].expand(p1);
          }
        });

      tbb::parallel_for(tbb::blocked_range<int>(0, nObsTri),
        [&](const tbb::blocked_range<int> &r) {
          for (int fi = r.begin(); fi < r.end(); ++fi) {
            V3d v0 = obsV(obs->triangles()(fi, 0));
            V3d v1 = obsV(obs->triangles()(fi, 1));
            V3d v2 = obsV(obs->triangles()(fi, 2));
            V3d d0 = obsDisp(obs->triangles()(fi, 0));
            V3d d1 = obsDisp(obs->triangles()(fi, 1));
            V3d d2 = obsDisp(obs->triangles()(fi, 2));
            obsTriBox[fi].init(v0, 0.0);
            obsTriBox[fi].expand(v1);
            obsTriBox[fi].expand(v2);
            obsTriBox[fi].expand(v0 + d0);
            obsTriBox[fi].expand(v1 + d1);
            obsTriBox[fi].expand(v2 + d2);
          }
        });

      tbb::parallel_for(tbb::blocked_range<int>(0, nObsEdge),
        [&](const tbb::blocked_range<int> &r) {
          for (int ei = r.begin(); ei < r.end(); ++ei) {
            V3d a0 = obsV(obs->uniqueEdges()(ei, 0));
            V3d a1 = obsV(obs->uniqueEdges()(ei, 1));
            V3d da0 = obsDisp(obs->uniqueEdges()(ei, 0));
            V3d da1 = obsDisp(obs->uniqueEdges()(ei, 1));
            obsEdgeBox[ei].init(a0, 0.0);
            obsEdgeBox[ei].expand(a1);
            obsEdgeBox[ei].expand(a0 + da0);
            obsEdgeBox[ei].expand(a1 + da1);
          }
        });

      // Cell size
      double avgBoxDiag = tbb::parallel_reduce(
        tbb::blocked_range<int>(0, nObsTri), 0.0,
        [&](const tbb::blocked_range<int> &r, double sum) {
          for (int fi = r.begin(); fi < r.end(); ++fi)
            sum += (obsTriBox[fi].hi - obsTriBox[fi].lo).norm();
          return sum;
        },
        std::plus<double>());
      double cellSize = nObsTri > 0 ? std::max(avgBoxDiag / nObsTri, 1e-6) : std::max(1e-6, dhat_external);

      // --- External PT CCD: insert obs triangles, query with dyn vertices ---
      {
        SpatialHashGrid obsTriHash(nObsTri);
        obsTriHash.setCellSize(cellSize);
        for (int fi = 0; fi < nObsTri; ++fi)
          obsTriHash.insert(obsTriBox[fi], fi);

        tbb::enumerable_thread_specific<std::vector<int>> tls_visited(
          [nObsTri]() { return std::vector<int>(nObsTri, 0); });
        tbb::enumerable_thread_specific<std::vector<int>> tls_candidates;

        alpha = tbb::parallel_reduce(
          tbb::blocked_range<int>(0, topology_.numVerts), alpha,
          [&](const tbb::blocked_range<int> &range, double localAlpha) {
            auto &visited = tls_visited.local();
            auto &candidates = tls_candidates.local();

            for (int vi = range.begin(); vi < range.end(); ++vi) {
              candidates.clear();
              obsTriHash.query(dynVertBox[vi], -1, visited, vi + 1, candidates);

              for (int fi : candidates) {
                if (!dynVertBox[vi].overlaps(obsTriBox[fi]))
                  continue;

                V3d p = getV(vi), dp = getdV(vi);
                V3d t0 = obsV(obs->triangles()(fi, 0));
                V3d t1 = obsV(obs->triangles()(fi, 1));
                V3d t2 = obsV(obs->triangles()(fi, 2));
                V3d dt0 = obsDisp(obs->triangles()(fi, 0));
                V3d dt1 = obsDisp(obs->triangles()(fi, 1));
                V3d dt2 = obsDisp(obs->triangles()(fi, 2));

                double toi = ccd::pointTriangleCCD(p, t0, t1, t2,
                  dp, dt0, dt1, dt2, 0.0, localAlpha);
                if (toi < localAlpha)
                  localAlpha = toi * slackness;
              }
            }
            return localAlpha;
          },
          [](double a, double b) { return std::min(a, b); });
      }

      // --- External TP CCD: insert dyn triangles, query with obs vertices ---
      {
        SpatialHashGrid dynTriHash(nDynTri);
        dynTriHash.setCellSize(cellSize);
        for (int fi = 0; fi < nDynTri; ++fi)
          dynTriHash.insert(dynTriBox[fi], fi);

        tbb::enumerable_thread_specific<std::vector<int>> tls_visited(
          [nDynTri]() { return std::vector<int>(nDynTri, 0); });
        tbb::enumerable_thread_specific<std::vector<int>> tls_candidates;

        alpha = tbb::parallel_reduce(
          tbb::blocked_range<int>(0, nObsVert), alpha,
          [&](const tbb::blocked_range<int> &range, double localAlpha) {
            auto &visited = tls_visited.local();
            auto &candidates = tls_candidates.local();

            for (int ovi = range.begin(); ovi < range.end(); ++ovi) {
              candidates.clear();
              dynTriHash.query(obsVertBox[ovi], -1, visited, ovi + 1, candidates);

              for (int fi : candidates) {
                if (!obsVertBox[ovi].overlaps(dynTriBox[fi]))
                  continue;

                auto &tri = topology_.triangles[fi];
                V3d p = obsV(ovi), dp = obsDisp(ovi);
                V3d t0 = getV(tri[0]);
                V3d t1 = getV(tri[1]);
                V3d t2 = getV(tri[2]);
                V3d dt0 = getdV(tri[0]);
                V3d dt1 = getdV(tri[1]);
                V3d dt2 = getdV(tri[2]);

                double toi = ccd::pointTriangleCCD(p, t0, t1, t2,
                  dp, dt0, dt1, dt2, 0.0, localAlpha);
                if (toi < localAlpha)
                  localAlpha = toi * slackness;
              }
            }
            return localAlpha;
          },
          [](double a, double b) { return std::min(a, b); });
      }

      // --- External EE CCD ---
      {
        SpatialHashGrid obsEdgeHash(nObsEdge);
        obsEdgeHash.setCellSize(cellSize);
        for (int ei = 0; ei < nObsEdge; ++ei)
          obsEdgeHash.insert(obsEdgeBox[ei], ei);

        tbb::enumerable_thread_specific<std::vector<int>> tls_visited(
          [nObsEdge]() { return std::vector<int>(nObsEdge, 0); });
        tbb::enumerable_thread_specific<std::vector<int>> tls_candidates;

        alpha = tbb::parallel_reduce(
          tbb::blocked_range<int>(0, nDynEdge), alpha,
          [&](const tbb::blocked_range<int> &range, double localAlpha) {
            auto &visited = tls_visited.local();
            auto &candidates = tls_candidates.local();

            for (int ei = range.begin(); ei < range.end(); ++ei) {
              candidates.clear();
              obsEdgeHash.query(dynEdgeBox[ei], -1, visited, ei + 1, candidates);

              int a0 = topology_.edges[ei][0], a1 = topology_.edges[ei][1];
              for (int ej : candidates) {
                if (!dynEdgeBox[ei].overlaps(obsEdgeBox[ej]))
                  continue;

                int b0 = obs->uniqueEdges()(ej, 0);
                int b1 = obs->uniqueEdges()(ej, 1);
                V3d va0 = getV(a0), va1 = getV(a1);
                V3d vb0 = obsV(b0), vb1 = obsV(b1);
                V3d da0 = getdV(a0), da1 = getdV(a1);
                V3d db0 = obsDisp(b0), db1 = obsDisp(b1);

                double toi = ccd::edgeEdgeCCD(va0, va1, vb0, vb1,
                  da0, da1, db0, db1, 0.0, localAlpha);
                if (toi < localAlpha)
                  localAlpha = toi * slackness;
              }
            }
            return localAlpha;
          },
          [](double a, double b) { return std::min(a, b); });
      }
    }
  }

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
  double e = SurfaceIPCBarrierAssembler().computeEnergy(preparedState_.positions, preparedState_.ptPairs, preparedState_.eePairs, topology_.numVerts, dhat, kappa, eps_ee);
  if (!obstacles_.empty()) {
    e += SurfaceIPCBarrierAssembler().computeExternalEnergy(
      preparedState_.positions, obstacles_, preparedState_.externalPTPairs, preparedState_.externalTPPairs, preparedState_.externalEEPairs, dhat_external, kappa, eps_ee);
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
  SurfaceIPCBarrierAssembler().computeGradient(preparedState_.positions, preparedState_.ptPairs, preparedState_.eePairs, topology_.numVerts, dhat, kappa, eps_ee, grad);
  if (!obstacles_.empty()) {
    SurfaceIPCBarrierAssembler().computeExternalGradient(
      preparedState_.positions, obstacles_, preparedState_.externalPTPairs, preparedState_.externalTPPairs, preparedState_.externalEEPairs, topology_.numVerts, dhat_external, kappa, eps_ee, grad);
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
  SurfaceIPCBarrierAssembler().computeHessian(preparedState_.positions, preparedState_.ptPairs, preparedState_.eePairs, topology_.numVerts, dhat, kappa, eps_ee, hess);
  if (!obstacles_.empty()) {
    SurfaceIPCBarrierAssembler().computeExternalHessian(
      preparedState_.positions, obstacles_, preparedState_.externalPTPairs, preparedState_.externalTPPairs, preparedState_.externalEEPairs, topology_.numVerts, dhat_external, kappa, eps_ee, hess);
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
  SurfaceIPCBarrierAssembler().computeAll(preparedState_.positions, preparedState_.ptPairs, preparedState_.eePairs, topology_.numVerts, dhat, kappa, eps_ee, energy, grad, hess);
  if (!obstacles_.empty()) {
    double extEnergy = 0.0;
    SurfaceIPCBarrierAssembler().computeExternalAll(
      preparedState_.positions, obstacles_, preparedState_.externalPTPairs, preparedState_.externalTPPairs, preparedState_.externalEEPairs, topology_.numVerts, dhat_external, kappa, eps_ee, extEnergy, grad, hess);
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
