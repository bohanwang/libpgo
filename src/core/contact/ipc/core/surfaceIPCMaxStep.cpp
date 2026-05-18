#include "ipc/core/surfaceIPCMaxStep.h"

#include "ipc/broadPhase/spatialHashGrid.h"
#include "ipc/geometry/ipcCCD.h"
#include "scopedProfileSection.h"
#include "ipc/profiling/surfaceIPCProfiling.h"

#include <tbb/blocked_range.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>

#include <algorithm>
#include <functional>
#include <vector>

namespace pgo {
namespace Contact {
namespace CIPC {
using namespace pgo::EigenSupport;

// =========================================================================
//  Self max-step (CCD)
// =========================================================================

double computeSelfMaxStep(
  const SurfaceIPCTopology &topology,
  EigenSupport::ConstRefVecXd x,
  EigenSupport::ConstRefVecXd dx,
  double dhat,
  double slackness,
  double thickness)
{
  const VXd pos = VXd(x);

  auto getV = [&](int i) -> V3d {
    return pos.segment<3>(3 * i);
  };
  auto getdV = [&](int i) -> V3d {
    return dx.segment<3>(3 * i);
  };

  int nTri = (int)topology.triangles.size();
  int nEdge = (int)topology.edges.size();

  std::vector<SpatialHashGrid::AABB> vertBox(topology.numVerts);
  std::vector<SpatialHashGrid::AABB> triBox(nTri);
  std::vector<SpatialHashGrid::AABB> edgeBox(nEdge);
  double cellSize = 0.0;

  {
    Profiling::ScopedProfileSection scopedProfile(SurfaceIPCProfileSections::kPairBuildSwept);

    // Inflate swept AABBs by `thickness` on every side so the broad-phase
    // prune stays sound for min-separation CCD (contact at distance == thickness).
    tbb::parallel_for(tbb::blocked_range<int>(0, topology.numVerts),
      [&](const tbb::blocked_range<int> &r) {
        for (int vi = r.begin(); vi < r.end(); ++vi) {
          V3d p0 = getV(vi), p1 = p0 + getdV(vi);
          vertBox[vi].init(p0, thickness);
          vertBox[vi].expand(p1, thickness);
        }
      });

    tbb::parallel_for(tbb::blocked_range<int>(0, nTri),
      [&](const tbb::blocked_range<int> &r) {
        for (int fi = r.begin(); fi < r.end(); ++fi) {
          auto &tri = topology.triangles[fi];
          V3d v0 = getV(tri[0]), v1 = getV(tri[1]), v2 = getV(tri[2]);
          V3d d0 = getdV(tri[0]), d1 = getdV(tri[1]), d2 = getdV(tri[2]);
          triBox[fi].init(v0, thickness);
          triBox[fi].expand(v1, thickness);
          triBox[fi].expand(v2, thickness);
          triBox[fi].expand(v0 + d0, thickness);
          triBox[fi].expand(v1 + d1, thickness);
          triBox[fi].expand(v2 + d2, thickness);
        }
      });

    tbb::parallel_for(tbb::blocked_range<int>(0, nEdge),
      [&](const tbb::blocked_range<int> &r) {
        for (int ei = r.begin(); ei < r.end(); ++ei) {
          V3d a0 = getV(topology.edges[ei][0]), a1 = getV(topology.edges[ei][1]);
          V3d da0 = getdV(topology.edges[ei][0]), da1 = getdV(topology.edges[ei][1]);
          edgeBox[ei].init(a0, thickness);
          edgeBox[ei].expand(a1, thickness);
          edgeBox[ei].expand(a0 + da0, thickness);
          edgeBox[ei].expand(a1 + da1, thickness);
        }
      });

    double avgBoxDiag = tbb::parallel_reduce(
      tbb::blocked_range<int>(0, nTri), 0.0,
      [&](const tbb::blocked_range<int> &r, double sum) {
        for (int fi = r.begin(); fi < r.end(); ++fi)
          sum += (triBox[fi].hi - triBox[fi].lo).norm();
        return sum;
      },
      std::plus<double>());
    cellSize = nTri > 0 ? std::max(avgBoxDiag / nTri, 1e-6) : std::max(1e-6, dhat);
  }

  double alpha = 1.0;

  // --- PT CCD ---
  {
    Profiling::ScopedProfileSection scopedProfile(SurfaceIPCProfileSections::kMaxStepPT);

    SpatialHashGrid triHash(nTri);
    triHash.setCellSize(cellSize);
    for (int fi = 0; fi < nTri; ++fi)
      triHash.insert(triBox[fi], fi);

    tbb::enumerable_thread_specific<std::vector<int>> tls_visited(
      [nTri]() { return std::vector<int>(nTri, 0); });
    tbb::enumerable_thread_specific<std::vector<int>> tls_candidates;

    alpha = tbb::parallel_reduce(
      tbb::blocked_range<int>(0, topology.numVerts), 1.0,
      [&](const tbb::blocked_range<int> &range, double localAlpha) {
        auto &visited = tls_visited.local();
        auto &candidates = tls_candidates.local();

        for (int vi = range.begin(); vi < range.end(); ++vi) {
          candidates.clear();
          triHash.query(vertBox[vi], -1, visited, vi + 1, candidates);

          for (int fi : candidates) {
            auto &tri = topology.triangles[fi];
            if (vi == tri[0] || vi == tri[1] || vi == tri[2])
              continue;
            if (!vertBox[vi].overlaps(triBox[fi]))
              continue;

            V3d p = getV(vi), dp = getdV(vi);
            V3d t0 = getV(tri[0]), dt0 = getdV(tri[0]);
            V3d t1 = getV(tri[1]), dt1 = getdV(tri[1]);
            V3d t2 = getV(tri[2]), dt2 = getdV(tri[2]);

            double toi = ccd::pointTriangleCCD(p, t0, t1, t2,
              dp, dt0, dt1, dt2,
              thickness, localAlpha);
            if (toi < localAlpha)
              localAlpha = toi * slackness;
          }
        }
        return localAlpha;
      },
      [](double a, double b) { return std::min(a, b); });
  }

  // --- EE CCD ---
  {
    Profiling::ScopedProfileSection scopedProfile(SurfaceIPCProfileSections::kMaxStepEE);

    SpatialHashGrid edgeHash(nEdge);
    edgeHash.setCellSize(cellSize);
    for (int ei = 0; ei < nEdge; ++ei)
      edgeHash.insert(edgeBox[ei], ei);

    tbb::enumerable_thread_specific<std::vector<int>> tls_visited(
      [nEdge]() { return std::vector<int>(nEdge, 0); });
    tbb::enumerable_thread_specific<std::vector<int>> tls_candidates;

    alpha = tbb::parallel_reduce(
      tbb::blocked_range<int>(0, nEdge), alpha,
      [&](const tbb::blocked_range<int> &range, double localAlpha) {
        auto &visited = tls_visited.local();
        auto &candidates = tls_candidates.local();

        for (int ei = range.begin(); ei < range.end(); ++ei) {
          candidates.clear();
          edgeHash.query(edgeBox[ei], ei, visited, ei + 1, candidates);

          int a0 = topology.edges[ei][0], a1 = topology.edges[ei][1];
          for (int ej : candidates) {
            if (ej <= ei)
              continue;

            int b0 = topology.edges[ej][0], b1 = topology.edges[ej][1];
            if (a0 == b0 || a0 == b1 || a1 == b0 || a1 == b1)
              continue;
            if (!edgeBox[ei].overlaps(edgeBox[ej]))
              continue;

            V3d va0 = getV(a0), va1 = getV(a1);
            V3d vb0 = getV(b0), vb1 = getV(b1);
            V3d da0 = getdV(a0), da1 = getdV(a1);
            V3d db0 = getdV(b0), db1 = getdV(b1);

            double toi = ccd::edgeEdgeCCD(va0, va1, vb0, vb1,
              da0, da1, db0, db1,
              thickness, localAlpha);
            if (toi < localAlpha)
              localAlpha = toi * slackness;
          }
        }
        return localAlpha;
      },
      [](double a, double b) { return std::min(a, b); });
  }

  return alpha;
}

// =========================================================================
//  External max-step (CCD)
// =========================================================================

double computeExternalMaxStep(
  const SurfaceIPCTopology &topology,
  EigenSupport::ConstRefVecXd x,
  EigenSupport::ConstRefVecXd dx,
  const std::vector<ObstacleSurface> &obstacles,
  double dhatExternal,
  double slackness,
  double thickness)
{
  if (obstacles.empty())
    return 1.0;

  double alpha = 1.0;
  const VXd pos = VXd(x);

  auto getV = [&](int i) -> V3d {
    return pos.segment<3>(3 * i);
  };
  auto getdV = [&](int i) -> V3d {
    return dx.segment<3>(3 * i);
  };

  int nDynTri = (int)topology.triangles.size();
  int nDynEdge = (int)topology.edges.size();

  // Build dynamic swept AABBs, inflated by `thickness` for min-separation CCD.
  std::vector<SpatialHashGrid::AABB> dynVertBox(topology.numVerts);
  std::vector<SpatialHashGrid::AABB> dynTriBox(nDynTri);
  std::vector<SpatialHashGrid::AABB> dynEdgeBox(nDynEdge);

  tbb::parallel_for(tbb::blocked_range<int>(0, topology.numVerts),
    [&](const tbb::blocked_range<int> &r) {
      for (int vi = r.begin(); vi < r.end(); ++vi) {
        V3d p0 = getV(vi), p1 = p0 + getdV(vi);
        dynVertBox[vi].init(p0, thickness);
        dynVertBox[vi].expand(p1, thickness);
      }
    });

  tbb::parallel_for(tbb::blocked_range<int>(0, nDynTri),
    [&](const tbb::blocked_range<int> &r) {
      for (int fi = r.begin(); fi < r.end(); ++fi) {
        auto &tri = topology.triangles[fi];
        V3d v0 = getV(tri[0]), v1 = getV(tri[1]), v2 = getV(tri[2]);
        V3d d0 = getdV(tri[0]), d1 = getdV(tri[1]), d2 = getdV(tri[2]);
        dynTriBox[fi].init(v0, thickness);
        dynTriBox[fi].expand(v1, thickness);
        dynTriBox[fi].expand(v2, thickness);
        dynTriBox[fi].expand(v0 + d0, thickness);
        dynTriBox[fi].expand(v1 + d1, thickness);
        dynTriBox[fi].expand(v2 + d2, thickness);
      }
    });

  tbb::parallel_for(tbb::blocked_range<int>(0, nDynEdge),
    [&](const tbb::blocked_range<int> &r) {
      for (int ei = r.begin(); ei < r.end(); ++ei) {
        V3d a0 = getV(topology.edges[ei][0]), a1 = getV(topology.edges[ei][1]);
        V3d da0 = getdV(topology.edges[ei][0]), da1 = getdV(topology.edges[ei][1]);
        dynEdgeBox[ei].init(a0, thickness);
        dynEdgeBox[ei].expand(a1, thickness);
        dynEdgeBox[ei].expand(a0 + da0, thickness);
        dynEdgeBox[ei].expand(a1 + da1, thickness);
      }
    });

  for (const auto &obs : obstacles) {
    const VXd &obsCur = obs.currentPositions();
    int nObsVert = (int)obsCur.size() / 3;
    int nObsTri = (int)obs.triangles().rows();
    int nObsEdge = (int)obs.uniqueEdges().rows();

    // The obstacle is treated as fixed at its sampled pose during the
    // Newton line-search; intra-frame obstacle motion is not swept here.
    auto obsV = [&](int i) -> V3d { return obsCur.segment<3>(3 * i); };
    auto obsDisp = [](int) -> V3d { return V3d::Zero(); };

    // Obstacle AABBs, hashes, and cell size are cached on ObstacleSurface
    // (rebuilt by update()); they are stored un-inflated. The dyn-side swept
    // AABBs above already include the `thickness` margin, which by Minkowski
    // equivalence is enough to keep this broad-phase sound under min-separation
    // CCD (and is exact when thickness == 0).
    const ObstaclePoseCache &poseCache = obs.cache();
    const std::vector<SpatialHashGrid::AABB> &obsVertBox = poseCache.vertBoxes;
    const std::vector<SpatialHashGrid::AABB> &obsTriBox = poseCache.triBoxes;
    const std::vector<SpatialHashGrid::AABB> &obsEdgeBox = poseCache.edgeBoxes;
    const double cellSize = poseCache.cellSize > 0.0 ? poseCache.cellSize : std::max(1e-6, dhatExternal);

    // --- External PT CCD ---
    {
      const SpatialHashGrid &obsTriHash = poseCache.triHash;

      tbb::enumerable_thread_specific<std::vector<int>> tls_visited(
        [nObsTri]() { return std::vector<int>(nObsTri, 0); });
      tbb::enumerable_thread_specific<std::vector<int>> tls_candidates;

      alpha = tbb::parallel_reduce(
        tbb::blocked_range<int>(0, topology.numVerts), alpha,
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
              V3d t0 = obsV(obs.triangles()(fi, 0));
              V3d t1 = obsV(obs.triangles()(fi, 1));
              V3d t2 = obsV(obs.triangles()(fi, 2));
              V3d dt0 = obsDisp(obs.triangles()(fi, 0));
              V3d dt1 = obsDisp(obs.triangles()(fi, 1));
              V3d dt2 = obsDisp(obs.triangles()(fi, 2));

              double toi = ccd::pointTriangleCCD(p, t0, t1, t2,
                dp, dt0, dt1, dt2, thickness, localAlpha);
              if (toi < localAlpha)
                localAlpha = toi * slackness;
            }
          }
          return localAlpha;
        },
        [](double a, double b) { return std::min(a, b); });
    }

    // --- External TP CCD ---
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

              auto &tri = topology.triangles[fi];
              V3d p = obsV(ovi), dp = obsDisp(ovi);
              V3d t0 = getV(tri[0]);
              V3d t1 = getV(tri[1]);
              V3d t2 = getV(tri[2]);
              V3d dt0 = getdV(tri[0]);
              V3d dt1 = getdV(tri[1]);
              V3d dt2 = getdV(tri[2]);

              double toi = ccd::pointTriangleCCD(p, t0, t1, t2,
                dp, dt0, dt1, dt2, thickness, localAlpha);
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
      const SpatialHashGrid &obsEdgeHash = poseCache.edgeHash;

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

            int a0 = topology.edges[ei][0], a1 = topology.edges[ei][1];
            for (int ej : candidates) {
              if (!dynEdgeBox[ei].overlaps(obsEdgeBox[ej]))
                continue;

              int b0 = obs.uniqueEdges()(ej, 0);
              int b1 = obs.uniqueEdges()(ej, 1);
              V3d va0 = getV(a0), va1 = getV(a1);
              V3d vb0 = obsV(b0), vb1 = obsV(b1);
              V3d da0 = getdV(a0), da1 = getdV(a1);
              V3d db0 = obsDisp(b0), db1 = obsDisp(b1);

              double toi = ccd::edgeEdgeCCD(va0, va1, vb0, vb1,
                da0, da1, db0, db1, thickness, localAlpha);
              if (toi < localAlpha)
                localAlpha = toi * slackness;
            }
          }
          return localAlpha;
        },
        [](double a, double b) { return std::min(a, b); });
    }
  }

  return alpha;
}

}  // namespace CIPC
}  // namespace Contact
}  // namespace pgo
