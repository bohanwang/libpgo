#include "ipc/broadPhase/surfaceIPCBroadPhase.h"

#include "scopedProfileSection.h"
#include "ipc/profiling/surfaceIPCProfiling.h"
#include "ipc/broadPhase/spatialHashGrid.h"
#include "ipc/geometry/ipcDistancePrimitives.h"

#include <tbb/blocked_range.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>

#include <algorithm>
#include <functional>
#include <vector>

namespace pgo
{
namespace Contact
{
namespace CIPC
{

// =========================================================================
//  Self broad phase: buildSelfPairs
// =========================================================================

void buildSelfPairs(
  const SurfaceIPCTopology &topology,
  EigenSupport::ConstRefVecXd positions,
  double dhat,
  SelfPairSet &pairs)
{
  using namespace pgo::EigenSupport;

  Profiling::ScopedProfileSection scopedProfile(SurfaceIPCProfileSections::kPairBuildStatic);

  pairs.clear();

  const double inflate = dhat;

  auto getV = [&](int i) -> V3d {
    return positions.segment<3>(3 * i);
  };

  int nTri = (int)topology.triangles.size();
  int nEdge = (int)topology.edges.size();

  // --- Build AABBs (parallel) ---
  std::vector<SpatialHashGrid::AABB> vertBox(topology.numVerts);
  tbb::parallel_for(tbb::blocked_range<int>(0, topology.numVerts),
    [&](const tbb::blocked_range<int> &r) {
      for (int vi = r.begin(); vi < r.end(); ++vi)
        vertBox[vi].init(getV(vi), inflate);
    });

  std::vector<SpatialHashGrid::AABB> triBox(nTri);
  tbb::parallel_for(tbb::blocked_range<int>(0, nTri),
    [&](const tbb::blocked_range<int> &r) {
      for (int fi = r.begin(); fi < r.end(); ++fi) {
        auto &tri = topology.triangles[fi];
        triBox[fi].init(getV(tri[0]), inflate);
        triBox[fi].expand(getV(tri[1]), inflate);
        triBox[fi].expand(getV(tri[2]), inflate);
      }
    });

  std::vector<SpatialHashGrid::AABB> edgeBox(nEdge);
  tbb::parallel_for(tbb::blocked_range<int>(0, nEdge),
    [&](const tbb::blocked_range<int> &r) {
      for (int ei = r.begin(); ei < r.end(); ++ei) {
        edgeBox[ei].init(getV(topology.edges[ei][0]), inflate);
        edgeBox[ei].expand(getV(topology.edges[ei][1]), inflate);
      }
    });

  double avgBoxDiag = 0;
  for (const auto &aabb : triBox) {
    avgBoxDiag += (aabb.hi - aabb.lo).norm();
  }

  double cellSize = nTri > 0 ? std::max(avgBoxDiag / nTri, 1e-6) : std::max(1e-6, dhat);
  double dhat2 = dhat * dhat;

  // --- PT pairs: insert triangles (serial), query with vertices (parallel) ---
  {
    SpatialHashGrid triHash(nTri);
    triHash.setCellSize(cellSize);
    for (int fi = 0; fi < nTri; ++fi)
      triHash.insert(triBox[fi], fi);

    tbb::enumerable_thread_specific<std::vector<int>> tls_visited(
      [nTri]() { return std::vector<int>(nTri, 0); });
    tbb::enumerable_thread_specific<std::vector<int>> tls_candidates;
    tbb::enumerable_thread_specific<std::vector<PTPair>> tls_pairs;

    tbb::parallel_for(tbb::blocked_range<int>(0, topology.numVerts),
      [&](const tbb::blocked_range<int> &range) {
        auto &visited = tls_visited.local();
        auto &candidates = tls_candidates.local();
        auto &localPairs = tls_pairs.local();

        for (int vi = range.begin(); vi < range.end(); ++vi) {
          candidates.clear();
          triHash.query(vertBox[vi], -1, visited, vi + 1, candidates);

          for (int fi : candidates) {
            auto &tri = topology.triangles[fi];
            if (vi == tri[0] || vi == tri[1] || vi == tri[2])
              continue;
            if (!vertBox[vi].overlaps(triBox[fi]))
              continue;

            V3d vp = getV(vi);
            V3d vt0 = getV(tri[0]), vt1 = getV(tri[1]), vt2 = getV(tri[2]);
            double d2 = distance::computePTSqDist(vp, vt0, vt1, vt2);
            if (d2 < dhat2)
              localPairs.push_back({ vi, tri[0], tri[1], tri[2],
                topology.vertexArea[vi] * topology.triArea[fi] });
          }
        }
      });

    for (auto &lp : tls_pairs)
      pairs.ptPairs.insert(pairs.ptPairs.end(), lp.begin(), lp.end());
  }

  // --- EE pairs: insert edges (serial), query with edges (parallel) ---
  {
    SpatialHashGrid edgeHash(nEdge);
    edgeHash.setCellSize(cellSize);
    for (int ei = 0; ei < nEdge; ++ei)
      edgeHash.insert(edgeBox[ei], ei);

    tbb::enumerable_thread_specific<std::vector<int>> tls_visited(
      [nEdge]() { return std::vector<int>(nEdge, 0); });
    tbb::enumerable_thread_specific<std::vector<int>> tls_candidates;
    tbb::enumerable_thread_specific<std::vector<EEPair>> tls_pairs;

    tbb::parallel_for(tbb::blocked_range<int>(0, nEdge),
      [&](const tbb::blocked_range<int> &range) {
        auto &visited = tls_visited.local();
        auto &candidates = tls_candidates.local();
        auto &localPairs = tls_pairs.local();

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
            double d2 = distance::computeEESqDist(va0, va1, vb0, vb1);
            if (d2 < dhat2)
              localPairs.push_back({ a0, a1, b0, b1,
                topology.edgeLength[ei] * topology.edgeLength[ej] });
          }
        }
      });

    for (auto &lp : tls_pairs)
      pairs.eePairs.insert(pairs.eePairs.end(), lp.begin(), lp.end());
  }
}

// =========================================================================
//  External broad phase: buildExternalPairs
// =========================================================================

static EigenSupport::V3d obsVtx(const EigenSupport::VXd &pos, int i)
{
  return pos.segment<3>(3 * i);
}

void buildExternalPairs(
  const SurfaceIPCTopology &topology,
  EigenSupport::ConstRefVecXd positions,
  const std::vector<std::shared_ptr<ObstacleSurface>> &obstacles,
  double dhatExternal,
  ExternalPairSet &pairs)
{
  pairs.clear();

  if (obstacles.empty())
    return;

  const double inflate = dhatExternal;
  const double dhat2 = dhatExternal * dhatExternal;

  auto getV = [&](int i) -> EigenSupport::V3d {
    return positions.segment<3>(3 * i);
  };

  int nDynTri = (int)topology.triangles.size();
  int nDynEdge = (int)topology.edges.size();

  // Build dynamic-side AABBs
  std::vector<SpatialHashGrid::AABB> dynVertBox(topology.numVerts);
  tbb::parallel_for(tbb::blocked_range<int>(0, topology.numVerts),
    [&](const tbb::blocked_range<int> &r) {
      for (int vi = r.begin(); vi < r.end(); ++vi)
        dynVertBox[vi].init(getV(vi), inflate);
    });

  std::vector<SpatialHashGrid::AABB> dynTriBox(nDynTri);
  tbb::parallel_for(tbb::blocked_range<int>(0, nDynTri),
    [&](const tbb::blocked_range<int> &r) {
      for (int fi = r.begin(); fi < r.end(); ++fi) {
        auto &tri = topology.triangles[fi];
        dynTriBox[fi].init(getV(tri[0]), inflate);
        dynTriBox[fi].expand(getV(tri[1]), inflate);
        dynTriBox[fi].expand(getV(tri[2]), inflate);
      }
    });

  std::vector<SpatialHashGrid::AABB> dynEdgeBox(nDynEdge);
  tbb::parallel_for(tbb::blocked_range<int>(0, nDynEdge),
    [&](const tbb::blocked_range<int> &r) {
      for (int ei = r.begin(); ei < r.end(); ++ei) {
        dynEdgeBox[ei].init(getV(topology.edges[ei][0]), inflate);
        dynEdgeBox[ei].expand(getV(topology.edges[ei][1]), inflate);
      }
    });

  for (const auto &obs : obstacles) {
    const EigenSupport::VXd &obsPos = obs->currentPositions();
    int nObsVert = (int)obsPos.size() / 3;
    int nObsTri = (int)obs->triangles().rows();
    int nObsEdge = (int)obs->uniqueEdges().rows();

    int32_t obsId = obs->objectId();

    // Precompute obstacle triangle areas and edge lengths from current positions
    std::vector<double> obsTriArea(nObsTri, 0.0);
    for (int fi = 0; fi < nObsTri; ++fi) {
      EigenSupport::V3d v0 = obsVtx(obsPos, obs->triangles()(fi, 0));
      EigenSupport::V3d v1 = obsVtx(obsPos, obs->triangles()(fi, 1));
      EigenSupport::V3d v2 = obsVtx(obsPos, obs->triangles()(fi, 2));
      obsTriArea[fi] = 0.5 * (v1 - v0).cross(v2 - v0).norm();
    }
    std::vector<double> obsEdgeLen(nObsEdge, 0.0);
    for (int ei = 0; ei < nObsEdge; ++ei) {
      EigenSupport::V3d e0 = obsVtx(obsPos, obs->uniqueEdges()(ei, 0));
      EigenSupport::V3d e1 = obsVtx(obsPos, obs->uniqueEdges()(ei, 1));
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
    double cellSize = nObsTri > 0 ? std::max(avgBoxDiag / nObsTri, 1e-6) : std::max(1e-6, dhatExternal);

    // ---- External PT: dyn vertex x obs triangle ----
    {
      SpatialHashGrid obsTriHash(nObsTri);
      obsTriHash.setCellSize(cellSize);
      for (int fi = 0; fi < nObsTri; ++fi)
        obsTriHash.insert(obsTriBox[fi], fi);

      tbb::enumerable_thread_specific<std::vector<int>> tls_visited(
        [nObsTri]() { return std::vector<int>(nObsTri, 0); });
      tbb::enumerable_thread_specific<std::vector<int>> tls_candidates;
      tbb::enumerable_thread_specific<std::vector<ExternalPTPair>> tls_pairs;

      tbb::parallel_for(tbb::blocked_range<int>(0, topology.numVerts),
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

              EigenSupport::V3d vp = getV(vi);
              EigenSupport::V3d vt0 = obsVtx(obsPos, obs->triangles()(fi, 0));
              EigenSupport::V3d vt1 = obsVtx(obsPos, obs->triangles()(fi, 1));
              EigenSupport::V3d vt2 = obsVtx(obsPos, obs->triangles()(fi, 2));
              double d2 = distance::computePTSqDist(vp, vt0, vt1, vt2);
              if (d2 < dhat2 && d2 > 0.0) {
                double w = topology.vertexArea[vi] * obsTriArea[fi];
                localPairs.push_back({ obsId, vi,
                  {{ obs->triangles()(fi, 0), obs->triangles()(fi, 1), obs->triangles()(fi, 2) }},
                  w });
              }
            }
          }
        });

      for (auto &lp : tls_pairs)
        pairs.ptPairs.insert(pairs.ptPairs.end(), lp.begin(), lp.end());
    }

    // ---- External TP: obs vertex x dyn triangle ----
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

              auto &tri = topology.triangles[fi];
              EigenSupport::V3d vp = obsVtx(obsPos, ovi);
              EigenSupport::V3d vt0 = getV(tri[0]);
              EigenSupport::V3d vt1 = getV(tri[1]);
              EigenSupport::V3d vt2 = getV(tri[2]);
              double d2 = distance::computePTSqDist(vp, vt0, vt1, vt2);
              if (d2 < dhat2 && d2 > 0.0) {
                double w = 1.0 * topology.triArea[fi];
                localPairs.push_back({ obsId,
                  {{ tri[0], tri[1], tri[2] }}, ovi, w });
              }
            }
          }
        });

      for (auto &lp : tls_pairs)
        pairs.tpPairs.insert(pairs.tpPairs.end(), lp.begin(), lp.end());
    }

    // ---- External EE: dyn edge x obs edge ----
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

            int a0 = topology.edges[ei][0], a1 = topology.edges[ei][1];
            for (int ej : candidates) {
              if (!dynEdgeBox[ei].overlaps(obsEdgeBox[ej]))
                continue;

              int b0 = obs->uniqueEdges()(ej, 0);
              int b1 = obs->uniqueEdges()(ej, 1);
              EigenSupport::V3d va0 = getV(a0), va1 = getV(a1);
              EigenSupport::V3d vb0 = obsVtx(obsPos, b0), vb1 = obsVtx(obsPos, b1);
              double d2 = distance::computeEESqDist(va0, va1, vb0, vb1);
              if (d2 < dhat2 && d2 > 0.0) {
                double w = topology.edgeLength[ei] * obsEdgeLen[ej];
                localPairs.push_back({ obsId,
                  {{ a0, a1 }}, {{ b0, b1 }}, w });
              }
            }
          }
        });

      for (auto &lp : tls_pairs)
        pairs.eePairs.insert(pairs.eePairs.end(), lp.begin(), lp.end());
    }
  }
}

}  // namespace CIPC
}  // namespace Contact
}  // namespace pgo
