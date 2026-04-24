/*
copyright to Bohan Wang
*/

#include "surfaceIPCSelfBroadPhase.h"

#include "scopedProfileSection.h"
#include "ipc/profiling/surfaceIPCProfiling.h"
#include "../geometry/ipcDistancePrimitives.h"
#include "ipc/broadPhase/spatialHashGrid.h"

#include <tbb/blocked_range.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>

#include <algorithm>
#include <functional>
#include <vector>

namespace pgo {
namespace Contact {
namespace CIPC {
using namespace pgo::EigenSupport;

void SurfaceIPCSelfBroadPhase::buildPairs(
  const SurfaceIPCTopology &topology,
  EigenSupport::ConstRefVecXd positions,
  double dhat,
  std::vector<PTPair> &ptPairs,
  std::vector<EEPair> &eePairs) const
{
  Profiling::ScopedProfileSection scopedProfile(SurfaceIPCProfileSections::kPairBuildStatic);

  ptPairs.clear();
  eePairs.clear();

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

    // Thread-local buffers: each thread gets its own visited array,
    // candidates vector, and result pairs. The visited stamp dedup is
    // per-query (prevents same triangle appearing from multiple cells),
    // not cross-query, so thread-local is correct.
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
      ptPairs.insert(ptPairs.end(), lp.begin(), lp.end());
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
      eePairs.insert(eePairs.end(), lp.begin(), lp.end());
  }
}


}  // namespace CIPC
}  // namespace Contact
}  // namespace pgo
