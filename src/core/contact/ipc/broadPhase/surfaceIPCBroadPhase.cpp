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
//  Helper: build AABBs (parallel)
// =========================================================================

template <typename GetV>
static void buildVertexAABBs(
  std::vector<SpatialHashGrid::AABB> &boxes, int n,
  GetV &&getV, double inflate)
{
  tbb::parallel_for(tbb::blocked_range<int>(0, n),
    [&](const tbb::blocked_range<int> &r) {
      for (int i = r.begin(); i < r.end(); ++i)
        boxes[i].init(getV(i), inflate);
    });
}

template <typename GetV, typename Triangles>
static void buildTriangleAABBs(
  std::vector<SpatialHashGrid::AABB> &boxes, int n,
  const Triangles &triangles, GetV &&getV, double inflate)
{
  tbb::parallel_for(tbb::blocked_range<int>(0, n),
    [&](const tbb::blocked_range<int> &r) {
      for (int fi = r.begin(); fi < r.end(); ++fi) {
        auto &tri = triangles[fi];
        boxes[fi].init(getV(tri[0]), inflate);
        boxes[fi].expand(getV(tri[1]), inflate);
        boxes[fi].expand(getV(tri[2]), inflate);
      }
    });
}

template <typename GetV, typename Edges>
static void buildEdgeAABBs(
  std::vector<SpatialHashGrid::AABB> &boxes, int n,
  const Edges &edges, GetV &&getV, double inflate)
{
  tbb::parallel_for(tbb::blocked_range<int>(0, n),
    [&](const tbb::blocked_range<int> &r) {
      for (int ei = r.begin(); ei < r.end(); ++ei) {
        boxes[ei].init(getV(edges[ei][0]), inflate);
        boxes[ei].expand(getV(edges[ei][1]), inflate);
      }
    });
}

// =========================================================================
//  Helper: TLS parallel query + collect + merge
// =========================================================================

template <typename PairType, typename Body>
static void collectPairsParallel(
  int nTarget, int queryBegin, int queryEnd,
  Body &&body, std::vector<PairType> &outputPairs)
{
  tbb::enumerable_thread_specific<std::vector<int>> tls_visited(
    [nTarget]() { return std::vector<int>(nTarget, 0); });
  tbb::enumerable_thread_specific<std::vector<int>> tls_candidates;
  tbb::enumerable_thread_specific<std::vector<PairType>> tls_pairs;

  tbb::parallel_for(tbb::blocked_range<int>(queryBegin, queryEnd),
    [&](const tbb::blocked_range<int> &range) {
      auto &visited = tls_visited.local();
      auto &candidates = tls_candidates.local();
      auto &localPairs = tls_pairs.local();
      body(range, visited, candidates, localPairs);
    });

  for (auto &lp : tls_pairs)
    outputPairs.insert(outputPairs.end(), lp.begin(), lp.end());
}

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

  // --- Build AABBs ---
  std::vector<SpatialHashGrid::AABB> vertBox(topology.numVerts);
  std::vector<SpatialHashGrid::AABB> triBox(nTri);
  std::vector<SpatialHashGrid::AABB> edgeBox(nEdge);
  {
    Profiling::ScopedProfileSection scopedProfile(SurfaceIPCProfileSections::kPairBuildSelfAABB);
    buildVertexAABBs(vertBox, topology.numVerts, getV, inflate);
    buildTriangleAABBs(triBox, nTri, topology.triangles, getV, inflate);
    buildEdgeAABBs(edgeBox, nEdge, topology.edges, getV, inflate);
  }

  double avgBoxDiag = 0;
  for (const auto &aabb : triBox) {
    avgBoxDiag += (aabb.hi - aabb.lo).norm();
  }

  double cellSize = nTri > 0 ? std::max(avgBoxDiag / nTri, 1e-6) : std::max(1e-6, dhat);
  double dhat2 = dhat * dhat;

  // --- PT pairs ---
  {
    Profiling::ScopedProfileSection scopedProfile(SurfaceIPCProfileSections::kPairBuildSelfPTHashQuery);
    SpatialHashGrid triHash(nTri);
    triHash.setCellSize(cellSize);
    for (int fi = 0; fi < nTri; ++fi)
      triHash.insert(triBox[fi], fi);

    collectPairsParallel<PTPair>(nTri, 0, topology.numVerts,
      [&](const tbb::blocked_range<int> &range,
          std::vector<int> &visited,
          std::vector<int> &candidates,
          std::vector<PTPair> &localPairs) {
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
      }, pairs.ptPairs);
  }

  // --- EE pairs ---
  {
    Profiling::ScopedProfileSection scopedProfile(SurfaceIPCProfileSections::kPairBuildSelfEEHashQuery);
    SpatialHashGrid edgeHash(nEdge);
    edgeHash.setCellSize(cellSize);
    for (int ei = 0; ei < nEdge; ++ei)
      edgeHash.insert(edgeBox[ei], ei);

    collectPairsParallel<EEPair>(nEdge, 0, nEdge,
      [&](const tbb::blocked_range<int> &range,
          std::vector<int> &visited,
          std::vector<int> &candidates,
          std::vector<EEPair> &localPairs) {
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
      }, pairs.eePairs);
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
  const std::vector<ObstacleSurface> &obstacles,
  double dhatExternal,
  ExternalPairSet &pairs)
{
  pairs.clear();

  if (obstacles.empty())
    return;

  Profiling::ScopedProfileSection scopedExternalProfile(SurfaceIPCProfileSections::kPairBuildExternal);

  const double inflate = dhatExternal;
  const double dhat2 = dhatExternal * dhatExternal;

  auto getV = [&](int i) -> EigenSupport::V3d {
    return positions.segment<3>(3 * i);
  };

  int nDynTri = (int)topology.triangles.size();
  int nDynEdge = (int)topology.edges.size();

  // Build dynamic-side AABBs
  std::vector<SpatialHashGrid::AABB> dynVertBox(topology.numVerts);
  std::vector<SpatialHashGrid::AABB> dynTriBox(nDynTri);
  std::vector<SpatialHashGrid::AABB> dynEdgeBox(nDynEdge);
  {
    Profiling::ScopedProfileSection scopedProfile(SurfaceIPCProfileSections::kPairBuildExternalAABB);
    buildVertexAABBs(dynVertBox, topology.numVerts, getV, inflate);
    buildTriangleAABBs(dynTriBox, nDynTri, topology.triangles, getV, inflate);
    buildEdgeAABBs(dynEdgeBox, nDynEdge, topology.edges, getV, inflate);
  }

  for (const auto &obs : obstacles) {
    const EigenSupport::VXd &obsPos = obs.currentPositions();
    int nObsVert = (int)obsPos.size() / 3;
    int nObsTri = (int)obs.triangles().rows();
    int nObsEdge = (int)obs.uniqueEdges().rows();

    int32_t obsId = obs.objectId();

    const ObstaclePoseCache &poseCache = obs.cache();
    const std::vector<double> &obsTriArea = poseCache.triAreas;
    const std::vector<double> &obsEdgeLen = poseCache.edgeLengths;

    const std::vector<SpatialHashGrid::AABB> &obsVertBox = poseCache.vertBoxes;
    const std::vector<SpatialHashGrid::AABB> &obsTriBox = poseCache.triBoxes;
    const std::vector<SpatialHashGrid::AABB> &obsEdgeBox = poseCache.edgeBoxes;
    const double cellSize = poseCache.cellSize > 0.0 ? poseCache.cellSize : std::max(1e-6, dhatExternal);

    // ---- External PT: dyn vertex x obs triangle ----
    {
      Profiling::ScopedProfileSection scopedProfile(SurfaceIPCProfileSections::kPairBuildExternalPT);
      const SpatialHashGrid &obsTriHash = poseCache.triHash;

      collectPairsParallel<ExternalPTPair>(nObsTri, 0, topology.numVerts,
        [&](const tbb::blocked_range<int> &range,
            std::vector<int> &visited,
            std::vector<int> &candidates,
            std::vector<ExternalPTPair> &localPairs) {
          for (int vi = range.begin(); vi < range.end(); ++vi) {
            candidates.clear();
            obsTriHash.query(dynVertBox[vi], -1, visited, vi + 1, candidates);

            for (int fi : candidates) {
              if (!dynVertBox[vi].overlaps(obsTriBox[fi]))
                continue;

              EigenSupport::V3d vp = getV(vi);
              EigenSupport::V3d vt0 = obsVtx(obsPos, obs.triangles()(fi, 0));
              EigenSupport::V3d vt1 = obsVtx(obsPos, obs.triangles()(fi, 1));
              EigenSupport::V3d vt2 = obsVtx(obsPos, obs.triangles()(fi, 2));
              double d2 = distance::computePTSqDist(vp, vt0, vt1, vt2);
              if (d2 < dhat2 && d2 > 0.0) {
                double w = topology.vertexArea[vi] * obsTriArea[fi];
                localPairs.push_back({ obsId, vi,
                  {{ obs.triangles()(fi, 0), obs.triangles()(fi, 1), obs.triangles()(fi, 2) }},
                  w });
              }
            }
          }
        }, pairs.ptPairs);
    }

    // ---- External TP: obs vertex x dyn triangle ----
    {
      Profiling::ScopedProfileSection scopedProfile(SurfaceIPCProfileSections::kPairBuildExternalTP);
      SpatialHashGrid dynTriHash(nDynTri);
      dynTriHash.setCellSize(cellSize);
      for (int fi = 0; fi < nDynTri; ++fi)
        dynTriHash.insert(dynTriBox[fi], fi);

      collectPairsParallel<ExternalTPPair>(nDynTri, 0, nObsVert,
        [&](const tbb::blocked_range<int> &range,
            std::vector<int> &visited,
            std::vector<int> &candidates,
            std::vector<ExternalTPPair> &localPairs) {
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
        }, pairs.tpPairs);
    }

    // ---- External EE: dyn edge x obs edge ----
    {
      Profiling::ScopedProfileSection scopedProfile(SurfaceIPCProfileSections::kPairBuildExternalEE);
      const SpatialHashGrid &obsEdgeHash = poseCache.edgeHash;

      collectPairsParallel<ExternalEEPair>(nObsEdge, 0, nDynEdge,
        [&](const tbb::blocked_range<int> &range,
            std::vector<int> &visited,
            std::vector<int> &candidates,
            std::vector<ExternalEEPair> &localPairs) {
          for (int ei = range.begin(); ei < range.end(); ++ei) {
            candidates.clear();
            obsEdgeHash.query(dynEdgeBox[ei], -1, visited, ei + 1, candidates);

            int a0 = topology.edges[ei][0], a1 = topology.edges[ei][1];
            for (int ej : candidates) {
              if (!dynEdgeBox[ei].overlaps(obsEdgeBox[ej]))
                continue;

              int b0 = obs.uniqueEdges()(ej, 0);
              int b1 = obs.uniqueEdges()(ej, 1);
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
        }, pairs.eePairs);
    }
  }
}

}  // namespace CIPC
}  // namespace Contact
}  // namespace pgo
