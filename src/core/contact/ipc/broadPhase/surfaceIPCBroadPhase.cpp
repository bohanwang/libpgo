#include "ipc/broadPhase/surfaceIPCBroadPhase.h"

#include "scopedProfileSection.h"
#include "ipc/profiling/surfaceIPCProfiling.h"
#include "ipc/broadPhase/spatialHashGrid.h"
#include "ipc/geometry/ipcDistancePrimitives.h"

#include <tbb/blocked_range.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>

#include <algorithm>
#include <cstdint>
#include <functional>
#include <string_view>
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

template<typename GetV>
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

template<typename GetV, typename Triangles>
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

template<typename GetV, typename Edges>
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

template<typename GetV, typename GetDV>
static void buildSweptVertexAABBs(
  std::vector<SpatialHashGrid::AABB> &boxes, int n,
  GetV &&getV, GetDV &&getDV, double inflate)
{
  tbb::parallel_for(tbb::blocked_range<int>(0, n),
    [&](const tbb::blocked_range<int> &r) {
      for (int i = r.begin(); i < r.end(); ++i) {
        const EigenSupport::V3d v0 = getV(i);
        boxes[i].init(v0, inflate);
        boxes[i].expand(v0 + getDV(i), inflate);
      }
    });
}

template<typename GetV, typename GetDV, typename Triangles>
static void buildSweptTriangleAABBs(
  std::vector<SpatialHashGrid::AABB> &boxes, int n,
  const Triangles &triangles, GetV &&getV, GetDV &&getDV, double inflate)
{
  tbb::parallel_for(tbb::blocked_range<int>(0, n),
    [&](const tbb::blocked_range<int> &r) {
      for (int fi = r.begin(); fi < r.end(); ++fi) {
        auto &tri = triangles[fi];
        const EigenSupport::V3d v0 = getV(tri[0]);
        const EigenSupport::V3d v1 = getV(tri[1]);
        const EigenSupport::V3d v2 = getV(tri[2]);
        boxes[fi].init(v0, inflate);
        boxes[fi].expand(v1, inflate);
        boxes[fi].expand(v2, inflate);
        boxes[fi].expand(v0 + getDV(tri[0]), inflate);
        boxes[fi].expand(v1 + getDV(tri[1]), inflate);
        boxes[fi].expand(v2 + getDV(tri[2]), inflate);
      }
    });
}

template<typename GetV, typename GetDV, typename Edges>
static void buildSweptEdgeAABBs(
  std::vector<SpatialHashGrid::AABB> &boxes, int n,
  const Edges &edges, GetV &&getV, GetDV &&getDV, double inflate)
{
  tbb::parallel_for(tbb::blocked_range<int>(0, n),
    [&](const tbb::blocked_range<int> &r) {
      for (int ei = r.begin(); ei < r.end(); ++ei) {
        const EigenSupport::V3d v0 = getV(edges[ei][0]);
        const EigenSupport::V3d v1 = getV(edges[ei][1]);
        boxes[ei].init(v0, inflate);
        boxes[ei].expand(v1, inflate);
        boxes[ei].expand(v0 + getDV(edges[ei][0]), inflate);
        boxes[ei].expand(v1 + getDV(edges[ei][1]), inflate);
      }
    });
}

static bool computeUnionAABB(
  const std::vector<SpatialHashGrid::AABB> &boxes,
  SpatialHashGrid::AABB &out)
{
  if (boxes.empty())
    return false;

  out = boxes.front();
  for (std::size_t i = 1; i < boxes.size(); ++i) {
    out.lo = out.lo.cwiseMin(boxes[i].lo);
    out.hi = out.hi.cwiseMax(boxes[i].hi);
  }
  return true;
}

struct PairQueryCounts
{
  std::uint64_t hashCandidates = 0;
  std::uint64_t exactTests = 0;
};

static void addCounts(PairQueryCounts &dst, const PairQueryCounts &src)
{
  dst.hashCandidates += src.hashCandidates;
  dst.exactTests += src.exactTests;
}

static void recordPairQueryCounters(
  std::string_view hashCandidateName,
  std::string_view exactTestName,
  std::string_view acceptedPairName,
  const PairQueryCounts &counts,
  std::size_t acceptedPairCount)
{
  if (!Profiling::isProfilingEnabled())
    return;

  Profiling::recordProfileCounter(hashCandidateName, counts.hashCandidates);
  Profiling::recordProfileCounter(exactTestName, counts.exactTests);
  Profiling::recordProfileCounter(acceptedPairName, static_cast<std::uint64_t>(acceptedPairCount));
}

// =========================================================================
//  Helper: TLS parallel query + collect + merge
// =========================================================================

template<typename PairType, typename Body>
static PairQueryCounts collectPairsParallel(
  int nTarget, int queryBegin, int queryEnd,
  Body &&body, std::vector<PairType> &outputPairs)
{
  tbb::enumerable_thread_specific<std::vector<int>> tls_visited(
    [nTarget]() { return std::vector<int>(nTarget, 0); });
  tbb::enumerable_thread_specific<std::vector<int>> tls_candidates;
  tbb::enumerable_thread_specific<std::vector<PairType>> tls_pairs;
  tbb::enumerable_thread_specific<PairQueryCounts> tls_counts;

  tbb::parallel_for(tbb::blocked_range<int>(queryBegin, queryEnd),
    [&](const tbb::blocked_range<int> &range) {
      auto &visited = tls_visited.local();
      auto &candidates = tls_candidates.local();
      auto &localPairs = tls_pairs.local();
      auto &localCounts = tls_counts.local();
      body(range, visited, candidates, localPairs, localCounts);
    });

  for (auto &lp : tls_pairs)
    outputPairs.insert(outputPairs.end(), lp.begin(), lp.end());

  PairQueryCounts counts;
  for (const auto &localCounts : tls_counts)
    addCounts(counts, localCounts);
  return counts;
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
  const bool profilingEnabled = Profiling::isProfilingEnabled();

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
    {
      Profiling::ScopedProfileSection hashProfile(SurfaceIPCProfileSections::kPairBuildSelfPTHashInsert);
      triHash.build(triBox);
    }
    {
      Profiling::ScopedProfileSection queryProfile(SurfaceIPCProfileSections::kPairBuildSelfPTQuery);
      const PairQueryCounts counts = collectPairsParallel<PTPair>(nTri, 0, topology.numVerts,
        [&](const tbb::blocked_range<int> &range,
          std::vector<int> &visited,
          std::vector<int> &candidates,
          std::vector<PTPair> &localPairs,
          PairQueryCounts &localCounts) {
          for (int vi = range.begin(); vi < range.end(); ++vi) {
            candidates.clear();
            const std::uint64_t hashCandidates =
              triHash.queryOverlapping(vertBox[vi], triBox, -1, visited, vi + 1, candidates);
            if (profilingEnabled)
              localCounts.hashCandidates += hashCandidates;

            V3d vp = getV(vi);
            for (int fi : candidates) {
              auto &tri = topology.triangles[fi];
              if (vi == tri[0] || vi == tri[1] || vi == tri[2])
                continue;

              if (profilingEnabled)
                localCounts.exactTests += 1;
              V3d vt0 = getV(tri[0]), vt1 = getV(tri[1]), vt2 = getV(tri[2]);
              double d2 = distance::computePTSqDist(vp, vt0, vt1, vt2);
              if (d2 < dhat2)
                localPairs.push_back({ vi, tri[0], tri[1], tri[2],
                  topology.vertexArea[vi] * topology.triArea[fi] });
            }
          }
        },
        pairs.ptPairs);
      recordPairQueryCounters(
        SurfaceIPCProfileSections::kPairBuildSelfPTHashCandidates,
        SurfaceIPCProfileSections::kPairBuildSelfPTDistanceTests,
        SurfaceIPCProfileSections::kPairBuildSelfPTAcceptedPairs,
        counts, pairs.ptPairs.size());
    }
  }

  // --- EE pairs ---
  {
    Profiling::ScopedProfileSection scopedProfile(SurfaceIPCProfileSections::kPairBuildSelfEEHashQuery);
    SpatialHashGrid edgeHash(nEdge);
    edgeHash.setCellSize(cellSize);
    {
      Profiling::ScopedProfileSection hashProfile(SurfaceIPCProfileSections::kPairBuildSelfEEHashInsert);
      edgeHash.build(edgeBox);
    }
    {
      Profiling::ScopedProfileSection queryProfile(SurfaceIPCProfileSections::kPairBuildSelfEEQuery);
      const PairQueryCounts counts = collectPairsParallel<EEPair>(nEdge, 0, nEdge,
        [&](const tbb::blocked_range<int> &range,
          std::vector<int> &visited,
          std::vector<int> &candidates,
          std::vector<EEPair> &localPairs,
          PairQueryCounts &localCounts) {
          for (int ei = range.begin(); ei < range.end(); ++ei) {
            candidates.clear();
            const std::uint64_t hashCandidates =
              edgeHash.queryOverlappingAfter(edgeBox[ei], edgeBox, ei, visited, ei + 1, candidates);
            if (profilingEnabled)
              localCounts.hashCandidates += hashCandidates;

            int a0 = topology.edges[ei][0], a1 = topology.edges[ei][1];
            V3d va0 = getV(a0), va1 = getV(a1);
            for (int ej : candidates) {
              int b0 = topology.edges[ej][0], b1 = topology.edges[ej][1];
              if (a0 == b0 || a0 == b1 || a1 == b0 || a1 == b1)
                continue;

              if (profilingEnabled)
                localCounts.exactTests += 1;
              V3d vb0 = getV(b0), vb1 = getV(b1);
              double d2 = distance::computeEESqDist(va0, va1, vb0, vb1);
              if (d2 < dhat2)
                localPairs.push_back({ a0, a1, b0, b1,
                  topology.edgeLength[ei] * topology.edgeLength[ej] });
            }
          }
        },
        pairs.eePairs);
      recordPairQueryCounters(
        SurfaceIPCProfileSections::kPairBuildSelfEEHashCandidates,
        SurfaceIPCProfileSections::kPairBuildSelfEEDistanceTests,
        SurfaceIPCProfileSections::kPairBuildSelfEEAcceptedPairs,
        counts, pairs.eePairs.size());
    }
  }
}

void buildSelfPairsLineSearchSuperset(
  const SurfaceIPCTopology &topology,
  EigenSupport::ConstRefVecXd positions,
  EigenSupport::ConstRefVecXd displacements,
  double dhat,
  SelfPairSet &pairs)
{
  using namespace pgo::EigenSupport;

  Profiling::ScopedProfileSection scopedProfile(SurfaceIPCProfileSections::kPairBuildStatic);
  const bool profilingEnabled = Profiling::isProfilingEnabled();

  pairs.clear();

  auto getV = [&](int i) -> V3d {
    return positions.segment<3>(3 * i);
  };
  auto getDV = [&](int i) -> V3d {
    return displacements.segment<3>(3 * i);
  };

  const int nTri = static_cast<int>(topology.triangles.size());
  const int nEdge = static_cast<int>(topology.edges.size());

  std::vector<SpatialHashGrid::AABB> vertBox(topology.numVerts);
  std::vector<SpatialHashGrid::AABB> triBox(nTri);
  std::vector<SpatialHashGrid::AABB> edgeBox(nEdge);
  {
    Profiling::ScopedProfileSection scopedProfile(SurfaceIPCProfileSections::kPairBuildSelfAABB);
    buildSweptVertexAABBs(vertBox, topology.numVerts, getV, getDV, dhat);
    buildSweptTriangleAABBs(triBox, nTri, topology.triangles, getV, getDV, dhat);
    buildSweptEdgeAABBs(edgeBox, nEdge, topology.edges, getV, getDV, dhat);
  }

  double avgBoxDiag = 0.0;
  for (const auto &aabb : triBox)
    avgBoxDiag += (aabb.hi - aabb.lo).norm();
  const double cellSize = nTri > 0 ? std::max(avgBoxDiag / nTri, 1e-6) : std::max(1e-6, dhat);

  {
    Profiling::ScopedProfileSection scopedProfile(SurfaceIPCProfileSections::kPairBuildSelfPTHashQuery);
    SpatialHashGrid triHash(nTri);
    triHash.setCellSize(cellSize);
    {
      Profiling::ScopedProfileSection hashProfile(SurfaceIPCProfileSections::kPairBuildSelfPTHashInsert);
      triHash.build(triBox);
    }
    {
      Profiling::ScopedProfileSection queryProfile(SurfaceIPCProfileSections::kPairBuildSelfPTQuery);
      const PairQueryCounts counts = collectPairsParallel<PTPair>(nTri, 0, topology.numVerts,
        [&](const tbb::blocked_range<int> &range,
          std::vector<int> &visited,
          std::vector<int> &candidates,
          std::vector<PTPair> &localPairs,
          PairQueryCounts &localCounts) {
          for (int vi = range.begin(); vi < range.end(); ++vi) {
            candidates.clear();
            const std::uint64_t hashCandidates =
              triHash.queryOverlapping(vertBox[vi], triBox, -1, visited, vi + 1, candidates);
            if (profilingEnabled)
              localCounts.hashCandidates += hashCandidates;

            for (int fi : candidates) {
              auto &tri = topology.triangles[fi];
              if (vi == tri[0] || vi == tri[1] || vi == tri[2])
                continue;
              localPairs.push_back({ vi, tri[0], tri[1], tri[2],
                topology.vertexArea[vi] * topology.triArea[fi] });
            }
          }
        },
        pairs.ptPairs);
      recordPairQueryCounters(
        SurfaceIPCProfileSections::kPairBuildSelfPTHashCandidates,
        SurfaceIPCProfileSections::kPairBuildSelfPTDistanceTests,
        SurfaceIPCProfileSections::kPairBuildSelfPTAcceptedPairs,
        counts, pairs.ptPairs.size());
    }
  }

  {
    Profiling::ScopedProfileSection scopedProfile(SurfaceIPCProfileSections::kPairBuildSelfEEHashQuery);
    SpatialHashGrid edgeHash(nEdge);
    edgeHash.setCellSize(cellSize);
    {
      Profiling::ScopedProfileSection hashProfile(SurfaceIPCProfileSections::kPairBuildSelfEEHashInsert);
      edgeHash.build(edgeBox);
    }
    {
      Profiling::ScopedProfileSection queryProfile(SurfaceIPCProfileSections::kPairBuildSelfEEQuery);
      const PairQueryCounts counts = collectPairsParallel<EEPair>(nEdge, 0, nEdge,
        [&](const tbb::blocked_range<int> &range,
          std::vector<int> &visited,
          std::vector<int> &candidates,
          std::vector<EEPair> &localPairs,
          PairQueryCounts &localCounts) {
          for (int ei = range.begin(); ei < range.end(); ++ei) {
            candidates.clear();
            const std::uint64_t hashCandidates =
              edgeHash.queryOverlappingAfter(edgeBox[ei], edgeBox, ei, visited, ei + 1, candidates);
            if (profilingEnabled)
              localCounts.hashCandidates += hashCandidates;

            const int a0 = topology.edges[ei][0];
            const int a1 = topology.edges[ei][1];
            for (int ej : candidates) {
              const int b0 = topology.edges[ej][0];
              const int b1 = topology.edges[ej][1];
              if (a0 == b0 || a0 == b1 || a1 == b0 || a1 == b1)
                continue;
              localPairs.push_back({ a0, a1, b0, b1,
                topology.edgeLength[ei] * topology.edgeLength[ej] });
            }
          }
        },
        pairs.eePairs);
      recordPairQueryCounters(
        SurfaceIPCProfileSections::kPairBuildSelfEEHashCandidates,
        SurfaceIPCProfileSections::kPairBuildSelfEEDistanceTests,
        SurfaceIPCProfileSections::kPairBuildSelfEEAcceptedPairs,
        counts, pairs.eePairs.size());
    }
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
  const bool profilingEnabled = Profiling::isProfilingEnabled();

  const double inflate = dhatExternal;
  const double dhat2 = dhatExternal * dhatExternal;

  auto getV = [&](int i) -> EigenSupport::V3d {
    return positions.segment<3>(3 * i);
  };

  int nDynTri = (int)topology.triangles.size();
  int nDynEdge = (int)topology.edges.size();

  std::vector<SpatialHashGrid::AABB> dynVertBox(topology.numVerts);
  {
    Profiling::ScopedProfileSection scopedProfile(SurfaceIPCProfileSections::kPairBuildExternalAABB);
    buildVertexAABBs(dynVertBox, topology.numVerts, getV, inflate);
  }

  SpatialHashGrid::AABB dynSurfaceBox;
  const bool hasDynSurfaceBox = computeUnionAABB(dynVertBox, dynSurfaceBox);

  std::vector<const ObstacleSurface *> overlappingObstacles;
  overlappingObstacles.reserve(obstacles.size());
  for (const auto &obs : obstacles) {
    const ObstaclePoseCache &poseCache = obs.cache();
    if (hasDynSurfaceBox && poseCache.hasSurfaceBox && !dynSurfaceBox.overlaps(poseCache.surfaceBox))
      continue;
    overlappingObstacles.push_back(&obs);
  }

  if (profilingEnabled)
    Profiling::recordProfileCounter(SurfaceIPCProfileSections::kPairBuildExternalOverlappingObstacles,
      static_cast<std::uint64_t>(overlappingObstacles.size()));

  if (overlappingObstacles.empty()) {
    recordPairQueryCounters(
      SurfaceIPCProfileSections::kPairBuildExternalPTHashCandidates,
      SurfaceIPCProfileSections::kPairBuildExternalPTDistanceTests,
      SurfaceIPCProfileSections::kPairBuildExternalPTAcceptedPairs,
      PairQueryCounts{}, 0);
    recordPairQueryCounters(
      SurfaceIPCProfileSections::kPairBuildExternalTPHashCandidates,
      SurfaceIPCProfileSections::kPairBuildExternalTPDistanceTests,
      SurfaceIPCProfileSections::kPairBuildExternalTPAcceptedPairs,
      PairQueryCounts{}, 0);
    recordPairQueryCounters(
      SurfaceIPCProfileSections::kPairBuildExternalEEHashCandidates,
      SurfaceIPCProfileSections::kPairBuildExternalEEDistanceTests,
      SurfaceIPCProfileSections::kPairBuildExternalEEAcceptedPairs,
      PairQueryCounts{}, 0);
    return;
  }

  std::vector<SpatialHashGrid::AABB> dynTriBox(nDynTri);
  std::vector<SpatialHashGrid::AABB> dynEdgeBox(nDynEdge);
  {
    Profiling::ScopedProfileSection scopedProfile(SurfaceIPCProfileSections::kPairBuildExternalAABB);
    buildTriangleAABBs(dynTriBox, nDynTri, topology.triangles, getV, inflate);
    buildEdgeAABBs(dynEdgeBox, nDynEdge, topology.edges, getV, inflate);
  }

  for (const ObstacleSurface *obsPtr : overlappingObstacles) {
    const ObstacleSurface &obs = *obsPtr;
    const EigenSupport::VXd &obsPos = obs.currentPositions();
    int nObsVert = (int)obsPos.size() / 3;
    int nObsTri = (int)obs.triangles().rows();
    const EigenSupport::MXi &obsContactEdges = obs.contactEdges();
    int nObsEdge = (int)obsContactEdges.rows();

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

      const std::size_t acceptedBefore = pairs.ptPairs.size();
      const PairQueryCounts counts = collectPairsParallel<ExternalPTPair>(nObsTri, 0, topology.numVerts,
        [&](const tbb::blocked_range<int> &range,
          std::vector<int> &visited,
          std::vector<int> &candidates,
          std::vector<ExternalPTPair> &localPairs,
          PairQueryCounts &localCounts) {
          for (int vi = range.begin(); vi < range.end(); ++vi) {
            candidates.clear();
            const std::uint64_t hashCandidates =
              obsTriHash.queryOverlapping(dynVertBox[vi], obsTriBox, -1, visited, vi + 1, candidates);
            if (profilingEnabled)
              localCounts.hashCandidates += hashCandidates;

            EigenSupport::V3d vp = getV(vi);
            for (int fi : candidates) {
              if (profilingEnabled)
                localCounts.exactTests += 1;
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
        },
        pairs.ptPairs);
      recordPairQueryCounters(
        SurfaceIPCProfileSections::kPairBuildExternalPTHashCandidates,
        SurfaceIPCProfileSections::kPairBuildExternalPTDistanceTests,
        SurfaceIPCProfileSections::kPairBuildExternalPTAcceptedPairs,
        counts, pairs.ptPairs.size() - acceptedBefore);
    }

    // ---- External TP: obs vertex x dyn triangle ----
    {
      Profiling::ScopedProfileSection scopedProfile(SurfaceIPCProfileSections::kPairBuildExternalTP);
      SpatialHashGrid dynTriHash(nDynTri);
      dynTriHash.setCellSize(cellSize);
      dynTriHash.build(dynTriBox);

      const std::size_t acceptedBefore = pairs.tpPairs.size();
      const PairQueryCounts counts = collectPairsParallel<ExternalTPPair>(nDynTri, 0, nObsVert,
        [&](const tbb::blocked_range<int> &range,
          std::vector<int> &visited,
          std::vector<int> &candidates,
          std::vector<ExternalTPPair> &localPairs,
          PairQueryCounts &localCounts) {
          for (int ovi = range.begin(); ovi < range.end(); ++ovi) {
            candidates.clear();
            const std::uint64_t hashCandidates =
              dynTriHash.queryOverlapping(obsVertBox[ovi], dynTriBox, -1, visited, ovi + 1, candidates);
            if (profilingEnabled)
              localCounts.hashCandidates += hashCandidates;

            EigenSupport::V3d vp = obsVtx(obsPos, ovi);
            for (int fi : candidates) {
              if (profilingEnabled)
                localCounts.exactTests += 1;
              auto &tri = topology.triangles[fi];
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
        },
        pairs.tpPairs);
      recordPairQueryCounters(
        SurfaceIPCProfileSections::kPairBuildExternalTPHashCandidates,
        SurfaceIPCProfileSections::kPairBuildExternalTPDistanceTests,
        SurfaceIPCProfileSections::kPairBuildExternalTPAcceptedPairs,
        counts, pairs.tpPairs.size() - acceptedBefore);
    }

    // ---- External EE: dyn edge x obs edge ----
    {
      Profiling::ScopedProfileSection scopedProfile(SurfaceIPCProfileSections::kPairBuildExternalEE);
      const SpatialHashGrid &obsEdgeHash = poseCache.edgeHash;

      const std::size_t acceptedBefore = pairs.eePairs.size();
      const PairQueryCounts counts = collectPairsParallel<ExternalEEPair>(nObsEdge, 0, nDynEdge,
        [&](const tbb::blocked_range<int> &range,
          std::vector<int> &visited,
          std::vector<int> &candidates,
          std::vector<ExternalEEPair> &localPairs,
          PairQueryCounts &localCounts) {
          for (int ei = range.begin(); ei < range.end(); ++ei) {
            candidates.clear();
            const std::uint64_t hashCandidates =
              obsEdgeHash.queryOverlapping(dynEdgeBox[ei], obsEdgeBox, -1, visited, ei + 1, candidates);
            if (profilingEnabled)
              localCounts.hashCandidates += hashCandidates;

            int a0 = topology.edges[ei][0], a1 = topology.edges[ei][1];
            EigenSupport::V3d va0 = getV(a0), va1 = getV(a1);
            for (int ej : candidates) {
              if (profilingEnabled)
                localCounts.exactTests += 1;
              int b0 = obsContactEdges(ej, 0);
              int b1 = obsContactEdges(ej, 1);
              EigenSupport::V3d vb0 = obsVtx(obsPos, b0), vb1 = obsVtx(obsPos, b1);
              double d2 = distance::computeEESqDist(va0, va1, vb0, vb1);
              if (d2 < dhat2 && d2 > 0.0) {
                double w = topology.edgeLength[ei] * obsEdgeLen[ej];
                localPairs.push_back({ obsId,
                  {{ a0, a1 }}, {{ b0, b1 }}, w });
              }
            }
          }
        },
        pairs.eePairs);
      recordPairQueryCounters(
        SurfaceIPCProfileSections::kPairBuildExternalEEHashCandidates,
        SurfaceIPCProfileSections::kPairBuildExternalEEDistanceTests,
        SurfaceIPCProfileSections::kPairBuildExternalEEAcceptedPairs,
        counts, pairs.eePairs.size() - acceptedBefore);
    }
  }
}

void buildExternalPairsLineSearchSuperset(
  const SurfaceIPCTopology &topology,
  EigenSupport::ConstRefVecXd positions,
  EigenSupport::ConstRefVecXd displacements,
  const std::vector<ObstacleSurface> &obstacles,
  double dhatExternal,
  ExternalPairSet &pairs)
{
  pairs.clear();

  if (obstacles.empty())
    return;

  Profiling::ScopedProfileSection scopedExternalProfile(SurfaceIPCProfileSections::kPairBuildExternal);
  const bool profilingEnabled = Profiling::isProfilingEnabled();

  auto getV = [&](int i) -> EigenSupport::V3d {
    return positions.segment<3>(3 * i);
  };
  auto getDV = [&](int i) -> EigenSupport::V3d {
    return displacements.segment<3>(3 * i);
  };

  const int nDynTri = static_cast<int>(topology.triangles.size());
  const int nDynEdge = static_cast<int>(topology.edges.size());

  std::vector<SpatialHashGrid::AABB> dynVertBox(topology.numVerts);
  {
    Profiling::ScopedProfileSection scopedProfile(SurfaceIPCProfileSections::kPairBuildExternalAABB);
    buildSweptVertexAABBs(dynVertBox, topology.numVerts, getV, getDV, dhatExternal);
  }

  SpatialHashGrid::AABB dynSurfaceBox;
  const bool hasDynSurfaceBox = computeUnionAABB(dynVertBox, dynSurfaceBox);

  std::vector<const ObstacleSurface *> overlappingObstacles;
  overlappingObstacles.reserve(obstacles.size());
  for (const auto &obs : obstacles) {
    const ObstaclePoseCache &poseCache = obs.cache();
    if (hasDynSurfaceBox && poseCache.hasSurfaceBox && !dynSurfaceBox.overlaps(poseCache.surfaceBox))
      continue;
    overlappingObstacles.push_back(&obs);
  }

  if (profilingEnabled)
    Profiling::recordProfileCounter(SurfaceIPCProfileSections::kPairBuildExternalOverlappingObstacles,
      static_cast<std::uint64_t>(overlappingObstacles.size()));

  if (overlappingObstacles.empty()) {
    recordPairQueryCounters(
      SurfaceIPCProfileSections::kPairBuildExternalPTHashCandidates,
      SurfaceIPCProfileSections::kPairBuildExternalPTDistanceTests,
      SurfaceIPCProfileSections::kPairBuildExternalPTAcceptedPairs,
      PairQueryCounts{}, 0);
    recordPairQueryCounters(
      SurfaceIPCProfileSections::kPairBuildExternalTPHashCandidates,
      SurfaceIPCProfileSections::kPairBuildExternalTPDistanceTests,
      SurfaceIPCProfileSections::kPairBuildExternalTPAcceptedPairs,
      PairQueryCounts{}, 0);
    recordPairQueryCounters(
      SurfaceIPCProfileSections::kPairBuildExternalEEHashCandidates,
      SurfaceIPCProfileSections::kPairBuildExternalEEDistanceTests,
      SurfaceIPCProfileSections::kPairBuildExternalEEAcceptedPairs,
      PairQueryCounts{}, 0);
    return;
  }

  std::vector<SpatialHashGrid::AABB> dynTriBox(nDynTri);
  std::vector<SpatialHashGrid::AABB> dynEdgeBox(nDynEdge);
  {
    Profiling::ScopedProfileSection scopedProfile(SurfaceIPCProfileSections::kPairBuildExternalAABB);
    buildSweptTriangleAABBs(dynTriBox, nDynTri, topology.triangles, getV, getDV, dhatExternal);
    buildSweptEdgeAABBs(dynEdgeBox, nDynEdge, topology.edges, getV, getDV, dhatExternal);
  }

  for (const ObstacleSurface *obsPtr : overlappingObstacles) {
    const ObstacleSurface &obs = *obsPtr;
    const EigenSupport::VXd &obsPos = obs.currentPositions();
    const int nObsVert = static_cast<int>(obsPos.size()) / 3;
    const int nObsTri = static_cast<int>(obs.triangles().rows());
    const EigenSupport::MXi &obsContactEdges = obs.contactEdges();
    const int nObsEdge = static_cast<int>(obsContactEdges.rows());
    const int32_t obsId = obs.objectId();

    const ObstaclePoseCache &poseCache = obs.cache();
    const std::vector<double> &obsTriArea = poseCache.triAreas;
    const std::vector<double> &obsEdgeLen = poseCache.edgeLengths;

    const std::vector<SpatialHashGrid::AABB> &obsVertBox = poseCache.vertBoxes;
    const std::vector<SpatialHashGrid::AABB> &obsTriBox = poseCache.triBoxes;
    const std::vector<SpatialHashGrid::AABB> &obsEdgeBox = poseCache.edgeBoxes;
    const double cellSize = poseCache.cellSize > 0.0 ? poseCache.cellSize : std::max(1e-6, dhatExternal);

    {
      Profiling::ScopedProfileSection scopedProfile(SurfaceIPCProfileSections::kPairBuildExternalPT);
      const SpatialHashGrid &obsTriHash = poseCache.triHash;

      const std::size_t acceptedBefore = pairs.ptPairs.size();
      const PairQueryCounts counts = collectPairsParallel<ExternalPTPair>(nObsTri, 0, topology.numVerts,
        [&](const tbb::blocked_range<int> &range,
          std::vector<int> &visited,
          std::vector<int> &candidates,
          std::vector<ExternalPTPair> &localPairs,
          PairQueryCounts &localCounts) {
          for (int vi = range.begin(); vi < range.end(); ++vi) {
            candidates.clear();
            const std::uint64_t hashCandidates =
              obsTriHash.queryOverlapping(dynVertBox[vi], obsTriBox, -1, visited, vi + 1, candidates);
            if (profilingEnabled)
              localCounts.hashCandidates += hashCandidates;

            for (int fi : candidates) {
              const double w = topology.vertexArea[vi] * obsTriArea[fi];
              localPairs.push_back({ obsId, vi,
                {{ obs.triangles()(fi, 0), obs.triangles()(fi, 1), obs.triangles()(fi, 2) }},
                w });
            }
          }
        },
        pairs.ptPairs);
      recordPairQueryCounters(
        SurfaceIPCProfileSections::kPairBuildExternalPTHashCandidates,
        SurfaceIPCProfileSections::kPairBuildExternalPTDistanceTests,
        SurfaceIPCProfileSections::kPairBuildExternalPTAcceptedPairs,
        counts, pairs.ptPairs.size() - acceptedBefore);
    }

    {
      Profiling::ScopedProfileSection scopedProfile(SurfaceIPCProfileSections::kPairBuildExternalTP);
      SpatialHashGrid dynTriHash(nDynTri);
      dynTriHash.setCellSize(cellSize);
      dynTriHash.build(dynTriBox);

      const std::size_t acceptedBefore = pairs.tpPairs.size();
      const PairQueryCounts counts = collectPairsParallel<ExternalTPPair>(nDynTri, 0, nObsVert,
        [&](const tbb::blocked_range<int> &range,
          std::vector<int> &visited,
          std::vector<int> &candidates,
          std::vector<ExternalTPPair> &localPairs,
          PairQueryCounts &localCounts) {
          for (int ovi = range.begin(); ovi < range.end(); ++ovi) {
            candidates.clear();
            const std::uint64_t hashCandidates =
              dynTriHash.queryOverlapping(obsVertBox[ovi], dynTriBox, -1, visited, ovi + 1, candidates);
            if (profilingEnabled)
              localCounts.hashCandidates += hashCandidates;

            for (int fi : candidates) {
              auto &tri = topology.triangles[fi];
              const double w = topology.triArea[fi];
              localPairs.push_back({ obsId,
                {{ tri[0], tri[1], tri[2] }}, ovi, w });
            }
          }
        },
        pairs.tpPairs);
      recordPairQueryCounters(
        SurfaceIPCProfileSections::kPairBuildExternalTPHashCandidates,
        SurfaceIPCProfileSections::kPairBuildExternalTPDistanceTests,
        SurfaceIPCProfileSections::kPairBuildExternalTPAcceptedPairs,
        counts, pairs.tpPairs.size() - acceptedBefore);
    }

    {
      Profiling::ScopedProfileSection scopedProfile(SurfaceIPCProfileSections::kPairBuildExternalEE);
      const SpatialHashGrid &obsEdgeHash = poseCache.edgeHash;

      const std::size_t acceptedBefore = pairs.eePairs.size();
      const PairQueryCounts counts = collectPairsParallel<ExternalEEPair>(nObsEdge, 0, nDynEdge,
        [&](const tbb::blocked_range<int> &range,
          std::vector<int> &visited,
          std::vector<int> &candidates,
          std::vector<ExternalEEPair> &localPairs,
          PairQueryCounts &localCounts) {
          for (int ei = range.begin(); ei < range.end(); ++ei) {
            candidates.clear();
            const std::uint64_t hashCandidates =
              obsEdgeHash.queryOverlapping(dynEdgeBox[ei], obsEdgeBox, -1, visited, ei + 1, candidates);
            if (profilingEnabled)
              localCounts.hashCandidates += hashCandidates;

            const int a0 = topology.edges[ei][0];
            const int a1 = topology.edges[ei][1];
            for (int ej : candidates) {
              const int b0 = obsContactEdges(ej, 0);
              const int b1 = obsContactEdges(ej, 1);
              const double w = topology.edgeLength[ei] * obsEdgeLen[ej];
              localPairs.push_back({ obsId,
                {{ a0, a1 }}, {{ b0, b1 }}, w });
            }
          }
        },
        pairs.eePairs);
      recordPairQueryCounters(
        SurfaceIPCProfileSections::kPairBuildExternalEEHashCandidates,
        SurfaceIPCProfileSections::kPairBuildExternalEEDistanceTests,
        SurfaceIPCProfileSections::kPairBuildExternalEEAcceptedPairs,
        counts, pairs.eePairs.size() - acceptedBefore);
    }
  }
}

}  // namespace CIPC
}  // namespace Contact
}  // namespace pgo
