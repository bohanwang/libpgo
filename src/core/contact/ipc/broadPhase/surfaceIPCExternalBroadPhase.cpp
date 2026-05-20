#include "ipc/broadPhase/surfaceIPCBroadPhase.h"

#include "ipc/broadPhase/surfaceIPCBroadPhaseInternal.h"
#include "ipc/geometry/ipcDistancePrimitives.h"
#include "ipc/profiling/surfaceIPCProfiling.h"
#include "scopedProfileSection.h"

#include <tbb/blocked_range.h>

#include <algorithm>
#include <cstdint>
#include <vector>

namespace pgo
{
namespace Contact
{
namespace CIPC
{

using namespace broad_phase_detail;

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
