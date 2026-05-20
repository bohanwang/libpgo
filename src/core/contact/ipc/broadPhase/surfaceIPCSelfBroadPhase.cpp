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

}  // namespace CIPC
}  // namespace Contact
}  // namespace pgo
