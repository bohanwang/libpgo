#pragma once

#include "EigenDef.h"
#include "ipc/broadPhase/spatialHashGrid.h"
#include "scopedProfileSection.h"

#include <tbb/blocked_range.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>

#include <cstddef>
#include <cstdint>
#include <string_view>
#include <vector>

namespace pgo
{
namespace Contact
{
namespace CIPC
{
namespace broad_phase_detail
{

template<typename GetV>
void buildVertexAABBs(
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
void buildTriangleAABBs(
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
void buildEdgeAABBs(
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
void buildSweptVertexAABBs(
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
void buildSweptTriangleAABBs(
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
void buildSweptEdgeAABBs(
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

inline bool computeUnionAABB(
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

inline void addCounts(PairQueryCounts &dst, const PairQueryCounts &src)
{
  dst.hashCandidates += src.hashCandidates;
  dst.exactTests += src.exactTests;
}

inline void recordPairQueryCounters(
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

template<typename PairType, typename Body>
PairQueryCounts collectPairsParallel(
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

inline EigenSupport::V3d obsVtx(const EigenSupport::VXd &pos, int i)
{
  return pos.segment<3>(3 * i);
}

}  // namespace broad_phase_detail
}  // namespace CIPC
}  // namespace Contact
}  // namespace pgo
