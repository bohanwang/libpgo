#include "obstaclePoseCache.h"

#include <algorithm>

namespace pgo
{
namespace Contact
{
namespace CIPC
{

void buildObstaclePoseCache(
  const EigenSupport::VXd &positions,
  const EigenSupport::MXi &triangles,
  const EigenSupport::MXi &uniqueEdges,
  ObstaclePoseCache &cache)
{
  const int nVerts = static_cast<int>(positions.size() / 3);
  const int nTri = static_cast<int>(triangles.rows());
  const int nEdge = static_cast<int>(uniqueEdges.rows());

  cache.triAreas.assign(nTri, 0.0);
  cache.edgeLengths.assign(nEdge, 0.0);
  cache.vertBoxes.resize(nVerts);
  cache.triBoxes.resize(nTri);
  cache.edgeBoxes.resize(nEdge);

  // Per-vertex degenerate point AABBs.
  for (int vi = 0; vi < nVerts; ++vi) {
    const EigenSupport::V3d p = positions.segment<3>(3 * vi);
    cache.vertBoxes[vi].init(p, 0.0);
  }

  // Triangle AABBs + cached areas in a single pass; sum diagonals for cellSize.
  double diagSum = 0.0;
  for (int fi = 0; fi < nTri; ++fi) {
    const EigenSupport::V3d v0 = positions.segment<3>(3 * triangles(fi, 0));
    const EigenSupport::V3d v1 = positions.segment<3>(3 * triangles(fi, 1));
    const EigenSupport::V3d v2 = positions.segment<3>(3 * triangles(fi, 2));
    cache.triAreas[fi] = 0.5 * (v1 - v0).cross(v2 - v0).norm();
    cache.triBoxes[fi].init(v0, 0.0);
    cache.triBoxes[fi].expand(v1);
    cache.triBoxes[fi].expand(v2);
    diagSum += (cache.triBoxes[fi].hi - cache.triBoxes[fi].lo).norm();
  }

  // Edge AABBs + cached lengths.
  for (int ei = 0; ei < nEdge; ++ei) {
    const EigenSupport::V3d e0 = positions.segment<3>(3 * uniqueEdges(ei, 0));
    const EigenSupport::V3d e1 = positions.segment<3>(3 * uniqueEdges(ei, 1));
    cache.edgeLengths[ei] = (e1 - e0).norm();
    cache.edgeBoxes[ei].init(e0, 0.0);
    cache.edgeBoxes[ei].expand(e1);
  }

  cache.cellSize = nTri > 0 ? std::max(diagSum / nTri, 1e-6) : 1e-6;

  // Rebuild spatial hashes against the new pose. clear() retains bucket
  // capacity (libstdc++/libc++ behavior) so steady-state refresh is rehash-free.
  cache.triHash.clear();
  cache.triHash.setCellSize(cache.cellSize);
  for (int fi = 0; fi < nTri; ++fi)
    cache.triHash.insert(cache.triBoxes[fi], fi);

  cache.edgeHash.clear();
  cache.edgeHash.setCellSize(cache.cellSize);
  for (int ei = 0; ei < nEdge; ++ei)
    cache.edgeHash.insert(cache.edgeBoxes[ei], ei);
}

}  // namespace CIPC
}  // namespace Contact
}  // namespace pgo
