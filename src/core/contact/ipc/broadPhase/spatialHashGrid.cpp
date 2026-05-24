/*
copyright to Bohan Wang
*/

#include "ipc/broadPhase/spatialHashGrid.h"

#include <algorithm>
#include <cmath>

namespace pgo
{
namespace Contact
{
namespace IPC
{

void SpatialHashGrid::AABB::init(const EigenSupport::V3d &v, double pad)
{
  lo = v.array() - pad;
  hi = v.array() + pad;
}

void SpatialHashGrid::AABB::expand(const EigenSupport::V3d &v, double pad)
{
  lo = (v.array() - pad).cwiseMin(lo.array());
  hi = (v.array() + pad).cwiseMax(hi.array());
}

void SpatialHashGrid::AABB::expand(const EigenSupport::V3d &v)
{
  lo = lo.cwiseMin(v);
  hi = hi.cwiseMax(v);
}

bool SpatialHashGrid::AABB::overlaps(const AABB &other) const
{
  return (lo.array() <= other.hi.array()).all() &&
    (other.lo.array() <= hi.array()).all();
}

SpatialHashGrid::SpatialHashGrid(int capacity)
{
  cells_.reserve(capacity);
}

void SpatialHashGrid::setCellSize(double cellSize)
{
  cellSize_ = cellSize;
}

void SpatialHashGrid::clear()
{
  cells_.clear();
}

void SpatialHashGrid::build(const std::vector<AABB> &boxes)
{
  clear();

  struct GridRange
  {
    int loIx;
    int loIy;
    int loIz;
    int hiIx;
    int hiIy;
    int hiIz;
  };

  std::vector<GridRange> ranges(boxes.size());
  std::size_t cellReferenceCount = 0;
  for (int primitiveId = 0; primitiveId < static_cast<int>(boxes.size()); ++primitiveId) {
    GridRange &range = ranges[primitiveId];
    toGrid(boxes[primitiveId].lo, range.loIx, range.loIy, range.loIz);
    toGrid(boxes[primitiveId].hi, range.hiIx, range.hiIy, range.hiIz);

    const std::size_t nx = static_cast<std::size_t>(range.hiIx - range.loIx + 1);
    const std::size_t ny = static_cast<std::size_t>(range.hiIy - range.loIy + 1);
    const std::size_t nz = static_cast<std::size_t>(range.hiIz - range.loIz + 1);
    cellReferenceCount += nx * ny * nz;
  }

  cells_.reserve(std::max(boxes.size(), cellReferenceCount));

  for (int primitiveId = 0; primitiveId < static_cast<int>(ranges.size()); ++primitiveId) {
    const GridRange &range = ranges[primitiveId];
    insertGridRange(range.loIx, range.loIy, range.loIz,
      range.hiIx, range.hiIy, range.hiIz, primitiveId);
  }
}

std::uint64_t SpatialHashGrid::hashCoord(int ix, int iy, int iz)
{
  // Pack low coordinate bits into one key. Broad-phase false positives remain
  // acceptable beyond this range, but typical IPC grids avoid the dense
  // low-coordinate collisions from XOR-prime spatial hashing.
  constexpr std::uint64_t bitsPerAxis = 21;
  constexpr std::uint64_t axisMask = (1ULL << bitsPerAxis) - 1ULL;
  return ((static_cast<std::uint64_t>(ix) & axisMask) << (2 * bitsPerAxis)) |
    ((static_cast<std::uint64_t>(iy) & axisMask) << bitsPerAxis) |
    (static_cast<std::uint64_t>(iz) & axisMask);
}

void SpatialHashGrid::toGrid(const EigenSupport::V3d &p, int &ix, int &iy, int &iz) const
{
  ix = static_cast<int>(std::floor(p.x() / cellSize_));
  iy = static_cast<int>(std::floor(p.y() / cellSize_));
  iz = static_cast<int>(std::floor(p.z() / cellSize_));
}

void SpatialHashGrid::insert(const AABB &box, int primitiveId)
{
  int lo_ix, lo_iy, lo_iz, hi_ix, hi_iy, hi_iz;
  toGrid(box.lo, lo_ix, lo_iy, lo_iz);
  toGrid(box.hi, hi_ix, hi_iy, hi_iz);

  insertGridRange(lo_ix, lo_iy, lo_iz, hi_ix, hi_iy, hi_iz, primitiveId);
}

void SpatialHashGrid::insertGridRange(int loIx, int loIy, int loIz,
  int hiIx, int hiIy, int hiIz, int primitiveId)
{
  for (int iz = loIz; iz <= hiIz; ++iz)
    for (int iy = loIy; iy <= hiIy; ++iy)
      for (int ix = loIx; ix <= hiIx; ++ix)
        cells_[hashCoord(ix, iy, iz)].push_back(primitiveId);
}

void SpatialHashGrid::query(const AABB &box, int selfPrimitiveId,
  std::vector<int> &visitedStamp, int stamp,
  std::vector<int> &result) const
{
  queryFiltered(box, selfPrimitiveId, -1, visitedStamp, stamp, result);
}

void SpatialHashGrid::queryAfter(const AABB &box, int minPrimitiveId,
  std::vector<int> &visitedStamp, int stamp,
  std::vector<int> &result) const
{
  queryFiltered(box, -1, minPrimitiveId, visitedStamp, stamp, result);
}

std::uint64_t SpatialHashGrid::queryOverlapping(const AABB &box, const std::vector<AABB> &candidateBoxes, int selfPrimitiveId,
  std::vector<int> &visitedStamp, int stamp,
  std::vector<int> &result) const
{
  return queryFilteredOverlapping(box, candidateBoxes, selfPrimitiveId, -1, visitedStamp, stamp, result);
}

std::uint64_t SpatialHashGrid::queryOverlappingAfter(
  const AABB &box, const std::vector<AABB> &candidateBoxes, int minPrimitiveId,
  std::vector<int> &visitedStamp, int stamp,
  std::vector<int> &result) const
{
  return queryFilteredOverlapping(box, candidateBoxes, -1, minPrimitiveId, visitedStamp, stamp, result);
}

void SpatialHashGrid::queryFiltered(const AABB &box, int selfPrimitiveId, int minPrimitiveId,
  std::vector<int> &visitedStamp, int stamp,
  std::vector<int> &result) const
{
  int lo_ix, lo_iy, lo_iz, hi_ix, hi_iy, hi_iz;
  toGrid(box.lo, lo_ix, lo_iy, lo_iz);
  toGrid(box.hi, hi_ix, hi_iy, hi_iz);

  for (int iz = lo_iz; iz <= hi_iz; ++iz)
    for (int iy = lo_iy; iy <= hi_iy; ++iy)
      for (int ix = lo_ix; ix <= hi_ix; ++ix) {
        auto it = cells_.find(hashCoord(ix, iy, iz));
        if (it == cells_.end())
          continue;
        for (int idx : it->second) {
          if (idx <= minPrimitiveId)
            continue;
          if (idx == selfPrimitiveId)
            continue;
          if (visitedStamp[idx] == stamp)
            continue;
          visitedStamp[idx] = stamp;
          result.push_back(idx);
        }
      }
}

std::uint64_t SpatialHashGrid::queryFilteredOverlapping(
  const AABB &box, const std::vector<AABB> &candidateBoxes, int selfPrimitiveId, int minPrimitiveId,
  std::vector<int> &visitedStamp, int stamp,
  std::vector<int> &result) const
{
  int lo_ix, lo_iy, lo_iz, hi_ix, hi_iy, hi_iz;
  toGrid(box.lo, lo_ix, lo_iy, lo_iz);
  toGrid(box.hi, hi_ix, hi_iy, hi_iz);

  std::uint64_t hashCandidates = 0;
  for (int iz = lo_iz; iz <= hi_iz; ++iz)
    for (int iy = lo_iy; iy <= hi_iy; ++iy)
      for (int ix = lo_ix; ix <= hi_ix; ++ix) {
        auto it = cells_.find(hashCoord(ix, iy, iz));
        if (it == cells_.end())
          continue;
        for (int idx : it->second) {
          if (idx <= minPrimitiveId)
            continue;
          if (idx == selfPrimitiveId)
            continue;
          if (visitedStamp[idx] == stamp)
            continue;
          visitedStamp[idx] = stamp;
          ++hashCandidates;
          if (box.overlaps(candidateBoxes[idx]))
            result.push_back(idx);
        }
      }

  return hashCandidates;
}

}  // namespace IPC
}  // namespace Contact
}  // namespace pgo
