/*
copyright to Bohan Wang
*/

#include "ipc/broadPhase/spatialHashGrid.h"

#include <cmath>

namespace pgo
{
namespace Contact
{
namespace CIPC
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

std::int64_t SpatialHashGrid::hashCoord(int ix, int iy, int iz)
{
  // Large primes for spatial hashing (from Teschner et al. 2003).
  constexpr std::int64_t p1 = 73856093LL;
  constexpr std::int64_t p2 = 19349663LL;
  constexpr std::int64_t p3 = 83492791LL;
  return (static_cast<std::int64_t>(ix) * p1) ^
    (static_cast<std::int64_t>(iy) * p2) ^
    (static_cast<std::int64_t>(iz) * p3);
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

  for (int iz = lo_iz; iz <= hi_iz; ++iz)
    for (int iy = lo_iy; iy <= hi_iy; ++iy)
      for (int ix = lo_ix; ix <= hi_ix; ++ix)
        cells_[hashCoord(ix, iy, iz)].push_back(primitiveId);
}

void SpatialHashGrid::query(const AABB &box, int selfPrimitiveId,
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
          if (idx == selfPrimitiveId)
            continue;
          if (visitedStamp[idx] == stamp)
            continue;
          visitedStamp[idx] = stamp;
          result.push_back(idx);
        }
      }
}

}  // namespace CIPC
}  // namespace Contact
}  // namespace pgo
