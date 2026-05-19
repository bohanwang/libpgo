/*
copyright to Bohan Wang
*/

#pragma once

#include "EigenDef.h"

#include <cstdint>
#include <unordered_map>
#include <vector>

namespace pgo
{
namespace Contact
{
namespace CIPC
{

class SpatialHashGrid
{
public:
  struct AABB
  {
    EigenSupport::V3d lo, hi;

    void init(const EigenSupport::V3d &v, double pad);
    void expand(const EigenSupport::V3d &v, double pad);
    void expand(const EigenSupport::V3d &v);
    bool overlaps(const AABB &other) const;
  };

  explicit SpatialHashGrid(int capacity);

  void setCellSize(double cellSize);
  void clear();

  void build(const std::vector<AABB> &boxes);
  void insert(const AABB &box, int primitiveId);
  void query(const AABB &box, int selfPrimitiveId,
    std::vector<int> &visitedStamp, int stamp,
    std::vector<int> &result) const;
  void queryAfter(const AABB &box, int minPrimitiveId,
    std::vector<int> &visitedStamp, int stamp,
    std::vector<int> &result) const;
  std::uint64_t queryOverlapping(const AABB &box, const std::vector<AABB> &candidateBoxes, int selfPrimitiveId,
    std::vector<int> &visitedStamp, int stamp,
    std::vector<int> &result) const;
  std::uint64_t queryOverlappingAfter(const AABB &box, const std::vector<AABB> &candidateBoxes, int minPrimitiveId,
    std::vector<int> &visitedStamp, int stamp,
    std::vector<int> &result) const;

private:
  static std::uint64_t hashCoord(int ix, int iy, int iz);
  void toGrid(const EigenSupport::V3d &p, int &ix, int &iy, int &iz) const;
  void insertGridRange(int loIx, int loIy, int loIz,
    int hiIx, int hiIy, int hiIz, int primitiveId);
  void queryFiltered(const AABB &box, int selfPrimitiveId, int minPrimitiveId,
    std::vector<int> &visitedStamp, int stamp,
    std::vector<int> &result) const;
  std::uint64_t queryFilteredOverlapping(
    const AABB &box, const std::vector<AABB> &candidateBoxes, int selfPrimitiveId, int minPrimitiveId,
    std::vector<int> &visitedStamp, int stamp,
    std::vector<int> &result) const;

  double cellSize_ = 1.0;
  std::unordered_map<std::uint64_t, std::vector<int>> cells_;
};

}  // namespace CIPC
}  // namespace Contact
}  // namespace pgo
