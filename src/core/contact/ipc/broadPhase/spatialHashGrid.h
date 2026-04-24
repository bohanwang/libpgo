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

  void insert(const AABB &box, int primitiveId);
  void query(const AABB &box, int selfPrimitiveId,
    std::vector<int> &visitedStamp, int stamp,
    std::vector<int> &result) const;

private:
  static std::int64_t hashCoord(int ix, int iy, int iz);
  void toGrid(const EigenSupport::V3d &p, int &ix, int &iy, int &iz) const;

  double cellSize_ = 1.0;
  std::unordered_map<std::int64_t, std::vector<int>> cells_;
};

}  // namespace CIPC
}  // namespace Contact
}  // namespace pgo
