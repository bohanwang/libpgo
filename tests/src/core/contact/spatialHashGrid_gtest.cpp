#include <gtest/gtest.h>

#include "ipc/broadPhase/spatialHashGrid.h"

namespace
{
namespace ES = pgo::EigenSupport;
using pgo::Contact::CIPC::SpatialHashGrid;

SpatialHashGrid::AABB makeBox(const ES::V3d &lo, const ES::V3d &hi)
{
  SpatialHashGrid::AABB box;
  box.lo = lo;
  box.hi = hi;
  return box;
}

std::vector<int> queryBox(const SpatialHashGrid &grid, const SpatialHashGrid::AABB &box, int selfId, int candidateCount)
{
  std::vector<int> visitedStamp(candidateCount, 0);
  std::vector<int> result;
  grid.query(box, selfId, visitedStamp, 1, result);
  return result;
}
}  // namespace

TEST(SpatialHashGridGTest, QueryFindsInsertedAABB)
{
  SpatialHashGrid grid(1);
  grid.setCellSize(1.0);
  grid.insert(makeBox(ES::V3d(0.1, 0.1, 0.1), ES::V3d(0.2, 0.2, 0.2)), 0);

  const std::vector<int> result = queryBox(
    grid,
    makeBox(ES::V3d(0.0, 0.0, 0.0), ES::V3d(0.3, 0.3, 0.3)),
    -1,
    1);

  ASSERT_EQ(result.size(), 1u);
  EXPECT_EQ(result[0], 0);
}

TEST(SpatialHashGridGTest, QueryDeduplicatesPrimitiveAcrossMultipleCells)
{
  SpatialHashGrid grid(1);
  grid.setCellSize(1.0);
  grid.insert(makeBox(ES::V3d(0.25, 0.25, 0.25), ES::V3d(1.25, 1.25, 1.25)), 0);

  const std::vector<int> result = queryBox(
    grid,
    makeBox(ES::V3d(0.0, 0.0, 0.0), ES::V3d(1.5, 1.5, 1.5)),
    -1,
    1);

  ASSERT_EQ(result.size(), 1u);
  EXPECT_EQ(result[0], 0);
}

TEST(SpatialHashGridGTest, QueryExcludesSelfPrimitiveId)
{
  SpatialHashGrid grid(1);
  grid.setCellSize(1.0);
  grid.insert(makeBox(ES::V3d(0.1, 0.1, 0.1), ES::V3d(0.2, 0.2, 0.2)), 7);

  const std::vector<int> result = queryBox(
    grid,
    makeBox(ES::V3d(0.0, 0.0, 0.0), ES::V3d(0.3, 0.3, 0.3)),
    7,
    8);

  EXPECT_TRUE(result.empty());
}

TEST(SpatialHashGridGTest, QueryEmptyForUnoccupiedCells)
{
  SpatialHashGrid grid(1);
  grid.setCellSize(1.0);
  grid.insert(makeBox(ES::V3d(0.1, 0.1, 0.1), ES::V3d(0.2, 0.2, 0.2)), 0);

  const std::vector<int> result = queryBox(
    grid,
    makeBox(ES::V3d(5.0, 5.0, 5.0), ES::V3d(5.3, 5.3, 5.3)),
    -1,
    1);

  EXPECT_TRUE(result.empty());
}

TEST(SpatialHashGridGTest, ClearRemovesInsertedPrimitives)
{
  SpatialHashGrid grid(1);
  grid.setCellSize(1.0);
  grid.insert(makeBox(ES::V3d(0.1, 0.1, 0.1), ES::V3d(0.2, 0.2, 0.2)), 0);
  grid.clear();

  const std::vector<int> result = queryBox(
    grid,
    makeBox(ES::V3d(0.0, 0.0, 0.0), ES::V3d(0.3, 0.3, 0.3)),
    -1,
    1);

  EXPECT_TRUE(result.empty());
}
