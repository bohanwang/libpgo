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

TEST(SpatialHashGridGTest, QueryDoesNotMergeLegacySpatialHashCollisionCells)
{
  SpatialHashGrid grid(2);
  grid.setCellSize(1.0);
  grid.insert(makeBox(ES::V3d(0.1, 0.1, 0.1), ES::V3d(0.2, 0.2, 0.2)), 0);
  grid.insert(makeBox(
                ES::V3d(19349663.1, 73856093.1, 0.1),
                ES::V3d(19349663.2, 73856093.2, 0.2)),
    1);

  const std::vector<int> result = queryBox(
    grid,
    makeBox(ES::V3d(0.0, 0.0, 0.0), ES::V3d(0.3, 0.3, 0.3)),
    -1,
    2);

  ASSERT_EQ(result.size(), 1u);
  EXPECT_EQ(result[0], 0);
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

TEST(SpatialHashGridGTest, BuildFromAABBsMatchesRepeatedInsert)
{
  std::vector<SpatialHashGrid::AABB> boxes;
  boxes.push_back(makeBox(ES::V3d(0.1, 0.1, 0.1), ES::V3d(0.2, 0.2, 0.2)));
  boxes.push_back(makeBox(ES::V3d(0.8, 0.8, 0.8), ES::V3d(1.2, 1.2, 1.2)));
  boxes.push_back(makeBox(ES::V3d(2.1, 0.1, 0.1), ES::V3d(2.2, 0.2, 0.2)));

  SpatialHashGrid repeated(3);
  repeated.setCellSize(1.0);
  for (int i = 0; i < static_cast<int>(boxes.size()); ++i)
    repeated.insert(boxes[i], i);

  SpatialHashGrid bulk(3);
  bulk.setCellSize(1.0);
  bulk.build(boxes);

  const SpatialHashGrid::AABB query = makeBox(
    ES::V3d(0.0, 0.0, 0.0),
    ES::V3d(1.5, 1.5, 1.5));

  std::vector<int> repeatedVisited(boxes.size(), 0);
  std::vector<int> bulkVisited(boxes.size(), 0);
  std::vector<int> repeatedResult;
  std::vector<int> bulkResult;
  repeated.query(query, -1, repeatedVisited, 1, repeatedResult);
  bulk.query(query, -1, bulkVisited, 1, bulkResult);

  std::sort(repeatedResult.begin(), repeatedResult.end());
  std::sort(bulkResult.begin(), bulkResult.end());

  EXPECT_EQ(bulkResult, repeatedResult);
}

TEST(SpatialHashGridGTest, QueryAfterOnlyReturnsLargerPrimitiveIds)
{
  SpatialHashGrid grid(3);
  grid.setCellSize(1.0);
  grid.insert(makeBox(ES::V3d(0.1, 0.1, 0.1), ES::V3d(0.2, 0.2, 0.2)), 0);
  grid.insert(makeBox(ES::V3d(0.1, 0.1, 0.1), ES::V3d(0.2, 0.2, 0.2)), 1);
  grid.insert(makeBox(ES::V3d(0.1, 0.1, 0.1), ES::V3d(0.2, 0.2, 0.2)), 2);

  std::vector<int> visitedStamp(3, 0);
  std::vector<int> result;
  grid.queryAfter(
    makeBox(ES::V3d(0.0, 0.0, 0.0), ES::V3d(0.3, 0.3, 0.3)),
    1,
    visitedStamp,
    1,
    result);

  ASSERT_EQ(result.size(), 1u);
  EXPECT_EQ(result[0], 2);
}

TEST(SpatialHashGridGTest, QueryOverlappingFiltersSameCellNonOverlappingAABBs)
{
  std::vector<SpatialHashGrid::AABB> boxes;
  boxes.push_back(makeBox(ES::V3d(0.1, 0.1, 0.1), ES::V3d(0.2, 0.2, 0.2)));
  boxes.push_back(makeBox(ES::V3d(8.0, 8.0, 8.0), ES::V3d(9.0, 9.0, 9.0)));

  SpatialHashGrid grid(2);
  grid.setCellSize(10.0);
  grid.build(boxes);

  std::vector<int> visitedStamp(boxes.size(), 0);
  std::vector<int> result;
  const std::uint64_t hashCandidates = grid.queryOverlapping(
    makeBox(ES::V3d(0.0, 0.0, 0.0), ES::V3d(0.3, 0.3, 0.3)),
    boxes,
    -1,
    visitedStamp,
    1,
    result);

  ASSERT_EQ(result.size(), 1u);
  EXPECT_EQ(result[0], 0);
  EXPECT_EQ(hashCandidates, 2u);
}

TEST(SpatialHashGridGTest, QueryOverlappingAfterCountsOnlyEligiblePrimitiveIds)
{
  std::vector<SpatialHashGrid::AABB> boxes;
  boxes.push_back(makeBox(ES::V3d(0.1, 0.1, 0.1), ES::V3d(0.2, 0.2, 0.2)));
  boxes.push_back(makeBox(ES::V3d(0.1, 0.1, 0.1), ES::V3d(0.2, 0.2, 0.2)));
  boxes.push_back(makeBox(ES::V3d(8.0, 8.0, 8.0), ES::V3d(9.0, 9.0, 9.0)));

  SpatialHashGrid grid(3);
  grid.setCellSize(10.0);
  grid.build(boxes);

  std::vector<int> visitedStamp(boxes.size(), 0);
  std::vector<int> result;
  const std::uint64_t hashCandidates = grid.queryOverlappingAfter(
    makeBox(ES::V3d(0.0, 0.0, 0.0), ES::V3d(0.3, 0.3, 0.3)),
    boxes,
    1,
    visitedStamp,
    1,
    result);

  EXPECT_TRUE(result.empty());
  EXPECT_EQ(hashCandidates, 1u);
}
