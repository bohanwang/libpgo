#include <gtest/gtest.h>

#include "barycentricCoordinates.h"

#include <chrono>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
namespace ES = pgo::EigenSupport;
namespace fs = std::filesystem;
using pgo::InterpolationCoordinates::BarycentricCoordinates;

const std::vector<int> kIndices = { 0, 1, 2, 3, 1, 2, 3, 4 };
const std::vector<double> kWeights = { 0.25, 0.25, 0.25, 0.25, 0.1, 0.2, 0.3, 0.4 };
}  // namespace

TEST(BarycentricCoordinatesGTest, WeightOnlyConstructorsKnowTheirCageSize)
{
  // inferred: the largest referenced vertex is 4, so the cage has 5 vertices
  BarycentricCoordinates inferred(2, 4, kIndices.data(), kWeights.data());
  EXPECT_EQ(inferred.getNumCageVertices(), 5);
  const ES::SpMatD inferredMatrix = inferred.generateInterpolationMatrix();
  EXPECT_EQ(inferredMatrix.rows(), 6);
  EXPECT_EQ(inferredMatrix.cols(), 15);
  EXPECT_DOUBLE_EQ(inferredMatrix.coeff(3, 4 * 3), 0.4);

  // explicit: unreferenced trailing cage vertices still get columns
  BarycentricCoordinates explicitCage(2, 4, kIndices.data(), kWeights.data(), nullptr, 8);
  EXPECT_EQ(explicitCage.getNumCageVertices(), 8);
  EXPECT_EQ(explicitCage.generateInterpolationMatrix().cols(), 24);

  const std::vector<ES::V4i> tets = { ES::V4i(0, 1, 2, 3), ES::V4i(1, 2, 3, 4) };
  const std::vector<ES::V4d> weights = { ES::V4d(0.25, 0.25, 0.25, 0.25), ES::V4d(0.1, 0.2, 0.3, 0.4) };
  BarycentricCoordinates fromTets(tets, weights, nullptr, 6);
  EXPECT_EQ(fromTets.generateInterpolationMatrix().cols(), 18);

  EXPECT_THROW(BarycentricCoordinates(2, 4, kIndices.data(), kWeights.data(), nullptr, 4), std::invalid_argument);
  EXPECT_THROW(BarycentricCoordinates().generateInterpolationMatrix(), std::invalid_argument);
}

TEST(BarycentricCoordinatesGTest, SavedWeightsCanBeLoadedAndTurnedIntoAMatrix)
{
  const fs::path path = fs::temp_directory_path() /
    ("libpgo-barycentric-gtest-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) + ".interp");
  BarycentricCoordinates original(2, 4, kIndices.data(), kWeights.data(), nullptr, 8);
  ASSERT_EQ(original.saveInterpolationWeights(path.string()), 0);

  BarycentricCoordinates loaded(path.string(), 8);
  EXPECT_EQ(loaded.getNumLocations(), 2);
  EXPECT_EQ(loaded.getNumCageVertices(), 8);
  const ES::MXd expected = ES::MXd(original.generateInterpolationMatrix());
  const ES::MXd actual = ES::MXd(loaded.generateInterpolationMatrix());
  EXPECT_EQ(actual.rows(), expected.rows());
  EXPECT_EQ(actual.cols(), expected.cols());
  EXPECT_LT((actual - expected).norm(), 1e-12);

  BarycentricCoordinates loadedInferred(path.string());
  EXPECT_EQ(loadedInferred.getNumCageVertices(), 5);
  EXPECT_EQ(loadedInferred.generateInterpolationMatrix().cols(), 15);

  std::error_code ec;
  fs::remove(path, ec);
}
