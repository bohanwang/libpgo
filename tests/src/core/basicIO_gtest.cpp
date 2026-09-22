#include <gtest/gtest.h>

#include "basicIO.h"

#include <chrono>
#include <filesystem>
#include <iterator>
#include <string>
#include <vector>

namespace
{
namespace fs = std::filesystem;

fs::path uniqueTempFile(const char *stem)
{
  return fs::temp_directory_path() /
    (std::string("libpgo-basicIO-gtest-") + stem + "-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
}
}  // namespace

TEST(BasicIOGTest, BinaryScalarRoundTripReadsThroughTheScalarAddress)
{
  const std::vector<double> doubles = { 1.5, -2.25, 3e-7, 4e12 };
  const fs::path doublePath = uniqueTempFile("double.bin");
  ASSERT_EQ(pgo::BasicIO::write1DBinary(doublePath.string().c_str(), doubles.begin(), doubles.end()), 0);
  std::vector<double> readDoubles;
  ASSERT_EQ(pgo::BasicIO::read1DBinary(doublePath.string().c_str(), std::back_inserter(readDoubles)), 0);
  EXPECT_EQ(readDoubles, doubles);

  const std::vector<int> ints = { 7, -3, 0, 123456 };
  const fs::path intPath = uniqueTempFile("int.bin");
  ASSERT_EQ(pgo::BasicIO::write1DBinary(intPath.string().c_str(), ints.begin(), ints.end()), 0);
  std::vector<int> readInts;
  ASSERT_EQ(pgo::BasicIO::read1DBinary(intPath.string().c_str(), std::back_inserter(readInts), 3), 0);
  EXPECT_EQ(readInts, std::vector<int>(ints.begin(), ints.begin() + 3));

  std::error_code ec;
  fs::remove(doublePath, ec);
  fs::remove(intPath, ec);
}
