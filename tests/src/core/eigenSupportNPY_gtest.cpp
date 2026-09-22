#include <gtest/gtest.h>

#include "EigenSupport.h"

#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

namespace
{
namespace ES = pgo::EigenSupport;
namespace fs = std::filesystem;

fs::path uniqueTempFile(const char *stem)
{
  return fs::temp_directory_path() /
    (std::string("libpgo-npy-gtest-") + stem + "-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) + ".npy");
}

// Write a version 1.0 .npy file by hand so the test controls the header text.
void writeNpy(const fs::path &path, bool fortranOrder, int rows, int cols, const std::vector<double> &data)
{
  std::string header = std::string("{'descr': '<f8', 'fortran_order': ") + (fortranOrder ? "True" : "False") +
    ", 'shape': (" + std::to_string(rows) + ", " + std::to_string(cols) + "), }";
  const size_t pad = 64 - ((10 + header.size()) % 64);
  header += std::string(pad, ' ');
  header.back() = '\n';

  std::ofstream out(path, std::ios::binary);
  ASSERT_TRUE(out.is_open());
  out.write("\x93NUMPY", 6);
  const uint8_t version[2] = { 1, 0 };
  out.write(reinterpret_cast<const char *>(version), 2);
  const uint16_t headerLen = static_cast<uint16_t>(header.size());
  out.write(reinterpret_cast<const char *>(&headerLen), 2);
  out.write(header.data(), header.size());
  out.write(reinterpret_cast<const char *>(data.data()), data.size() * sizeof(double));
}
}  // namespace

TEST(EigenSupportNPYGTest, FortranOrderFlagIsParsedAndHonored)
{
  ES::MXd expected(2, 3);
  expected << 1, 2, 3, 4, 5, 6;

  // numpy stores a Fortran-order array column by column
  const fs::path fortranPath = uniqueTempFile("fortran");
  writeNpy(fortranPath, true, 2, 3, { 1, 4, 2, 5, 3, 6 });
  ES::MXd fromFortran;
  ASSERT_EQ(ES::readNPY(fortranPath.string().c_str(), fromFortran), 0);
  EXPECT_EQ((fromFortran - expected).norm(), 0.0);

  const fs::path cPath = uniqueTempFile("c");
  writeNpy(cPath, false, 2, 3, { 1, 2, 3, 4, 5, 6 });
  ES::MXd fromC;
  ASSERT_EQ(ES::readNPY(cPath.string().c_str(), fromC), 0);
  EXPECT_EQ((fromC - expected).norm(), 0.0);

  // the library writer round-trips through the reader
  const fs::path roundTripPath = uniqueTempFile("roundtrip");
  ASSERT_EQ(ES::writeNPY(roundTripPath.string().c_str(), expected, 0), 0);
  ES::MXd roundTrip;
  ASSERT_EQ(ES::readNPY(roundTripPath.string().c_str(), roundTrip), 0);
  EXPECT_EQ((roundTrip - expected).norm(), 0.0);

  std::error_code ec;
  fs::remove(fortranPath, ec);
  fs::remove(cPath, ec);
  fs::remove(roundTripPath, ec);
}
