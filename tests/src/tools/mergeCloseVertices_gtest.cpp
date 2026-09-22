#include <gtest/gtest.h>

#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

namespace
{
namespace fs = std::filesystem;

std::string quotePath(const fs::path &path)
{
  return "\"" + path.generic_string() + "\"";
}

int countObjLines(const fs::path &path, const std::string &prefix)
{
  std::ifstream in(path);
  int count = 0;
  std::string line;
  while (std::getline(in, line))
    if (line.rfind(prefix, 0) == 0)
      ++count;
  return count;
}
}  // namespace

// The input mixes a quad with triangles and has an unstitched seam whose
// duplicate vertices lie within the merge distance. Rebuilding faces used to
// write the quad's fourth corner past a three-entry array.
TEST(MergeCloseVerticesCli, PolygonFacesDoNotOverflowTheTriangleBuffer)
{
  const fs::path binary(PGO_TEST_MERGE_CLOSE_VERTICES_BIN);
  ASSERT_TRUE(fs::exists(binary));

  const fs::path dir = fs::temp_directory_path() /
    ("libpgo-mergeCloseVertices-gtest-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
  fs::create_directories(dir);
  const fs::path input = dir / "seam.obj";
  const fs::path output = dir / "merged.obj";
  {
    std::ofstream out(input);
    ASSERT_TRUE(out.is_open());
    // quad 1-2-3-4 on the left, triangles on the right sharing the seam x = 1
    // through duplicated vertices 5/6 that sit 1e-6 away from 2/3
    out << "v 0 0 0\nv 1 0 0\nv 1 1 0\nv 0 1 0\n"
        << "v 1.000001 0 0\nv 1.000001 1 0\nv 2 0 0\nv 2 1 0\n"
        << "f 1 2 3 4\n"
        << "f 5 7 8\nf 5 8 6\n";
  }

  const int status = std::system((quotePath(binary) + " " + quotePath(input) + " " + quotePath(output) + " 1e-3").c_str());
  EXPECT_EQ(status, 0);
  ASSERT_TRUE(fs::exists(output));
  // The two seam duplicates are merged away and the triangles survive. The
  // tool rebuilds only triangles, so the quad is dropped and its two private
  // corners are removed as isolated vertices.
  EXPECT_EQ(countObjLines(output, "v "), 4);
  EXPECT_EQ(countObjLines(output, "f "), 2);

  std::error_code ec;
  fs::remove_all(dir, ec);
}
