#include <gtest/gtest.h>

#include "loadMshFile.h"
#include "tetMesh.h"

#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <string>

namespace
{
namespace fs = std::filesystem;

// MSH 4.1 ASCII: one tetrahedron whose node tags are 10, 20, 30, 40 rather
// than 1..4. Treating tags as indices would produce out-of-range connectivity.
const char *kSparseTagMesh =
  "$MeshFormat\n4.1 0 8\n$EndMeshFormat\n"
  "$Entities\n0 0 0 1\n1 0 0 0 1 1 1 0 0\n$EndEntities\n"
  "$Nodes\n1 4 10 40\n3 1 0 4\n10\n20\n30\n40\n0 0 0\n1 0 0\n0 1 0\n0 0 1\n$EndNodes\n"
  "$Elements\n1 1 1 1\n3 1 4 1\n1 10 20 30 40\n$EndElements\n";
}  // namespace

TEST(LoadMshFileGTest, MapsSparseNodeTagsToArrayIndices)
{
  const fs::path path = fs::temp_directory_path() /
    ("libpgo-loadMshFile-gtest-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) + ".msh");
  {
    std::ofstream out(path);
    ASSERT_TRUE(out.is_open());
    out << kSparseTagMesh;
  }

  const pgo::VolumetricMeshes::TetMesh mesh = pgo::VolumetricMeshes::loadMshFile(path.string().c_str());
  ASSERT_EQ(mesh.getNumVertices(), 4);
  ASSERT_EQ(mesh.getNumElements(), 1);
  const int *v = mesh.getVertexIndices(0);
  for (int j = 0; j < 4; ++j) {
    EXPECT_GE(v[j], 0);
    EXPECT_LT(v[j], 4);
  }
  const pgo::Vec3d a = mesh.getVertex(v[0]), b = mesh.getVertex(v[1]), c = mesh.getVertex(v[2]), d = mesh.getVertex(v[3]);
  const pgo::Vec3d ab = b - a, ac = c - a, ad = d - a;
  const double volume = (ab[0] * (ac[1] * ad[2] - ac[2] * ad[1]) - ab[1] * (ac[0] * ad[2] - ac[2] * ad[0]) +
                          ab[2] * (ac[0] * ad[1] - ac[1] * ad[0])) / 6.0;
  EXPECT_NEAR(volume, 1.0 / 6.0, 1e-12);

  std::error_code ec;
  fs::remove(path, ec);
}
