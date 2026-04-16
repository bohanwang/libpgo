#include <gtest/gtest.h>

#include "animationLoader.h"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <string>

namespace
{
namespace fs = std::filesystem;

class ScopedTempDir
{
public:
  ScopedTempDir()
  {
    const auto uniqueSuffix = std::to_string(std::chrono::steady_clock::now().time_since_epoch().count());
    dir = fs::temp_directory_path() / ("libpgo-animation-loader-gtest-" + uniqueSuffix);
    fs::create_directories(dir);
  }

  ~ScopedTempDir()
  {
    std::error_code ec;
    fs::remove_all(dir, ec);
  }

  const fs::path &path() const { return dir; }

private:
  fs::path dir;
};

void writeTextFile(const fs::path &path, const std::string &contents)
{
  fs::create_directories(path.parent_path());
  std::ofstream out(path);
  ASSERT_TRUE(out.is_open());
  out << contents;
}
}  // namespace

TEST(AnimationLoaderGTest, LoadsConfigRelativeDrivingMeshAndSequence)
{
  ScopedTempDir tempDir;
  const fs::path configPath = tempDir.path() / "anim" / "anim.json";
  const fs::path meshPath = configPath.parent_path() / "mesh.obj";
  const fs::path frame0Path = configPath.parent_path() / "frames" / "frame0000.obj";
  const fs::path frame1Path = configPath.parent_path() / "frames" / "frame0001.obj";

  writeTextFile(meshPath, "v 0 0 0\nv 1 0 0\nv 0 1 0\nf 1 2 3\n");
  writeTextFile(frame0Path, "v 0 0 0\nv 1 0 0\nv 0 1 0\nf 1 2 3\n");
  writeTextFile(frame1Path, "v 0 0 0\nv 1.2 0 0\nv 0 1 0\nf 1 2 3\n");
  writeTextFile(configPath, R"({
  "save-cache": 0,
  "meshes": [
    {
      "name": "tri",
      "driving-mesh": "mesh.obj",
      "sequence": "frames/frame{:04d}.obj",
      "sequence-type": "objmesh",
      "sequence-range": [0, 2]
    }
  ]
})");

  pgo::AnimationIO::AnimationLoader loader;
  EXPECT_EQ(loader.load(configPath.string().c_str()), 0);
}
