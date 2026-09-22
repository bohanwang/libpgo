#include <gtest/gtest.h>

#include "animationLoader.h"
#include "EigenSupport.h"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <algorithm>
#include <string>
#include <vector>

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

namespace
{
namespace ES = pgo::EigenSupport;

const char *kTriangle = "v 0 0 0\nv 1 0 0\nv 0 1 0\nf 1 2 3\n";

void writeTriangleSequence(const fs::path &dir, int numFrames, const std::vector<int> &skip = {})
{
  writeTextFile(dir / "mesh.obj", kTriangle);
  for (int i = 0; i < numFrames; ++i) {
    if (std::find(skip.begin(), skip.end(), i) != skip.end())
      continue;
    const double dx = 0.1 * i;
    writeTextFile(dir / "frames" / ("frame000" + std::to_string(i) + ".obj"),
      "v 0 0 0\nv " + std::to_string(1.0 + dx) + " 0 0\nv 0 1 0\nf 1 2 3\n");
  }
}

std::string triangleConfig(int numFrames, const std::string &extra = "", const std::string &sequenceType = "objmesh",
  const std::string &sequence = "frames/frame{:04d}.obj")
{
  return "{\n  \"save-cache\": 0,\n  \"meshes\": [ { \"name\": \"tri\", \"driving-mesh\": \"mesh.obj\", \"sequence\": \"" + sequence +
    "\", \"sequence-type\": \"" + sequenceType + "\", \"sequence-range\": [0, " + std::to_string(numFrames) + "]" + extra + " } ]\n}\n";
}
}  // namespace

TEST(AnimationLoaderGTest, MissingOrEmptyFrameSequencesFailInsteadOfShiftingTime)
{
  ScopedTempDir tempDir;
  const fs::path dir = tempDir.path() / "anim";
  const fs::path configPath = dir / "anim.json";

  // frame 1 of 3 missing: previously frames 2.. were silently written one sample early
  writeTriangleSequence(dir, 3, { 1 });
  writeTextFile(configPath, triangleConfig(3));
  pgo::AnimationIO::AnimationLoader missing;
  EXPECT_NE(missing.load(configPath.string().c_str()), 0);

  // nothing loads at all: previously indexed disp[0] on an empty vector
  writeTextFile(configPath, triangleConfig(3, "", "objmesh", "frames/nope{:04d}.obj"));
  pgo::AnimationIO::AnimationLoader none;
  EXPECT_NE(none.load(configPath.string().c_str()), 0);

  // complete sequence loads, with one column per frame in order
  writeTriangleSequence(dir, 3);
  writeTextFile(configPath, triangleConfig(3));
  pgo::AnimationIO::AnimationLoader ok;
  ASSERT_EQ(ok.load(configPath.string().c_str()), 0);
  ASSERT_EQ(ok.sequences().size(), 1u);
  const ES::MXd &u = ok.sequences()[0].drivingDisplacements;
  ASSERT_EQ(u.rows(), 9);
  ASSERT_EQ(u.cols(), 3);
  EXPECT_NEAR(u(3, 1), 0.1, 1e-12);
  EXPECT_NEAR(u(3, 2), 0.2, 1e-12);
  EXPECT_EQ(ok.saveABC(dir.string().c_str()), 0);
  EXPECT_GT(fs::file_size(dir / "tri.abc"), 0u);
}

TEST(AnimationLoaderGTest, DisplayMeshIsLoadedAndEmbeddedIntoTheDrivingMesh)
{
  ScopedTempDir tempDir;
  const fs::path dir = tempDir.path() / "anim";
  const fs::path configPath = dir / "anim.json";
  writeTriangleSequence(dir, 2);
  // a smaller triangle whose first vertex is the midpoint of the driving edge 0-1
  writeTextFile(dir / "skin.obj", "v 0.5 0 0\nv 0.25 0.25 0\nv 0.1 0.1 0\nf 1 2 3\n");
  writeTextFile(configPath, triangleConfig(2, ", \"display-mesh\": \"skin.obj\""));

  pgo::AnimationIO::AnimationLoader loader;
  ASSERT_EQ(loader.load(configPath.string().c_str()), 0);
  ASSERT_EQ(loader.sequences().size(), 1u);
  const auto &seq = loader.sequences()[0];
  EXPECT_EQ(seq.displayMeshV.size(), 9);
  EXPECT_EQ(seq.displayMeshF.size(), 1u);
  ASSERT_EQ(seq.displayDisplacements.rows(), 9);
  ASSERT_EQ(seq.displayDisplacements.cols(), 2);
  // vertex 1 of the driving mesh moves by 0.1 in frame 1, so its edge midpoint moves by 0.05
  EXPECT_NEAR(seq.displayDisplacements(0, 1), 0.05, 1e-9);
  EXPECT_NEAR(seq.displayDisplacements(1, 1), 0.0, 1e-9);
  EXPECT_EQ(loader.saveABC(dir.string().c_str()), 0);
  EXPECT_GT(fs::file_size(dir / "tri.abc"), 0u);

  // a display mesh that cannot be loaded is an error, not an empty archive
  writeTextFile(configPath, triangleConfig(2, ", \"display-mesh\": \"missing.obj\""));
  pgo::AnimationIO::AnimationLoader missing;
  EXPECT_NE(missing.load(configPath.string().c_str()), 0);
}

TEST(AnimationLoaderGTest, DisplacementCacheIsWrittenReusedAndValidated)
{
  ScopedTempDir tempDir;
  const fs::path dir = tempDir.path() / "anim";
  const fs::path configPath = dir / "anim.json";
  const fs::path cachePath = dir / "tri-uAll.u";
  writeTriangleSequence(dir, 3);
  std::string config = triangleConfig(3);
  config.replace(config.find("\"save-cache\": 0"), 15, "\"save-cache\": 1");
  writeTextFile(configPath, config);

  pgo::AnimationIO::AnimationLoader first;
  ASSERT_EQ(first.load(configPath.string().c_str()), 0);
  ASSERT_TRUE(fs::exists(cachePath));
  const ES::MXd fromFrames = first.sequences()[0].drivingDisplacements;

  // remove the frames: the second load must come from the cache with identical data
  fs::remove_all(dir / "frames");
  pgo::AnimationIO::AnimationLoader second;
  ASSERT_EQ(second.load(configPath.string().c_str()), 0);
  const ES::MXd &fromCache = second.sequences()[0].drivingDisplacements;
  ASSERT_EQ(fromCache.rows(), fromFrames.rows());
  ASSERT_EQ(fromCache.cols(), fromFrames.cols());
  EXPECT_EQ((fromCache - fromFrames).norm(), 0.0);
  EXPECT_EQ(second.saveABC(dir.string().c_str()), 0);

  // a cache of the wrong size is rejected instead of producing an empty archive
  ES::MXd wrong = ES::MXd::Zero(9, 2);
  ASSERT_EQ(ES::writeMatrix(cachePath.string().c_str(), wrong), 0);
  pgo::AnimationIO::AnimationLoader stale;
  EXPECT_NE(stale.load(configPath.string().c_str()), 0);
}

TEST(AnimationLoaderGTest, DisplacementFileSequencesLoadAndUnknownTypesFail)
{
  ScopedTempDir tempDir;
  const fs::path dir = tempDir.path() / "anim";
  const fs::path configPath = dir / "anim.json";
  writeTextFile(dir / "mesh.obj", kTriangle);
  fs::create_directories(dir / "u");

  ES::MXd uAll = ES::MXd::Zero(9, 3);
  for (int i = 0; i < 3; ++i) {
    uAll(3, i) = 0.1 * i;
    ES::MXd frame = uAll.col(i);
    ASSERT_EQ(ES::writeMatrix((dir / "u" / ("u000" + std::to_string(i) + ".u")).string().c_str(), frame), 0);
  }
  ASSERT_EQ(ES::writeMatrix((dir / "uall.u").string().c_str(), uAll), 0);

  writeTextFile(configPath, triangleConfig(3, "", "u", "u/u{:04d}.u"));
  pgo::AnimationIO::AnimationLoader perFrame;
  ASSERT_EQ(perFrame.load(configPath.string().c_str()), 0);
  EXPECT_EQ((perFrame.sequences()[0].drivingDisplacements - uAll).norm(), 0.0);

  writeTextFile(configPath, triangleConfig(3, "", "uall", "uall.u"));
  pgo::AnimationIO::AnimationLoader all;
  ASSERT_EQ(all.load(configPath.string().c_str()), 0);
  EXPECT_EQ((all.sequences()[0].drivingDisplacements - uAll).norm(), 0.0);
  EXPECT_EQ(all.saveABC(dir.string().c_str()), 0);

  writeTextFile(configPath, triangleConfig(3, "", "objMesh"));
  pgo::AnimationIO::AnimationLoader unknown;
  EXPECT_NE(unknown.load(configPath.string().c_str()), 0);
}
