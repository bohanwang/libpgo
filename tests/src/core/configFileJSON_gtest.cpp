#include <gtest/gtest.h>

#include "configFileJSON.h"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <string>

namespace
{
namespace fs = std::filesystem;

constexpr const char *kBoxJsonPath = LIBPGO_TEST_BOX_JSON;
constexpr const char *kCubicBoxJsonPath = LIBPGO_TEST_CUBIC_BOX_JSON;
constexpr const char *kShellJsonPath = LIBPGO_TEST_SHELL_JSON;

class ScopedTempDir
{
public:
  ScopedTempDir()
  {
    const auto uniqueSuffix = std::to_string(std::chrono::steady_clock::now().time_since_epoch().count());
    dir = fs::temp_directory_path() / ("libpgo-config-json-gtest-" + uniqueSuffix);
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

TEST(ConfigFileJSONGTest, ResolvesExamplePathsAgainstConfigDirectory)
{
  pgo::ConfigFileJSON config;
  ASSERT_TRUE(config.open(kBoxJsonPath));

  const fs::path configPath = fs::path(kBoxJsonPath).lexically_normal();
  const fs::path configDir = configPath.parent_path();
  EXPECT_EQ(fs::path(config.getConfigFilename()), configPath);
  EXPECT_EQ(fs::path(config.getConfigDirectory()), configDir);
  EXPECT_EQ(fs::path(config.getResolvedPath("tet-mesh", 1)), (configDir / "box.veg").lexically_normal());
  EXPECT_EQ(fs::path(config.getResolvedPath("surface-mesh", 1)), (configDir / "box.obj").lexically_normal());
  EXPECT_EQ(fs::path(config.getResolvedPath("output", 1)), (configDir / "ret-box").lexically_normal());
}

TEST(ConfigFileJSONGTest, ResolvesCubicExampleSurfaceToOriginalObj)
{
  pgo::ConfigFileJSON config;
  ASSERT_TRUE(config.open(kCubicBoxJsonPath));

  const fs::path configDir = fs::path(kCubicBoxJsonPath).lexically_normal().parent_path();
  EXPECT_EQ(fs::path(config.getResolvedPath("cubic-mesh", 1)), (configDir / "box.veg").lexically_normal());
  EXPECT_EQ(fs::path(config.getResolvedPath("surface-mesh", 1)), (configDir / "box.obj").lexically_normal());
  EXPECT_EQ(fs::path(config.getResolvedPath("output", 1)), (configDir / "ret-cubic-box").lexically_normal());
}

TEST(ConfigFileJSONGTest, ResolvesShellNestedPathsAgainstConfigDirectory)
{
  pgo::ConfigFileJSON config;
  ASSERT_TRUE(config.open(kShellJsonPath));

  const fs::path configDir = fs::path(kShellJsonPath).lexically_normal().parent_path();
  ASSERT_TRUE(config.exist("fixed-vertices"));
  ASSERT_TRUE(config.exist("external-objects"));
  ASSERT_FALSE(config.handle()["fixed-vertices"].empty());
  ASSERT_FALSE(config.handle()["external-objects"].empty());

  const std::string fixedFilename = config.handle()["fixed-vertices"][0]["filename"].get<std::string>();
  const std::string externalFilename = config.handle()["external-objects"][0]["filename"].get<std::string>();
  EXPECT_EQ(fs::path(config.resolvePath(fixedFilename)), (configDir / fixedFilename).lexically_normal());
  EXPECT_EQ(fs::path(config.resolvePath(externalFilename)), (configDir / externalFilename).lexically_normal());
}

TEST(ConfigFileJSONGTest, ResolvesRawRelativeAndAbsolutePaths)
{
  ScopedTempDir tempDir;
  const fs::path configPath = tempDir.path() / "configs" / "scene.json";
  writeTextFile(configPath, R"({"paths":["a.txt","sub/../b.txt"]})");

  pgo::ConfigFileJSON config;
  ASSERT_TRUE(config.open(configPath.string().c_str()));

  const fs::path configDir = configPath.parent_path();
  const fs::path absolutePath = tempDir.path() / "absolute" / "mesh.obj";
  EXPECT_EQ(fs::path(config.resolvePath("./foo/../bar.obj")), (configDir / "bar.obj").lexically_normal());
  EXPECT_EQ(fs::path(config.resolvePath("../shared/mesh.obj")), (configDir / "../shared/mesh.obj").lexically_normal());
  EXPECT_EQ(fs::path(config.resolvePath(absolutePath.string())), absolutePath.lexically_normal());
}

TEST(ConfigFileJSONGTest, ResolvesVectorPathEntriesAgainstConfigDirectory)
{
  ScopedTempDir tempDir;
  const fs::path configPath = tempDir.path() / "nested" / "scene.json";
  writeTextFile(configPath, R"({"paths":["a.txt","sub/../b.txt"]})");

  pgo::ConfigFileJSON config;
  ASSERT_TRUE(config.open(configPath.string().c_str()));

  const auto paths = config.getVectorPath("paths", 1);
  ASSERT_EQ(paths.size(), 2u);
  EXPECT_EQ(fs::path(paths[0]), (configPath.parent_path() / "a.txt").lexically_normal());
  EXPECT_EQ(fs::path(paths[1]), (configPath.parent_path() / "b.txt").lexically_normal());
}

TEST(ConfigFileJSONGTest, PreservesLegacyTokenHelpers)
{
  pgo::ConfigFileJSON config;
  config.handle()["work-path"] = "{work}/out";
  config.handle()["asset-path"] = "{asset}/mesh.obj";

  EXPECT_EQ(config.getPathFromKey("work-path", "/tmp/work"), "/tmp/work/out");
  EXPECT_EQ(config.getPathFromKey("asset-path", "asset", "/tmp/assets"), "/tmp/assets/mesh.obj");
  EXPECT_EQ(pgo::ConfigFileJSON::getPathFromInput("{asset}/mesh.obj", "asset", "/tmp/assets"), "/tmp/assets/mesh.obj");
}
