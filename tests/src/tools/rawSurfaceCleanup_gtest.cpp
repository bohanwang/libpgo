#include <gtest/gtest.h>
#include <nlohmann/json.hpp>

#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#ifndef PGO_TEST_RAW_SURFACE_CLEANUP_BIN
#define PGO_TEST_RAW_SURFACE_CLEANUP_BIN ""
#endif

namespace
{

struct Vec3
{
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;
};

struct Tri
{
  int a = 0;
  int b = 0;
  int c = 0;
};

std::string quotePath(const std::filesystem::path &path)
{
  return "\"" + path.string() + "\"";
}

std::string shellExecutable(const std::filesystem::path &path)
{
#ifdef _WIN32
  return "call " + quotePath(path);
#else
  return quotePath(path);
#endif
}

std::filesystem::path getRawSurfaceCleanupBinaryPath()
{
  const char *envBin = std::getenv("RAW_SURFACE_CLEANUP_BIN");
  if (envBin != nullptr && *envBin != '\0')
    return std::filesystem::path(envBin);

  if (std::string(PGO_TEST_RAW_SURFACE_CLEANUP_BIN).empty())
    return {};

  return std::filesystem::path(PGO_TEST_RAW_SURFACE_CLEANUP_BIN);
}

int runCommand(const std::string &command)
{
  return std::system(command.c_str());
}

class ScopedTempDir
{
public:
  ScopedTempDir()
  {
    const auto base = std::filesystem::temp_directory_path();
    for (int attempt = 0; attempt < 32; ++attempt) {
      const auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
      path_ = base / ("libpgo-rawSurfaceCleanup-test-" + std::to_string(stamp) + "-" + std::to_string(std::rand()) + "-" + std::to_string(attempt));

      std::error_code ec;
      if (std::filesystem::create_directories(path_, ec))
        return;
    }

    throw std::runtime_error("Failed to create a unique temporary directory");
  }

  ~ScopedTempDir()
  {
    std::error_code ec;
    std::filesystem::remove_all(path_, ec);
  }

  const std::filesystem::path &path() const { return path_; }

private:
  std::filesystem::path path_;
};

void writeObj(const std::filesystem::path &path, const std::vector<Vec3> &vertices, const std::vector<Tri> &triangles)
{
  std::ofstream out(path);
  ASSERT_TRUE(out.is_open());

  for (const auto &v : vertices)
    out << "v " << v.x << " " << v.y << " " << v.z << "\n";

  for (const auto &tri : triangles)
    out << "f " << (tri.a + 1) << " " << (tri.b + 1) << " " << (tri.c + 1) << "\n";
}

void writeCubeWithDuplicateDegenerateFace(const std::filesystem::path &path)
{
  writeObj(path,
    {
      { 0.0, 0.0, 0.0 }, { 1.0, 0.0, 0.0 }, { 1.0, 1.0, 0.0 }, { 0.0, 1.0, 0.0 },
      { 0.0, 0.0, 1.0 }, { 1.0, 0.0, 1.0 }, { 1.0, 1.0, 1.0 }, { 0.0, 1.0, 1.0 }
    },
    {
      { 0, 2, 1 }, { 0, 3, 2 },
      { 4, 5, 6 }, { 4, 6, 7 },
      { 0, 1, 5 }, { 0, 5, 4 },
      { 3, 7, 6 }, { 3, 6, 2 },
      { 0, 4, 7 }, { 0, 7, 3 },
      { 1, 2, 6 }, { 1, 6, 5 },
      { 0, 1, 1 },
    });
}

void writeTwoComponentsWithDegenerateTriangle(const std::filesystem::path &path)
{
  writeObj(path,
    {
      { 0.0, 0.0, 0.0 }, { 1.0, 0.0, 0.0 }, { 1.0, 1.0, 0.0 }, { 0.0, 1.0, 0.0 },
      { 0.0, 0.0, 1.0 }, { 1.0, 0.0, 1.0 }, { 1.0, 1.0, 1.0 }, { 0.0, 1.0, 1.0 },
      { 3.0, 0.0, 0.0 }, { 3.5, 0.0, 0.0 }, { 3.5, 0.0, 0.0 },
    },
    {
      { 0, 2, 1 }, { 0, 3, 2 },
      { 4, 5, 6 }, { 4, 6, 7 },
      { 0, 1, 5 }, { 0, 5, 4 },
      { 3, 7, 6 }, { 3, 6, 2 },
      { 0, 4, 7 }, { 0, 7, 3 },
      { 1, 2, 6 }, { 1, 6, 5 },
      { 8, 9, 10 },
    });
}

nlohmann::json readJson(const std::filesystem::path &path)
{
  std::ifstream in(path);
  EXPECT_TRUE(in.is_open());
  nlohmann::json json;
  in >> json;
  return json;
}

std::filesystem::path requireBinary()
{
  const std::filesystem::path binary = getRawSurfaceCleanupBinaryPath();
  if (binary.empty())
    throw std::runtime_error("rawSurfaceCleanup binary path is not configured");
  return binary;
}

}  // namespace

TEST(rawSurfaceCleanupTool, DropsDegenerateFaceOnlyWhenTopologyGatePasses)
{
  ScopedTempDir tempDir;
  const std::filesystem::path inputObj = tempDir.path() / "cube_degenerate.obj";
  const std::filesystem::path outputObj = tempDir.path() / "cube_clean.obj";
  const std::filesystem::path reportPath = tempDir.path() / "cube_clean.json";
  writeCubeWithDuplicateDegenerateFace(inputObj);

  const std::string command = shellExecutable(requireBinary()) +
    " --input " + quotePath(inputObj) +
    " --output " + quotePath(outputObj) +
    " --json " + quotePath(reportPath) +
    " --expected-components 1";
  EXPECT_EQ(runCommand(command), 0);

  const nlohmann::json report = readJson(reportPath);
  EXPECT_TRUE(std::filesystem::exists(outputObj));
  EXPECT_EQ(report["invalid_triangles_before"], 1);
  EXPECT_EQ(report["invalid_triangles_after"], 0);
  EXPECT_EQ(report["accepted_deletions"], 1);
  EXPECT_EQ(report["components_after"], 1);
  EXPECT_EQ(report["boundary_or_nonmanifold_edges_after"], 0);
  EXPECT_TRUE(report["is_manifold_after"]);
  EXPECT_TRUE(report["topology_preserved"]);
}

TEST(rawSurfaceCleanupTool, PreservesExpectedComponentCountWhenCleanupWouldDropAComponent)
{
  ScopedTempDir tempDir;
  const std::filesystem::path inputObj = tempDir.path() / "two_components.obj";
  const std::filesystem::path outputObj = tempDir.path() / "two_components_clean.obj";
  const std::filesystem::path reportPath = tempDir.path() / "two_components_clean.json";
  writeTwoComponentsWithDegenerateTriangle(inputObj);

  const std::string command = shellExecutable(requireBinary()) +
    " --input " + quotePath(inputObj) +
    " --output " + quotePath(outputObj) +
    " --json " + quotePath(reportPath) +
    " --expected-components 2";
  EXPECT_EQ(runCommand(command), 0);

  const nlohmann::json report = readJson(reportPath);
  EXPECT_EQ(report["invalid_triangles_before"], 1);
  EXPECT_EQ(report["invalid_triangles_after"], 1);
  EXPECT_EQ(report["accepted_deletions"], 0);
  EXPECT_EQ(report["components_after"], 2);
  EXPECT_FALSE(report["topology_preserved"]);
}

TEST(rawSurfaceCleanupTool, DryRunReportsCleanupWithoutWritingMesh)
{
  ScopedTempDir tempDir;
  const std::filesystem::path inputObj = tempDir.path() / "cube_degenerate.obj";
  const std::filesystem::path outputObj = tempDir.path() / "cube_clean.obj";
  const std::filesystem::path reportPath = tempDir.path() / "cube_clean.json";
  writeCubeWithDuplicateDegenerateFace(inputObj);

  const std::string command = shellExecutable(requireBinary()) +
    " --input " + quotePath(inputObj) +
    " --output " + quotePath(outputObj) +
    " --json " + quotePath(reportPath) +
    " --expected-components 1 --dry-run";
  EXPECT_EQ(runCommand(command), 0);

  const nlohmann::json report = readJson(reportPath);
  EXPECT_TRUE(report["dry_run"]);
  EXPECT_FALSE(std::filesystem::exists(outputObj));
  EXPECT_EQ(report["invalid_triangles_before"], 1);
  EXPECT_EQ(report["invalid_triangles_after"], 0);
  EXPECT_EQ(report["accepted_deletions"], 1);
}
