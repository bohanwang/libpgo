#include <gtest/gtest.h>
#include <nlohmann/json.hpp>

#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#ifndef PGO_TEST_REMESH_SURFACE_BIN
#define PGO_TEST_REMESH_SURFACE_BIN ""
#endif

#ifndef PGO_TEST_MESH_QUALITY_CHECK_BIN
#define PGO_TEST_MESH_QUALITY_CHECK_BIN ""
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

std::filesystem::path getBinaryPath(const char *envName, const char *compiledPath)
{
  const char *envBin = std::getenv(envName);
  if (envBin != nullptr && *envBin != '\0')
    return std::filesystem::path(envBin);

  if (std::string(compiledPath).empty())
    return {};

  return std::filesystem::path(compiledPath);
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
      path_ = base / ("libpgo-remeshSurface-test-" + std::to_string(stamp) + "-" + std::to_string(std::rand()) + "-" + std::to_string(attempt));

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

  for (const Vec3 &v : vertices)
    out << "v " << v.x << ' ' << v.y << ' ' << v.z << '\n';

  for (const Tri &t : triangles)
    out << "f " << (t.a + 1) << ' ' << (t.b + 1) << ' ' << (t.c + 1) << '\n';
}

void writeIntersectingTetrahedraObj(const std::filesystem::path &path)
{
  writeObj(path,
    {
      { 0.0, 0.0, 0.0 },
      { 1.0, 0.0, 0.0 },
      { 0.0, 1.0, 0.0 },
      { 0.0, 0.0, 1.0 },
      { 0.35, 0.35, -0.35 },
      { 1.35, 0.35, -0.35 },
      { 0.35, 1.35, -0.35 },
      { 0.35, 0.35, 0.65 },
    },
    {
      { 0, 2, 1 },
      { 0, 1, 3 },
      { 1, 2, 3 },
      { 2, 0, 3 },
      { 4, 6, 5 },
      { 4, 5, 7 },
      { 5, 6, 7 },
      { 6, 4, 7 },
    });
}

nlohmann::json readJson(const std::filesystem::path &path)
{
  std::ifstream in(path);
  if (in.is_open() == false)
    throw std::runtime_error("Failed to open JSON file: " + path.string());

  nlohmann::json j;
  in >> j;
  return j;
}

}  // namespace

TEST(RemeshSurface, CgalRepairSelfIntersectionsRemovesIntersectingFaces)
{
  const std::filesystem::path remeshSurface = getBinaryPath("REMESH_SURFACE_BIN", PGO_TEST_REMESH_SURFACE_BIN);
  const std::filesystem::path meshQualityCheck = getBinaryPath("MESH_QUALITY_CHECK_BIN", PGO_TEST_MESH_QUALITY_CHECK_BIN);
  ASSERT_FALSE(remeshSurface.empty());
  ASSERT_FALSE(meshQualityCheck.empty());

  ScopedTempDir tempDir;
  const std::filesystem::path inputObj = tempDir.path() / "intersecting.obj";
  const std::filesystem::path repairedObj = tempDir.path() / "repaired.obj";
  const std::filesystem::path reportPath = tempDir.path() / "quality.json";
  writeIntersectingTetrahedraObj(inputObj);

  const std::string repairCommand =
    shellExecutable(remeshSurface) + " cgal_repair_self_intersections --input-mesh " + quotePath(inputObj) +
    " --output-mesh " + quotePath(repairedObj);
  EXPECT_EQ(runCommand(repairCommand), 0);
  ASSERT_TRUE(std::filesystem::exists(repairedObj));

  const std::string qualityCommand =
    shellExecutable(meshQualityCheck) +
    " surface --check-level full --expected-components 1 --invalid-triangles-policy fail --self-intersection-backend exact-count"
    " --self-intersection-triangle-limit 10000 --input " +
    quotePath(repairedObj) + " --json " + quotePath(reportPath);
  EXPECT_EQ(runCommand(qualityCommand), 0);

  const nlohmann::json report = readJson(reportPath);
  EXPECT_EQ(report["self_intersections"], 0);
  EXPECT_TRUE(report["is_manifold"]);
  EXPECT_EQ(report["components_by_edge"], 1);
}
