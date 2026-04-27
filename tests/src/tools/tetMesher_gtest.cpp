#include "tetMesh.h"
#include "triMeshGeo.h"

#include <gtest/gtest.h>

#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#ifndef PGO_TEST_TET_MESHER_BIN
#define PGO_TEST_TET_MESHER_BIN ""
#endif

#ifndef PGO_TEST_TET_MESHER_HAS_TET_WILD
#define PGO_TEST_TET_MESHER_HAS_TET_WILD 0
#endif

namespace
{

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

std::filesystem::path getTetMesherBinaryPath()
{
  const char *envBin = std::getenv("TET_MESHER_BIN");
  if (envBin != nullptr && *envBin != '\0')
    return std::filesystem::path(envBin);

  if (std::string(PGO_TEST_TET_MESHER_BIN).empty())
    return {};

  return std::filesystem::path(PGO_TEST_TET_MESHER_BIN);
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
      path_ = base / ("libpgo-tetMesher-test-" + std::to_string(stamp) + "-" + std::to_string(std::rand()) + "-" + std::to_string(attempt));

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

void writeUnitCubeObj(const std::filesystem::path &path)
{
  std::ofstream out(path);
  ASSERT_TRUE(out.is_open());

  out << "v 0 0 0\n";
  out << "v 1 0 0\n";
  out << "v 1 1 0\n";
  out << "v 0 1 0\n";
  out << "v 0 0 1\n";
  out << "v 1 0 1\n";
  out << "v 1 1 1\n";
  out << "v 0 1 1\n";
  out << "f 1 3 2\n";
  out << "f 1 4 3\n";
  out << "f 5 6 7\n";
  out << "f 5 7 8\n";
  out << "f 1 2 6\n";
  out << "f 1 6 5\n";
  out << "f 4 8 7\n";
  out << "f 4 7 3\n";
  out << "f 1 5 8\n";
  out << "f 1 8 4\n";
  out << "f 2 3 7\n";
  out << "f 2 7 6\n";
}

void writeHollowCubeObj(const std::filesystem::path &path)
{
  std::ofstream out(path);
  ASSERT_TRUE(out.is_open());

  out << "v 0 0 0\n";
  out << "v 1 0 0\n";
  out << "v 1 1 0\n";
  out << "v 0 1 0\n";
  out << "v 0 0 1\n";
  out << "v 1 0 1\n";
  out << "v 1 1 1\n";
  out << "v 0 1 1\n";
  out << "v 0.35 0.35 0.35\n";
  out << "v 0.65 0.35 0.35\n";
  out << "v 0.65 0.65 0.35\n";
  out << "v 0.35 0.65 0.35\n";
  out << "v 0.35 0.35 0.65\n";
  out << "v 0.65 0.35 0.65\n";
  out << "v 0.65 0.65 0.65\n";
  out << "v 0.35 0.65 0.65\n";

  const int outerFaces[12][3] = {
    { 1, 3, 2 }, { 1, 4, 3 }, { 5, 6, 7 }, { 5, 7, 8 },
    { 1, 2, 6 }, { 1, 6, 5 }, { 4, 8, 7 }, { 4, 7, 3 },
    { 1, 5, 8 }, { 1, 8, 4 }, { 2, 3, 7 }, { 2, 7, 6 },
  };

  for (const auto &tri : outerFaces)
    out << "f " << tri[0] << " " << tri[1] << " " << tri[2] << "\n";

  // Reverse the inner cube orientation so TetGen treats it as a cavity.
  for (const auto &tri : outerFaces)
    out << "f " << (tri[0] + 8) << " " << (tri[2] + 8) << " " << (tri[1] + 8) << "\n";
}

std::string readTextFile(const std::filesystem::path &path)
{
  std::ifstream in(path);
  if (!in.is_open())
    return {};

  return std::string(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
}

void writeTetgenConfig(const std::filesystem::path &path)
{
  std::ofstream out(path);
  ASSERT_TRUE(out.is_open());

  out << "{\n";
  out << "  \"version\": 1,\n";
  out << "  \"backend\": \"tetgen\",\n";
  out << "  \"input_mesh\": \"cube.obj\",\n";
  out << "  \"output_mesh\": \"cube.veg\",\n";
  out << "  \"output_surface\": \"cube_surface.obj\",\n";
  out << "  \"print_stats\": true,\n";
  out << "  \"tetgen\": {\n";
  out << "    \"command\": \"pq1.414a0.01\"\n";
  out << "  }\n";
  out << "}\n";
}

void writeTetwildConfig(const std::filesystem::path &path)
{
  std::ofstream out(path);
  ASSERT_TRUE(out.is_open());

  out << "{\n";
  out << "  \"version\": 1,\n";
  out << "  \"backend\": \"tetwild\",\n";
  out << "  \"input_mesh\": \"cube.obj\",\n";
  out << "  \"output_mesh\": \"cube_tetwild.veg\",\n";
  out << "  \"output_surface\": \"cube_tetwild_surface.obj\",\n";
  out << "  \"quiet\": true,\n";
  out << "  \"tetwild\": {\n";
  out << "    \"lr\": 0.2,\n";
  out << "    \"epsr\": 0.001,\n";
  out << "    \"stop_energy\": 10,\n";
  out << "    \"max_threads\": 1\n";
  out << "  }\n";
  out << "}\n";
}

void writeTetgenHollowConfig(const std::filesystem::path &path)
{
  std::ofstream out(path);
  ASSERT_TRUE(out.is_open());

  out << "{\n";
  out << "  \"version\": 1,\n";
  out << "  \"backend\": \"tetgen\",\n";
  out << "  \"input_mesh\": \"hollow_cube.obj\",\n";
  out << "  \"output_mesh\": \"hollow_cube.veg\",\n";
  out << "  \"output_surface\": \"hollow_cube_surface.obj\",\n";
  out << "  \"tetgen\": {\n";
  out << "    \"command\": \"pq1.414a0.005\"\n";
  out << "  }\n";
  out << "}\n";
}

int countSurfaceComponents(const pgo::Mesh::TriMeshGeo &surface)
{
  std::unordered_map<int, std::vector<int>> vertexFaces;
  for (int triID = 0; triID < surface.numTriangles(); ++triID) {
    const pgo::Vec3i &tri = surface.tri(triID);
    for (int i = 0; i < 3; ++i)
      vertexFaces[tri[i]].push_back(triID);
  }

  std::vector<char> visited(surface.numTriangles(), 0);
  int componentCount = 0;
  std::vector<int> stack;
  for (int triID = 0; triID < surface.numTriangles(); ++triID) {
    if (visited[triID])
      continue;

    componentCount++;
    visited[triID] = 1;
    stack.push_back(triID);
    while (stack.empty() == false) {
      const int currentTriID = stack.back();
      stack.pop_back();
      const pgo::Vec3i &tri = surface.tri(currentTriID);
      for (int i = 0; i < 3; ++i) {
        for (int adjacentTriID : vertexFaces[tri[i]]) {
          if (visited[adjacentTriID])
            continue;

          visited[adjacentTriID] = 1;
          stack.push_back(adjacentTriID);
        }
      }
    }
  }

  return componentCount;
}

}  // namespace

TEST(TetMesherCli, JsonTetgenGeneratesVegAndSurfaceAssets)
{
  const std::filesystem::path tetMesherBin = getTetMesherBinaryPath();
  ASSERT_FALSE(tetMesherBin.empty());
  ASSERT_TRUE(std::filesystem::exists(tetMesherBin));

  ScopedTempDir tempDir;
  const std::filesystem::path inputObj = tempDir.path() / "cube.obj";
  const std::filesystem::path outputVeg = tempDir.path() / "cube.veg";
  const std::filesystem::path outputSurface = tempDir.path() / "cube_surface.obj";
  const std::filesystem::path configPath = tempDir.path() / "tetgen.json";
  writeUnitCubeObj(inputObj);
  writeTetgenConfig(configPath);

  const std::string cmd = shellExecutable(tetMesherBin) +
    " --config " + quotePath(configPath);

  ASSERT_EQ(runCommand(cmd), 0);
  ASSERT_TRUE(std::filesystem::exists(outputVeg));
  ASSERT_TRUE(std::filesystem::exists(outputSurface));

  pgo::VolumetricMeshes::TetMesh tetMesh(outputVeg.string().c_str());
  EXPECT_GT(tetMesh.getNumVertices(), 0);
  EXPECT_GT(tetMesh.getNumElements(), 0);

  pgo::Mesh::TriMeshGeo surface;
  ASSERT_TRUE(surface.load(outputSurface.string()));
  EXPECT_GT(surface.numVertices(), 0);
  EXPECT_GT(surface.numTriangles(), 0);
}

TEST(TetMesherCli, TetwildDisabledReportsClearError)
{
  if constexpr (PGO_TEST_TET_MESHER_HAS_TET_WILD) {
    GTEST_SKIP() << "TetWild backend is enabled in this build.";
  }

  const std::filesystem::path tetMesherBin = getTetMesherBinaryPath();
  ASSERT_FALSE(tetMesherBin.empty());
  ASSERT_TRUE(std::filesystem::exists(tetMesherBin));

  ScopedTempDir tempDir;
  const std::filesystem::path inputObj = tempDir.path() / "cube.obj";
  const std::filesystem::path outputVeg = tempDir.path() / "cube_tetwild.veg";
  const std::filesystem::path configPath = tempDir.path() / "tetwild.json";
  const std::filesystem::path stderrPath = tempDir.path() / "tetwild.stderr";
  writeUnitCubeObj(inputObj);
  writeTetwildConfig(configPath);

  const std::string cmd = shellExecutable(tetMesherBin) +
    " --config " + quotePath(configPath) +
    " 2> " + quotePath(stderrPath);

  EXPECT_NE(runCommand(cmd), 0);
  EXPECT_FALSE(std::filesystem::exists(outputVeg));

  const std::string stderrText = readTextFile(stderrPath);
  EXPECT_NE(stderrText.find("tetwild backend is not enabled"), std::string::npos);
  EXPECT_NE(stderrText.find("PGO_TET_MESHER_USE_TET_WILD=ON"), std::string::npos);
}

TEST(TetMesherCli, JsonTetgenPreservesNegativeInnerBoundaryComponents)
{
  const std::filesystem::path tetMesherBin = getTetMesherBinaryPath();
  ASSERT_FALSE(tetMesherBin.empty());
  ASSERT_TRUE(std::filesystem::exists(tetMesherBin));

  ScopedTempDir tempDir;
  const std::filesystem::path inputObj = tempDir.path() / "hollow_cube.obj";
  const std::filesystem::path outputSurface = tempDir.path() / "hollow_cube_surface.obj";
  const std::filesystem::path configPath = tempDir.path() / "tetgen_hollow.json";
  writeHollowCubeObj(inputObj);
  writeTetgenHollowConfig(configPath);

  const std::string cmd = shellExecutable(tetMesherBin) +
    " --config " + quotePath(configPath);

  ASSERT_EQ(runCommand(cmd), 0);
  ASSERT_TRUE(std::filesystem::exists(outputSurface));

  pgo::Mesh::TriMeshGeo surface;
  ASSERT_TRUE(surface.load(outputSurface.string()));
  EXPECT_EQ(countSurfaceComponents(surface), 2);
}

TEST(TetMesherCli, JsonTetwildGeneratesVegAndSurfaceAssetsWhenEnabled)
{
  if constexpr (!PGO_TEST_TET_MESHER_HAS_TET_WILD) {
    GTEST_SKIP() << "TetWild backend is not enabled in this build.";
  }

  const std::filesystem::path tetMesherBin = getTetMesherBinaryPath();
  ASSERT_FALSE(tetMesherBin.empty());
  ASSERT_TRUE(std::filesystem::exists(tetMesherBin));

  ScopedTempDir tempDir;
  const std::filesystem::path inputObj = tempDir.path() / "cube.obj";
  const std::filesystem::path outputVeg = tempDir.path() / "cube_tetwild.veg";
  const std::filesystem::path outputSurface = tempDir.path() / "cube_tetwild_surface.obj";
  const std::filesystem::path configPath = tempDir.path() / "tetwild.json";
  writeUnitCubeObj(inputObj);
  writeTetwildConfig(configPath);

  const std::string cmd = shellExecutable(tetMesherBin) +
    " --config " + quotePath(configPath);

  ASSERT_EQ(runCommand(cmd), 0);
  ASSERT_TRUE(std::filesystem::exists(outputVeg));
  ASSERT_TRUE(std::filesystem::exists(outputSurface));

  pgo::VolumetricMeshes::TetMesh tetMesh(outputVeg.string().c_str());
  EXPECT_GT(tetMesh.getNumVertices(), 0);
  EXPECT_GT(tetMesh.getNumElements(), 0);

  pgo::Mesh::TriMeshGeo surface;
  ASSERT_TRUE(surface.load(outputSurface.string()));
  EXPECT_GT(surface.numVertices(), 0);
  EXPECT_GT(surface.numTriangles(), 0);
}

TEST(TetMesherCli, OldSubcommandEntryIsRejected)
{
  const std::filesystem::path tetMesherBin = getTetMesherBinaryPath();
  ASSERT_FALSE(tetMesherBin.empty());
  ASSERT_TRUE(std::filesystem::exists(tetMesherBin));

  ScopedTempDir tempDir;
  const std::filesystem::path stderrPath = tempDir.path() / "old-entry.stderr";

  const std::string cmd = shellExecutable(tetMesherBin) +
    " tetgen --help 2> " + quotePath(stderrPath);

  EXPECT_NE(runCommand(cmd), 0);

  const std::string stderrText = readTextFile(stderrPath);
  EXPECT_NE(stderrText.find("--config"), std::string::npos);
}
