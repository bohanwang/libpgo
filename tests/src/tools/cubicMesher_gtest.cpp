#include "cubicMesherIO.h"
#include "triangleMeshVoxelizer.h"
#include "triMeshGeo.h"
#include "volumetricMeshENuMaterial.h"

#include <gtest/gtest.h>

#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#ifndef PGO_TEST_CUBIC_MESHER_BIN
#define PGO_TEST_CUBIC_MESHER_BIN ""
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

std::filesystem::path getCubicMesherBinaryPath()
{
  const char *envBin = std::getenv("CUBIC_MESHER_BIN");
  if (envBin != nullptr && *envBin != '\0')
    return std::filesystem::path(envBin);

  if (std::string(PGO_TEST_CUBIC_MESHER_BIN).empty())
    return {};

  return std::filesystem::path(PGO_TEST_CUBIC_MESHER_BIN);
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
      path_ = base / ("libpgo-cubicMesher-test-" + std::to_string(stamp) + "-" + std::to_string(std::rand()) + "-" + std::to_string(attempt));

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

void writeObj(const std::filesystem::path &path, const std::vector<pgo::Vec3d> &vertices, const std::vector<pgo::Vec3i> &triangles)
{
  std::ofstream out(path);
  ASSERT_TRUE(out.is_open());

  for (const auto &v : vertices)
    out << "v " << v[0] << " " << v[1] << " " << v[2] << "\n";

  for (const auto &tri : triangles)
    out << "f " << (tri[0] + 1) << " " << (tri[1] + 1) << " " << (tri[2] + 1) << "\n";
}

void writeUnitCubeObj(const std::filesystem::path &path)
{
  writeObj(path,
    {
      pgo::Vec3d(0.0, 0.0, 0.0), pgo::Vec3d(1.0, 0.0, 0.0), pgo::Vec3d(1.0, 1.0, 0.0), pgo::Vec3d(0.0, 1.0, 0.0),
      pgo::Vec3d(0.0, 0.0, 1.0), pgo::Vec3d(1.0, 0.0, 1.0), pgo::Vec3d(1.0, 1.0, 1.0), pgo::Vec3d(0.0, 1.0, 1.0)
    },
    {
      pgo::Vec3i(0, 2, 1), pgo::Vec3i(0, 3, 2),
      pgo::Vec3i(4, 5, 6), pgo::Vec3i(4, 6, 7),
      pgo::Vec3i(0, 1, 5), pgo::Vec3i(0, 5, 4),
      pgo::Vec3i(3, 7, 6), pgo::Vec3i(3, 6, 2),
      pgo::Vec3i(0, 4, 7), pgo::Vec3i(0, 7, 3),
      pgo::Vec3i(1, 2, 6), pgo::Vec3i(1, 6, 5)
    });
}

void writeOpenMeshObj(const std::filesystem::path &path)
{
  writeObj(path,
    {
      pgo::Vec3d(0.0, 0.0, 0.0), pgo::Vec3d(1.0, 0.0, 0.0), pgo::Vec3d(1.0, 1.0, 0.0), pgo::Vec3d(0.0, 1.0, 0.0)
    },
    {
      pgo::Vec3i(0, 1, 2), pgo::Vec3i(0, 2, 3)
    });
}

void writeSelfIntersectingMeshObj(const std::filesystem::path &path)
{
  writeObj(path,
    {
      pgo::Vec3d(0.0, 0.0, 0.0), pgo::Vec3d(1.0, 0.0, 0.0), pgo::Vec3d(0.5, 1.0, 0.0),
      pgo::Vec3d(0.5, 0.25, -1.0), pgo::Vec3d(0.5, 0.25, 1.0), pgo::Vec3d(0.5, 0.75, 0.0)
    },
    {
      pgo::Vec3i(0, 1, 2), pgo::Vec3i(3, 4, 5)
    });
}

void writeNonManifoldMeshObj(const std::filesystem::path &path)
{
  writeObj(path,
    {
      pgo::Vec3d(0.0, 0.0, 0.0), pgo::Vec3d(1.0, 0.0, 0.0), pgo::Vec3d(0.5, 1.0, 0.0), pgo::Vec3d(0.5, -1.0, 0.0), pgo::Vec3d(0.5, 0.0, 1.0)
    },
    {
      pgo::Vec3i(0, 1, 2), pgo::Vec3i(1, 0, 3), pgo::Vec3i(0, 1, 4)
    });
}

}  // namespace

TEST(CubicMesherHelper, CreatesResolution2CubeMeshFromUnitCubeObj)
{
  ScopedTempDir tempDir;
  const std::filesystem::path inputObj = tempDir.path() / "cube.obj";
  writeUnitCubeObj(inputObj);

  cubic_mesher::TriangleMeshVoxelizerOptions options;
  options.inputMesh = inputObj.string();
  options.resolution = 2;
  options.E = 6000.0;
  options.nu = 0.4;
  options.density = 1000.0;

  std::unique_ptr<pgo::VolumetricMeshes::CubicMesh> cubicMesh = cubic_mesher::createTriangleMeshCubicMesh(options);
  ASSERT_NE(cubicMesh, nullptr);
  EXPECT_EQ(cubicMesh->getNumVertices(), 27);
  EXPECT_EQ(cubicMesh->getNumElements(), 8);

  const pgo::Mesh::BoundingBox bb = cubicMesh->getBoundingBox();
  EXPECT_DOUBLE_EQ(bb.bmin()[0], 0.0);
  EXPECT_DOUBLE_EQ(bb.bmin()[1], 0.0);
  EXPECT_DOUBLE_EQ(bb.bmin()[2], 0.0);
  EXPECT_DOUBLE_EQ(bb.bmax()[0], 1.0);
  EXPECT_DOUBLE_EQ(bb.bmax()[1], 1.0);
  EXPECT_DOUBLE_EQ(bb.bmax()[2], 1.0);
}

TEST(CubicMesherHelper, SavesReloadsMaterialAndSurfaceMesh)
{
  ScopedTempDir tempDir;
  const std::filesystem::path inputObj = tempDir.path() / "cube.obj";
  const std::filesystem::path outputVeg = tempDir.path() / "cube.veg";
  const std::filesystem::path outputObj = tempDir.path() / "cube_surface.obj";
  writeUnitCubeObj(inputObj);

  cubic_mesher::TriangleMeshVoxelizerOptions options;
  options.inputMesh = inputObj.string();
  options.resolution = 2;
  options.E = 6000.0;
  options.nu = 0.4;
  options.density = 1000.0;

  std::unique_ptr<pgo::VolumetricMeshes::CubicMesh> cubicMesh = cubic_mesher::createTriangleMeshCubicMesh(options);
  ASSERT_NE(cubicMesh, nullptr);

  cubic_mesher::saveCubicMesh(*cubicMesh, outputVeg.string());
  cubic_mesher::writeSurfaceMesh(*cubicMesh, outputObj.string());

  ASSERT_TRUE(std::filesystem::exists(outputVeg));
  ASSERT_TRUE(std::filesystem::exists(outputObj));

  pgo::VolumetricMeshes::CubicMesh reloadedMesh(outputVeg.string().c_str());
  EXPECT_EQ(reloadedMesh.getNumVertices(), 27);
  EXPECT_EQ(reloadedMesh.getNumElements(), 8);

  const auto *material = pgo::VolumetricMeshes::downcastENuMaterial(reloadedMesh.getMaterial(0));
  ASSERT_NE(material, nullptr);
  EXPECT_DOUBLE_EQ(material->getE(), 6000.0);
  EXPECT_DOUBLE_EQ(material->getNu(), 0.4);
  EXPECT_DOUBLE_EQ(material->getDensity(), 1000.0);

  pgo::Mesh::TriMeshGeo surface;
  ASSERT_TRUE(surface.load(outputObj.string()));
  EXPECT_GT(surface.numVertices(), 0);
  EXPECT_GT(surface.numTriangles(), 0);
}

TEST(CubicMesherHelper, RejectsOpenMesh)
{
  ScopedTempDir tempDir;
  const std::filesystem::path inputObj = tempDir.path() / "open.obj";
  writeOpenMeshObj(inputObj);

  cubic_mesher::TriangleMeshVoxelizerOptions options;
  options.inputMesh = inputObj.string();
  options.resolution = 2;
  options.E = 6000.0;
  options.nu = 0.4;
  options.density = 1000.0;

  EXPECT_THROW(
    {
      try {
        (void)cubic_mesher::createTriangleMeshCubicMesh(options);
      }
      catch (const std::runtime_error &err) {
        EXPECT_NE(std::string(err.what()).find("closed"), std::string::npos);
        throw;
      }
    },
    std::runtime_error);
}

TEST(CubicMesherHelper, RejectsSelfIntersectingMesh)
{
  ScopedTempDir tempDir;
  const std::filesystem::path inputObj = tempDir.path() / "self_intersect.obj";
  writeSelfIntersectingMeshObj(inputObj);

  cubic_mesher::TriangleMeshVoxelizerOptions options;
  options.inputMesh = inputObj.string();
  options.resolution = 2;
  options.E = 6000.0;
  options.nu = 0.4;
  options.density = 1000.0;

  EXPECT_THROW(
    {
      try {
        (void)cubic_mesher::createTriangleMeshCubicMesh(options);
      }
      catch (const std::runtime_error &err) {
        EXPECT_NE(std::string(err.what()).find("self intersections"), std::string::npos);
        throw;
      }
    },
    std::runtime_error);
}

TEST(CubicMesherHelper, RejectsNonManifoldMesh)
{
  ScopedTempDir tempDir;
  const std::filesystem::path inputObj = tempDir.path() / "non_manifold.obj";
  writeNonManifoldMeshObj(inputObj);

  cubic_mesher::TriangleMeshVoxelizerOptions options;
  options.inputMesh = inputObj.string();
  options.resolution = 2;
  options.E = 6000.0;
  options.nu = 0.4;
  options.density = 1000.0;

  EXPECT_THROW(
    {
      try {
        (void)cubic_mesher::createTriangleMeshCubicMesh(options);
      }
      catch (const std::runtime_error &err) {
        EXPECT_NE(std::string(err.what()).find("manifold"), std::string::npos);
        throw;
      }
    },
    std::runtime_error);
}

TEST(CubicMesherHelper, RejectsNonPositiveResolution)
{
  ScopedTempDir tempDir;
  const std::filesystem::path inputObj = tempDir.path() / "cube.obj";
  writeUnitCubeObj(inputObj);

  cubic_mesher::TriangleMeshVoxelizerOptions options;
  options.inputMesh = inputObj.string();
  options.resolution = 0;
  options.E = 6000.0;
  options.nu = 0.4;
  options.density = 1000.0;

  EXPECT_THROW(
    {
      try {
        (void)cubic_mesher::createTriangleMeshCubicMesh(options);
      }
      catch (const std::runtime_error &err) {
        EXPECT_NE(std::string(err.what()).find("Resolution"), std::string::npos);
        throw;
      }
    },
    std::runtime_error);
}

TEST(CubicMesherCli, GeneratesVegAndSurfaceAssets)
{
  const std::filesystem::path cubicMesherBin = getCubicMesherBinaryPath();
  ASSERT_FALSE(cubicMesherBin.empty());
  ASSERT_TRUE(std::filesystem::exists(cubicMesherBin));

  ScopedTempDir tempDir;
  const std::filesystem::path inputObj = tempDir.path() / "cube.obj";
  const std::filesystem::path outputVeg = tempDir.path() / "cube.veg";
  const std::filesystem::path outputObj = tempDir.path() / "cube_surface.obj";
  writeUnitCubeObj(inputObj);

  std::ostringstream cmd;
  cmd << shellExecutable(cubicMesherBin)
      << " --input-mesh " << quotePath(inputObj)
      << " --resolution 2"
      << " --output-mesh " << quotePath(outputVeg)
      << " --output-surface " << quotePath(outputObj)
      << " --E 6000 --nu 0.4 --density 1000";

  ASSERT_EQ(runCommand(cmd.str()), 0);
  ASSERT_TRUE(std::filesystem::exists(outputVeg));
  ASSERT_TRUE(std::filesystem::exists(outputObj));

  pgo::VolumetricMeshes::CubicMesh cubicMesh(outputVeg.string().c_str());
  EXPECT_EQ(cubicMesh.getNumVertices(), 27);
  EXPECT_EQ(cubicMesh.getNumElements(), 8);

  pgo::Mesh::TriMeshGeo surface;
  ASSERT_TRUE(surface.load(outputObj.string()));
  EXPECT_GT(surface.numTriangles(), 0);
}
