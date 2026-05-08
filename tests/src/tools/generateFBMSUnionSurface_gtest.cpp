#include "boundingBox.h"
#include "triMeshGeo.h"

#include <gtest/gtest.h>

#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <stdexcept>
#include <string>

#ifndef PGO_TEST_GENERATE_FBMS_UNION_SURFACE_BIN
#define PGO_TEST_GENERATE_FBMS_UNION_SURFACE_BIN ""
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

std::filesystem::path generatorBinaryPath()
{
  const char *envBin = std::getenv("GENERATE_FBMS_UNION_SURFACE_BIN");
  if (envBin != nullptr && *envBin != '\0')
    return std::filesystem::path(envBin);

  if (std::string(PGO_TEST_GENERATE_FBMS_UNION_SURFACE_BIN).empty())
    return {};

  return std::filesystem::path(PGO_TEST_GENERATE_FBMS_UNION_SURFACE_BIN);
}

class ScopedTempDir
{
public:
  ScopedTempDir()
  {
    const auto base = std::filesystem::temp_directory_path();
    for (int attempt = 0; attempt < 32; ++attempt) {
      const auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
      path_ = base / ("libpgo-fbms-union-test-" + std::to_string(stamp) + "-" + std::to_string(std::rand()) + "-" + std::to_string(attempt));

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

void writeOpenFbmsPatch(const std::filesystem::path &path)
{
  std::ofstream out(path);
  ASSERT_TRUE(out.is_open());

  out << "v -0.35 -0.35 0\n";
  out << "v 0.35 -0.35 0\n";
  out << "v 0.35 0.35 0\n";
  out << "v -0.35 0.35 0\n";
  out << "f 1 2 3\n";
  out << "f 1 3 4\n";
}

void writeHeaderSphereOctahedron(const std::filesystem::path &path)
{
  std::ofstream out(path);
  ASSERT_TRUE(out.is_open());

  out << "# center = (0, 0, 0)\n";
  out << "# radius = 1\n";
  out << "v 1 0 0\n";
  out << "v -1 0 0\n";
  out << "v 0 1 0\n";
  out << "v 0 -1 0\n";
  out << "v 0 0 1\n";
  out << "v 0 0 -1\n";
  out << "f 1 3 5\n";
  out << "f 3 2 5\n";
  out << "f 2 4 5\n";
  out << "f 4 1 5\n";
  out << "f 3 1 6\n";
  out << "f 2 3 6\n";
  out << "f 4 2 6\n";
  out << "f 1 4 6\n";
}

std::string readTextFile(const std::filesystem::path &path)
{
  std::ifstream in(path);
  EXPECT_TRUE(in.is_open());
  return std::string(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
}

}  // namespace

TEST(GenerateFBMSUnionSurfaceCli, CreatesNonEmptyUnionShellFromOpenPatchAndHeaderSphere)
{
  const std::filesystem::path binary = generatorBinaryPath();
  ASSERT_FALSE(binary.empty());
  ASSERT_TRUE(std::filesystem::exists(binary)) << binary;

  ScopedTempDir tempDir;
  const std::filesystem::path fbms = tempDir.path() / "patch.obj";
  const std::filesystem::path sphere = tempDir.path() / "sphere.obj";
  const std::filesystem::path output = tempDir.path() / "union.obj";
  const std::filesystem::path log = tempDir.path() / "generate.log";

  writeOpenFbmsPatch(fbms);
  writeHeaderSphereOctahedron(sphere);

  const std::string command =
    shellExecutable(binary) +
    " --fbms " + quotePath(fbms) +
    " --sphere " + quotePath(sphere) +
    " --fbms-thickness 0.20"
    " --sphere-thickness 0.20"
    " --resolution 28"
    " --padding-ratio 0.02"
    " --output-surface " + quotePath(output) +
    " > " + quotePath(log) + " 2>&1";

  ASSERT_EQ(std::system(command.c_str()), 0) << command;
  ASSERT_TRUE(std::filesystem::exists(output));

  pgo::Mesh::TriMeshGeo mesh;
  ASSERT_TRUE(mesh.load(output.string()));
  EXPECT_GT(mesh.numVertices(), 0);
  EXPECT_GT(mesh.numTriangles(), 0);

  const pgo::Mesh::BoundingBox bb(mesh.positions());
  EXPECT_LT(bb.bmin()[0], -0.95);
  EXPECT_GT(bb.bmax()[0], 0.95);
  EXPECT_GT(bb.bmax()[0], 1.02);
  EXPECT_LT(bb.bmax()[0], 1.30);
}

TEST(GenerateFBMSUnionSurfaceCli, SupportsFixedIsoOffsetExtraction)
{
  const std::filesystem::path binary = generatorBinaryPath();
  ASSERT_FALSE(binary.empty());
  ASSERT_TRUE(std::filesystem::exists(binary)) << binary;

  ScopedTempDir tempDir;
  const std::filesystem::path fbms = tempDir.path() / "patch.obj";
  const std::filesystem::path sphere = tempDir.path() / "sphere.obj";
  const std::filesystem::path output = tempDir.path() / "union-offset.obj";
  const std::filesystem::path log = tempDir.path() / "generate-offset.log";

  writeOpenFbmsPatch(fbms);
  writeHeaderSphereOctahedron(sphere);

  const std::string command =
    shellExecutable(binary) +
    " --fbms " + quotePath(fbms) +
    " --sphere " + quotePath(sphere) +
    " --fbms-thickness 0.20"
    " --sphere-thickness 0.20"
    " --resolution 28"
    " --padding-ratio 0.02"
    " --iso-offset-mode fixed"
    " --iso-offset 0.001"
    " --output-surface " + quotePath(output) +
    " > " + quotePath(log) + " 2>&1";

  ASSERT_EQ(std::system(command.c_str()), 0) << command;
  ASSERT_TRUE(std::filesystem::exists(output));

  pgo::Mesh::TriMeshGeo mesh;
  ASSERT_TRUE(mesh.load(output.string()));
  EXPECT_GT(mesh.numVertices(), 0);
  EXPECT_GT(mesh.numTriangles(), 0);
}

TEST(GenerateFBMSUnionSurfaceCli, SupportsMarchingCubesSubcommand)
{
  const std::filesystem::path binary = generatorBinaryPath();
  ASSERT_FALSE(binary.empty());
  ASSERT_TRUE(std::filesystem::exists(binary)) << binary;

  ScopedTempDir tempDir;
  const std::filesystem::path fbms = tempDir.path() / "patch.obj";
  const std::filesystem::path sphere = tempDir.path() / "sphere.obj";
  const std::filesystem::path output = tempDir.path() / "union-subcommand.obj";
  const std::filesystem::path log = tempDir.path() / "generate-subcommand.log";

  writeOpenFbmsPatch(fbms);
  writeHeaderSphereOctahedron(sphere);

  const std::string command =
    shellExecutable(binary) +
    " marching-cubes"
    " --fbms " + quotePath(fbms) +
    " --sphere " + quotePath(sphere) +
    " --fbms-thickness 0.20"
    " --sphere-thickness 0.20"
    " --resolution 28"
    " --padding-ratio 0.02"
    " --iso-offset-mode fixed"
    " --iso-offset 0.001"
    " --output-surface " + quotePath(output) +
    " > " + quotePath(log) + " 2>&1";

  ASSERT_EQ(std::system(command.c_str()), 0) << command;
  ASSERT_TRUE(std::filesystem::exists(output));

  pgo::Mesh::TriMeshGeo mesh;
  ASSERT_TRUE(mesh.load(output.string()));
  EXPECT_GT(mesh.numVertices(), 0);
  EXPECT_GT(mesh.numTriangles(), 0);
}

TEST(GenerateFBMSUnionSurfaceCli, SupportsSmallComponentFilterOptions)
{
  const std::filesystem::path binary = generatorBinaryPath();
  ASSERT_FALSE(binary.empty());
  ASSERT_TRUE(std::filesystem::exists(binary)) << binary;

  ScopedTempDir tempDir;
  const std::filesystem::path fbms = tempDir.path() / "patch.obj";
  const std::filesystem::path sphere = tempDir.path() / "sphere.obj";
  const std::filesystem::path output = tempDir.path() / "union-filtered.obj";
  const std::filesystem::path log = tempDir.path() / "generate-filtered.log";

  writeOpenFbmsPatch(fbms);
  writeHeaderSphereOctahedron(sphere);

  const std::string command =
    shellExecutable(binary) +
    " marching-cubes"
    " --fbms " + quotePath(fbms) +
    " --sphere " + quotePath(sphere) +
    " --fbms-thickness 0.20"
    " --sphere-thickness 0.20"
    " --resolution 28"
    " --padding-ratio 0.02"
    " --output-surface " + quotePath(output) +
    " --filter-small-components"
    " --min-component-triangles 2"
    " > " + quotePath(log) + " 2>&1";

  ASSERT_EQ(std::system(command.c_str()), 0) << command << "\n" << readTextFile(log);
  ASSERT_TRUE(std::filesystem::exists(output));

  const std::string logText = readTextFile(log);
  EXPECT_NE(logText.find("Small component filter = enabled"), std::string::npos) << logText;
  EXPECT_NE(logText.find("Min component triangles = 2"), std::string::npos) << logText;
}

TEST(GenerateFBMSUnionSurfaceCli, OpenVDBUnionKeepsFbmsShellByDefault)
{
  const std::filesystem::path binary = generatorBinaryPath();
  ASSERT_FALSE(binary.empty());
  ASSERT_TRUE(std::filesystem::exists(binary)) << binary;

  ScopedTempDir tempDir;
  const std::filesystem::path fbms = tempDir.path() / "patch.obj";
  const std::filesystem::path sphere = tempDir.path() / "sphere.obj";
  const std::filesystem::path unionOutput = tempDir.path() / "union-openvdb.obj";
  const std::filesystem::path sphereOutput = tempDir.path() / "sphere-openvdb.obj";
  const std::filesystem::path differenceOutput = tempDir.path() / "difference-openvdb.obj";
  const std::filesystem::path unionLog = tempDir.path() / "union-openvdb.log";
  const std::filesystem::path sphereLog = tempDir.path() / "sphere-openvdb.log";
  const std::filesystem::path differenceLog = tempDir.path() / "difference-openvdb.log";

  writeOpenFbmsPatch(fbms);
  writeHeaderSphereOctahedron(sphere);

  const std::string commonArgs =
    " --fbms " + quotePath(fbms) +
    " --sphere " + quotePath(sphere) +
    " --fbms-thickness 0.20"
    " --sphere-thickness 0.20"
    " --resolution 36"
    " --padding-ratio 0.02";

  const std::string unionCommand =
    shellExecutable(binary) +
    " openvdb" +
    commonArgs +
    " --output-surface " + quotePath(unionOutput) +
    " > " + quotePath(unionLog) + " 2>&1";

  const int unionStatus = std::system(unionCommand.c_str());
  if (unionStatus != 0 && readTextFile(unionLog).find("OpenVDB backend is unavailable") != std::string::npos)
    GTEST_SKIP() << "OpenVDB backend is unavailable in this build.";
  ASSERT_EQ(unionStatus, 0) << unionCommand << "\n" << readTextFile(unionLog);

  const std::string sphereCommand =
    shellExecutable(binary) +
    " openvdb" +
    commonArgs +
    " --debug-field-mode sphere"
    " --output-surface " + quotePath(sphereOutput) +
    " > " + quotePath(sphereLog) + " 2>&1";
  ASSERT_EQ(std::system(sphereCommand.c_str()), 0) << sphereCommand << "\n" << readTextFile(sphereLog);

  const std::string differenceCommand =
    shellExecutable(binary) +
    " openvdb" +
    commonArgs +
    " --surface-mode union-minus-sphere"
    " --output-surface " + quotePath(differenceOutput) +
    " > " + quotePath(differenceLog) + " 2>&1";
  ASSERT_EQ(std::system(differenceCommand.c_str()), 0) << differenceCommand << "\n" << readTextFile(differenceLog);

  pgo::Mesh::TriMeshGeo unionMesh;
  pgo::Mesh::TriMeshGeo sphereMesh;
  pgo::Mesh::TriMeshGeo differenceMesh;
  ASSERT_TRUE(unionMesh.load(unionOutput.string()));
  ASSERT_TRUE(sphereMesh.load(sphereOutput.string()));
  ASSERT_TRUE(differenceMesh.load(differenceOutput.string()));

  EXPECT_GT(unionMesh.numTriangles(), sphereMesh.numTriangles());
  EXPECT_GT(differenceMesh.numTriangles(), 0);
  EXPECT_LT(differenceMesh.numTriangles(), unionMesh.numTriangles());
}

TEST(GenerateFBMSUnionSurfaceCli, OpenVDBSupportsVolumeBudgetThickness)
{
  const std::filesystem::path binary = generatorBinaryPath();
  ASSERT_FALSE(binary.empty());
  ASSERT_TRUE(std::filesystem::exists(binary)) << binary;

  ScopedTempDir tempDir;
  const std::filesystem::path fbms = tempDir.path() / "patch.obj";
  const std::filesystem::path sphere = tempDir.path() / "sphere.obj";
  const std::filesystem::path output = tempDir.path() / "union-openvdb-budget.obj";
  const std::filesystem::path log = tempDir.path() / "union-openvdb-budget.log";

  writeOpenFbmsPatch(fbms);
  writeHeaderSphereOctahedron(sphere);

  const std::string command =
    shellExecutable(binary) +
    " openvdb"
    " --fbms " + quotePath(fbms) +
    " --sphere " + quotePath(sphere) +
    " --fbms-thickness -3.0"
    " --sphere-thickness 0.20"
    " --resolution 36"
    " --padding-ratio 0.02"
    " --output-surface " + quotePath(output) +
    " > " + quotePath(log) + " 2>&1";

  const int status = std::system(command.c_str());
  if (status != 0 && readTextFile(log).find("OpenVDB backend is unavailable") != std::string::npos)
    GTEST_SKIP() << "OpenVDB backend is unavailable in this build.";
  ASSERT_EQ(status, 0) << command << "\n" << readTextFile(log);
  ASSERT_TRUE(std::filesystem::exists(output));

  const std::string logText = readTextFile(log);
  EXPECT_NE(logText.find("[volume-search] selected thickness"), std::string::npos) << logText;

  pgo::Mesh::TriMeshGeo mesh;
  ASSERT_TRUE(mesh.load(output.string()));
  EXPECT_GT(mesh.numVertices(), 0);
  EXPECT_GT(mesh.numTriangles(), 0);
}
