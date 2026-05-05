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

std::filesystem::path getMeshQualityCheckBinaryPath()
{
  const char *envBin = std::getenv("MESH_QUALITY_CHECK_BIN");
  if (envBin != nullptr && *envBin != '\0')
    return std::filesystem::path(envBin);

  if (std::string(PGO_TEST_MESH_QUALITY_CHECK_BIN).empty())
    return {};

  return std::filesystem::path(PGO_TEST_MESH_QUALITY_CHECK_BIN);
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
      path_ = base / ("libpgo-meshQualityCheck-test-" + std::to_string(stamp) + "-" + std::to_string(std::rand()) + "-" + std::to_string(attempt));

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

void writeUnitCubeObj(const std::filesystem::path &path)
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
      { 1, 2, 6 }, { 1, 6, 5 }
    });
}

void writeOpenMeshObj(const std::filesystem::path &path)
{
  writeObj(path,
    {
      { 0.0, 0.0, 0.0 }, { 1.0, 0.0, 0.0 }, { 1.0, 1.0, 0.0 }, { 0.0, 1.0, 0.0 }
    },
    {
      { 0, 1, 2 }, { 0, 2, 3 }
    });
}

void writeSelfIntersectingMeshObj(const std::filesystem::path &path)
{
  writeObj(path,
    {
      { 0.0, 0.0, 0.0 }, { 1.0, 0.0, 0.0 }, { 0.5, 1.0, 0.0 },
      { 0.5, 0.25, -1.0 }, { 0.5, 0.25, 1.0 }, { 0.5, 0.75, 0.0 }
    },
    {
      { 0, 1, 2 }, { 3, 4, 5 }
    });
}

void writeNonManifoldMeshObj(const std::filesystem::path &path)
{
  writeObj(path,
    {
      { 0.0, 0.0, 0.0 }, { 1.0, 0.0, 0.0 }, { 0.5, 1.0, 0.0 }, { 0.5, -1.0, 0.0 }, { 0.5, 0.0, 1.0 }
    },
    {
      { 0, 1, 2 }, { 1, 0, 3 }, { 0, 1, 4 }
    });
}

void writeDegenerateTriangleMeshObj(const std::filesystem::path &path)
{
  writeObj(path,
    {
      { 0.0, 0.0, 0.0 }, { 1.0, 0.0, 0.0 }, { 2.0, 0.0, 0.0 }
    },
    {
      { 0, 1, 2 }
    });
}

void writeClosedCubeWithInconsistentOrientationObj(const std::filesystem::path &path)
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
      { 1, 6, 2 }, { 1, 6, 5 }
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
  const std::filesystem::path binary = getMeshQualityCheckBinaryPath();
  if (binary.empty())
    throw std::runtime_error("meshQualityCheck binary path is not configured");
  return binary;
}

}  // namespace

TEST(meshQualityCheckTool, AcceptsClosedCube)
{
  ScopedTempDir tempDir;
  const std::filesystem::path inputObj = tempDir.path() / "cube.obj";
  const std::filesystem::path reportPath = tempDir.path() / "cube.quality.json";
  writeUnitCubeObj(inputObj);

  const std::string command = shellExecutable(requireBinary()) + " surface --input " + quotePath(inputObj) + " --json " + quotePath(reportPath);
  EXPECT_EQ(runCommand(command), 0);

  const nlohmann::json report = readJson(reportPath);
  EXPECT_EQ(report["type"], "surface");
  EXPECT_EQ(report["input"], inputObj.string());
  EXPECT_TRUE(report["passed"]);
  EXPECT_EQ(report["vertices"], 8);
  EXPECT_EQ(report["triangles"], 12);
  EXPECT_EQ(report["boundary_or_exterior_edges"], 0);
  EXPECT_TRUE(report["is_manifold"]);
  EXPECT_EQ(report["self_intersections"], 0);
  EXPECT_EQ(report["self_intersection_backend"], "cgal-bool");
  EXPECT_FALSE(report["self_intersections_exact"]);
  EXPECT_EQ(report["check_status"]["basic_geometry"], "passed");
  EXPECT_EQ(report["check_status"]["topology"], "passed");
  EXPECT_EQ(report["check_status"]["winding"], "passed");
  EXPECT_EQ(report["check_status"]["self_intersection"], "passed");
  EXPECT_TRUE(report["is_winding_consistent"]);
  EXPECT_EQ(report["oriented_boundary_or_exterior_edges"], 0);
  ASSERT_TRUE(report.contains("enclosed_volume"));
  ASSERT_TRUE(report["enclosed_volume"].is_number());
  EXPECT_NEAR(report["enclosed_volume"].get<double>(), 1.0, 1e-12);
}

TEST(meshQualityCheckTool, ExpectedComponentsAcceptsMatchingTopology)
{
  ScopedTempDir tempDir;
  const std::filesystem::path inputObj = tempDir.path() / "cube.obj";
  const std::filesystem::path reportPath = tempDir.path() / "cube.expected.quality.json";
  writeUnitCubeObj(inputObj);

  const std::string command = shellExecutable(requireBinary()) +
    " surface --expected-components 1 --input " + quotePath(inputObj) +
    " --json " + quotePath(reportPath);
  EXPECT_EQ(runCommand(command), 0);

  const nlohmann::json report = readJson(reportPath);
  EXPECT_TRUE(report["passed"]);
  EXPECT_EQ(report["components_by_edge"], 1);
  EXPECT_TRUE(report["errors"].empty());
}

TEST(meshQualityCheckTool, ExpectedComponentsRejectsMismatchingTopology)
{
  ScopedTempDir tempDir;
  const std::filesystem::path inputObj = tempDir.path() / "cube.obj";
  const std::filesystem::path reportPath = tempDir.path() / "cube.expected_mismatch.quality.json";
  writeUnitCubeObj(inputObj);

  const std::string command = shellExecutable(requireBinary()) +
    " surface --expected-components 3 --input " + quotePath(inputObj) +
    " --json " + quotePath(reportPath);
  EXPECT_NE(runCommand(command), 0);

  const nlohmann::json report = readJson(reportPath);
  EXPECT_FALSE(report["passed"]);
  EXPECT_EQ(report["components_by_edge"], 1);
  EXPECT_EQ(report["check_status"]["topology"], "failed");

  bool foundExpectedComponentsError = false;
  for (const auto &error : report["errors"]) {
    if (error.get<std::string>().find("Expected 3 edge-connected components, found 1") != std::string::npos)
      foundExpectedComponentsError = true;
  }
  EXPECT_TRUE(foundExpectedComponentsError) << report.dump(2);
}

TEST(meshQualityCheckTool, RejectsOpenMeshAndWritesReport)
{
  ScopedTempDir tempDir;
  const std::filesystem::path inputObj = tempDir.path() / "open.obj";
  const std::filesystem::path reportPath = tempDir.path() / "open.quality.json";
  writeOpenMeshObj(inputObj);

  const std::string command = shellExecutable(requireBinary()) + " surface --input " + quotePath(inputObj) + " --json " + quotePath(reportPath);
  EXPECT_NE(runCommand(command), 0);

  const nlohmann::json report = readJson(reportPath);
  EXPECT_FALSE(report["passed"]);
  EXPECT_GT(report["boundary_or_exterior_edges"], 0);
  EXPECT_EQ(report["check_status"]["basic_geometry"], "passed");
  EXPECT_EQ(report["check_status"]["topology"], "failed");
  EXPECT_EQ(report["check_status"]["winding"], "failed");
  EXPECT_EQ(report["check_status"]["self_intersection"], "passed");
  EXPECT_FALSE(report["is_winding_consistent"]);
  EXPECT_GT(report["oriented_boundary_or_exterior_edges"], 0);
  ASSERT_TRUE(report.contains("enclosed_volume"));
  EXPECT_TRUE(report["enclosed_volume"].is_null());
}

TEST(meshQualityCheckTool, RejectsSelfIntersectingMesh)
{
  ScopedTempDir tempDir;
  const std::filesystem::path inputObj = tempDir.path() / "self_intersect.obj";
  const std::filesystem::path reportPath = tempDir.path() / "self_intersect.quality.json";
  writeSelfIntersectingMeshObj(inputObj);

  const std::string command = shellExecutable(requireBinary()) + " surface --input " + quotePath(inputObj) + " --json " + quotePath(reportPath);
  EXPECT_NE(runCommand(command), 0);

  const nlohmann::json report = readJson(reportPath);
  EXPECT_FALSE(report["passed"]);
  EXPECT_EQ(report["self_intersections"], 1);
  EXPECT_EQ(report["self_intersection_backend"], "cgal-bool");
  EXPECT_FALSE(report["self_intersections_exact"]);
  EXPECT_EQ(report["check_status"]["basic_geometry"], "passed");
  EXPECT_EQ(report["check_status"]["topology"], "failed");
  EXPECT_EQ(report["check_status"]["winding"], "failed");
  EXPECT_EQ(report["check_status"]["self_intersection"], "failed");
}

TEST(meshQualityCheckTool, ExactCountBackendReportsPreciseSelfIntersectionCount)
{
  ScopedTempDir tempDir;
  const std::filesystem::path inputObj = tempDir.path() / "self_intersect.obj";
  const std::filesystem::path reportPath = tempDir.path() / "self_intersect.exact.quality.json";
  writeSelfIntersectingMeshObj(inputObj);

  const std::string command = shellExecutable(requireBinary()) +
    " surface --self-intersection-backend exact-count --input " + quotePath(inputObj) +
    " --json " + quotePath(reportPath);
  EXPECT_NE(runCommand(command), 0);

  const nlohmann::json report = readJson(reportPath);
  EXPECT_FALSE(report["passed"]);
  EXPECT_GT(report["self_intersections"], 0);
  EXPECT_EQ(report["self_intersection_backend"], "exact-count");
  EXPECT_TRUE(report["self_intersections_exact"]);
  EXPECT_EQ(report["check_status"]["self_intersection"], "failed");
}

TEST(meshQualityCheckTool, RejectsNonManifoldMesh)
{
  ScopedTempDir tempDir;
  const std::filesystem::path inputObj = tempDir.path() / "non_manifold.obj";
  const std::filesystem::path reportPath = tempDir.path() / "non_manifold.quality.json";
  writeNonManifoldMeshObj(inputObj);

  const std::string command = shellExecutable(requireBinary()) + " surface --input " + quotePath(inputObj) + " --json " + quotePath(reportPath);
  EXPECT_NE(runCommand(command), 0);

  const nlohmann::json report = readJson(reportPath);
  EXPECT_FALSE(report["passed"]);
  EXPECT_FALSE(report["is_manifold"]);
  EXPECT_EQ(report["check_status"]["basic_geometry"], "passed");
  EXPECT_EQ(report["check_status"]["topology"], "failed");
  EXPECT_EQ(report["check_status"]["winding"], "failed");
  EXPECT_EQ(report["check_status"]["self_intersection"], "passed");
}

TEST(meshQualityCheckTool, MarksTopologyAndSelfIntersectionUnknownForDegenerateMesh)
{
  ScopedTempDir tempDir;
  const std::filesystem::path inputObj = tempDir.path() / "degenerate.obj";
  const std::filesystem::path reportPath = tempDir.path() / "degenerate.quality.json";
  writeDegenerateTriangleMeshObj(inputObj);

  const std::string command = shellExecutable(requireBinary()) + " surface --input " + quotePath(inputObj) + " --json " + quotePath(reportPath);
  EXPECT_NE(runCommand(command), 0);

  const nlohmann::json report = readJson(reportPath);
  EXPECT_FALSE(report["passed"]);
  EXPECT_GT(report["invalid_triangles"], 0);
  EXPECT_FALSE(report["components_by_edge"].is_null());
  EXPECT_FALSE(report["boundary_or_exterior_edges"].is_null());
  EXPECT_FALSE(report["is_manifold"].is_null());
  EXPECT_FALSE(report["self_intersections"].is_null());
  EXPECT_FALSE(report["is_winding_consistent"].is_null());
  EXPECT_FALSE(report["oriented_boundary_or_exterior_edges"].is_null());
  EXPECT_EQ(report["check_status"]["basic_geometry"], "failed");
  EXPECT_NE(report["check_status"]["topology"], "skipped");
  EXPECT_NE(report["check_status"]["winding"], "skipped");
  EXPECT_NE(report["check_status"]["self_intersection"], "skipped");
}

TEST(meshQualityCheckTool, DegenerateOnlyModeSkipsTopologyWindingAndSelfIntersection)
{
  ScopedTempDir tempDir;
  const std::filesystem::path inputObj = tempDir.path() / "open.obj";
  const std::filesystem::path reportPath = tempDir.path() / "open.quality.json";
  writeOpenMeshObj(inputObj);

  const std::string command = shellExecutable(requireBinary()) + " surface --check-level degenerate-only --input " + quotePath(inputObj) + " --json " + quotePath(reportPath);
  EXPECT_EQ(runCommand(command), 0);

  const nlohmann::json report = readJson(reportPath);
  EXPECT_TRUE(report["passed"]);
  EXPECT_EQ(report["invalid_triangles"], 0);
  EXPECT_EQ(report["check_status"]["basic_geometry"], "passed");
  EXPECT_EQ(report["check_status"]["topology"], "skipped");
  EXPECT_EQ(report["check_status"]["winding"], "skipped");
  EXPECT_EQ(report["check_status"]["self_intersection"], "skipped");
  EXPECT_TRUE(report["components_by_edge"].is_null());
  EXPECT_TRUE(report["boundary_or_exterior_edges"].is_null());
  EXPECT_TRUE(report["is_manifold"].is_null());
  EXPECT_TRUE(report["is_winding_consistent"].is_null());
  EXPECT_TRUE(report["oriented_boundary_or_exterior_edges"].is_null());
  EXPECT_TRUE(report["self_intersections"].is_null());
  EXPECT_TRUE(report["self_intersection_backend"].is_null());
  EXPECT_TRUE(report["self_intersections_exact"].is_null());
  ASSERT_TRUE(report.contains("enclosed_volume"));
  EXPECT_TRUE(report["enclosed_volume"].is_null());
}

TEST(meshQualityCheckTool, SkipsSelfIntersectionWhenTriangleLimitIsExceeded)
{
  ScopedTempDir tempDir;
  const std::filesystem::path inputObj = tempDir.path() / "cube.obj";
  const std::filesystem::path reportPath = tempDir.path() / "cube.limit.quality.json";
  writeUnitCubeObj(inputObj);

  const std::string command = shellExecutable(requireBinary()) +
    " surface --self-intersection-triangle-limit 1 --input " + quotePath(inputObj) +
    " --json " + quotePath(reportPath);
  EXPECT_NE(runCommand(command), 0);

  const nlohmann::json report = readJson(reportPath);
  EXPECT_FALSE(report["passed"]);
  EXPECT_EQ(report["check_status"]["self_intersection"], "skipped");
  EXPECT_TRUE(report["self_intersections"].is_null());
  EXPECT_TRUE(report["self_intersection_backend"].is_null());
  EXPECT_TRUE(report["self_intersections_exact"].is_null());
}

TEST(meshQualityCheckTool, RawModeReportsTopologyComponentsAndAllowsDegenerateTriangles)
{
  ScopedTempDir tempDir;
  const std::filesystem::path inputObj = tempDir.path() / "degenerate.obj";
  const std::filesystem::path reportPath = tempDir.path() / "degenerate.raw.quality.json";
  writeDegenerateTriangleMeshObj(inputObj);

  const std::string command = shellExecutable(requireBinary()) + " surface --check-level raw --input " + quotePath(inputObj) + " --json " + quotePath(reportPath);
  EXPECT_NE(runCommand(command), 0);

  const nlohmann::json report = readJson(reportPath);
  EXPECT_FALSE(report["passed"]);
  EXPECT_GT(report["invalid_triangles"], 0);
  EXPECT_EQ(report["check_status"]["basic_geometry"], "passed");
  EXPECT_NE(report["check_status"]["topology"], "skipped");
  EXPECT_EQ(report["check_status"]["winding"], "skipped");
  EXPECT_EQ(report["check_status"]["self_intersection"], "skipped");
  EXPECT_FALSE(report["components_by_edge"].is_null());
  EXPECT_FALSE(report["component_triangle_counts_by_edge"].is_null());
  EXPECT_EQ(report["components_by_edge"], report["component_triangle_counts_by_edge"].size());
  EXPECT_FALSE(report["is_manifold"].is_null());
  EXPECT_TRUE(report["is_winding_consistent"].is_null());
  EXPECT_TRUE(report["self_intersections"].is_null());
  EXPECT_TRUE(report["self_intersection_backend"].is_null());
  EXPECT_TRUE(report["self_intersections_exact"].is_null());
}

TEST(meshQualityCheckTool, RejectsInvalidSelfIntersectionBackend)
{
  ScopedTempDir tempDir;
  const std::filesystem::path inputObj = tempDir.path() / "cube.obj";
  const std::filesystem::path reportPath = tempDir.path() / "cube.invalid_backend.quality.json";
  writeUnitCubeObj(inputObj);

  const std::string command = shellExecutable(requireBinary()) +
    " surface --self-intersection-backend not-a-backend --input " + quotePath(inputObj) +
    " --json " + quotePath(reportPath);
  EXPECT_NE(runCommand(command), 0);
  EXPECT_FALSE(std::filesystem::exists(reportPath));
}

TEST(meshQualityCheckTool, RawModeReportsEnclosedVolumeForClosedMesh)
{
  ScopedTempDir tempDir;
  const std::filesystem::path inputObj = tempDir.path() / "cube.obj";
  const std::filesystem::path reportPath = tempDir.path() / "cube.raw.quality.json";
  writeUnitCubeObj(inputObj);

  const std::string command = shellExecutable(requireBinary()) + " surface --check-level raw --input " + quotePath(inputObj) + " --json " + quotePath(reportPath);
  EXPECT_EQ(runCommand(command), 0);

  const nlohmann::json report = readJson(reportPath);
  EXPECT_TRUE(report["passed"]);
  EXPECT_EQ(report["check_status"]["topology"], "passed");
  EXPECT_EQ(report["check_status"]["winding"], "skipped");
  ASSERT_TRUE(report.contains("enclosed_volume"));
  ASSERT_TRUE(report["enclosed_volume"].is_number());
  EXPECT_NEAR(report["enclosed_volume"].get<double>(), 1.0, 1e-12);
}

TEST(meshQualityCheckTool, RejectsClosedMeshWithInconsistentOrientation)
{
  ScopedTempDir tempDir;
  const std::filesystem::path inputObj = tempDir.path() / "inconsistent_orientation.obj";
  const std::filesystem::path reportPath = tempDir.path() / "inconsistent_orientation.quality.json";
  writeClosedCubeWithInconsistentOrientationObj(inputObj);

  const std::string command = shellExecutable(requireBinary()) + " surface --input " + quotePath(inputObj) + " --json " + quotePath(reportPath);
  EXPECT_NE(runCommand(command), 0);

  const nlohmann::json report = readJson(reportPath);
  EXPECT_FALSE(report["passed"]);
  EXPECT_TRUE(report["is_manifold"]);
  EXPECT_EQ(report["boundary_or_exterior_edges"], 0);
  EXPECT_FALSE(report["is_winding_consistent"]);
  EXPECT_GT(report["oriented_boundary_or_exterior_edges"], 0);
  EXPECT_EQ(report["check_status"]["basic_geometry"], "passed");
  EXPECT_EQ(report["check_status"]["topology"], "passed");
  EXPECT_EQ(report["check_status"]["winding"], "failed");
  EXPECT_EQ(report["check_status"]["self_intersection"], "passed");
}

TEST(meshQualityCheckTool, InvalidTrianglesWarnPolicyWritesWarningAndDoesNotFailBasicGeometry)
{
  ScopedTempDir tempDir;
  const std::filesystem::path inputObj = tempDir.path() / "degenerate.obj";
  const std::filesystem::path reportPath = tempDir.path() / "degenerate.warn.quality.json";
  writeDegenerateTriangleMeshObj(inputObj);

  const std::string command = shellExecutable(requireBinary()) +
    " surface --check-level degenerate-only --invalid-triangles-policy warn --input " + quotePath(inputObj) +
    " --json " + quotePath(reportPath);
  EXPECT_EQ(runCommand(command), 0);

  const nlohmann::json report = readJson(reportPath);
  EXPECT_TRUE(report["passed"]);
  EXPECT_GT(report["invalid_triangles"], 0);
  EXPECT_EQ(report["check_status"]["basic_geometry"], "passed");
  EXPECT_TRUE(report["errors"].empty());

  bool foundInvalidWarning = false;
  for (const auto &warning : report["warnings"]) {
    if (warning.get<std::string>().find("Input triangle mesh contains invalid or degenerate triangles") != std::string::npos)
      foundInvalidWarning = true;
  }
  EXPECT_TRUE(foundInvalidWarning) << report.dump(2);
}

TEST(meshQualityCheckTool, RejectsMissingInput)
{
  ScopedTempDir tempDir;
  const std::filesystem::path inputObj = tempDir.path() / "missing.obj";
  const std::filesystem::path reportPath = tempDir.path() / "missing.quality.json";

  const std::string command = shellExecutable(requireBinary()) + " surface --input " + quotePath(inputObj) + " --json " + quotePath(reportPath);
  EXPECT_NE(runCommand(command), 0);
  EXPECT_TRUE(std::filesystem::exists(reportPath));

  const nlohmann::json report = readJson(reportPath);
  EXPECT_FALSE(report["passed"]);
  EXPECT_EQ(report["input"], inputObj.string());
  EXPECT_EQ(report["check_status"]["basic_geometry"], "failed");
  EXPECT_EQ(report["check_status"]["topology"], "skipped");
  EXPECT_EQ(report["check_status"]["winding"], "skipped");
  EXPECT_EQ(report["check_status"]["self_intersection"], "skipped");
  EXPECT_TRUE(report["components_by_edge"].is_null());
  EXPECT_TRUE(report["boundary_or_exterior_edges"].is_null());
  EXPECT_TRUE(report["is_manifold"].is_null());
  EXPECT_TRUE(report["self_intersections"].is_null());
  EXPECT_TRUE(report["self_intersection_backend"].is_null());
  EXPECT_TRUE(report["self_intersections_exact"].is_null());
  EXPECT_TRUE(report["is_winding_consistent"].is_null());
  EXPECT_TRUE(report["oriented_boundary_or_exterior_edges"].is_null());
}
