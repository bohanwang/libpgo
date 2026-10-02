#include <gtest/gtest.h>

#include "EigenSupport.h"

#include <nlohmann/json.hpp>

#include <array>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <tuple>
#include <vector>

namespace
{
namespace fs = std::filesystem;
namespace ES = pgo::EigenSupport;
using Json = nlohmann::json;
using Row = std::array<std::string, 13>;

std::string quote(const fs::path &path)
{
  std::string result = "'";
  for (char c : path.string())
    result += c == '\'' ? "'\\''" : std::string(1, c);
  return result + "'";
}

std::string scalar(double value)
{
  std::ostringstream stream;
  stream << std::setprecision(17) << value;
  return stream.str();
}

class IPCTrajectoryValidatorTest : public ::testing::Test
{
protected:
  void SetUp() override
  {
    directory = fs::temp_directory_path() / ("pgo-ipc-trajectory-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    ASSERT_TRUE(fs::create_directory(directory));
    output = directory / "states";
    fs::create_directories(output / "substeps");
    std::ofstream(directory / "plane.obj") << "v -2 0 -2\nv 2 0 -2\nv 2 0 2\nv -2 0 2\nf 1 3 2\nf 1 4 3\n";
    const fs::path assets(PGO_TEST_VALIDATOR_BOX_DIR);
    std::ifstream mesh(assets / "box-tet.veg");
    std::string line;
    while (std::getline(mesh, line) && line != "*VERTICES") {
    }
    int count, dimension, attributes, markers;
    ASSERT_TRUE(bool(mesh >> count >> dimension >> attributes >> markers));
    ASSERT_EQ(dimension, 3);
    rest.resize(3 * count);
    for (int i = 0; i < count; ++i) {
      int index;
      ASSERT_TRUE(bool(mesh >> index >> rest[3 * i] >> rest[3 * i + 1] >> rest[3 * i + 2]));
    }
    config = {
      { "tet-mesh", (assets / "box-tet.veg").string() }, { "surface-mesh", (assets / "box.obj").string() },
      { "fixed-vertices", Json::array() },
      { "ipc-external-objects", Json::array({ { { "filename", (directory / "plane.obj").string() }, { "scale", 1.0 }, { "initial-translation", { 0, 0, 0 } }, { "movement", { 0, 0, 0 } } } }) },
      { "g", { 0, -9.81, 0 } }, { "init-vel", { 0, 0, 0 } }, { "init-disp", { 0, 0, 0 } },
      { "scale", 1.0 }, { "timestep", 0.001 }, { "num-timestep", 1 }, { "damping-params", { 0, 0 } },
      { "sim-type", "dynamic" }, { "solver-eps", 1e-6 }, { "solver-max-iter", 200 },
      { "elastic-material", "stable-neo" }, { "dump-interval", 1 }, { "output", output.string() },
      { "ipc-dhat", 0.002 }, { "ipc-kappa", 3000.0 }, { "contact-model", "ipc" },
      { "enable-material-max-step", true }
    };
    zero = ES::VXd::Zero(rest.size());
    writeState("deform0000.u", zero);
    writeState("substeps/frame0000-step0000.u", zero);
  }

  void TearDown() override
  {
    std::error_code error;
    fs::remove_all(directory, error);
  }

  void writeState(const std::string &filename, const ES::VXd &u)
  {
    ES::MXd state = ES::MXd::Zero(u.size(), 3);
    state.col(0) = u;
    ASSERT_EQ(ES::writeMatrix((output / filename).string().c_str(), state), 0);
  }

  Row row(int frame, int index, int count)
  {
    std::ostringstream filename;
    if (index + 1 == count)
      filename << "deform" << std::setfill('0') << std::setw(4) << frame << ".u";
    else
      filename << "substeps/frame" << std::setfill('0') << std::setw(4) << frame << "-step" << std::setw(4) << index << ".u";
    const double fraction = double(index + 1) / count;
    return { std::to_string(frame), std::to_string(index), std::to_string(count), scalar(fraction),
      scalar((frame + fraction) * 0.001), filename.str(), "2", "5", "1e-8", "1e-6", "0", "gradient_tolerance", "0" };
  }

  void writeManifest(const std::vector<Row> &rows)
  {
    std::ofstream stream(output / "newton-substeps.csv");
    stream << "frame,substep,substeps,fraction,physical_time,state_file,newton_iterations,pcg_iterations,gradient_inf,tolerance,status,stop,active_pairs\n";
    for (const auto &row : rows) {
      for (std::size_t i = 0; i < row.size(); ++i)
        stream << (i ? "," : "") << row[i];
      stream << '\n';
    }
  }

  void makeShell(const std::vector<double> &vertices, const std::vector<int> &triangles)
  {
    std::ofstream stream(directory / "shell.obj");
    stream << std::setprecision(17);
    for (std::size_t v = 0; v < vertices.size(); v += 3)
      stream << "v " << vertices[v] << ' ' << vertices[v + 1] << ' ' << vertices[v + 2] << '\n';
    for (std::size_t t = 0; t < triangles.size(); t += 3)
      stream << "f " << triangles[t] + 1 << ' ' << triangles[t + 1] + 1 << ' ' << triangles[t + 2] + 1 << '\n';
    stream.close();
    config.erase("tet-mesh");
    config["surface-mesh"] = (directory / "shell.obj").string();
    config["elastic-material"] = "koiter-stvk";
    rest = Eigen::Map<const ES::VXd>(vertices.data(), vertices.size());
    zero = ES::VXd::Zero(rest.size());
    writeState("deform0000.u", zero);
    writeState("substeps/frame0000-step0000.u", zero);
  }

  void makeCubic(bool invertedRestCorner = false)
  {
    const std::array<double, 24> vertices{ 0, 2, 0, 1, 2, 0, 1, 3, 0, 0, 3, 0,
      0, 2, 1, 1, 2, 1, 1, 3, 1, 0, 3, 1 };
    rest = Eigen::Map<const ES::VXd>(vertices.data(), vertices.size());
    if (invertedRestCorner)
      rest.segment<3>(18).array() -= 0.4;
    std::ofstream mesh(directory / "cube.veg");
    mesh << std::setprecision(17) << "*VERTICES\n8 3 0 0\n";
    for (int vertex = 0; vertex < 8; ++vertex)
      mesh << vertex + 1 << ' ' << rest[3 * vertex] << ' ' << rest[3 * vertex + 1] << ' ' << rest[3 * vertex + 2] << '\n';
    mesh << "*ELEMENTS\nCUBIC\n1 8 0\n1 1 2 3 4 5 6 7 8\n"
         << "*MATERIAL defaultMaterial\nENU, 1000, 100000, 0.3\n*REGION\nallElements, defaultMaterial\n";
    mesh.close();
    // An embedded interior triangle isolates the material criterion from
    // surface folding while still exercising the independent surface checker.
    std::ofstream(directory / "embedded.obj") << "v .1 2.1 .2\nv .3 2.1 .2\nv .1 2.3 .2\nf 1 2 3\n";
    config.erase("tet-mesh");
    config["cubic-mesh"] = (directory / "cube.veg").string();
    config["surface-mesh"] = (directory / "embedded.obj").string();
    zero = ES::VXd::Zero(24);
    writeState("deform0000.u", zero);
    writeState("substeps/frame0000-step0000.u", zero);
  }

  Json writeStaticPath(const std::vector<std::vector<ES::VXd>> &solves)
  {
    config["sim-type"] = "static";
    fs::create_directories(output / "newton-path");
    Json manifest = { { "schema_version", 1 }, { "path_type", "static_newton_iterates" },
      { "parameter", "accepted_newton_iteration" }, { "physical_time", false }, { "complete", true },
      { "num_dofs", zero.size() }, { "states", Json::array() }, { "solves", Json::array() } };
    for (std::size_t solve = 0; solve < solves.size(); ++solve) {
      const auto first = manifest["states"].size();
      for (std::size_t iteration = 0; iteration < solves[solve].size(); ++iteration) {
        const std::string filename = "newton-path/solve" + std::to_string(solve) + "-iteration" + std::to_string(iteration) + ".u";
        writeState(filename, solves[solve][iteration]);
        manifest["states"].push_back({ { "solve", solve }, { "iteration", iteration }, { "state_file", filename },
          { "energy", 100.0 - iteration }, { "gradient_inf", iteration + 1 == solves[solve].size() ? 1e-8 : 100.0 } });
      }
      manifest["solves"].push_back({ { "solve", solve }, { "first_state", first },
        { "last_state", manifest["states"].size() - 1 }, { "newton_iterations", solves[solve].size() - 1 },
        { "status", 0 }, { "stop", "gradient_tolerance" }, { "tolerance", 1e-6 }, { "gradient_inf", 1e-8 } });
    }
    writeState("deform0000.u", solves.back().back());
    writeStaticManifest(manifest);
    return manifest;
  }

  void writeStaticManifest(const Json &manifest)
  {
    std::ofstream(output / "static-newton-path.json") << manifest.dump(2);
  }

  std::pair<int, Json> check(const std::string &arguments = "", bool swept = true, bool expectReport = true)
  {
    fs::remove(directory / "report.json");
    std::ofstream(directory / "config.json") << config.dump(2);
    const std::string command = "OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 " +
      quote(PGO_TEST_VALIDATE_IPC_TRAJECTORY_BIN) + " " + quote(directory / "config.json") +
      (swept ? " --swept --report " : " --report ") + quote(directory / "report.json") + " --csv " + quote(directory / "frames.csv") +
      " " + arguments + " > " + quote(directory / "checker.log") + " 2>&1";
    const int status = std::system(command.c_str());
    std::ifstream report(directory / "report.json");
    if (!report && expectReport)
      throw std::runtime_error("Checker did not write a report; see " + (directory / "checker.log").string());
    if (!report)
      return { status, Json::object() };
    Json result;
    report >> result;
    return { status, std::move(result) };
  }

  fs::path directory, output;
  Json config;
  ES::VXd rest, zero;
};
}  // namespace

TEST_F(IPCTrajectoryValidatorTest, OldStatesWithoutManifestRemainSupported)
{
  const auto [status, report] = check();
  EXPECT_EQ(status, 0);
  EXPECT_TRUE(report["passed"].get<bool>());
  EXPECT_FALSE(report["physical_substeps"]["manifest_present"].get<bool>());
}

TEST_F(IPCTrajectoryValidatorTest, CubicGaussAcceptsHiddenCornerInversionAndRetainsDiagnostics)
{
  makeCubic();
  ES::VXd folded = zero;
  folded.segment<3>(18).setConstant(-0.4);
  writeState("deform0000.u", folded);
  writeManifest({ row(0, 0, 1) });
  auto [status, report] = check();
  EXPECT_NE(status, 0);
  EXPECT_EQ(report["element_validity_criterion"], "whole-cell");
  std::tie(status, report) = check("--cubic-validity gauss");
  ASSERT_EQ(status, 0) << report.dump(2);
  EXPECT_TRUE(report["passed"].get<bool>());
  EXPECT_EQ(report["element_validity_criterion"], "gauss");
  EXPECT_GT(report["minimum_element_determinant"].get<double>(), 1e-8);
  EXPECT_TRUE(report["minimum_is_certified_lower_bound"].get<bool>());
  EXPECT_FALSE(report["whole_cell_diagnostics"]["passed"].get<bool>());
  EXPECT_FALSE(report["whole_cell_diagnostics"]["required_for_validation"].get<bool>());
  EXPECT_GT(report["whole_cell_diagnostics"]["nonpositive_endpoint_element_evaluations"].get<int>(), 0);
  EXPECT_GT(report["whole_cell_diagnostics"]["failed_segment_checks"].get<int>(), 0);
  EXPECT_EQ(report["whole_cell_diagnostics"]["confirmed_inverted_element_count"], 1);
  EXPECT_EQ(report["whole_cell_diagnostics"]["confirmed_inverted_element_ids"], Json::array({ 0 }));
  EXPECT_EQ(report["whole_cell_diagnostics"]["frames_with_confirmed_inversion"], Json::array({ 0 }));
  EXPECT_EQ(report["frames"][0]["cubic_determinants"]["whole_cell_diagnostics"]["negative_corner_element_ids"], Json::array({ 0 }));
  EXPECT_TRUE(report["physical_substeps"]["passed"].get<bool>());
}

TEST_F(IPCTrajectoryValidatorTest, CubicGaussRejectsInteriorCollapseBetweenPositiveEndpoints)
{
  makeCubic();
  ES::VXd rotated = zero;
  for (int vertex = 0; vertex < 8; ++vertex) {
    rotated[3 * vertex] = 1 - 2 * rest[3 * vertex];
    rotated[3 * vertex + 1] = 5 - 2 * rest[3 * vertex + 1];
  }
  writeState("deform0000.u", rotated);
  writeManifest({ row(0, 0, 1) });
  const auto [status, report] = check("--cubic-validity gauss", false);
  EXPECT_NE(status, 0);
  const auto &geometry = report["frames"][0]["cubic_determinants"];
  EXPECT_TRUE(geometry["passed"].get<bool>());
  EXPECT_FALSE(geometry["linear_segment_passed"].get<bool>());
  EXPECT_GT(geometry["linear_segment_uncertified_count"].get<int>(), 0);
  EXPECT_FALSE(report["physical_substeps"]["passed"].get<bool>());
  EXPECT_FALSE(report["minimum_is_certified_lower_bound"].get<bool>());
  EXPECT_EQ(report["whole_cell_diagnostics"]["confirmed_inverted_element_count"], 0);
  EXPECT_EQ(report["whole_cell_diagnostics"]["uncertified_without_negative_witness_element_ids"], Json::array({ 0 }));
}

TEST_F(IPCTrajectoryValidatorTest, CubicInversionCountsFollowActualPhysicalPathWithoutDuplicateElementIDs)
{
  makeCubic();
  ES::VXd folded = zero;
  folded.segment<3>(18).setConstant(-0.4);
  writeState("substeps/frame0000-step0000.u", folded);
  writeManifest({ row(0, 0, 2), row(0, 1, 2) });
  const auto [status, report] = check("--cubic-validity gauss");
  ASSERT_EQ(status, 0) << report.dump(2);
  EXPECT_TRUE(report["coarse_chord_diagnostics"]["passed"].get<bool>());
  EXPECT_EQ(report["frames"][0]["cubic_determinants"]["whole_cell_diagnostics"]["confirmed_inverted_element_count"], 0);
  const auto &whole = report["whole_cell_diagnostics"];
  EXPECT_EQ(whole["confirmed_inverted_element_count"], 1);
  EXPECT_EQ(whole["confirmed_inverted_element_ids"], Json::array({ 0 }));
  EXPECT_EQ(whole["confirmed_inverted_endpoint_evaluations"], 1);
  EXPECT_EQ(whole["confirmed_inverted_segment_evaluations"], 2);
  EXPECT_EQ(whole["frames_with_confirmed_inversion"], Json::array({ 0 }));
}

TEST_F(IPCTrajectoryValidatorTest, CubicGaussChecksActualIntermediateStatesAndOriginalFloor)
{
  makeCubic();
  ES::VXd compressed = zero;
  for (int vertex = 0; vertex < 8; ++vertex)
    compressed[3 * vertex + 2] = (0.5e-8 - 1) * rest[3 * vertex + 2];
  writeState("substeps/frame0000-step0000.u", compressed);
  writeManifest({ row(0, 0, 2), row(0, 1, 2) });
  const auto [status, report] = check("--cubic-validity gauss");
  EXPECT_NE(status, 0);
  EXPECT_TRUE(report["coarse_chord_diagnostics"]["passed"].get<bool>());
  const auto &geometry = report["frames"][0]["physical_substeps"]["steps"][0]["cubic_determinants"];
  EXPECT_FALSE(geometry["passed"].get<bool>());
  EXPECT_GT(geometry["at_or_below_floor_or_nonfinite_count"].get<int>(), 0);
  EXPECT_EQ(geometry["nonpositive_or_nonfinite_count"], 0);
  EXPECT_FALSE(report["physical_substeps"]["passed"].get<bool>());
}

TEST_F(IPCTrajectoryValidatorTest, CubicGaussDoesNotWeakenSurfaceCCD)
{
  makeCubic();
  ES::VXd translated = zero;
  for (int vertex = 0; vertex < 8; ++vertex)
    translated[3 * vertex + 1] = -4;
  writeState("deform0000.u", translated);
  const auto [status, report] = check("--cubic-validity gauss");
  EXPECT_NE(status, 0);
  EXPECT_TRUE(report["frames"][0]["cubic_determinants"]["linear_segment_passed"].get<bool>());
  EXPECT_TRUE(report["frames"][0]["collision"]["passed"].get<bool>());
  EXPECT_FALSE(report["frames"][0]["deformable_motion_sweep"]["certified"].get<bool>());
}

TEST_F(IPCTrajectoryValidatorTest, CubicValidityMetadataConfigAndExplicitOverride)
{
  makeCubic();
  ES::VXd folded = zero;
  folded.segment<3>(18).setConstant(-0.4);
  writeState("deform0000.u", folded);
  Json metadata = { { "schema_version", 1 }, { "cubic_validity", "gauss" },
    { "material_max_step_enabled", true }, { "determinant_floor", 1e-8 }, { "element_type", "cubic" } };
  std::ofstream(output / "material-validity.json") << metadata.dump();
  auto [status, report] = check();
  ASSERT_EQ(status, 0);
  EXPECT_EQ(report["cubic_validity_source"], "material-validity.json");
  std::tie(status, report) = check("--cubic-validity whole-cell");
  EXPECT_NE(status, 0);
  EXPECT_EQ(report["element_validity_criterion"], "whole-cell");
  config["cuda-cubic-validity"] = "whole-cell";
  std::tie(status, report) = check("", true, false);
  EXPECT_NE(status, 0);
  EXPECT_TRUE(report.empty());
  std::tie(status, report) = check("--cubic-validity gauss");
  EXPECT_EQ(status, 0);
  EXPECT_TRUE(report["config_metadata_criterion_conflict_overridden"].get<bool>());
  std::tie(status, report) = check("--cubic-validity invalid", true, false);
  EXPECT_NE(status, 0);
  EXPECT_TRUE(report.empty());
  fs::remove(output / "material-validity.json");
  config["cuda-cubic-validity"] = "";
  std::tie(status, report) = check("", true, false);
  EXPECT_NE(status, 0);
  EXPECT_TRUE(report.empty());
  config["cuda-cubic-validity"] = "gauss";
  std::tie(status, report) = check();
  EXPECT_EQ(status, 0);
  EXPECT_EQ(report["cubic_validity_source"], "config");
}

TEST_F(IPCTrajectoryValidatorTest, CubicGaussRestValidityDoesNotRequireWholeCellRestValidity)
{
  makeCubic(true);
  const auto [status, report] = check("--cubic-validity gauss");
  ASSERT_EQ(status, 0) << report.dump(2);
  EXPECT_TRUE(report["passed"].get<bool>());
  EXPECT_EQ(report["whole_cell_diagnostics"]["uncertified_rest_elements"], 1);
  EXPECT_FALSE(report["whole_cell_diagnostics"]["passed"].get<bool>());
  EXPECT_TRUE(report["whole_cell_diagnostics"]["minimum_endpoint_determinant"].is_null());
  EXPECT_EQ(report["whole_cell_diagnostics"]["confirmed_inverted_element_count"], 0);
  const auto [strictStatus, strictReport] = check("--cubic-validity whole-cell", true, false);
  EXPECT_NE(strictStatus, 0);
  EXPECT_TRUE(strictReport.empty());
}

TEST_F(IPCTrajectoryValidatorTest, StaticNewtonDetourIsAuthoritativeWithoutPhysicalTime)
{
  auto translate = [&](double x, double y) {
    ES::VXd u = zero;
    for (int vertex = 0; vertex < u.size() / 3; ++vertex) {
      u[3 * vertex] = x;
      u[3 * vertex + 1] = y;
    }
    return u;
  };
  writeStaticPath({ { zero, translate(3, 0), translate(3, -1.3), translate(0, -1.3) } });
  auto terminal = row(0, 0, 1);
  terminal[4] = "0";
  writeManifest({ terminal });
  const auto [status, report] = check("", false);
  EXPECT_EQ(status, 0);
  EXPECT_EQ(report["trajectory_path"], "recorded_static_newton_iterates");
  EXPECT_TRUE(report["swept_checks_performed"].get<bool>());
  EXPECT_TRUE(report["static_newton_path"]["passed"].get<bool>());
  EXPECT_EQ(report["static_newton_path"]["states_checked"], 4);
  EXPECT_EQ(report["static_newton_path"]["segments_checked"], 3);
  EXPECT_FALSE(report["static_newton_path"]["physical_time"].get<bool>());
  EXPECT_FALSE(report["coarse_chord_diagnostics"]["passed"].get<bool>());
  EXPECT_FALSE(report["coarse_chord_diagnostics"]["required_for_validation"].get<bool>());
  EXPECT_TRUE(report["physical_substeps"]["passed"].get<bool>());
  EXPECT_FALSE(report["physical_substeps"]["authoritative_for_motion"].get<bool>());
  EXPECT_EQ(report["uncertified_sweep_segments"], 0);
  EXPECT_NEAR(report["minimum_element_determinant"].get<double>(), 1, 1e-10);
}

TEST_F(IPCTrajectoryValidatorTest, StaticNewtonShellRotationAllowsChangingNormals)
{
  makeShell({ 0, 1, 0, 1, 1, 0, 0, 2, 0 }, { 0, 1, 2 });
  config["ipc-external-objects"] = Json::array();
  auto rotate = [&](double angle) {
    ES::VXd u = zero;
    for (int v = 0; v < 3; ++v) {
      const double x = rest[3 * v];
      u[3 * v] = (std::cos(angle) - 1) * x;
      u[3 * v + 2] = std::sin(angle) * x;
    }
    return u;
  };
  writeStaticPath({ { zero, rotate(std::acos(-1.0) / 2), rotate(std::acos(-1.0)) } });
  const auto [status, report] = check();
  EXPECT_EQ(status, 0);
  EXPECT_TRUE(report["static_newton_path"]["passed"].get<bool>());
  EXPECT_FALSE(report["coarse_chord_diagnostics"]["passed"].get<bool>());
  EXPECT_GT(report["minimum_triangle_area_squared_ratio"].get<double>(), 0);
  EXPECT_TRUE(report["minimum_is_certified_lower_bound"].get<bool>());
}

TEST_F(IPCTrajectoryValidatorTest, StaticNewtonIntermediateCollisionCannotHideBehindClearEquilibrium)
{
  ES::VXd crossing = zero;
  for (int vertex = 0; vertex < crossing.size() / 3; ++vertex)
    crossing[3 * vertex + 1] = -0.6;
  writeStaticPath({ { zero, crossing, zero } });
  const auto [status, report] = check();
  EXPECT_NE(status, 0);
  EXPECT_TRUE(report["coarse_chord_diagnostics"]["passed"].get<bool>());
  EXPECT_FALSE(report["static_newton_path"]["passed"].get<bool>());
  EXPECT_GT(report["static_newton_path"]["external_edge_triangle_intersections"].get<int>(), 0);
}

TEST_F(IPCTrajectoryValidatorTest, StaticNewtonInteriorTriangleCollapseIsRejected)
{
  makeShell({ 0, 1, 0, 1, 1, 0, 0, 2, 0 }, { 0, 1, 2 });
  config["ipc-external-objects"] = Json::array();
  ES::VXd reversed = zero;
  reversed[3] = -2;
  writeStaticPath({ { zero, reversed, zero } });
  const auto [status, report] = check();
  EXPECT_NE(status, 0);
  EXPECT_FALSE(report["static_newton_path"]["passed"].get<bool>());
  EXPECT_TRUE(report["frames"][0]["static_newton_path"]["states"][1]["triangle_areas"]["passed"].get<bool>());
  EXPECT_FALSE(report["frames"][0]["static_newton_path"]["states"][1]["triangle_areas"]["linear_segment_passed"].get<bool>());
  EXPECT_FALSE(report["minimum_is_certified_lower_bound"].get<bool>());
}

TEST_F(IPCTrajectoryValidatorTest, StaticNewtonRejectsMalformedOrConflictingManifests)
{
  for (int error = 0; error < 15; ++error) {
    auto manifest = writeStaticPath({ { zero, zero } });
    if (error == 0)
      manifest["complete"] = false;
    if (error == 1)
      manifest["states"][1]["iteration"] = 0;
    if (error == 2)
      manifest["states"][1]["state_file"] = "missing.u";
    if (error == 3)
      manifest["states"][1]["state_file"] = "../outside.u";
    if (error == 4)
      manifest["solves"][0]["last_state"] = 0;
    if (error == 5)
      manifest["solves"][0]["gradient_inf"] = 1e-3;
    if (error == 6)
      manifest["solves"][0]["tolerance"] = 1e-3;
    if (error == 7)
      manifest["states"][1]["gradient_inf"] = 2e-8;
    if (error == 8)
      manifest["states"][0]["energy"] = nullptr;
    if (error == 9)
      manifest["physical_time"] = true;
    if (error == 10)
      manifest["states"][1]["iteration"] = 1.5;
    if (error == 11)
      writeState("deform0000.u", ES::VXd(zero.array() + 1e-8));
    if (error == 12)
      writeState("newton-path/solve0-iteration0.u", ES::VXd(zero.array() + 1e-8));
    if (error == 13)
      manifest["solves"][0]["newton_iterations"] = 201;
    if (error == 14) {
      auto terminal = row(0, 0, 1);
      terminal[4] = "0";
      terminal[8] = "2e-8";
      writeManifest({ terminal });
    }
    writeStaticManifest(manifest);
    const auto [status, report] = check();
    EXPECT_NE(status, 0) << error;
    EXPECT_FALSE(report["static_newton_path"]["passed"].get<bool>()) << error;
    EXPECT_FALSE(report["static_newton_path"]["errors"].empty()) << error;
    EXPECT_EQ(report["trajectory_path"], "recorded_static_newton_iterates");
  }
}

TEST_F(IPCTrajectoryValidatorTest, StaticNewtonLaggedSolvesRequireConvergedContinuousChains)
{
  config["contact-friction-coeff"] = 0.5;
  config["ipc-friction-iterations"] = 2;
  auto manifest = writeStaticPath({ { zero, zero }, { zero, zero } });
  auto [status, report] = check();
  EXPECT_EQ(status, 0);
  EXPECT_EQ(report["static_newton_path"]["recorded_solves"], 2);
  manifest["solves"][0]["status"] = 1;
  writeStaticManifest(manifest);
  std::tie(status, report) = check();
  EXPECT_NE(status, 0);
  manifest = writeStaticPath({ { zero, zero }, { ES::VXd(zero.array() + 1e-8), zero } });
  std::tie(status, report) = check();
  EXPECT_NE(status, 0);
  EXPECT_FALSE(report["static_newton_path"]["passed"].get<bool>());
}

TEST_F(IPCTrajectoryValidatorTest, ValidPhysicalSubstepsAndOriginalTimeInPartialRange)
{
  writeManifest({ row(0, 0, 2), row(0, 1, 2) });
  auto [status, report] = check();
  EXPECT_EQ(status, 0);
  EXPECT_TRUE(report["physical_substeps"]["passed"].get<bool>());
  EXPECT_EQ(report["physical_substeps"]["checked"], 2);
  writeManifest({ row(0, 0, 1) });
  std::tie(status, report) = check();
  EXPECT_EQ(status, 0);
  EXPECT_TRUE(report["frames"][0]["physical_substeps"]["steps"][0]["motion_sweep"]["required_for_validation"].get<bool>());
  EXPECT_EQ(report["frames"][0]["physical_substeps"]["steps"][0]["motion_sweep"]["scope"], "recorded physical substep");

  config["num-timestep"] = 3;
  writeState("deform0001.u", zero);
  writeState("deform0002.u", zero);
  writeState("substeps/frame0002-step0000.u", zero);
  writeManifest({ row(2, 0, 2), row(2, 1, 2) });
  std::tie(status, report) = check("--first-frame 2 --last-frame 2");
  EXPECT_EQ(status, 0);
  EXPECT_TRUE(report["passed"].get<bool>());
  EXPECT_DOUBLE_EQ(report["frames"][0]["physical_substeps"]["steps"][1]["physical_time"].get<double>(), 0.003);
  // The same manifest cannot claim complete coverage of the full trajectory.
  std::tie(status, report) = check();
  EXPECT_NE(status, 0);
  EXPECT_FALSE(report["physical_substeps"]["passed"].get<bool>());
}

TEST_F(IPCTrajectoryValidatorTest, SafeRecordedDetourPassesDespiteCollidingCoarseShortcut)
{
  // Translate the whole box around the finite floor's edge and underneath it.
  // Every recorded segment is clear and preserves all tetrahedron volumes;
  // the direct start/end shortcut crosses the floor.
  auto translate = [&](double x, double y) {
    ES::VXd u = zero;
    for (int vertex = 0; vertex < u.size() / 3; ++vertex) {
      u[3 * vertex] = x;
      u[3 * vertex + 1] = y;
    }
    return u;
  };
  writeState("substeps/frame0000-step0000.u", translate(3.0, 0.0));
  writeState("substeps/frame0000-step0001.u", translate(3.0, -1.3));
  writeState("deform0000.u", translate(0.0, -1.3));
  writeManifest({ row(0, 0, 3), row(0, 1, 3), row(0, 2, 3) });
  for (bool requestSweeps : { true, false }) {
    const auto [status, report] = check("", requestSweeps);
    EXPECT_EQ(status, 0);
    EXPECT_TRUE(report["passed"].get<bool>());
    EXPECT_EQ(report["trajectory_path"], "recorded_physical_substeps");
    EXPECT_EQ(report["swept_checks_requested"], requestSweeps);
    EXPECT_TRUE(report["swept_checks_performed"].get<bool>());
    EXPECT_TRUE(report["physical_substeps"]["passed"].get<bool>());
    EXPECT_EQ(report["physical_substeps"]["checked"], 3);
    EXPECT_EQ(report["uncertified_sweep_segments"], 0);
    EXPECT_EQ(report["sweep_midpoint_external_edge_triangle_intersections"], 0);
    EXPECT_EQ(report["actual_sweep_segments_checked"], 4);  // Prescribed motion plus three physical segments.
    EXPECT_FALSE(report["coarse_chord_diagnostics"]["passed"].get<bool>());
    EXPECT_FALSE(report["coarse_chord_diagnostics"]["required_for_validation"].get<bool>());
    EXPECT_GT(report["coarse_chord_diagnostics"]["uncertified_sweep_segments"].get<int>(), 0);
    EXPECT_FALSE(report["frames"][0]["deformable_motion_sweep"]["required_for_validation"].get<bool>());
    EXPECT_NEAR(report["minimum_tet_determinant"].get<double>(), 1.0, 1e-10);
  }
  // With no physical manifest, the same unsafe shortcut is the only recorded
  // path and must fail. An incomplete manifest must fail without fallback.
  fs::remove(output / "newton-substeps.csv");
  // Independently establish a real crossing on the direct shortcut. CCD's
  // first unresolved interval can lie just before the first touch, so its
  // localized midpoint alone need not show an intersection.
  writeState("deform0000.u", translate(0.0, -0.65));
  auto [status, report] = check();
  EXPECT_NE(status, 0);
  EXPECT_GT(report["frames"][0]["collision"]["external_edge_triangle_intersections"].get<int>(), 0);
  writeState("deform0000.u", translate(0.0, -1.3));
  std::tie(status, report) = check();
  EXPECT_NE(status, 0);
  EXPECT_FALSE(report["passed"].get<bool>());
  EXPECT_EQ(report["trajectory_path"], "legacy_macro_chords");
  EXPECT_TRUE(report["coarse_chord_diagnostics"]["required_for_validation"].get<bool>());
  writeManifest({ row(0, 0, 3), row(0, 2, 3) });
  std::tie(status, report) = check();
  EXPECT_NE(status, 0);
  EXPECT_FALSE(report["physical_substeps"]["passed"].get<bool>());
}

TEST_F(IPCTrajectoryValidatorTest, PrescribedMotionRemainsRequiredWithPhysicalManifest)
{
  config["num-timestep"] = 2;
  config["ipc-external-objects"][0]["movement"] = { 0.0, 2.0, 0.0 };
  writeState("deform0001.u", zero);
  writeManifest({ row(0, 0, 1), row(1, 0, 1) });
  const auto [status, report] = check("", false);
  EXPECT_NE(status, 0);
  EXPECT_FALSE(report["passed"].get<bool>());
  EXPECT_TRUE(report["physical_substeps"]["passed"].get<bool>());
  EXPECT_TRUE(report["frames"][1]["collision"]["passed"].get<bool>());
  EXPECT_FALSE(report["frames"][1]["prescribed_motion_sweep"]["certified"].get<bool>());
  EXPECT_GT(report["uncertified_sweep_segments"].get<int>(), 0);
}

TEST_F(IPCTrajectoryValidatorTest, IntermediateCollisionFailsWithClearMacroChord)
{
  ES::VXd crossing = zero;
  for (int i = 0; i < crossing.size() / 3; ++i)
    crossing[3 * i + 1] = -0.6;
  writeState("substeps/frame0000-step0000.u", crossing);
  writeManifest({ row(0, 0, 2), row(0, 1, 2) });
  const auto [status, report] = check();
  EXPECT_NE(status, 0);
  EXPECT_TRUE(report["frames"][0]["collision"]["passed"].get<bool>());
  EXPECT_TRUE(report["frames"][0]["deformable_motion_sweep"]["certified"].get<bool>());
  EXPECT_FALSE(report["physical_substeps"]["passed"].get<bool>());
  EXPECT_GT(report["physical_substeps"]["external_edge_triangle_intersections"].get<int>(), 0);
}

TEST_F(IPCTrajectoryValidatorTest, IntermediateAndInteriorInversionsAreRejected)
{
  for (bool rotation : { false, true }) {
    ES::VXd displacement = zero;
    for (int i = 0; i < rest.size() / 3; ++i) {
      displacement[3 * i] = (rotation ? -3 : -2) * rest[3 * i];
      if (rotation)
        displacement[3 * i + 1] = -1.5 * (rest[3 * i + 1] - 0.628135);
    }
    writeState("substeps/frame0000-step0000.u", displacement);
    writeManifest({ row(0, 0, 2), row(0, 1, 2) });
    const auto [status, report] = check();
    EXPECT_NE(status, 0);
    EXPECT_TRUE(report["frames"][0]["tet_determinants"]["passed"].get<bool>());
    EXPECT_FALSE(report["physical_substeps"]["passed"].get<bool>());
    const auto &determinants = report["frames"][0]["physical_substeps"]["steps"][0]["tet_determinants"];
    if (rotation)
      EXPECT_TRUE(determinants["passed"].get<bool>());
    else
      EXPECT_FALSE(determinants["passed"].get<bool>());
    EXPECT_GT(determinants["linear_segment_nonpositive_or_nonfinite_count"].get<int>(), 0);
  }
}

TEST_F(IPCTrajectoryValidatorTest, RejectsMalformedCoverageTimingAndSolverMetadata)
{
  for (int error = 0; error < 9; ++error) {
    std::vector<Row> rows{ row(0, 0, 2), row(0, 1, 2) };
    if (error == 0)
      rows.pop_back();
    if (error == 1)
      rows[1][1] = "0";
    if (error == 2)
      rows[0][3] = "0.4";
    if (error == 3)
      rows[1][4] = "0.002";
    if (error == 4)
      rows[0][8] = "1e-3";
    if (error == 5)
      rows[0][9] = "1e-5";
    if (error == 6)
      rows[0][10] = "1";
    if (error == 7)
      rows[0][8] = "nan";
    if (error == 8)
      rows[1][5] = "substeps/frame0000-step0000.u";
    writeManifest(rows);
    const auto [status, report] = check();
    EXPECT_NE(status, 0) << error;
    EXPECT_FALSE(report["passed"].get<bool>()) << error;
    EXPECT_FALSE(report["physical_substeps"]["passed"].get<bool>()) << error;
  }
}

TEST_F(IPCTrajectoryValidatorTest, RejectsMissingAndNonfinitePhysicalStates)
{
  writeManifest({ row(0, 0, 2), row(0, 1, 2) });
  fs::remove(output / "substeps/frame0000-step0000.u");
  auto [status, report] = check();
  EXPECT_NE(status, 0);
  EXPECT_FALSE(report["physical_substeps"]["passed"].get<bool>());
  ES::VXd invalid = zero;
  invalid[0] = std::numeric_limits<double>::quiet_NaN();
  writeState("substeps/frame0000-step0000.u", invalid);
  std::tie(status, report) = check();
  EXPECT_NE(status, 0);
  EXPECT_FALSE(report["physical_substeps"]["passed"].get<bool>());
}

TEST_F(IPCTrajectoryValidatorTest, ShellNormalRotatesThroughSafeRecordedSubsteps)
{
  makeShell({ 0, 2, 0, 1, 2, 0, 0, 3, 0 }, { 0, 1, 2 });
  ES::VXd quarter = zero, reversed = zero;
  quarter[7] = -1;
  quarter[8] = 1;
  reversed[7] = -2;
  writeState("substeps/frame0000-step0000.u", quarter);
  writeState("deform0000.u", reversed);
  writeManifest({ row(0, 0, 2), row(0, 1, 2) });
  const auto [status, report] = check("", false);
  EXPECT_EQ(status, 0);
  EXPECT_TRUE(report["passed"].get<bool>());
  EXPECT_EQ(report["element_type"], "shell");
  EXPECT_FALSE(report.contains("minimum_element_determinant"));
  EXPECT_FALSE(report.contains("minimum_tet_determinant"));
  EXPECT_TRUE(report["minimum_is_certified_lower_bound"].get<bool>());
  const double minimum = report["minimum_triangle_area_squared_ratio"].get<double>();
  EXPECT_GT(minimum, 0);
  EXPECT_LE(minimum, 0.5);
  EXPECT_TRUE(report["physical_substeps"]["passed"].get<bool>());
  EXPECT_EQ(report["physical_substeps"]["checked"], 2);
  EXPECT_FALSE(report["coarse_chord_diagnostics"]["passed"].get<bool>());
  EXPECT_FALSE(report["coarse_chord_diagnostics"]["required_for_validation"].get<bool>());
  EXPECT_FALSE(report["frames"][0]["triangle_areas"]["linear_segment_passed"].get<bool>());
  EXPECT_TRUE(report["frames"][0]["triangle_areas"]["passed"].get<bool>());
  std::ifstream actualCsv(directory / "frames.csv");
  std::string header;
  std::getline(actualCsv, header);
  EXPECT_NE(header.find("min_area_squared_ratio"), std::string::npos);
  EXPECT_EQ(header.find("det_F"), std::string::npos);
}

TEST_F(IPCTrajectoryValidatorTest, ShellInteriorCollapseFailsDespiteClearEndpoints)
{
  makeShell({ 0, 2, 0, 1, 2, 0, 0, 3, 0 }, { 0, 1, 2 });
  ES::VXd reversed = zero;
  reversed[7] = -2;
  writeState("deform0000.u", reversed);
  const auto [status, report] = check();
  EXPECT_NE(status, 0);
  EXPECT_FALSE(report["passed"].get<bool>());
  const auto &areas = report["frames"][0]["triangle_areas"];
  EXPECT_TRUE(areas["passed"].get<bool>());
  EXPECT_FALSE(areas["linear_segment_passed"].get<bool>());
  EXPECT_GT(areas["linear_segment_uncertified_count"].get<int>(), 0);
  EXPECT_FALSE(report["minimum_is_certified_lower_bound"].get<bool>());
  EXPECT_TRUE(report["frames"][0]["collision"]["passed"].get<bool>());
}

TEST_F(IPCTrajectoryValidatorTest, ShellCollapsedPhysicalEndpointFailsWithClearMacroChord)
{
  makeShell({ 0, 2, 0, 1, 2, 0, 0, 3, 0 }, { 0, 1, 2 });
  ES::VXd collapsed = zero;
  collapsed[7] = -1;
  writeState("substeps/frame0000-step0000.u", collapsed);
  writeManifest({ row(0, 0, 2), row(0, 1, 2) });
  const auto [status, report] = check();
  EXPECT_NE(status, 0);
  EXPECT_TRUE(report["coarse_chord_diagnostics"]["passed"].get<bool>());
  EXPECT_FALSE(report["physical_substeps"]["passed"].get<bool>());
  EXPECT_FALSE(report["frames"][0]["physical_substeps"]["steps"][0]["triangle_areas"]["passed"].get<bool>());
  EXPECT_FALSE(report["physical_substeps"]["endpoint_minimum_is_certified_lower_bound"].get<bool>());
}

TEST_F(IPCTrajectoryValidatorTest, ShellSelfCrossingBetweenClearEndpointsFails)
{
  makeShell({ 0, 1, 0, 1, 1, 0, 0, 2, 0, 0.2, 1.2, 1, 0.8, 1.2, 1, 0.2, 1.8, 1 }, { 0, 1, 2, 3, 4, 5 });
  ES::VXd u = zero;
  for (int v = 3; v < 6; ++v)
    u[3 * v + 2] = -2;
  writeState("deform0000.u", u);
  writeManifest({ row(0, 0, 1) });
  const auto [status, report] = check();
  EXPECT_NE(status, 0);
  EXPECT_TRUE(report["frames"][0]["collision"]["passed"].get<bool>());
  EXPECT_TRUE(report["frames"][0]["triangle_areas"]["linear_segment_passed"].get<bool>());
  EXPECT_FALSE(report["physical_substeps"]["passed"].get<bool>());
  EXPECT_GT(report["physical_substeps"]["uncertified_sweep_segments"].get<int>(), 0);
}

TEST_F(IPCTrajectoryValidatorTest, ShellObstacleCrossingBetweenClearEndpointsFails)
{
  makeShell({ 0, 1, 0, 1, 1, 0, 0, 1, 1 }, { 0, 1, 2 });
  ES::VXd u = zero;
  for (int v = 0; v < 3; ++v)
    u[3 * v + 1] = -2;
  writeState("deform0000.u", u);
  const auto [status, report] = check();
  EXPECT_NE(status, 0);
  EXPECT_TRUE(report["frames"][0]["collision"]["passed"].get<bool>());
  EXPECT_TRUE(report["frames"][0]["triangle_areas"]["linear_segment_passed"].get<bool>());
  EXPECT_FALSE(report["frames"][0]["deformable_motion_sweep"]["certified"].get<bool>());
}

TEST_F(IPCTrajectoryValidatorTest, ShellManifestCoverageAndPartialRangeRetainOriginalTime)
{
  makeShell({ 0, 1, 0, 1, 1, 0, 0, 2, 0 }, { 0, 1, 2 });
  config["num-timestep"] = 2;
  writeState("deform0001.u", zero);
  writeManifest({ row(1, 0, 1) });
  auto [status, report] = check("--first-frame 1 --last-frame 1");
  EXPECT_EQ(status, 0);
  EXPECT_TRUE(report["minimum_is_certified_lower_bound"].get<bool>());
  EXPECT_DOUBLE_EQ(report["frames"][0]["physical_substeps"]["steps"][0]["physical_time"].get<double>(), 0.002);
  std::tie(status, report) = check();
  EXPECT_NE(status, 0);
  EXPECT_FALSE(report["physical_substeps"]["passed"].get<bool>());
  auto bad = row(1, 0, 1);
  bad[8] = "0.001";
  writeManifest({ bad });
  std::tie(status, report) = check("--first-frame 1 --last-frame 1");
  EXPECT_NE(status, 0);
  EXPECT_FALSE(report["physical_substeps"]["passed"].get<bool>());
}
