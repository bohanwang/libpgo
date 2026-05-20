#include <gtest/gtest.h>

#include "EigenSupport.h"
#include "barycentricCoordinates.h"
#include "configFileJSON.h"
#include "deformationModelEnergy.h"
#include "initPredicates.h"
#include "pgoLogging.h"
#include "runIPCSimApp.h"
#include "runIPCSimCli.h"
#include "runIPCSimConfig.h"
#include "runIPCSimOutput.h"
#include "runIPCSimSetup.h"
#include "runSimCliLogging.h"
#include "runSimVolumeMeshIO.h"
#include "tetMesh.h"
#include "triMeshGeo.h"
#include "volumetricMesh.h"

#include <nlohmann/json.hpp>

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
namespace ES = pgo::EigenSupport;
namespace fs = std::filesystem;

constexpr const char *kShellExampleDir = LIBPGO_TEST_SHELL_EXAMPLE_DIR;
constexpr const char *kLegacyTetBoxDir = LIBPGO_TEST_LEGACY_TET_BOX_DIR;
constexpr const char *kLegacyCubicBoxDir = LIBPGO_TEST_LEGACY_CUBIC_BOX_DIR;
constexpr const char *kTetIPCExampleDir = LIBPGO_TEST_IPC_TET_EXAMPLE_DIR;
constexpr const char *kCubicIPCExampleDir = LIBPGO_TEST_IPC_CUBIC_EXAMPLE_DIR;

std::string quotePath(const fs::path &path)
{
  return "\"" + path.string() + "\"";
}

std::string shellExecutable(const fs::path &path)
{
#ifdef _WIN32
  return "call " + quotePath(path);
#else
  return quotePath(path);
#endif
}

std::string frameFilename(const char *prefix, int frame, const char *extension)
{
  std::ostringstream filename;
  filename << prefix << std::setfill('0') << std::setw(4) << frame << extension;
  return filename.str();
}

fs::path statePath(const fs::path &outputDir, int frame)
{
  return outputDir / "states" / frameFilename("deform", frame, ".u");
}

fs::path surfacePath(const fs::path &outputDir, int frame)
{
  return outputDir / "surface" / frameFilename("ret", frame, ".obj");
}

fs::path stressPath(const fs::path &outputDir, int frame)
{
  return outputDir / "stress" / frameFilename("von_mises", frame, ".json");
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
    const auto base = fs::temp_directory_path();
    for (int attempt = 0; attempt < 32; ++attempt) {
      const auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
      path_ = base / ("libpgo-runIPCSim-gtest-" + std::to_string(stamp) + "-" + std::to_string(std::rand()) + "-" + std::to_string(attempt));

      std::error_code ec;
      if (fs::create_directories(path_, ec))
        return;
    }

    throw std::runtime_error("Failed to create a unique temporary directory");
  }

  ~ScopedTempDir()
  {
    std::error_code ec;
    fs::remove_all(path_, ec);
  }

  const fs::path &path() const { return path_; }

private:
  fs::path path_;
};

void writeTextFile(const fs::path &path, const std::string &contents)
{
  fs::create_directories(path.parent_path());
  std::ofstream out(path);
  ASSERT_TRUE(out.is_open());
  out << contents;
}

std::string readTextFile(const fs::path &path)
{
  std::ifstream in(path);
  EXPECT_TRUE(in.is_open());
  return std::string((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
}

nlohmann::json readJsonFile(const fs::path &path)
{
  std::ifstream in(path);
  EXPECT_TRUE(in.is_open());
  nlohmann::json json;
  in >> json;
  return json;
}

std::string addBoolConfigField(std::string json, const std::string &name, bool value)
{
  const std::string marker = "\n}\n";
  const std::size_t pos = json.rfind(marker);
  if (pos == std::string::npos)
    throw std::runtime_error("Failed to add bool config field to test JSON.");

  json.insert(pos, ",\n  \"" + name + "\": " + (value ? "true" : "false"));
  return json;
}

std::string addSurfacePressureForceConfig(std::string json, bool enabled, double pressure = 1000.0, std::optional<int> rampSteps = 20,
  const std::string &centerField = "[0, 0, 0]")
{
  const std::string marker = "\n}\n";
  const std::size_t pos = json.rfind(marker);
  if (pos == std::string::npos)
    throw std::runtime_error("Failed to add surface pressure force config to test JSON.");

  std::ostringstream field;
  field << ",\n"
        << "  \"surface-pressure-force\": {\n"
        << "    \"enabled\": " << (enabled ? "true" : "false") << ",\n"
        << "    \"center\": " << centerField << ",\n"
        << "    \"pressure\": " << pressure;
  if (rampSteps.has_value()) {
    field << ",\n"
          << "    \"ramp-steps\": " << *rampSteps;
  }
  field << "\n"
        << "  }";
  json.insert(pos, field.str());
  return json;
}

std::string addSurfacePressureForceConfigWithoutCenter(std::string json, bool enabled, double pressure = 1000.0, std::optional<int> rampSteps = 20)
{
  const std::string marker = "\n}\n";
  const std::size_t pos = json.rfind(marker);
  if (pos == std::string::npos)
    throw std::runtime_error("Failed to add surface pressure force config to test JSON.");

  std::ostringstream field;
  field << ",\n"
        << "  \"surface-pressure-force\": {\n"
        << "    \"enabled\": " << (enabled ? "true" : "false") << ",\n"
        << "    \"pressure\": " << pressure;
  if (rampSteps.has_value()) {
    field << ",\n"
          << "    \"ramp-steps\": " << *rampSteps;
  }
  field << "\n"
        << "  }";
  json.insert(pos, field.str());
  return json;
}

void writeZeroShellRestartState(const fs::path &outputDir, int frame)
{
  pgo::Mesh::TriMeshGeo mesh;
  ASSERT_TRUE(mesh.load((fs::path(kShellExampleDir) / "shell.obj").string()));

  fs::create_directories(outputDir / "states");
  ES::MXd restartState = ES::MXd::Zero(mesh.numVertices() * 3, 3);
  ASSERT_EQ(ES::writeMatrix(statePath(outputDir, frame).string().c_str(), restartState), 0);
}

void appendFloorFields(std::ostringstream &json, bool useFloor,
  std::optional<std::string> floorAxis = std::nullopt,
  std::optional<double> floorHeight = std::nullopt,
  std::optional<double> floorKappa = std::nullopt)
{
  if (!useFloor && !floorAxis.has_value() && !floorHeight.has_value() && !floorKappa.has_value())
    return;

  json << ",\n"
       << "  \"floors\": [\n"
       << "    {\n";
  if (floorAxis.has_value())
    json << "      \"axis\": \"" << *floorAxis << "\"";
  if (floorHeight.has_value())
    json << (floorAxis.has_value() ? ",\n" : "")
         << "      \"height\": " << *floorHeight;
  if (floorKappa.has_value())
    json << (floorAxis.has_value() || floorHeight.has_value() ? ",\n" : "")
         << "      \"kappa\": " << *floorKappa;
  json << "\n"
       << "    }\n"
       << "  ]";
}

std::string addTopLevelJsonField(std::string json, const std::string &field)
{
  const std::string marker = "\n}\n";
  const std::size_t pos = json.rfind(marker);
  if (pos == std::string::npos)
    throw std::runtime_error("Failed to add config field to test JSON.");

  json.insert(pos, ",\n" + field);
  return json;
}

std::string makeStaticConfig(std::string json)
{
  const std::string dynamicToken = "\"sim-type\": \"dynamic\"";
  const std::string staticToken = "\"sim-type\": \"static\"";
  const std::size_t pos = json.find(dynamicToken);
  if (pos == std::string::npos)
    throw std::runtime_error("test config does not contain dynamic sim-type");

  json.replace(pos, dynamicToken.size(), staticToken);
  return json;
}

std::string addEmptyFloorsConfig(std::string json)
{
  return addTopLevelJsonField(std::move(json), "  \"floors\": []");
}

std::string addMovingUpperFloorConfig(std::string json)
{
  return addTopLevelJsonField(std::move(json),
    "  \"floors\": [\n"
    "    {\n"
    "      \"axis\": \"y\",\n"
    "      \"side\": \"upper\",\n"
    "      \"kappa\": 4000.0,\n"
    "      \"motion\": {\n"
    "        \"height-start\": 1.0,\n"
    "        \"height-end\": 0.75,\n"
    "        \"frame-start\": 0,\n"
    "        \"frame-end\": 1\n"
    "      }\n"
    "    }\n"
    "  ]");
}

fs::path runIPCSimBinaryPath()
{
  if (std::string(PGO_TEST_RUN_IPC_SIM_BIN).empty())
    return {};

  return fs::path(PGO_TEST_RUN_IPC_SIM_BIN);
}

fs::path tetIPCExampleDir()
{
  return fs::path(kTetIPCExampleDir);
}

fs::path cubicIPCExampleDir()
{
  return fs::path(kCubicIPCExampleDir);
}

fs::path cubicBoxSquashIPCExampleDir()
{
  return fs::path(__FILE__).parent_path().parent_path().parent_path().parent_path() / "examples" / "ipc" / "cubic" / "box-squash";
}

fs::path tetIPCConfigPath()
{
  return tetIPCExampleDir() / "box-ipc.json";
}

fs::path cubicIPCConfigPath()
{
  return cubicIPCExampleDir() / "box-ipc.json";
}

void initializeRunIPCSimTestEnvironment()
{
  static const bool initialized = []() {
    pgo::Logging::init();
    pgo::Mesh::initPredicates();
    return true;
  }();
  (void)initialized;
}

std::string makeShellIPCConfig(const fs::path &tempDir, int numTimesteps,
  bool includeIPCFields = true, bool ipcHeuristic = false,
  double ipcDhat = 0.002, double ipcKappa = 3000.0, int dumpInterval = 1,
  const std::string &logLevel = "info",
  bool useFloor = false,
  std::optional<std::string> floorAxis = std::nullopt,
  std::optional<double> floorHeight = std::nullopt,
  std::optional<double> floorKappa = std::nullopt,
  int solverMaxIter = 5)
{
  const fs::path shellDir = fs::path(kShellExampleDir);
  const fs::path outputDir = tempDir / "shell-output";

  std::ostringstream json;
  json << "{\n"
       << "  \"surface-mesh\": " << quotePath(shellDir / "shell.obj") << ",\n"
       << "  \"fixed-vertices\": [\n"
       << "    {\n"
       << "      \"filename\": " << quotePath(shellDir / "shell-fixed.txt") << ",\n"
       << "      \"movement\": [0, 0, 0],\n"
       << "      \"coeff\": 1e5\n"
       << "    }\n"
       << "  ],\n"
       << "  \"g\": [0, -9.81, 0],\n"
       << "  \"init-vel\": [0, 0, 0],\n"
       << "  \"init-disp\": [0, 0, 0],\n"
       << "  \"scale\": 1.0,\n"
       << "  \"timestep\": 0.001,\n"
       << "  \"num-timestep\": " << numTimesteps << ",\n"
       << "  \"damping-params\": [0, 0],\n"
       << "  \"sim-type\": \"dynamic\",\n"
       << "  \"solver-eps\": 1e-4,\n"
      << "  \"solver-max-iter\": " << solverMaxIter << ",\n"
       << "  \"elastic-material\": \"koiter-stvk\",\n"
       << "  \"loglevel\": \"" << logLevel << "\",\n"
       << "  \"dump-interval\": " << dumpInterval << ",\n"
       << "  \"output\": " << quotePath(outputDir);

  appendFloorFields(json, useFloor, floorAxis, floorHeight, floorKappa);

  if (ipcHeuristic) {
    json << ",\n"
         << "  \"ipc-heuristic\": true";
  }

  if (includeIPCFields) {
    json << ",\n"
         << "  \"ipc-dhat\": " << ipcDhat << ",\n"
         << "  \"ipc-kappa\": " << ipcKappa << "\n";
  }
  else {
    json << "\n";
  }

  json << "}\n";
  return json.str();
}

std::string makeShellIPCConfigWithIgnoredLegacyContactFields(const fs::path &tempDir, int numTimesteps)
{
  const fs::path shellDir = fs::path(kShellExampleDir);
  const fs::path outputDir = tempDir / "shell-output";

  std::ostringstream json;
  json << "{\n"
       << "  \"surface-mesh\": " << quotePath(shellDir / "shell.obj") << ",\n"
       << "  \"fixed-vertices\": [\n"
       << "    {\n"
       << "      \"filename\": " << quotePath(shellDir / "shell-fixed.txt") << ",\n"
       << "      \"movement\": [0, 0, 0],\n"
       << "      \"coeff\": 1e5\n"
       << "    }\n"
       << "  ],\n"
       << "  \"g\": [0, -9.81, 0],\n"
       << "  \"init-vel\": [0, 0, 0],\n"
       << "  \"init-disp\": [0, 0, 0],\n"
       << "  \"scale\": 1.0,\n"
       << "  \"timestep\": 0.001,\n"
       << "  \"num-timestep\": " << numTimesteps << ",\n"
       << "  \"damping-params\": [0, 0],\n"
       << "  \"sim-type\": \"dynamic\",\n"
       << "  \"solver-eps\": 1e-4,\n"
       << "  \"solver-max-iter\": 5,\n"
       << "  \"elastic-material\": \"koiter-stvk\",\n"
       << "  \"dump-interval\": 1,\n"
       << "  \"output\": " << quotePath(outputDir) << ",\n"
       << "  \"ipc-dhat\": 0.002,\n"
       << "  \"ipc-kappa\": 3000.0,\n"
       << "  \"contact-stiffness\": 123.0,\n"
       << "  \"contact-samples\": 7,\n"
       << "  \"contact-friction-coeff\": 0.4,\n"
       << "  \"contact-vel-eps\": 1e-5\n"
       << "}\n";
  return json.str();
}

std::string makeVolumeIPCConfig(const fs::path &exampleDir, const char *meshKey, const fs::path &outputDir, int numTimesteps,
  double scale = 1.0, bool includeIPCFields = true, bool ipcHeuristic = false, const std::string &material = "stable-neo",
  double ipcDhat = 0.002, double ipcKappa = 3000.0, int dumpInterval = 1,
  bool enableMaterialMaxStep = true, const std::string &logLevel = "info",
  bool useFloor = false,
  std::optional<std::string> floorAxis = std::nullopt,
  std::optional<double> floorHeight = std::nullopt,
  std::optional<double> floorKappa = std::nullopt,
  int solverMaxIter = 5)
{
  std::ostringstream json;
  json << "{\n"
       << "  \"" << meshKey << "\": " << quotePath(exampleDir / "box.veg") << ",\n"
       << "  \"surface-mesh\": " << quotePath(exampleDir / "box.obj") << ",\n"
       << "  \"fixed-vertices\": [\n"
       << "    {\n"
       << "      \"filename\": " << quotePath(exampleDir / "box-fixed.txt") << ",\n"
       << "      \"movement\": [0, 0, 0],\n"
       << "      \"coeff\": 1e5\n"
       << "    }\n"
       << "  ],\n"
       << "  \"g\": [0, -9.81, 0],\n"
       << "  \"init-vel\": [0, 0, 0],\n"
       << "  \"init-disp\": [0, 0, 0],\n"
       << "  \"scale\": " << scale << ",\n"
       << "  \"timestep\": 0.001,\n"
       << "  \"num-timestep\": " << numTimesteps << ",\n"
       << "  \"damping-params\": [0, 0],\n"
       << "  \"sim-type\": \"dynamic\",\n"
       << "  \"solver-eps\": 1e-4,\n"
       << "  \"solver-max-iter\": " << solverMaxIter << ",\n"
       << "  \"elastic-material\": \"" << material << "\",\n"
       << "  \"loglevel\": \"" << logLevel << "\",\n"
       << "  \"dump-interval\": " << dumpInterval << ",\n"
       << "  \"output\": " << quotePath(outputDir) << ",\n"
       << "  \"enable-material-max-step\": " << (enableMaterialMaxStep ? "true" : "false");

  appendFloorFields(json, useFloor, floorAxis, floorHeight, floorKappa);

  if (ipcHeuristic) {
    json << ",\n"
         << "  \"ipc-heuristic\": true";
  }

  if (includeIPCFields) {
    json << ",\n"
         << "  \"ipc-dhat\": " << ipcDhat << ",\n"
         << "  \"ipc-kappa\": " << ipcKappa << "\n";
  }
  else {
    json << "\n";
  }

  json << "}\n";
  return json.str();
}

std::string makeTetIPCConfig(const fs::path &tempDir, int numTimesteps,
  double scale = 1.0, bool includeIPCFields = true, bool ipcHeuristic = false, const std::string &material = "stable-neo",
  int dumpInterval = 1, const std::string &logLevel = "info",
  bool useFloor = false,
  std::optional<std::string> floorAxis = std::nullopt,
  std::optional<double> floorHeight = std::nullopt,
  std::optional<double> floorKappa = std::nullopt,
  int solverMaxIter = 5)
{
  return makeVolumeIPCConfig(tetIPCExampleDir(), "tet-mesh", tempDir / "tet-output", numTimesteps,
    scale, includeIPCFields, ipcHeuristic, material, 0.002, 3000.0, dumpInterval, true, logLevel, useFloor, floorAxis, floorHeight, floorKappa, solverMaxIter);
}

std::string makeCubicIPCConfig(const fs::path &tempDir, int numTimesteps,
  double scale = 1.0, bool includeIPCFields = true, bool ipcHeuristic = false, const std::string &material = "stable-neo",
  int dumpInterval = 1, const std::string &logLevel = "info",
  bool useFloor = false,
  std::optional<std::string> floorAxis = std::nullopt,
  std::optional<double> floorHeight = std::nullopt,
  std::optional<double> floorKappa = std::nullopt,
  int solverMaxIter = 5)
{
  return makeVolumeIPCConfig(cubicIPCExampleDir(), "cubic-mesh", tempDir / "cubic-output", numTimesteps,
    scale, includeIPCFields, ipcHeuristic, material, 0.002, 3000.0, dumpInterval, true, logLevel, useFloor, floorAxis, floorHeight, floorKappa, solverMaxIter);
}

std::string makeLegacyVolumeConfig(const fs::path &volumeMesh, const fs::path &surfaceMesh,
  const fs::path &outputDir, const char *meshKey, const char *material, int numTimesteps)
{
  std::ostringstream json;
  json << "{\n"
       << "  \"" << meshKey << "\": " << quotePath(volumeMesh) << ",\n"
       << "  \"surface-mesh\": " << quotePath(surfaceMesh) << ",\n"
       << "  \"fixed-vertices\": [],\n"
       << "  \"g\": [0, -9.81, 0],\n"
       << "  \"init-vel\": [0, 0, 0],\n"
       << "  \"init-disp\": [0, 0, 0],\n"
       << "  \"scale\": 1.0,\n"
       << "  \"timestep\": 0.001,\n"
       << "  \"num-timestep\": " << numTimesteps << ",\n"
       << "  \"damping-params\": [0, 0],\n"
       << "  \"sim-type\": \"dynamic\",\n"
       << "  \"contact-stiffness\": 1000,\n"
       << "  \"contact-sample\": 2,\n"
       << "  \"contact-friction-coeff\": 0.0,\n"
       << "  \"contact-vel-eps\": 1e-5,\n"
       << "  \"solver-eps\": 1e-4,\n"
       << "  \"solver-max-iter\": 5,\n"
       << "  \"elastic-material\": \"" << material << "\",\n"
       << "  \"dump-interval\": 1,\n"
       << "  \"output\": " << quotePath(outputDir) << "\n"
       << "}\n";
  return json.str();
}

std::string makeCubicSquashIPCConfig(const fs::path &tempDir, int numTimesteps, const std::string &logLevel = "info")
{
  const fs::path exampleDir = cubicBoxSquashIPCExampleDir();
  const fs::path outputDir = tempDir / "cubic-squash-output";

  std::ostringstream json;
  json << "{\n"
       << "  \"cubic-mesh\": " << quotePath(exampleDir / "box.veg") << ",\n"
       << "  \"surface-mesh\": " << quotePath(exampleDir / "box.obj") << ",\n"
       << "  \"fixed-vertices\": [\n"
       << "    {\n"
       << "      \"filename\": " << quotePath(exampleDir / "box-zmin-fixed.txt") << ",\n"
       << "      \"movement\": [0, 0, 0],\n"
       << "      \"coeff\": 1e5\n"
       << "    },\n"
       << "    {\n"
       << "      \"filename\": " << quotePath(exampleDir / "box-zmax-push.txt") << ",\n"
       << "      \"movement\": [0, 0, 500.55],\n"
       << "      \"coeff\": 1e5\n"
       << "    }\n"
       << "  ],\n"
       << "  \"g\": [0, 0, 0],\n"
       << "  \"init-vel\": [0, 0, 0],\n"
       << "  \"init-disp\": [0, 0, 0],\n"
       << "  \"scale\": 1.0,\n"
       << "  \"timestep\": 0.001,\n"
       << "  \"num-timestep\": " << numTimesteps << ",\n"
       << "  \"damping-params\": [0, 0],\n"
       << "  \"sim-type\": \"dynamic\",\n"
       << "  \"solver-eps\": 1e-5,\n"
       << "  \"solver-max-iter\": 20,\n"
       << "  \"elastic-material\": \"stable-neo\",\n"
       << "  \"loglevel\": \"" << logLevel << "\",\n"
       << "  \"dump-interval\": 1,\n"
       << "  \"output\": " << quotePath(outputDir) << ",\n"
       << "  \"ipc-dhat\": 0.002,\n"
       << "  \"ipc-kappa\": 3000.0,\n"
       << "  \"enable-material-max-step\": true\n"
       << "}\n";
  return json.str();
}

ES::SpMatD computeExpectedEmbeddingMatrix(const fs::path &configPath)
{
  initializeRunIPCSimTestEnvironment();

  pgo::ConfigFileJSON config;
  if (!config.open(configPath.string().c_str())) {
    throw std::runtime_error("Failed to open IPC example config: " + configPath.string());
  }

  const double scale = config.getDouble("scale", 1);
  const pgo::RunSim::VolumeMeshInputConfig meshConfig = pgo::RunSim::parseVolumeMeshInputConfig(config);
  std::unique_ptr<pgo::VolumetricMeshes::VolumetricMesh> volumetricMesh =
    pgo::RunSim::loadValidatedVolumeMesh(meshConfig, scale);

  const pgo::RunSim::ResolvedRunSimPaths paths = pgo::RunSim::resolveRunSimPaths(config);
  pgo::Mesh::TriMeshGeo surfaceMesh;
  if (!surfaceMesh.load(paths.surfaceMeshFilename)) {
    throw std::runtime_error("Failed to load surface mesh for IPC example config: " + configPath.string());
  }
  for (int vi = 0; vi < surfaceMesh.numVertices(); ++vi)
    surfaceMesh.pos(vi) *= scale;

  ES::VXd surfaceRestPositions(surfaceMesh.numVertices() * 3);
  for (int vi = 0; vi < surfaceMesh.numVertices(); ++vi)
    surfaceRestPositions.segment<3>(vi * 3) = surfaceMesh.pos(vi);

  pgo::InterpolationCoordinates::BarycentricCoordinates bc(
    surfaceMesh.numVertices(), surfaceRestPositions.data(), volumetricMesh.get());
  return bc.generateInterpolationMatrix();
}

void expectSparseMatrixNear(const ES::SpMatD &actual, const ES::SpMatD &expected, double tol = 1e-12)
{
  ASSERT_EQ(actual.rows(), expected.rows());
  ASSERT_EQ(actual.cols(), expected.cols());

  const ES::MXd actualDense(actual);
  const ES::MXd expectedDense(expected);
  const double maxDiff = (actualDense - expectedDense).cwiseAbs().maxCoeff();
  EXPECT_LE(maxDiff, tol);
}

std::string vec3Json(const ES::V3d &v)
{
  std::ostringstream out;
  out << "[" << std::setprecision(17) << v[0] << ", " << v[1] << ", " << v[2] << "]";
  return out.str();
}

}  // namespace

TEST(RunIPCSimConfigGTest, RuntimeConfigParsesRequiredAndOptionalFields)
{
  initializeRunIPCSimTestEnvironment();

  ScopedTempDir tempDir;
  const fs::path configPath = tempDir.path() / "shell-runtime-config.json";
  writeTextFile(configPath,
    addBoolConfigField(
      addBoolConfigField(makeShellIPCConfig(tempDir.path(), 3, true, false, 0.002, 3000.0, 2),
        "restart-from-u", true),
      "dump_deform_every_frame", true));

  pgo::ConfigFileJSON config;
  ASSERT_TRUE(config.open(configPath.string().c_str()));
  config.handle()["profiling"] = true;
  config.handle()["output-von-mises"] = false;

  const pgo::RunIPCSim::RunIPCSimRuntimeConfig runtime =
    pgo::RunIPCSim::parseRunIPCSimRuntimeConfig(config);

  EXPECT_EQ(runtime.numSimSteps, 3);
  EXPECT_EQ(runtime.frameGap, 2);
  EXPECT_TRUE(runtime.restartFromU);
  EXPECT_TRUE(runtime.dumpDeformEveryFrame);
  EXPECT_FALSE(runtime.outputVonMises);
  EXPECT_TRUE(runtime.enableProfiling);
  EXPECT_DOUBLE_EQ(runtime.scale, 1.0);
  EXPECT_DOUBLE_EQ(runtime.timestep, 0.001);
  EXPECT_EQ(runtime.outputFolder.filename(), "shell-output");
  EXPECT_NEAR(runtime.gravity[1], -9.81, 1e-12);
}

TEST(RunIPCSimConfigGTest, RuntimeConfigAcceptsStaticMode)
{
  initializeRunIPCSimTestEnvironment();

  ScopedTempDir tempDir;
  const fs::path configPath = tempDir.path() / "static-runtime-config.json";
  writeTextFile(configPath, makeStaticConfig(makeShellIPCConfig(tempDir.path(), 1)));

  pgo::ConfigFileJSON config;
  ASSERT_TRUE(config.open(configPath.string().c_str()));

  const pgo::RunIPCSim::RunIPCSimRuntimeConfig runtime =
    pgo::RunIPCSim::parseRunIPCSimRuntimeConfig(config);

  EXPECT_EQ(runtime.simulationMode, pgo::RunIPCSim::RunIPCSimSimulationMode::Static);
}

TEST(RunIPCSimOutputGTest, OutputPathsPreserveCurrentLayout)
{
  ScopedTempDir tempDir;
  const fs::path outputDir = tempDir.path() / "ipc-output";
  const pgo::RunIPCSim::RunIPCSimOutput output(outputDir);

  EXPECT_EQ(output.directories().root, outputDir);
  EXPECT_EQ(output.directories().states, outputDir / "states");
  EXPECT_EQ(output.directories().surface, outputDir / "surface");
  EXPECT_EQ(output.directories().stress, outputDir / "stress");
  EXPECT_EQ(output.logPath(), outputDir / "runIPCSim.log");
  EXPECT_EQ(output.statePath(7), outputDir / "states" / "deform0007.u");
  EXPECT_EQ(output.surfacePath(3), outputDir / "surface" / "ret0003.obj");
  EXPECT_EQ(output.stressPath(11), outputDir / "stress" / "von_mises0011.json");
}

TEST(RunIPCSimAppGTest, RunFromConfigNoTimestepsMatchesCliSuccess)
{
  initializeRunIPCSimTestEnvironment();

  ScopedTempDir tempDir;
  const fs::path configPath = tempDir.path() / "shell-runner-zero-step.json";
  const fs::path logPath = tempDir.path() / "shell-output" / "runIPCSim.log";

  writeTextFile(configPath, makeShellIPCConfig(tempDir.path(), 0));

  pgo::RunIPCSim::RunIPCSimOptions options;
  options.enableCliLog = true;

  EXPECT_EQ(pgo::RunIPCSim::runFromConfig(configPath, options), 0);
  ASSERT_TRUE(fs::exists(logPath));
  const std::string contents = readTextFile(logPath);
  EXPECT_NE(contents.find("runIPCSim phase2 shell IPC parameters:"), std::string::npos);
  EXPECT_NE(contents.find("max-step summary"), std::string::npos);
}

TEST(RunIPCSimAppGTest, RunFromConfigReturnsFailureForMissingRequiredIPCFields)
{
  initializeRunIPCSimTestEnvironment();

  ScopedTempDir tempDir;
  const fs::path configPath = tempDir.path() / "shell-runner-missing-ipc.json";
  writeTextFile(configPath, makeShellIPCConfig(tempDir.path(), 0, false));

  pgo::RunIPCSim::RunIPCSimOptions options;
  EXPECT_NE(pgo::RunIPCSim::runFromConfig(configPath, options), 0);
}

TEST(RunIPCSimCliGTest, LegacyFlagSelectsLegacyPenaltyBackend)
{
  const char *argv[] = { "runIPCSim", "--legacy", "scene.json" };
  const auto options = pgo::RunIPCSim::parseRunIPCSimCli(3, const_cast<char **>(argv));
  EXPECT_EQ(options.configPath, fs::path("scene.json"));
  EXPECT_EQ(options.runOptions.contactBackendKind, pgo::RunIPCSim::ContactBackendKind::LegacyPenalty);
}

TEST(RunIPCSimCliGTest, DefaultBackendIsIpc)
{
  const char *argv[] = { "runIPCSim", "scene.json" };
  const auto options = pgo::RunIPCSim::parseRunIPCSimCli(2, const_cast<char **>(argv));
  EXPECT_EQ(options.runOptions.contactBackendKind, pgo::RunIPCSim::ContactBackendKind::Ipc);
}

TEST(RunIPCSimLegacyGTest, LegacyTetConfigRunsOneStepAndWritesUnifiedOutput)
{
  initializeRunIPCSimTestEnvironment();

  ScopedTempDir tempDir;
  const fs::path configPath = tempDir.path() / "legacy-tet.json";
  const fs::path outputDir = tempDir.path() / "legacy-tet-output";
  writeTextFile(configPath, makeLegacyVolumeConfig(
    fs::path(kLegacyTetBoxDir) / "box.veg",
    fs::path(kLegacyTetBoxDir) / "box.obj",
    outputDir, "tet-mesh", "stable-neo", 1));

  pgo::RunIPCSim::RunIPCSimOptions options;
  options.contactBackendKind = pgo::RunIPCSim::ContactBackendKind::LegacyPenalty;
  EXPECT_EQ(pgo::RunIPCSim::runFromConfig(configPath, options), 0);
  EXPECT_TRUE(fs::exists(outputDir / "states" / "deform0000.u"));
  EXPECT_TRUE(fs::exists(outputDir / "surface" / "ret0000.obj"));
}

TEST(RunIPCSimLegacyGTest, LegacyCubicConfigRunsOneStepAndWritesUnifiedOutput)
{
  initializeRunIPCSimTestEnvironment();

  ScopedTempDir tempDir;
  const fs::path configPath = tempDir.path() / "legacy-cubic.json";
  const fs::path outputDir = tempDir.path() / "legacy-cubic-output";
  writeTextFile(configPath, makeLegacyVolumeConfig(
    fs::path(kLegacyCubicBoxDir) / "box.veg",
    fs::path(kLegacyCubicBoxDir) / "box.obj",
    outputDir, "cubic-mesh", "stable-neo", 1));

  pgo::RunIPCSim::RunIPCSimOptions options;
  options.contactBackendKind = pgo::RunIPCSim::ContactBackendKind::LegacyPenalty;
  EXPECT_EQ(pgo::RunIPCSim::runFromConfig(configPath, options), 0);
  EXPECT_TRUE(fs::exists(outputDir / "states" / "deform0000.u"));
  EXPECT_TRUE(fs::exists(outputDir / "surface" / "ret0000.obj"));
}

TEST(RunIPCSimLegacyGTest, LegacyShellConfigIsRejected)
{
  initializeRunIPCSimTestEnvironment();

  ScopedTempDir tempDir;
  const fs::path configPath = tempDir.path() / "legacy-shell.json";
  writeTextFile(configPath, makeShellIPCConfig(tempDir.path(), 0));

  pgo::RunIPCSim::RunIPCSimOptions options;
  options.contactBackendKind = pgo::RunIPCSim::ContactBackendKind::LegacyPenalty;
  EXPECT_EQ(pgo::RunIPCSim::runFromConfig(configPath, options), 1);
}

TEST(RunIPCSimStaticGTest, StaticTetIpcWritesUnifiedSurfaceAndState)
{
  initializeRunIPCSimTestEnvironment();

  ScopedTempDir tempDir;
  const fs::path configPath = tempDir.path() / "tet-static-ipc.json";
  const fs::path outputDir = tempDir.path() / "tet-output";
  writeTextFile(configPath, makeStaticConfig(makeTetIPCConfig(
    tempDir.path(), 1, 1.0, true, false, "stable-neo", 1, "info",
    false, std::nullopt, std::nullopt, std::nullopt, 1)));

  EXPECT_EQ(pgo::RunIPCSim::runFromConfig(configPath, {}), 0);
  EXPECT_TRUE(fs::exists(surfacePath(outputDir, 0)));
  EXPECT_TRUE(fs::exists(statePath(outputDir, 0)));
}

TEST(RunIPCSimStaticGTest, StaticShellIpcWritesUnifiedSurfaceAndState)
{
  initializeRunIPCSimTestEnvironment();

  ScopedTempDir tempDir;
  const fs::path configPath = tempDir.path() / "shell-static-ipc.json";
  const fs::path outputDir = tempDir.path() / "shell-output";
  writeTextFile(configPath, makeStaticConfig(makeShellIPCConfig(
    tempDir.path(), 1, true, false, 0.002, 3000.0, 1, "info",
    false, std::nullopt, std::nullopt, std::nullopt, 1)));

  EXPECT_EQ(pgo::RunIPCSim::runFromConfig(configPath, {}), 0);
  EXPECT_TRUE(fs::exists(surfacePath(outputDir, 0)));
  EXPECT_TRUE(fs::exists(statePath(outputDir, 0)));
}

TEST(RunIPCSimStaticGTest, StaticCubicIpcWithFloorWritesUnifiedSurfaceAndState)
{
  initializeRunIPCSimTestEnvironment();

  ScopedTempDir tempDir;
  const fs::path configPath = tempDir.path() / "cubic-static-floor-ipc.json";
  const fs::path outputDir = tempDir.path() / "cubic-output";
  writeTextFile(configPath, makeStaticConfig(makeCubicIPCConfig(
    tempDir.path(), 1, 1.0, true, false, "stable-neo", 1, "info",
    true, std::string("y"), -0.15, 4000.0, 1)));

  EXPECT_EQ(pgo::RunIPCSim::runFromConfig(configPath, {}), 0);
  EXPECT_TRUE(fs::exists(surfacePath(outputDir, 0)));
  EXPECT_TRUE(fs::exists(statePath(outputDir, 0)));
}

TEST(RunIPCSimStaticGTest, StaticVolumeWritesVonMisesWhenRequested)
{
  initializeRunIPCSimTestEnvironment();

  ScopedTempDir tempDir;
  const fs::path configPath = tempDir.path() / "tet-static-von-mises.json";
  const fs::path outputDir = tempDir.path() / "tet-output";
  writeTextFile(configPath,
    addBoolConfigField(makeStaticConfig(makeTetIPCConfig(
      tempDir.path(), 1, 1.0, true, false, "stable-neo", 1, "info",
      false, std::nullopt, std::nullopt, std::nullopt, 1)),
      "output-von-mises", true));

  EXPECT_EQ(pgo::RunIPCSim::runFromConfig(configPath, {}), 0);
  EXPECT_TRUE(fs::exists(stressPath(outputDir, 0)));
}

TEST(RunIPCSimStaticGTest, StaticLegacyTetWritesUnifiedSurfaceAndState)
{
  initializeRunIPCSimTestEnvironment();

  ScopedTempDir tempDir;
  const fs::path configPath = tempDir.path() / "tet-static-legacy.json";
  const fs::path outputDir = tempDir.path() / "legacy-tet-static-output";
  writeTextFile(configPath, makeStaticConfig(makeLegacyVolumeConfig(
    fs::path(kLegacyTetBoxDir) / "box.veg",
    fs::path(kLegacyTetBoxDir) / "box.obj",
    outputDir, "tet-mesh", "stable-neo", 1)));

  pgo::RunIPCSim::RunIPCSimOptions options;
  options.contactBackendKind = pgo::RunIPCSim::ContactBackendKind::LegacyPenalty;
  EXPECT_EQ(pgo::RunIPCSim::runFromConfig(configPath, options), 0);
  EXPECT_TRUE(fs::exists(surfacePath(outputDir, 0)));
  EXPECT_TRUE(fs::exists(statePath(outputDir, 0)));
}

TEST(RunIPCSimStaticGTest, StaticLegacyShellIsRejected)
{
  initializeRunIPCSimTestEnvironment();

  ScopedTempDir tempDir;
  const fs::path configPath = tempDir.path() / "shell-static-legacy.json";
  writeTextFile(configPath, makeStaticConfig(makeShellIPCConfig(tempDir.path(), 1)));

  pgo::RunIPCSim::RunIPCSimOptions options;
  options.contactBackendKind = pgo::RunIPCSim::ContactBackendKind::LegacyPenalty;
  EXPECT_NE(pgo::RunIPCSim::runFromConfig(configPath, options), 0);
}

TEST(RunIPCSimStaticGTest, StaticRejectsRestartFromU)
{
  initializeRunIPCSimTestEnvironment();

  ScopedTempDir tempDir;
  const fs::path configPath = tempDir.path() / "tet-static-restart.json";
  writeTextFile(configPath,
    addBoolConfigField(makeStaticConfig(makeTetIPCConfig(tempDir.path(), 1)),
      "restart-from-u", true));

  EXPECT_NE(pgo::RunIPCSim::runFromConfig(configPath, {}), 0);
}

TEST(RunIPCSimCliGTest, VolumeSetupRespectsDisabledMaterialMaxStepFlag)
{
  initializeRunIPCSimTestEnvironment();

  ScopedTempDir tempDir;
  const fs::path configPath = tempDir.path() / "tet-ipc-disable-material-max-step.json";
  writeTextFile(configPath, makeTetIPCConfig(tempDir.path(), 0, 1.0, true, false, "stable-neo", 1));

  pgo::ConfigFileJSON config;
  ASSERT_TRUE(config.open(configPath.string().c_str()));
  config.handle()["enable-material-max-step"] = false;

  const pgo::RunIPCSim::IpcSimulationContext context = pgo::RunIPCSim::buildVolumeIpcSimulation(config);
  ASSERT_NE(context.elasticEnergy, nullptr);
  EXPECT_FALSE(context.elasticEnergy->isMaterialMaxStepEnabled());
}

TEST(RunIPCSimCliGTest, LogFlagWritesCliOutputIntoOutputDirectory)
{
  const fs::path binary = runIPCSimBinaryPath();
  ASSERT_FALSE(binary.empty());
  ASSERT_TRUE(fs::exists(binary));

  ScopedTempDir tempDir;
  const fs::path configPath = tempDir.path() / "shell-ipc.json";
  const fs::path logPath = tempDir.path() / "shell-output" / "runIPCSim.log";

  writeTextFile(configPath, makeShellIPCConfig(tempDir.path(), 0));

  std::ostringstream command;
  command << shellExecutable(binary)
          << " --log "
          << quotePath(configPath);

  ASSERT_EQ(runCommand(command.str()), 0);
  ASSERT_TRUE(fs::exists(logPath));

  const std::string contents = readTextFile(logPath);
  EXPECT_NE(contents.find("ipc-dhat"), std::string::npos);
  EXPECT_NE(contents.find("ipc-kappa"), std::string::npos);
  EXPECT_NE(contents.find("eps_ee"), std::string::npos);
  EXPECT_NE(contents.find("slackness"), std::string::npos);
  EXPECT_NE(contents.find("max-step summary"), std::string::npos);
  EXPECT_NE(contents.find("minFeasibleAlphaThisSolve"), std::string::npos);
  EXPECT_NE(contents.find("minLineSearchAlphaThisSolve"), std::string::npos);
  EXPECT_NE(contents.find("minEffectiveAlphaThisSolve"), std::string::npos);
  EXPECT_EQ(contents.find("finalAlpha"), std::string::npos);
  EXPECT_EQ(contents.find("lastMaterialAlpha"), std::string::npos);
  EXPECT_EQ(contents.find("lastContactAlpha"), std::string::npos);
  EXPECT_EQ(contents.find("runIPCSim profiling summary:"), std::string::npos);
}

TEST(RunIPCSimCliGTest, ProfilingConfigWritesSummaryToOutputLog)
{
  const fs::path binary = runIPCSimBinaryPath();
  ASSERT_FALSE(binary.empty());
  ASSERT_TRUE(fs::exists(binary));

  ScopedTempDir tempDir;
  const fs::path configPath = tempDir.path() / "shell-ipc-profile.json";
  const fs::path logPath = tempDir.path() / "shell-output" / "runIPCSim.log";

  writeTextFile(configPath, addBoolConfigField(makeShellIPCConfig(tempDir.path(), 1), "profiling", true));

  std::ostringstream command;
  command << shellExecutable(binary)
          << " --log "
          << quotePath(configPath);

  ASSERT_EQ(runCommand(command.str()), 0);
  ASSERT_TRUE(fs::exists(logPath));

  const std::string contents = readTextFile(logPath);
  EXPECT_NE(contents.find("runIPCSim profiling summary:"), std::string::npos);
  EXPECT_NE(contents.find("profile name=contact.surface.pair_build.static"), std::string::npos);
  EXPECT_NE(contents.find("profileCounter name=contact.surface.pair_build.self_pt.hash_candidates"), std::string::npos);
  EXPECT_NE(contents.find("callCount="), std::string::npos);
  EXPECT_EQ(contents.find("SurfaceIPCCore active pairs:"), std::string::npos);
  EXPECT_EQ(contents.find("# nonzeros in Hessian:"), std::string::npos);
}

TEST(RunIPCSimCliGTest, DefaultRunClearsOutputAndDoesNotRestartFromDeformState)
{
  const fs::path binary = runIPCSimBinaryPath();
  ASSERT_FALSE(binary.empty());
  ASSERT_TRUE(fs::exists(binary));

  ScopedTempDir tempDir;
  const fs::path outputDir = tempDir.path() / "shell-output";
  const fs::path sentinelPath = outputDir / "old-sentinel.txt";
  const fs::path configPath = tempDir.path() / "shell-ipc-restart-off.json";
  const fs::path logPath = outputDir / "runIPCSim.log";

  writeTextFile(sentinelPath, "old output");
  writeZeroShellRestartState(outputDir, 0);
  writeTextFile(configPath, makeShellIPCConfig(tempDir.path(), 1));

  std::ostringstream command;
  command << shellExecutable(binary)
          << " --log "
          << quotePath(configPath);

  ASSERT_EQ(runCommand(command.str()), 0);
  EXPECT_FALSE(fs::exists(sentinelPath));
  ASSERT_TRUE(fs::exists(logPath));

  const std::string contents = readTextFile(logPath);
  EXPECT_NE(contents.find("restart-from-u=false; clearing output folder"), std::string::npos);
  EXPECT_NE(contents.find("Starting from frame 0."), std::string::npos);
  EXPECT_EQ(contents.find("Restarting from frame"), std::string::npos);
}

TEST(RunIPCSimCliGTest, RestartFromUTrueKeepsExistingDeformState)
{
  const fs::path binary = runIPCSimBinaryPath();
  ASSERT_FALSE(binary.empty());
  ASSERT_TRUE(fs::exists(binary));

  ScopedTempDir tempDir;
  const fs::path outputDir = tempDir.path() / "shell-output";
  const fs::path configPath = tempDir.path() / "shell-ipc-restart-on.json";
  const fs::path logPath = outputDir / "runIPCSim.log";

  writeZeroShellRestartState(outputDir, 0);
  writeTextFile(configPath, addBoolConfigField(makeShellIPCConfig(tempDir.path(), 2), "restart-from-u", true));

  std::ostringstream command;
  command << shellExecutable(binary)
          << " --log "
          << quotePath(configPath);

  ASSERT_EQ(runCommand(command.str()), 0);
  ASSERT_TRUE(fs::exists(logPath));
  EXPECT_TRUE(fs::exists(statePath(outputDir, 1)));

  const std::string contents = readTextFile(logPath);
  EXPECT_NE(contents.find("Restarting from frame 0"), std::string::npos);
  EXPECT_EQ(contents.find("restart-from-u=false; clearing output folder"), std::string::npos);
}

TEST(RunIPCSimCliGTest, FloorEnabledLogPrintsFloorParameters)
{
  const fs::path binary = runIPCSimBinaryPath();
  ASSERT_FALSE(binary.empty());
  ASSERT_TRUE(fs::exists(binary));

  ScopedTempDir tempDir;
  const fs::path configPath = tempDir.path() / "shell-ipc-floor.json";
  const fs::path logPath = tempDir.path() / "shell-output" / "runIPCSim.log";

  writeTextFile(configPath, makeShellIPCConfig(tempDir.path(), 1, true, false, 0.002, 3000.0, 1, "info", true, "y", -0.15, 4321.0, 20));

  std::ostringstream command;
  command << shellExecutable(binary)
          << " --log "
          << quotePath(configPath);

  ASSERT_EQ(runCommand(command.str()), 0);
  ASSERT_TRUE(fs::exists(logPath));

  const std::string contents = readTextFile(logPath);
  EXPECT_NE(contents.find("floors=1"), std::string::npos);
  EXPECT_NE(contents.find("floor[0].axis=y"), std::string::npos);
  EXPECT_NE(contents.find("floor[0].side=lower"), std::string::npos);
  EXPECT_NE(contents.find("floor[0].height=-0.15"), std::string::npos);
  EXPECT_NE(contents.find("floor[0].kappa=4321"), std::string::npos);
}

TEST(RunIPCSimCliGTest, DebugLogLevelPrintsFullMaxStepSummary)
{
  const fs::path binary = runIPCSimBinaryPath();
  ASSERT_FALSE(binary.empty());
  ASSERT_TRUE(fs::exists(binary));

  ScopedTempDir tempDir;
  const fs::path configPath = tempDir.path() / "shell-ipc-debug.json";
  const fs::path logPath = tempDir.path() / "shell-output" / "runIPCSim.log";

  writeTextFile(configPath, makeShellIPCConfig(tempDir.path(), 1, true, false, 0.002, 3000.0, 1, "debug"));

  std::ostringstream command;
  command << shellExecutable(binary)
          << " --log "
          << quotePath(configPath);

  ASSERT_EQ(runCommand(command.str()), 0);
  ASSERT_TRUE(fs::exists(logPath));

  const std::string contents = readTextFile(logPath);
  EXPECT_NE(contents.find("max-step summary"), std::string::npos);
  EXPECT_NE(contents.find("minMaterialFeasibleAlphaThisSolve"), std::string::npos);
  EXPECT_NE(contents.find("minContactFeasibleAlphaThisSolve"), std::string::npos);
  EXPECT_NE(contents.find("minFeasibleAlphaThisSolve"), std::string::npos);
  EXPECT_NE(contents.find("minLineSearchAlphaThisSolve"), std::string::npos);
  EXPECT_NE(contents.find("minEffectiveAlphaThisSolve"), std::string::npos);
  EXPECT_NE(contents.find("rawStepMaxNorm="), std::string::npos);
  EXPECT_NE(contents.find("feasibleAlpha="), std::string::npos);
  EXPECT_NE(contents.find("lineSearchAlpha="), std::string::npos);
  EXPECT_NE(contents.find("effectiveAlpha="), std::string::npos);
  EXPECT_NE(contents.find("acceptedStepMaxNorm="), std::string::npos);
  EXPECT_NE(contents.find("acceptedEnergy="), std::string::npos);
  EXPECT_EQ(contents.find("finalAlpha"), std::string::npos);
  EXPECT_EQ(contents.find("lastMaterialAlpha"), std::string::npos);
  EXPECT_EQ(contents.find("lastContactAlpha"), std::string::npos);
}

TEST(RunIPCSimCliGTest, DebugLogLevelPrintsClampedFeasibleAlphaBreakdownForCubicIPC)
{
  const fs::path binary = runIPCSimBinaryPath();
  ASSERT_FALSE(binary.empty());
  ASSERT_TRUE(fs::exists(binary));

  ScopedTempDir tempDir;
  const fs::path configPath = tempDir.path() / "cubic-squash-ipc-debug.json";
  const fs::path logPath = tempDir.path() / "cubic-squash-output" / "runIPCSim.log";

  writeTextFile(configPath, makeCubicSquashIPCConfig(tempDir.path(), 3, "debug"));

  std::ostringstream command;
  command << shellExecutable(binary)
          << " --log "
          << quotePath(configPath);

  ASSERT_EQ(runCommand(command.str()), 0);
  ASSERT_TRUE(fs::exists(logPath));

  const std::string contents = readTextFile(logPath);
  EXPECT_NE(contents.find("feasible alpha clamped: material:"), std::string::npos);
  EXPECT_NE(contents.find("contact:"), std::string::npos);
  EXPECT_NE(contents.find("accepted=true"), std::string::npos);
}

TEST(RunIPCSimCliGTest, WarnLogLevelSuppressesMaxStepSummary)
{
  const fs::path binary = runIPCSimBinaryPath();
  ASSERT_FALSE(binary.empty());
  ASSERT_TRUE(fs::exists(binary));

  ScopedTempDir tempDir;
  const fs::path configPath = tempDir.path() / "shell-ipc-warn.json";
  const fs::path logPath = tempDir.path() / "shell-output" / "runIPCSim.log";

  writeTextFile(configPath, makeShellIPCConfig(tempDir.path(), 1, true, false, 0.002, 3000.0, 1, "warn"));

  std::ostringstream command;
  command << shellExecutable(binary)
          << " --log "
          << quotePath(configPath);

  ASSERT_EQ(runCommand(command.str()), 0);
  ASSERT_TRUE(fs::exists(logPath));

  const std::string contents = readTextFile(logPath);
  EXPECT_EQ(contents.find("max-step summary"), std::string::npos);
}

TEST(RunIPCSimCliGTest, MissingIPCDhatOrKappaFails)
{
  const fs::path binary = runIPCSimBinaryPath();
  ASSERT_FALSE(binary.empty());
  ASSERT_TRUE(fs::exists(binary));

  ScopedTempDir tempDir;
  const fs::path configPath = tempDir.path() / "shell-ipc-missing.json";

  writeTextFile(configPath, makeShellIPCConfig(tempDir.path(), 0, false));

  std::ostringstream command;
  command << shellExecutable(binary)
          << " "
          << quotePath(configPath);

  ASSERT_NE(runCommand(command.str()), 0);
}

TEST(RunIPCSimCliGTest, HeuristicAllowsMissingIPCDhatAndKappa)
{
  const fs::path binary = runIPCSimBinaryPath();
  ASSERT_FALSE(binary.empty());
  ASSERT_TRUE(fs::exists(binary));

  ScopedTempDir tempDir;
  const fs::path configPath = tempDir.path() / "shell-ipc-heuristic.json";

  writeTextFile(configPath, makeShellIPCConfig(tempDir.path(), 0, false, true));

  std::ostringstream command;
  command << shellExecutable(binary)
          << " "
          << quotePath(configPath);

  ASSERT_EQ(runCommand(command.str()), 0);
}

TEST(RunIPCSimCliGTest, OneTimestepShellSmokeSucceeds)
{
  const fs::path binary = runIPCSimBinaryPath();
  ASSERT_FALSE(binary.empty());
  ASSERT_TRUE(fs::exists(binary));

  ScopedTempDir tempDir;
  const fs::path configPath = tempDir.path() / "shell-ipc-step.json";

  writeTextFile(configPath, makeShellIPCConfig(tempDir.path(), 1));

  std::ostringstream command;
  command << shellExecutable(binary)
          << " "
          << quotePath(configPath);

  ASSERT_EQ(runCommand(command.str()), 0);
  ASSERT_TRUE(fs::exists(tempDir.path() / "shell-output"));
}

TEST(RunIPCSimCliGTest, OneTimestepShellFloorSmokeSucceeds)
{
  const fs::path binary = runIPCSimBinaryPath();
  ASSERT_FALSE(binary.empty());
  ASSERT_TRUE(fs::exists(binary));

  ScopedTempDir tempDir;
  const fs::path configPath = tempDir.path() / "shell-ipc-floor-step.json";

  writeTextFile(configPath, makeShellIPCConfig(tempDir.path(), 1, true, false, 0.002, 3000.0, 1, "info", true, "y", -0.15, 4000.0, 20));

  std::ostringstream command;
  command << shellExecutable(binary)
          << " "
          << quotePath(configPath);

  ASSERT_EQ(runCommand(command.str()), 0);
  const fs::path outputDir = tempDir.path() / "shell-output";
  EXPECT_TRUE(fs::exists(statePath(outputDir, 0)));
  EXPECT_TRUE(fs::exists(surfacePath(outputDir, 0)));
  EXPECT_FALSE(fs::exists(outputDir / "deform0000.u"));
  EXPECT_FALSE(fs::exists(outputDir / "ret0000.obj"));
}

TEST(RunIPCSimCliGTest, DeformStateDefaultsToDumpInterval)
{
  const fs::path binary = runIPCSimBinaryPath();
  ASSERT_FALSE(binary.empty());
  ASSERT_TRUE(fs::exists(binary));

  ScopedTempDir tempDir;
  const fs::path configPath = tempDir.path() / "shell-ipc-deform-dump-interval.json";

  writeTextFile(configPath, makeShellIPCConfig(tempDir.path(), 2, true, false, 0.002, 3000.0, 10));

  std::ostringstream command;
  command << shellExecutable(binary)
          << " "
          << quotePath(configPath);

  ASSERT_EQ(runCommand(command.str()), 0);
  const fs::path outputDir = tempDir.path() / "shell-output";
  EXPECT_TRUE(fs::exists(statePath(outputDir, 0)));
  EXPECT_FALSE(fs::exists(statePath(outputDir, 1)));
  EXPECT_TRUE(fs::exists(surfacePath(outputDir, 0)));
  EXPECT_FALSE(fs::exists(surfacePath(outputDir, 1)));
  EXPECT_FALSE(fs::exists(outputDir / "deform0000.u"));
  EXPECT_FALSE(fs::exists(outputDir / "ret0000.obj"));
}

TEST(RunIPCSimCliGTest, DumpDeformEveryFrameWritesEveryTimestep)
{
  const fs::path binary = runIPCSimBinaryPath();
  ASSERT_FALSE(binary.empty());
  ASSERT_TRUE(fs::exists(binary));

  ScopedTempDir tempDir;
  const fs::path configPath = tempDir.path() / "shell-ipc-deform-every-step.json";

  writeTextFile(configPath,
    addBoolConfigField(makeShellIPCConfig(tempDir.path(), 2, true, false, 0.002, 3000.0, 10),
      "dump_deform_every_frame", true));

  std::ostringstream command;
  command << shellExecutable(binary)
          << " "
          << quotePath(configPath);

  ASSERT_EQ(runCommand(command.str()), 0);
  const fs::path outputDir = tempDir.path() / "shell-output";
  EXPECT_TRUE(fs::exists(statePath(outputDir, 0)));
  EXPECT_TRUE(fs::exists(statePath(outputDir, 1)));
  EXPECT_TRUE(fs::exists(surfacePath(outputDir, 0)));
  EXPECT_FALSE(fs::exists(surfacePath(outputDir, 1)));
  EXPECT_FALSE(fs::exists(outputDir / "deform0000.u"));
  EXPECT_FALSE(fs::exists(outputDir / "ret0000.obj"));
}

TEST(RunIPCSimCliGTest, LegacyContactFieldsAreIgnoredWhenPresent)
{
  const fs::path binary = runIPCSimBinaryPath();
  ASSERT_FALSE(binary.empty());
  ASSERT_TRUE(fs::exists(binary));

  ScopedTempDir tempDir;
  const fs::path configPath = tempDir.path() / "shell-ipc-with-legacy-contact.json";

  writeTextFile(configPath, makeShellIPCConfigWithIgnoredLegacyContactFields(tempDir.path(), 0));

  std::ostringstream command;
  command << shellExecutable(binary)
          << " "
          << quotePath(configPath);

  ASSERT_EQ(runCommand(command.str()), 0);
}

TEST(RunIPCSimCliGTest, HeuristicOverridesExplicitIPCFieldsInLog)
{
  const fs::path binary = runIPCSimBinaryPath();
  ASSERT_FALSE(binary.empty());
  ASSERT_TRUE(fs::exists(binary));

  ScopedTempDir tempDir;
  const fs::path configPath = tempDir.path() / "shell-ipc-heuristic-override.json";
  const fs::path logPath = tempDir.path() / "shell-output" / "runIPCSim.log";

  writeTextFile(configPath, makeShellIPCConfig(tempDir.path(), 0, true, true, 9.0, 42.0));

  std::ostringstream command;
  command << shellExecutable(binary)
          << " --log "
          << quotePath(configPath);

  ASSERT_EQ(runCommand(command.str()), 0);
  ASSERT_TRUE(fs::exists(logPath));

  const std::string contents = readTextFile(logPath);
  EXPECT_NE(contents.find("ipc-heuristic=true"), std::string::npos);
  EXPECT_NE(contents.find("source=heuristic"), std::string::npos);
  EXPECT_NE(contents.find("ipc-dhat=0.00141421"), std::string::npos);
  EXPECT_NE(contents.find("ipc-kappa=3000"), std::string::npos);
  EXPECT_EQ(contents.find("ipc-dhat=9"), std::string::npos);
  EXPECT_EQ(contents.find("ipc-kappa=42"), std::string::npos);
}

TEST(RunIPCSimCliGTest, TetZeroTimestepSmokeCreatesNoOutputs)
{
  const fs::path binary = runIPCSimBinaryPath();
  ASSERT_FALSE(binary.empty());
  ASSERT_TRUE(fs::exists(binary));

  ScopedTempDir tempDir;
  const fs::path configPath = tempDir.path() / "tet-ipc-zero.json";

  writeTextFile(configPath, makeTetIPCConfig(tempDir.path(), 0));

  std::ostringstream command;
  command << shellExecutable(binary)
          << " "
          << quotePath(configPath);

  ASSERT_EQ(runCommand(command.str()), 0);
  EXPECT_TRUE(fs::exists(tempDir.path() / "tet-output"));
  EXPECT_TRUE(fs::exists(tempDir.path() / "tet-output" / "states"));
  EXPECT_TRUE(fs::exists(tempDir.path() / "tet-output" / "surface"));
  EXPECT_TRUE(fs::exists(tempDir.path() / "tet-output" / "stress"));
  EXPECT_FALSE(fs::exists(statePath(tempDir.path() / "tet-output", 0)));
  EXPECT_FALSE(fs::exists(surfacePath(tempDir.path() / "tet-output", 0)));
  EXPECT_FALSE(fs::exists(stressPath(tempDir.path() / "tet-output", 0)));
  EXPECT_FALSE(fs::exists(tempDir.path() / "tet-output" / "deform0000.u"));
  EXPECT_FALSE(fs::exists(tempDir.path() / "tet-output" / "ret0000.obj"));
}

TEST(RunIPCSimCliGTest, TetOneTimestepSmokeWritesDeformAndRet)
{
  const fs::path binary = runIPCSimBinaryPath();
  ASSERT_FALSE(binary.empty());
  ASSERT_TRUE(fs::exists(binary));

  ScopedTempDir tempDir;
  const fs::path configPath = tempDir.path() / "tet-ipc-one-step.json";

  writeTextFile(configPath, makeTetIPCConfig(tempDir.path(), 1));

  std::ostringstream command;
  command << shellExecutable(binary)
          << " "
          << quotePath(configPath);

  ASSERT_EQ(runCommand(command.str()), 0);
  const fs::path outputDir = tempDir.path() / "tet-output";
  EXPECT_TRUE(fs::exists(statePath(outputDir, 0)));
  EXPECT_TRUE(fs::exists(surfacePath(outputDir, 0)));
  EXPECT_FALSE(fs::exists(outputDir / "deform0000.u"));
  EXPECT_FALSE(fs::exists(outputDir / "ret0000.obj"));
}

TEST(RunIPCSimCliGTest, TetVonMisesOutputWritesElementStressJson)
{
  const fs::path binary = runIPCSimBinaryPath();
  ASSERT_FALSE(binary.empty());
  ASSERT_TRUE(fs::exists(binary));

  ScopedTempDir tempDir;
  const fs::path configPath = tempDir.path() / "tet-ipc-von-mises.json";
  const fs::path outputDir = tempDir.path() / "tet-output";

  writeTextFile(configPath, addBoolConfigField(makeTetIPCConfig(tempDir.path(), 1), "output-von-mises", true));

  std::ostringstream command;
  command << shellExecutable(binary)
          << " "
          << quotePath(configPath);

  ASSERT_EQ(runCommand(command.str()), 0);
  ASSERT_TRUE(fs::exists(stressPath(outputDir, 0)));

  const nlohmann::json stressJson = readJsonFile(stressPath(outputDir, 0));
  EXPECT_EQ(stressJson.at("frame").get<int>(), 0);
  EXPECT_DOUBLE_EQ(stressJson.at("time").get<double>(), 0.0);
  EXPECT_EQ(stressJson.at("stress_type").get<std::string>(), "von_mises");
  EXPECT_EQ(stressJson.at("location").get<std::string>(), "tet_element");

  const pgo::VolumetricMeshes::TetMesh tetMesh((tetIPCExampleDir() / "box.veg").string().c_str());
  const auto values = stressJson.at("values").get<std::vector<double>>();
  ASSERT_EQ(values.size(), static_cast<std::size_t>(tetMesh.getNumElements()));
  for (double value : values) {
    EXPECT_TRUE(std::isfinite(value));
    EXPECT_GE(value, 0.0);
  }
}

TEST(RunIPCSimCliGTest, TetSurfacePressureForceOneStepWritesOutputsAndStress)
{
  const fs::path binary = runIPCSimBinaryPath();
  ASSERT_FALSE(binary.empty());
  ASSERT_TRUE(fs::exists(binary));

  ScopedTempDir tempDir;
  const fs::path configPath = tempDir.path() / "tet-ipc-pressure.json";
  const fs::path outputDir = tempDir.path() / "tet-output";

  writeTextFile(configPath,
    addBoolConfigField(addSurfacePressureForceConfig(makeTetIPCConfig(tempDir.path(), 1), true), "output-von-mises", true));

  std::ostringstream command;
  command << shellExecutable(binary)
          << " "
          << quotePath(configPath);

  ASSERT_EQ(runCommand(command.str()), 0);
  EXPECT_TRUE(fs::exists(statePath(outputDir, 0)));
  EXPECT_TRUE(fs::exists(surfacePath(outputDir, 0)));
  ASSERT_TRUE(fs::exists(stressPath(outputDir, 0)));

  const nlohmann::json stressJson = readJsonFile(stressPath(outputDir, 0));
  const auto values = stressJson.at("values").get<std::vector<double>>();
  ASSERT_FALSE(values.empty());
  bool hasNonzeroStress = false;
  for (double value : values) {
    EXPECT_TRUE(std::isfinite(value));
    EXPECT_GE(value, 0.0);
    hasNonzeroStress = hasNonzeroStress || value > 0.0;
  }
  EXPECT_TRUE(hasNonzeroStress);
}

TEST(RunIPCSimCliGTest, CubicOneTimestepSmokeWritesDeformAndRet)
{
  const fs::path binary = runIPCSimBinaryPath();
  ASSERT_FALSE(binary.empty());
  ASSERT_TRUE(fs::exists(binary));

  ScopedTempDir tempDir;
  const fs::path configPath = tempDir.path() / "cubic-ipc-one-step.json";

  writeTextFile(configPath, makeCubicIPCConfig(tempDir.path(), 1));

  std::ostringstream command;
  command << shellExecutable(binary)
          << " "
          << quotePath(configPath);

  ASSERT_EQ(runCommand(command.str()), 0);
  const fs::path outputDir = tempDir.path() / "cubic-output";
  EXPECT_TRUE(fs::exists(statePath(outputDir, 0)));
  EXPECT_TRUE(fs::exists(surfacePath(outputDir, 0)));
  EXPECT_FALSE(fs::exists(outputDir / "deform0000.u"));
  EXPECT_FALSE(fs::exists(outputDir / "ret0000.obj"));
}

TEST(RunIPCSimCliGTest, CubicOneTimestepFloorSmokeWritesDeformAndRet)
{
  const fs::path binary = runIPCSimBinaryPath();
  ASSERT_FALSE(binary.empty());
  ASSERT_TRUE(fs::exists(binary));

  ScopedTempDir tempDir;
  const fs::path configPath = tempDir.path() / "cubic-ipc-floor-one-step.json";

  writeTextFile(configPath, makeCubicIPCConfig(tempDir.path(), 1, 1.0, true, false, "stable-neo", 1, "info", true, "y", -0.15, 4000.0));

  std::ostringstream command;
  command << shellExecutable(binary)
          << " "
          << quotePath(configPath);

  ASSERT_EQ(runCommand(command.str()), 0);
  const fs::path outputDir = tempDir.path() / "cubic-output";
  EXPECT_TRUE(fs::exists(statePath(outputDir, 0)));
  EXPECT_TRUE(fs::exists(surfacePath(outputDir, 0)));
  EXPECT_FALSE(fs::exists(outputDir / "deform0000.u"));
  EXPECT_FALSE(fs::exists(outputDir / "ret0000.obj"));
}

TEST(RunIPCSimCliGTest, TetMovingUpperFloorSmokeWritesStress)
{
  const fs::path binary = runIPCSimBinaryPath();
  ASSERT_FALSE(binary.empty());
  ASSERT_TRUE(fs::exists(binary));

  ScopedTempDir tempDir;
  const fs::path configPath = tempDir.path() / "tet-ipc-moving-upper-floor.json";

  writeTextFile(configPath, addBoolConfigField(addMovingUpperFloorConfig(makeTetIPCConfig(tempDir.path(), 2)), "output-von-mises", true));

  std::ostringstream command;
  command << shellExecutable(binary)
          << " "
          << quotePath(configPath);

  ASSERT_EQ(runCommand(command.str()), 0);
  const fs::path outputDir = tempDir.path() / "tet-output";
  EXPECT_TRUE(fs::exists(statePath(outputDir, 0)));
  EXPECT_TRUE(fs::exists(statePath(outputDir, 1)));
  EXPECT_TRUE(fs::exists(stressPath(outputDir, 1)));

  const nlohmann::json stress = readJsonFile(stressPath(outputDir, 1));
  ASSERT_TRUE(stress.contains("values"));
  bool hasNonzeroStress = false;
  for (const auto &value : stress["values"])
    hasNonzeroStress = hasNonzeroStress || value.get<double>() > 0.0;
  EXPECT_TRUE(hasNonzeroStress);
}

TEST(RunIPCSimCliGTest, TetRejectsIPCHeuristic)
{
  const fs::path binary = runIPCSimBinaryPath();
  ASSERT_FALSE(binary.empty());
  ASSERT_TRUE(fs::exists(binary));

  ScopedTempDir tempDir;
  const fs::path configPath = tempDir.path() / "tet-ipc-heuristic.json";

  writeTextFile(configPath, makeTetIPCConfig(tempDir.path(), 0, 1.0, true, true));

  std::ostringstream command;
  command << shellExecutable(binary)
          << " "
          << quotePath(configPath);

  ASSERT_NE(runCommand(command.str()), 0);
}

TEST(RunIPCSimCliGTest, CubicRejectsShellOnlyElasticMaterial)
{
  const fs::path binary = runIPCSimBinaryPath();
  ASSERT_FALSE(binary.empty());
  ASSERT_TRUE(fs::exists(binary));

  ScopedTempDir tempDir;
  const fs::path configPath = tempDir.path() / "cubic-ipc-koiter.json";

  writeTextFile(configPath, makeCubicIPCConfig(tempDir.path(), 0, 1.0, true, false, "koiter-stvk"));

  std::ostringstream command;
  command << shellExecutable(binary)
          << " "
          << quotePath(configPath);

  ASSERT_NE(runCommand(command.str()), 0);
}

TEST(RunIPCSimCliGTest, TetMissingIPCDhatOrKappaFails)
{
  const fs::path binary = runIPCSimBinaryPath();
  ASSERT_FALSE(binary.empty());
  ASSERT_TRUE(fs::exists(binary));

  ScopedTempDir tempDir;
  const fs::path configPath = tempDir.path() / "tet-ipc-missing.json";

  writeTextFile(configPath, makeTetIPCConfig(tempDir.path(), 0, 1.0, false));

  std::ostringstream command;
  command << shellExecutable(binary)
          << " "
          << quotePath(configPath);

  ASSERT_NE(runCommand(command.str()), 0);
}

TEST(RunIPCSimCliGTest, TetNonUnitScaleSmokeSucceeds)
{
  const fs::path binary = runIPCSimBinaryPath();
  ASSERT_FALSE(binary.empty());
  ASSERT_TRUE(fs::exists(binary));

  ScopedTempDir tempDir;
  const fs::path configPath = tempDir.path() / "tet-ipc-scale-two.json";

  writeTextFile(configPath, makeTetIPCConfig(tempDir.path(), 1, 2.0));

  std::ostringstream command;
  command << shellExecutable(binary)
          << " "
          << quotePath(configPath);

  ASSERT_EQ(runCommand(command.str()), 0);
  const fs::path outputDir = tempDir.path() / "tet-output";
  EXPECT_TRUE(fs::exists(statePath(outputDir, 0)));
  EXPECT_TRUE(fs::exists(surfacePath(outputDir, 0)));
  EXPECT_FALSE(fs::exists(outputDir / "deform0000.u"));
  EXPECT_FALSE(fs::exists(outputDir / "ret0000.obj"));
}

TEST(RunIPCSimCliGTest, CubicNonUnitScaleSmokeSucceeds)
{
  const fs::path binary = runIPCSimBinaryPath();
  ASSERT_FALSE(binary.empty());
  ASSERT_TRUE(fs::exists(binary));

  ScopedTempDir tempDir;
  const fs::path configPath = tempDir.path() / "cubic-ipc-scale-two.json";

  writeTextFile(configPath, makeCubicIPCConfig(tempDir.path(), 1, 2.0));

  std::ostringstream command;
  command << shellExecutable(binary)
          << " "
          << quotePath(configPath);

  ASSERT_EQ(runCommand(command.str()), 0);
  const fs::path outputDir = tempDir.path() / "cubic-output";
  EXPECT_TRUE(fs::exists(statePath(outputDir, 0)));
  EXPECT_TRUE(fs::exists(surfacePath(outputDir, 0)));
  EXPECT_FALSE(fs::exists(outputDir / "deform0000.u"));
  EXPECT_FALSE(fs::exists(outputDir / "ret0000.obj"));
}

TEST(RunIPCSimSetupGTest, TetEmbeddingMatrixMatchesBarycentricBaseline)
{
  initializeRunIPCSimTestEnvironment();

  pgo::ConfigFileJSON config;
  ASSERT_TRUE(config.open(tetIPCConfigPath().string().c_str()));

  const auto context = pgo::RunIPCSim::buildVolumeIpcSimulation(config);
  const ES::SpMatD expected = computeExpectedEmbeddingMatrix(tetIPCConfigPath());

  expectSparseMatrixNear(context.surfaceFromSimulationDispMap, expected);
}

TEST(RunIPCSimSetupGTest, CubicEmbeddingMatrixMatchesBarycentricBaseline)
{
  initializeRunIPCSimTestEnvironment();

  pgo::ConfigFileJSON config;
  ASSERT_TRUE(config.open(cubicIPCConfigPath().string().c_str()));

  const auto context = pgo::RunIPCSim::buildVolumeIpcSimulation(config);
  const ES::SpMatD expected = computeExpectedEmbeddingMatrix(cubicIPCConfigPath());

  expectSparseMatrixNear(context.surfaceFromSimulationDispMap, expected);
}

TEST(RunIPCSimSetupGTest, VolumeSurfacePressureForceProjectsToSimulationDofs)
{
  initializeRunIPCSimTestEnvironment();

  ScopedTempDir tempDir;
  const fs::path configPath = tempDir.path() / "tet-pressure-setup.json";
  writeTextFile(configPath, addSurfacePressureForceConfig(makeTetIPCConfig(tempDir.path(), 0), true, 1000.0, 20));

  pgo::ConfigFileJSON config;
  ASSERT_TRUE(config.open(configPath.string().c_str()));

  const auto context = pgo::RunIPCSim::buildVolumeIpcSimulation(config);
  EXPECT_TRUE(context.surfacePressureForceEnabled);
  EXPECT_EQ(context.surfacePressureRampSteps, 20);
  ASSERT_EQ(context.surfacePressureSimulationForce.size(), context.simulationRestPosition.size());
  EXPECT_GT(context.surfacePressureSimulationForce.norm(), 0.0);
}

TEST(RunIPCSimSetupGTest, VolumeSurfacePressureForceDefaultsRampStepsToOne)
{
  initializeRunIPCSimTestEnvironment();

  ScopedTempDir tempDir;
  const fs::path configPath = tempDir.path() / "tet-pressure-default-ramp.json";
  writeTextFile(configPath, addSurfacePressureForceConfig(makeTetIPCConfig(tempDir.path(), 0), true, 1000.0, std::nullopt));

  pgo::ConfigFileJSON config;
  ASSERT_TRUE(config.open(configPath.string().c_str()));

  const auto context = pgo::RunIPCSim::buildVolumeIpcSimulation(config);
  EXPECT_TRUE(context.surfacePressureForceEnabled);
  EXPECT_EQ(context.surfacePressureRampSteps, 1);
}

TEST(RunIPCSimSetupGTest, SurfacePressureForceDoesNotRequireCenter)
{
  initializeRunIPCSimTestEnvironment();

  ScopedTempDir tempDir;
  const fs::path configPath = tempDir.path() / "tet-pressure-no-center.json";
  writeTextFile(configPath, addSurfacePressureForceConfigWithoutCenter(makeTetIPCConfig(tempDir.path(), 0), true, 1000.0, 20));

  pgo::ConfigFileJSON config;
  ASSERT_TRUE(config.open(configPath.string().c_str()));

  const auto context = pgo::RunIPCSim::buildVolumeIpcSimulation(config);
  EXPECT_TRUE(context.surfacePressureForceEnabled);
  ASSERT_EQ(context.surfacePressureSimulationForce.size(), context.simulationRestPosition.size());
  EXPECT_GT(context.surfacePressureSimulationForce.norm(), 0.0);
}

TEST(RunIPCSimSetupGTest, SurfacePressureForceIgnoresDeprecatedCenterField)
{
  initializeRunIPCSimTestEnvironment();

  ScopedTempDir tempDir;
  const fs::path noCenterConfigPath = tempDir.path() / "tet-pressure-no-center-reference.json";
  const fs::path autoConfigPath = tempDir.path() / "tet-pressure-auto-center.json";
  const fs::path explicitConfigPath = tempDir.path() / "tet-pressure-explicit-center.json";
  writeTextFile(noCenterConfigPath, addSurfacePressureForceConfigWithoutCenter(makeTetIPCConfig(tempDir.path(), 0, 2.0), true, 1000.0, 20));
  writeTextFile(autoConfigPath, addSurfacePressureForceConfig(makeTetIPCConfig(tempDir.path(), 0, 2.0), true, 1000.0, 20, "\"auto\""));

  pgo::ConfigFileJSON noCenterConfig;
  ASSERT_TRUE(noCenterConfig.open(noCenterConfigPath.string().c_str()));
  const auto noCenterContext = pgo::RunIPCSim::buildVolumeIpcSimulation(noCenterConfig);

  pgo::ConfigFileJSON autoConfig;
  ASSERT_TRUE(autoConfig.open(autoConfigPath.string().c_str()));
  const auto autoContext = pgo::RunIPCSim::buildVolumeIpcSimulation(autoConfig);

  const ES::V3d explicitCenter(0.25, -0.5, 0.75);
  writeTextFile(explicitConfigPath,
    addSurfacePressureForceConfig(makeTetIPCConfig(tempDir.path(), 0, 2.0), true, 1000.0, 20, vec3Json(explicitCenter)));

  pgo::ConfigFileJSON explicitConfig;
  ASSERT_TRUE(explicitConfig.open(explicitConfigPath.string().c_str()));

  const auto explicitContext = pgo::RunIPCSim::buildVolumeIpcSimulation(explicitConfig);
  ASSERT_EQ(noCenterContext.surfacePressureSimulationForce.size(), autoContext.surfacePressureSimulationForce.size());
  ASSERT_EQ(noCenterContext.surfacePressureSimulationForce.size(), explicitContext.surfacePressureSimulationForce.size());
  EXPECT_GT(noCenterContext.surfacePressureSimulationForce.norm(), 0.0);
  EXPECT_NEAR((noCenterContext.surfacePressureSimulationForce - autoContext.surfacePressureSimulationForce).norm(), 0.0, 1e-8);
  EXPECT_NEAR((noCenterContext.surfacePressureSimulationForce - explicitContext.surfacePressureSimulationForce).norm(), 0.0, 1e-8);
  EXPECT_GT(autoContext.surfacePressureSimulationForce.norm(), 0.0);
}

TEST(RunIPCSimSetupGTest, SurfacePressureForceDisabledOrZeroPressureProducesNoContribution)
{
  initializeRunIPCSimTestEnvironment();

  ScopedTempDir tempDir;
  const fs::path disabledConfigPath = tempDir.path() / "tet-pressure-disabled.json";
  const fs::path zeroConfigPath = tempDir.path() / "tet-pressure-zero.json";
  writeTextFile(disabledConfigPath, addSurfacePressureForceConfig(makeTetIPCConfig(tempDir.path(), 0), false));
  writeTextFile(zeroConfigPath, addSurfacePressureForceConfig(makeTetIPCConfig(tempDir.path(), 0), true, 0.0, 20));

  pgo::ConfigFileJSON disabledConfig;
  ASSERT_TRUE(disabledConfig.open(disabledConfigPath.string().c_str()));
  const auto disabledContext = pgo::RunIPCSim::buildVolumeIpcSimulation(disabledConfig);
  EXPECT_FALSE(disabledContext.surfacePressureForceEnabled);
  EXPECT_EQ(disabledContext.surfacePressureSimulationForce.size(), 0);

  pgo::ConfigFileJSON zeroConfig;
  ASSERT_TRUE(zeroConfig.open(zeroConfigPath.string().c_str()));
  const auto zeroContext = pgo::RunIPCSim::buildVolumeIpcSimulation(zeroConfig);
  EXPECT_TRUE(zeroContext.surfacePressureForceEnabled);
  ASSERT_EQ(zeroContext.surfacePressureSimulationForce.size(), zeroContext.simulationRestPosition.size());
  EXPECT_DOUBLE_EQ(zeroContext.surfacePressureSimulationForce.norm(), 0.0);
}

TEST(RunIPCSimSetupGTest, ShellRejectsEnabledSurfacePressureForce)
{
  initializeRunIPCSimTestEnvironment();

  ScopedTempDir tempDir;
  const fs::path configPath = tempDir.path() / "shell-pressure-enabled.json";
  writeTextFile(configPath, addSurfacePressureForceConfig(makeShellIPCConfig(tempDir.path(), 0), true));

  pgo::ConfigFileJSON config;
  ASSERT_TRUE(config.open(configPath.string().c_str()));

  try {
    (void)pgo::RunIPCSim::buildShellIpcSimulation(config);
    FAIL() << "Expected enabled surface-pressure-force to be rejected on the shell path.";
  }
  catch (const std::invalid_argument &e) {
    EXPECT_NE(std::string(e.what()).find("surface-pressure-force"), std::string::npos);
    EXPECT_NE(std::string(e.what()).find("volume"), std::string::npos);
  }
}

TEST(RunIPCSimSetupGTest, FloorsArrayAcceptsEmptyArrayAndRequiresAxisHeightOrMotionAndKappa)
{
  initializeRunIPCSimTestEnvironment();

  ScopedTempDir tempDir;
  const fs::path emptyFloorsConfig = tempDir.path() / "shell-empty-floors.json";
  const fs::path missingAxisConfig = tempDir.path() / "shell-floor-missing-axis.json";
  const fs::path missingHeightConfig = tempDir.path() / "shell-floor-missing-height.json";
  const fs::path missingKappaConfig = tempDir.path() / "shell-floor-missing-kappa.json";

  writeTextFile(emptyFloorsConfig, addEmptyFloorsConfig(makeShellIPCConfig(tempDir.path(), 0)));
  writeTextFile(missingAxisConfig, makeShellIPCConfig(tempDir.path(), 0, true, false, 0.002, 3000.0, 1, "info", true, std::nullopt, -0.1, 4000.0));
  writeTextFile(missingHeightConfig, makeShellIPCConfig(tempDir.path(), 0, true, false, 0.002, 3000.0, 1, "info", true, "y", std::nullopt, 4000.0));
  writeTextFile(missingKappaConfig, makeShellIPCConfig(tempDir.path(), 0, true, false, 0.002, 3000.0, 1, "info", true, "y", -0.1, std::nullopt));

  pgo::ConfigFileJSON emptyFloors;
  ASSERT_TRUE(emptyFloors.open(emptyFloorsConfig.string().c_str()));
  EXPECT_EQ(pgo::RunIPCSim::buildShellIpcSimulation(emptyFloors).extraGeneralImplicitForceModels.size(), 0u);

  pgo::ConfigFileJSON missingAxis;
  ASSERT_TRUE(missingAxis.open(missingAxisConfig.string().c_str()));
  EXPECT_THROW(pgo::RunIPCSim::buildShellIpcSimulation(missingAxis), std::invalid_argument);

  pgo::ConfigFileJSON missingHeight;
  ASSERT_TRUE(missingHeight.open(missingHeightConfig.string().c_str()));
  EXPECT_THROW(pgo::RunIPCSim::buildShellIpcSimulation(missingHeight), std::invalid_argument);

  pgo::ConfigFileJSON missingKappa;
  ASSERT_TRUE(missingKappa.open(missingKappaConfig.string().c_str()));
  EXPECT_THROW(pgo::RunIPCSim::buildShellIpcSimulation(missingKappa), std::invalid_argument);
}

TEST(RunIPCSimSetupGTest, FloorsArrayRejectsLegacyFieldsAndInvalidHeightMotionCombinations)
{
  initializeRunIPCSimTestEnvironment();

  ScopedTempDir tempDir;
  const fs::path legacyConfigPath = tempDir.path() / "shell-legacy-floor.json";
  const fs::path heightAndMotionConfigPath = tempDir.path() / "shell-floor-height-and-motion.json";
  const fs::path missingHeightAndMotionConfigPath = tempDir.path() / "shell-floor-missing-height-and-motion.json";

  writeTextFile(legacyConfigPath, addBoolConfigField(makeShellIPCConfig(tempDir.path(), 0), "use-floor", true));
  writeTextFile(heightAndMotionConfigPath, addTopLevelJsonField(makeShellIPCConfig(tempDir.path(), 0),
    "  \"floors\": [\n"
    "    {\n"
    "      \"axis\": \"y\",\n"
    "      \"side\": \"upper\",\n"
    "      \"height\": -0.1,\n"
    "      \"kappa\": 4000.0,\n"
    "      \"motion\": {\n"
    "        \"height-start\": 1.0,\n"
    "        \"height-end\": 0.75,\n"
    "        \"frame-start\": 0,\n"
    "        \"frame-end\": 1\n"
    "      }\n"
    "    }\n"
    "  ]"));
  writeTextFile(missingHeightAndMotionConfigPath, addTopLevelJsonField(makeShellIPCConfig(tempDir.path(), 0),
    "  \"floors\": [\n"
    "    {\n"
    "      \"axis\": \"y\",\n"
    "      \"side\": \"upper\",\n"
    "      \"kappa\": 4000.0\n"
    "    }\n"
    "  ]"));

  pgo::ConfigFileJSON legacyConfig;
  ASSERT_TRUE(legacyConfig.open(legacyConfigPath.string().c_str()));
  EXPECT_THROW(pgo::RunIPCSim::buildShellIpcSimulation(legacyConfig), std::invalid_argument);

  pgo::ConfigFileJSON heightAndMotionConfig;
  ASSERT_TRUE(heightAndMotionConfig.open(heightAndMotionConfigPath.string().c_str()));
  EXPECT_THROW(pgo::RunIPCSim::buildShellIpcSimulation(heightAndMotionConfig), std::invalid_argument);

  pgo::ConfigFileJSON missingHeightAndMotionConfig;
  ASSERT_TRUE(missingHeightAndMotionConfig.open(missingHeightAndMotionConfigPath.string().c_str()));
  EXPECT_THROW(pgo::RunIPCSim::buildShellIpcSimulation(missingHeightAndMotionConfig), std::invalid_argument);
}

TEST(RunIPCSimSetupGTest, FloorEnabledSetupCreatesExtraGeneralImplicitForceModel)
{
  initializeRunIPCSimTestEnvironment();

  ScopedTempDir tempDir;
  const fs::path configPath = tempDir.path() / "cubic-floor-setup.json";
  writeTextFile(configPath, makeCubicIPCConfig(tempDir.path(), 0, 1.0, true, false, "stable-neo", 1, "info", true, "y", -0.15, 4000.0));

  pgo::ConfigFileJSON config;
  ASSERT_TRUE(config.open(configPath.string().c_str()));

  const auto context = pgo::RunIPCSim::buildVolumeIpcSimulation(config);
  EXPECT_EQ(context.extraGeneralImplicitForceModels.size(), 1u);
  EXPECT_EQ(context.floorPotentialEnergies.size(), 1u);
  EXPECT_EQ(context.floorMotionStates.size(), 1u);
}

TEST(RunIPCSimSetupGTest, MultipleFloorsCreateMultipleForceModels)
{
  initializeRunIPCSimTestEnvironment();

  ScopedTempDir tempDir;
  const fs::path configPath = tempDir.path() / "cubic-multi-floor-setup.json";
  writeTextFile(configPath, addTopLevelJsonField(makeCubicIPCConfig(tempDir.path(), 0),
    "  \"floors\": [\n"
    "    { \"axis\": \"y\", \"side\": \"lower\", \"height\": 0.0, \"kappa\": 4000.0 },\n"
    "    { \"axis\": \"y\", \"side\": \"upper\", \"height\": 0.8, \"kappa\": 4000.0 }\n"
    "  ]"));

  pgo::ConfigFileJSON config;
  ASSERT_TRUE(config.open(configPath.string().c_str()));

  const auto context = pgo::RunIPCSim::buildVolumeIpcSimulation(config);
  EXPECT_EQ(context.extraGeneralImplicitForceModels.size(), 2u);
  EXPECT_EQ(context.floorPotentialEnergies.size(), 2u);
  EXPECT_EQ(context.floorMotionStates.size(), 2u);
}
