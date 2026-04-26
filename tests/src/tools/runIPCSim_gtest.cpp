#include <gtest/gtest.h>

#include "EigenSupport.h"
#include "barycentricCoordinates.h"
#include "configFileJSON.h"
#include "deformationModelEnergy.h"
#include "initPredicates.h"
#include "pgoLogging.h"
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

std::string addSurfacePressureForceConfig(std::string json, bool enabled, double pressure = 1000.0, int rampSteps = 20)
{
  const std::string marker = "\n}\n";
  const std::size_t pos = json.rfind(marker);
  if (pos == std::string::npos)
    throw std::runtime_error("Failed to add surface pressure force config to test JSON.");

  std::ostringstream field;
  field << ",\n"
        << "  \"surface-pressure-force\": {\n"
        << "    \"enabled\": " << (enabled ? "true" : "false") << ",\n"
        << "    \"center\": [0, 0, 0],\n"
        << "    \"pressure\": " << pressure << ",\n"
        << "    \"ramp-steps\": " << rampSteps << "\n"
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
       << "  \"use-floor\": " << (useFloor ? "true" : "false");
  if (floorAxis.has_value())
    json << ",\n"
         << "  \"floor-axis\": \"" << *floorAxis << "\"";
  if (floorHeight.has_value())
    json << ",\n"
         << "  \"floor-height\": " << *floorHeight;
  if (floorKappa.has_value())
    json << ",\n"
         << "  \"floor-kappa\": " << *floorKappa;
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
  std::optional<double> floorKappa = std::nullopt)
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
       << "  \"solver-max-iter\": 5,\n"
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
  std::optional<double> floorKappa = std::nullopt)
{
  return makeVolumeIPCConfig(tetIPCExampleDir(), "tet-mesh", tempDir / "tet-output", numTimesteps,
    scale, includeIPCFields, ipcHeuristic, material, 0.002, 3000.0, dumpInterval, true, logLevel, useFloor, floorAxis, floorHeight, floorKappa);
}

std::string makeCubicIPCConfig(const fs::path &tempDir, int numTimesteps,
  double scale = 1.0, bool includeIPCFields = true, bool ipcHeuristic = false, const std::string &material = "stable-neo",
  int dumpInterval = 1, const std::string &logLevel = "info",
  bool useFloor = false,
  std::optional<std::string> floorAxis = std::nullopt,
  std::optional<double> floorHeight = std::nullopt,
  std::optional<double> floorKappa = std::nullopt)
{
  return makeVolumeIPCConfig(cubicIPCExampleDir(), "cubic-mesh", tempDir / "cubic-output", numTimesteps,
    scale, includeIPCFields, ipcHeuristic, material, 0.002, 3000.0, dumpInterval, true, logLevel, useFloor, floorAxis, floorHeight, floorKappa);
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
}  // namespace

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
  EXPECT_NE(contents.find("callCount="), std::string::npos);
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
  EXPECT_NE(contents.find("use-floor=true"), std::string::npos);
  EXPECT_NE(contents.find("floor-axis=y"), std::string::npos);
  EXPECT_NE(contents.find("floor-height=-0.15"), std::string::npos);
  EXPECT_NE(contents.find("floor-kappa=4321"), std::string::npos);
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

TEST(RunIPCSimCliGTest, DeformStateIsWrittenEveryTimestep)
{
  const fs::path binary = runIPCSimBinaryPath();
  ASSERT_FALSE(binary.empty());
  ASSERT_TRUE(fs::exists(binary));

  ScopedTempDir tempDir;
  const fs::path configPath = tempDir.path() / "shell-ipc-deform-every-step.json";

  writeTextFile(configPath, makeShellIPCConfig(tempDir.path(), 2, true, false, 0.002, 3000.0, 10));

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

TEST(RunIPCSimSetupGTest, UseFloorRequiresExplicitAxisHeightAndKappa)
{
  initializeRunIPCSimTestEnvironment();

  ScopedTempDir tempDir;
  const fs::path missingAxisConfig = tempDir.path() / "shell-floor-missing-axis.json";
  const fs::path missingHeightConfig = tempDir.path() / "shell-floor-missing-height.json";
  const fs::path missingKappaConfig = tempDir.path() / "shell-floor-missing-kappa.json";

  writeTextFile(missingAxisConfig, makeShellIPCConfig(tempDir.path(), 0, true, false, 0.002, 3000.0, 1, "info", true, std::nullopt, -0.1, 4000.0));
  writeTextFile(missingHeightConfig, makeShellIPCConfig(tempDir.path(), 0, true, false, 0.002, 3000.0, 1, "info", true, "y", std::nullopt, 4000.0));
  writeTextFile(missingKappaConfig, makeShellIPCConfig(tempDir.path(), 0, true, false, 0.002, 3000.0, 1, "info", true, "y", -0.1, std::nullopt));

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
}
