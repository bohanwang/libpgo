#include <gtest/gtest.h>

#include "EigenSupport.h"
#include "barycentricCoordinates.h"
#include "configFileJSON.h"
#include "deformationModelEnergy.h"
#include "embeddedSurfaceIPCPotentialEnergy.h"
#include "initPredicates.h"
#include "pgoLogging.h"
#include "runIPCSimSetup.h"
#include "runIPCExternalObjects.h"
#include "runSimCliLogging.h"
#include "runSimVolumeMeshIO.h"
#include "simulationRunner.h"
#include "triMeshGeo.h"
#include "volumetricMesh.h"

#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

namespace
{
namespace ES = pgo::EigenSupport;
namespace fs = std::filesystem;

constexpr const char *kShellExampleDir = LIBPGO_TEST_SHELL_EXAMPLE_DIR;
constexpr const char *kTetIPCExampleDir = LIBPGO_TEST_IPC_TET_EXAMPLE_DIR;
constexpr const char *kCubicIPCExampleDir = LIBPGO_TEST_IPC_CUBIC_EXAMPLE_DIR;

std::string quotePath(const fs::path &path)
{
  // Forward slashes are accepted by Windows command-line tools and do not
  // become accidental escape sequences when this helper is used in JSON.
  return "\"" + path.generic_string() + "\"";
}

std::string shellExecutable(const fs::path &path)
{
#ifdef _WIN32
  return "call " + quotePath(path);
#else
  return quotePath(path);
#endif
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

void writeExternalPlane(const fs::path &path, double y = 0.0, bool includeIsolatedVertex = false)
{
  std::ostringstream obj;
  obj << "v -1 " << y << " -1\n"
      << "v 1 " << y << " -1\n"
      << "v 1 " << y << " 1\n"
      << "v -1 " << y << " 1\n";
  if (includeIsolatedVertex)
    obj << "v 20 20 20\n";
  obj << "f 1 3 2\n"
      << "f 1 4 3\n";
  writeTextFile(path, obj.str());
}

void addIPCExternalObject(pgo::ConfigFileJSON &config, const fs::path &path,
  const nlohmann::json &movement = nlohmann::json::array({ 0.0, 0.0, 0.0 }),
  const nlohmann::json &scale = 1.0,
  const nlohmann::json &initialTranslation = nlohmann::json::array({ 0.0, 0.0, 0.0 }))
{
  config.handle()["ipc-external-objects"] = nlohmann::json::array(
    { { { "filename", path.generic_string() }, { "scale", scale }, { "initial-translation", initialTranslation }, { "movement", movement } } });
}

std::string readTextFile(const fs::path &path)
{
  std::ifstream in(path);
  EXPECT_TRUE(in.is_open());
  return std::string((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
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

std::string makeStaticConfig(std::string json)
{
  const std::string dynamicField = "\"sim-type\": \"dynamic\"";
  const std::size_t pos = json.find(dynamicField);
  if (pos == std::string::npos)
    throw std::runtime_error("Failed to find dynamic sim-type in test JSON.");

  json.replace(pos, dynamicField.size(), "\"sim-type\": \"static\"");
  return json;
}

std::string replaceTextOnce(std::string text, const std::string &from, const std::string &to)
{
  const std::size_t pos = text.find(from);
  if (pos == std::string::npos)
    throw std::runtime_error("Failed to find text in test configuration: " + from);

  text.replace(pos, from.size(), to);
  return text;
}

void expectFiniteStaticOutput(const fs::path &outputDir)
{
  const fs::path statePath = outputDir / "deform0000.u";
  const fs::path surfacePath = outputDir / "ret0000.obj";
  ASSERT_TRUE(fs::exists(statePath));
  ASSERT_TRUE(fs::exists(surfacePath));

  ES::MXd state;
  ASSERT_EQ(ES::readMatrix(statePath.string().c_str(), state), 0);
  ASSERT_EQ(state.cols(), 3);
  EXPECT_TRUE(state.allFinite());
  EXPECT_GT(state.col(0).norm(), 1e-12);

  pgo::Mesh::TriMeshGeo surface;
  ASSERT_TRUE(surface.load(surfacePath.string()));
  for (int vi = 0; vi < surface.numVertices(); ++vi)
    EXPECT_TRUE(surface.pos(vi).allFinite());
}

void writeZeroShellRestartState(const fs::path &outputDir, int frame)
{
  pgo::Mesh::TriMeshGeo mesh;
  ASSERT_TRUE(mesh.load((fs::path(kShellExampleDir) / "shell.obj").string()));

  fs::create_directories(outputDir);
  ES::MXd restartState = ES::MXd::Zero(mesh.numVertices() * 3, 3);
  std::ostringstream filename;
  filename << "deform" << std::setfill('0') << std::setw(4) << frame << ".u";
  ASSERT_EQ(ES::writeMatrix((outputDir / filename.str()).string().c_str(), restartState), 0);
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
  const std::string &logLevel = "info")
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
       << "  \"contact-model\": \"ipc\",\n"
       << "  \"solver-eps\": 1e-4,\n"
       << "  \"solver-max-iter\": 5,\n"
       << "  \"elastic-material\": \"koiter-stvk\",\n"
       << "  \"loglevel\": \"" << logLevel << "\",\n"
       << "  \"dump-interval\": " << dumpInterval << ",\n"
       << "  \"output\": " << quotePath(outputDir);

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
       << "  \"contact-model\": \"ipc\",\n"
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
  bool enableMaterialMaxStep = true, const std::string &logLevel = "info")
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
       << "  \"contact-model\": \"ipc\",\n"
       << "  \"solver-eps\": 1e-4,\n"
       << "  \"solver-max-iter\": 5,\n"
       << "  \"elastic-material\": \"" << material << "\",\n"
       << "  \"loglevel\": \"" << logLevel << "\",\n"
       << "  \"dump-interval\": " << dumpInterval << ",\n"
       << "  \"output\": " << quotePath(outputDir) << ",\n"
       << "  \"enable-material-max-step\": " << (enableMaterialMaxStep ? "true" : "false");

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
  int dumpInterval = 1, const std::string &logLevel = "info")
{
  return makeVolumeIPCConfig(tetIPCExampleDir(), "tet-mesh", tempDir / "tet-output", numTimesteps,
    scale, includeIPCFields, ipcHeuristic, material, 0.002, 3000.0, dumpInterval, true, logLevel);
}

std::string makeCubicIPCConfig(const fs::path &tempDir, int numTimesteps,
  double scale = 1.0, bool includeIPCFields = true, bool ipcHeuristic = false, const std::string &material = "stable-neo",
  int dumpInterval = 1, const std::string &logLevel = "info")
{
  return makeVolumeIPCConfig(cubicIPCExampleDir(), "cubic-mesh", tempDir / "cubic-output", numTimesteps,
    scale, includeIPCFields, ipcHeuristic, material, 0.002, 3000.0, dumpInterval, true, logLevel);
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
  EXPECT_EQ(contents.find("max-step summary"), std::string::npos);
  EXPECT_EQ(contents.find("minFeasibleAlphaThisSolve"), std::string::npos);
  EXPECT_EQ(contents.find("minLineSearchAlphaThisSolve"), std::string::npos);
  EXPECT_EQ(contents.find("minEffectiveAlphaThisSolve"), std::string::npos);
  EXPECT_EQ(contents.find("finalAlpha"), std::string::npos);
  EXPECT_EQ(contents.find("lastMaterialAlpha"), std::string::npos);
  EXPECT_EQ(contents.find("lastContactAlpha"), std::string::npos);
}

TEST(RunIPCSimCliGTest, DefaultRunPreservesUnrelatedOutputAndOverwritesFrames)
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
  EXPECT_TRUE(fs::exists(sentinelPath));
  ASSERT_TRUE(fs::exists(logPath));

  const std::string contents = readTextFile(logPath);
  EXPECT_NE(contents.find("existing output files will be overwritten as frames are written"), std::string::npos);
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
  std::string config = replaceTextOnce(makeShellIPCConfig(tempDir.path(), 2),
    "\"init-disp\": [0, 0, 0]", "\"init-disp\": [0.1, 0, 0]");
  config = replaceTextOnce(std::move(config),
    "\"init-vel\": [0, 0, 0]", "\"init-vel\": [100, 0, 0]");
  writeTextFile(configPath, addBoolConfigField(std::move(config), "restart-from-u", true));

  std::ostringstream command;
  command << shellExecutable(binary)
          << " --log "
          << quotePath(configPath);

  ASSERT_EQ(runCommand(command.str()), 0);
  ASSERT_TRUE(fs::exists(logPath));
  EXPECT_TRUE(fs::exists(outputDir / "deform0001.u"));

  const std::string contents = readTextFile(logPath);
  EXPECT_NE(contents.find("Restarting from frame 0"), std::string::npos);
  EXPECT_EQ(contents.find("existing output files will be overwritten as frames are written"), std::string::npos);

  ES::MXd state;
  ASSERT_EQ(ES::readMatrix((outputDir / "deform0001.u").string().c_str(), state), 0);
  ASSERT_EQ(state.cols(), 3);
  double maxAbsX = 0.0;
  for (Eigen::Index vi = 0; vi < state.rows() / 3; ++vi)
    maxAbsX = std::max(maxAbsX, std::abs(state(3 * vi, 0)));
  EXPECT_LT(maxAbsX, 0.05);
}

TEST(RunIPCSimCliGTest, MissingDynamicRestartFallsBackToInitialDisplacement)
{
  const fs::path binary = runIPCSimBinaryPath();
  ASSERT_FALSE(binary.empty());
  ASSERT_TRUE(fs::exists(binary));
  ScopedTempDir tempDir;
  const fs::path outputDir = tempDir.path() / "shell-output";
  const fs::path configPath = tempDir.path() / "shell-restart-fallback.json";
  std::string config = replaceTextOnce(makeShellIPCConfig(tempDir.path(), 1),
    "\"init-disp\": [0, 0, 0]", "\"init-disp\": [0.02, 0, 0]");
  writeTextFile(configPath, addBoolConfigField(std::move(config), "restart-from-u", true));

  ASSERT_EQ(runCommand(shellExecutable(binary) + " " + quotePath(configPath)), 0);
  ES::MXd state;
  ASSERT_EQ(ES::readMatrix((outputDir / "deform0000.u").string().c_str(), state), 0);
  double meanX = 0.0;
  for (Eigen::Index vi = 0; vi < state.rows() / 3; ++vi)
    meanX += state(3 * vi, 0);
  meanX /= static_cast<double>(state.rows() / 3);
  EXPECT_GT(meanX, 1e-3);
}

TEST(RunIPCSimCliGTest, RejectsInitiallyIntersectingExternalSurface)
{
  const fs::path binary = runIPCSimBinaryPath();
  ASSERT_FALSE(binary.empty());
  ASSERT_TRUE(fs::exists(binary));
  ScopedTempDir tempDir;
  const fs::path configPath = tempDir.path() / "shell-initial-intersection.json";
  writeTextFile(configPath, makeShellIPCConfig(tempDir.path(), 1));

  pgo::ConfigFileJSON config;
  ASSERT_TRUE(config.open(configPath.string().c_str()));
  config.handle()["ipc-external-objects"] = nlohmann::json::array(
    { { { "filename", (fs::path(kShellExampleDir) / "shell.obj").generic_string() } } });
  writeTextFile(configPath, config.handle().dump(2) + "\n");

  EXPECT_NE(runCommand(shellExecutable(binary) + " " + quotePath(configPath)), 0);
}

TEST(RunIPCSimCliGTest, DebugLogLevelPrintsNewtonStepDiagnosticsWithoutClassification)
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
  EXPECT_EQ(contents.find("max-step summary"), std::string::npos);
  EXPECT_EQ(contents.find("minMaterialFeasibleAlphaThisSolve"), std::string::npos);
  EXPECT_EQ(contents.find("minContactFeasibleAlphaThisSolve"), std::string::npos);
  EXPECT_EQ(contents.find("minFeasibleAlphaThisSolve"), std::string::npos);
  EXPECT_EQ(contents.find("minLineSearchAlphaThisSolve"), std::string::npos);
  EXPECT_EQ(contents.find("minEffectiveAlphaThisSolve"), std::string::npos);
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
  EXPECT_TRUE(fs::exists(tempDir.path() / "shell-output" / "deform0000.u"));
  EXPECT_TRUE(fs::exists(tempDir.path() / "shell-output" / "deform0001.u"));
  EXPECT_TRUE(fs::exists(tempDir.path() / "shell-output" / "ret0000.obj"));
  EXPECT_FALSE(fs::exists(tempDir.path() / "shell-output" / "ret0001.obj"));
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
  EXPECT_TRUE(fs::exists(tempDir.path() / "tet-output" / "deform0000.u"));
  EXPECT_TRUE(fs::exists(tempDir.path() / "tet-output" / "ret0000.obj"));
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
  EXPECT_TRUE(fs::exists(tempDir.path() / "cubic-output" / "deform0000.u"));
  EXPECT_TRUE(fs::exists(tempDir.path() / "cubic-output" / "ret0000.obj"));
}

TEST(RunIPCSimCliGTest, TetStaticIPCSmokeWritesFiniteStateAndSurface)
{
  const fs::path binary = runIPCSimBinaryPath();
  ASSERT_FALSE(binary.empty());
  ASSERT_TRUE(fs::exists(binary));

  ScopedTempDir tempDir;
  const fs::path configPath = tempDir.path() / "tet-ipc-static.json";
  const fs::path fixedVerticesPath = tempDir.path() / "tet-static-fixed.txt";
  writeTextFile(fixedVerticesPath, "0\n4\n20\n");
  std::string config = makeStaticConfig(makeTetIPCConfig(tempDir.path(), 1));
  config = replaceTextOnce(
    std::move(config),
    quotePath(tetIPCExampleDir() / "box-fixed.txt"),
    quotePath(fixedVerticesPath));
  writeTextFile(configPath, config);

  std::ostringstream command;
  command << shellExecutable(binary)
          << " "
          << quotePath(configPath);

  ASSERT_EQ(runCommand(command.str()), 0);
  expectFiniteStaticOutput(tempDir.path() / "tet-output");
}

TEST(RunIPCSimCliGTest, CubicStaticIPCSharedDispatcherWritesFiniteStateAndSurface)
{
  initializeRunIPCSimTestEnvironment();

  ScopedTempDir tempDir;
  const fs::path configPath = tempDir.path() / "cubic-ipc-static.json";
  const fs::path fixedVerticesPath = tempDir.path() / "cubic-static-fixed.txt";
  writeTextFile(fixedVerticesPath, "0\n4\n20\n");
  std::string config = makeStaticConfig(makeCubicIPCConfig(tempDir.path(), 1));
  config = replaceTextOnce(
    std::move(config),
    quotePath(cubicIPCExampleDir() / "box-fixed.txt"),
    quotePath(fixedVerticesPath));
  writeTextFile(configPath, config);

  ASSERT_EQ(pgo::SimulationRunner::runSimulationFromConfig(configPath), 0);
  expectFiniteStaticOutput(tempDir.path() / "cubic-output");
}

TEST(RunIPCSimCliGTest, ShellStaticIPCSmokeWritesFiniteStateAndSurface)
{
  const fs::path binary = runIPCSimBinaryPath();
  ASSERT_FALSE(binary.empty());
  ASSERT_TRUE(fs::exists(binary));

  ScopedTempDir tempDir;
  const fs::path configPath = tempDir.path() / "shell-ipc-static.json";
  std::string config = makeStaticConfig(makeShellIPCConfig(tempDir.path(), 1));
  config = replaceTextOnce(std::move(config),
    "\"init-disp\": [0, 0, 0]", "\"init-disp\": [0.01, 0, 0]");
  writeTextFile(configPath, config);

  std::ostringstream command;
  command << shellExecutable(binary)
          << " "
          << quotePath(configPath);

  ASSERT_EQ(runCommand(command.str()), 0);
  expectFiniteStaticOutput(tempDir.path() / "shell-output");
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
  EXPECT_TRUE(fs::exists(tempDir.path() / "tet-output" / "deform0000.u"));
  EXPECT_TRUE(fs::exists(tempDir.path() / "tet-output" / "ret0000.obj"));
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
  EXPECT_TRUE(fs::exists(tempDir.path() / "cubic-output" / "deform0000.u"));
  EXPECT_TRUE(fs::exists(tempDir.path() / "cubic-output" / "ret0000.obj"));
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

TEST(RunIPCSimSetupGTest, RemovedAnalyticFloorFieldsAreRejected)
{
  initializeRunIPCSimTestEnvironment();
  ScopedTempDir tempDir;
  const fs::path configPath = tempDir.path() / "shell-removed-floor.json";
  writeTextFile(configPath, makeShellIPCConfig(tempDir.path(), 0));

  pgo::ConfigFileJSON config;
  ASSERT_TRUE(config.open(configPath.string().c_str()));
  config.handle()["use-floor"] = true;
  EXPECT_THROW(pgo::RunIPCSim::buildShellIpcSimulation(config), std::invalid_argument);
}

TEST(RunIPCSimSetupGTest, ShellExternalObjectBuildsCombinedCollisionSurfaceOnly)
{
  initializeRunIPCSimTestEnvironment();
  ScopedTempDir tempDir;
  const fs::path configPath = tempDir.path() / "shell-external.json";
  const fs::path externalPath = tempDir.path() / "external-plane.obj";
  writeTextFile(configPath, makeShellIPCConfig(tempDir.path(), 0));
  writeExternalPlane(externalPath, -0.2, true);

  pgo::ConfigFileJSON config;
  ASSERT_TRUE(config.open(configPath.string().c_str()));
  addIPCExternalObject(config, externalPath);
  const auto context = pgo::RunIPCSim::buildShellIpcSimulation(config);

  ASSERT_NE(context.collisionHandler, nullptr);
  EXPECT_EQ(context.collisionHandler->getNumSurfaceVertices(), context.surfaceMesh.numVertices() + 4);
  EXPECT_EQ(context.collisionHandler->getNumSurfaceTriangles(), context.surfaceMesh.numTriangles() + 2);
  EXPECT_EQ(context.surfaceRestPositions.size(), context.surfaceMesh.numVertices() * 3);
  EXPECT_EQ(context.surfaceFromSimulationDispMap.rows(), context.surfaceMesh.numVertices() * 3);
  for (int vi = 0; vi < context.surfaceMesh.numVertices(); ++vi)
    EXPECT_TRUE(context.collisionHandler->isSurfaceVertexDeformable(vi));
  for (int vi = context.surfaceMesh.numVertices(); vi < context.collisionHandler->getNumSurfaceVertices(); ++vi)
    EXPECT_FALSE(context.collisionHandler->isSurfaceVertexDeformable(vi));
}

TEST(RunIPCSimSetupGTest, VolumeExternalObjectPreservesDeformableEmbeddingColumns)
{
  initializeRunIPCSimTestEnvironment();
  ScopedTempDir tempDir;
  const fs::path configPath = tempDir.path() / "tet-external.json";
  const fs::path externalPath = tempDir.path() / "external-plane.obj";
  writeTextFile(configPath, makeTetIPCConfig(tempDir.path(), 0));
  writeExternalPlane(externalPath, -0.2);

  pgo::ConfigFileJSON config;
  ASSERT_TRUE(config.open(configPath.string().c_str()));
  addIPCExternalObject(config, externalPath);
  const auto context = pgo::RunIPCSim::buildVolumeIpcSimulation(config);

  EXPECT_EQ(context.collisionHandler->getNumSurfaceVertices(), context.surfaceMesh.numVertices() + 4);
  EXPECT_EQ(context.collisionHandler->getNumSurfaceTriangles(), context.surfaceMesh.numTriangles() + 2);
  EXPECT_EQ(context.collisionHandler->getNumDOFs(), context.surfaceFromSimulationDispMap.cols());
  EXPECT_EQ(context.surfaceFromSimulationDispMap.rows(), context.surfaceMesh.numVertices() * 3);
}

TEST(RunIPCSimSetupGTest, ShellSupportsMultipleExternalObjects)
{
  initializeRunIPCSimTestEnvironment();
  ScopedTempDir tempDir;
  const fs::path configPath = tempDir.path() / "shell-multiple-external.json";
  const fs::path externalA = tempDir.path() / "external-a.obj";
  const fs::path externalB = tempDir.path() / "external-b.obj";
  writeTextFile(configPath, makeShellIPCConfig(tempDir.path(), 0));
  writeExternalPlane(externalA, -0.2);
  writeExternalPlane(externalB, -0.4);

  pgo::ConfigFileJSON config;
  ASSERT_TRUE(config.open(configPath.string().c_str()));
  config.handle()["ipc-external-objects"] = nlohmann::json::array(
    { { { "filename", externalA.generic_string() }, { "scale", 1.0 }, { "initial-translation", { 0.0, 0.0, 0.0 } }, { "movement", { 0.0, 0.0, 0.0 } } },
      { { "filename", externalB.generic_string() }, { "scale", 1.0 }, { "initial-translation", { 0.0, 0.0, 0.0 } }, { "movement", { 0.0, 0.0, 0.0 } } } });
  const auto context = pgo::RunIPCSim::buildShellIpcSimulation(config);
  EXPECT_EQ(context.collisionHandler->getNumSurfaceVertices(), context.surfaceMesh.numVertices() + 8);
  EXPECT_EQ(context.collisionHandler->getNumSurfaceTriangles(), context.surfaceMesh.numTriangles() + 4);
}

TEST(RunIPCSimSetupGTest, ShellSceneAndExternalScalesAreIndependentAndPrecedeInitialDisplacement)
{
  initializeRunIPCSimTestEnvironment();
  ScopedTempDir tempDir;
  const fs::path configPath = tempDir.path() / "shell-scaled-external.json";
  const fs::path externalPath = tempDir.path() / "external-plane.obj";
  writeTextFile(configPath, makeShellIPCConfig(tempDir.path(), 0));
  writeExternalPlane(externalPath, -0.25);

  pgo::ConfigFileJSON config;
  ASSERT_TRUE(config.open(configPath.string().c_str()));
  config.handle()["scale"] = 2.0;
  config.handle()["init-disp"] = { 0.1, -0.2, 0.3 };
  addIPCExternalObject(config, externalPath, nlohmann::json::array({ 0.0, 0.0, 0.0 }), 3.0,
    nlohmann::json::array({ 0.5, 1.0, -0.25 }));
  const auto context = pgo::RunIPCSim::buildShellIpcSimulation(config);

  pgo::Mesh::TriMeshGeo rawSurface;
  ASSERT_TRUE(rawSurface.load((fs::path(kShellExampleDir) / "shell.obj").string()));
  EXPECT_TRUE(context.surfaceMesh.pos(0).isApprox(rawSurface.pos(0) * 2.0, 1e-14));
  EXPECT_TRUE(context.initialDisplacement.segment<3>(0).isApprox(ES::V3d(0.1, -0.2, 0.3), 1e-14));
  EXPECT_TRUE(context.initialDisplacement.tail<3>().isApprox(ES::V3d(0.1, -0.2, 0.3), 1e-14));

  const ES::VXd &combinedRest = context.collisionHandler->getSurfaceRestPositions();
  const int externalOffset = context.surfaceMesh.numVertices() * 3;
  EXPECT_TRUE(combinedRest.segment<3>(externalOffset).isApprox(ES::V3d(-2.5, 0.25, -3.25), 1e-14));
}

TEST(RunIPCSimSetupGTest, ExternalObjectValidationRejectsUnsupportedInputs)
{
  initializeRunIPCSimTestEnvironment();
  ScopedTempDir tempDir;
  const fs::path configPath = tempDir.path() / "shell-external-errors.json";
  const fs::path validPath = tempDir.path() / "valid.obj";
  const fs::path emptyPath = tempDir.path() / "empty.obj";
  const fs::path repeatedPath = tempDir.path() / "repeated.obj";
  const fs::path zeroAreaPath = tempDir.path() / "zero-area.obj";
  const fs::path nonFiniteAreaPath = tempDir.path() / "non-finite-area.obj";
  writeTextFile(configPath, makeShellIPCConfig(tempDir.path(), 0));
  writeExternalPlane(validPath);
  writeTextFile(emptyPath, "v 0 0 0\n");
  writeTextFile(repeatedPath, "v 0 0 0\nv 1 0 0\nf 1 2 2\n");
  writeTextFile(zeroAreaPath, "v 0 0 0\nv 1 0 0\nv 2 0 0\nf 1 2 3\n");
  writeTextFile(nonFiniteAreaPath,
    "v 1e308 0 0\nv 0 1e308 0\nv 0 0 1e308\nf 1 2 3\n");

  pgo::ConfigFileJSON config;
  ASSERT_TRUE(config.open(configPath.string().c_str()));

  config.handle()["ipc-external-objects"] = nlohmann::json::array(
    { { { "filename", validPath.generic_string() }, { "movement", { 0.0, 0.0, 0.0 } } } });
  EXPECT_NO_THROW(pgo::RunIPCSim::buildShellIpcSimulation(config));
  addIPCExternalObject(config, validPath, nlohmann::json::array({ 0.0, 0.0, 0.0 }), 0.0);
  EXPECT_THROW(pgo::RunIPCSim::buildShellIpcSimulation(config), std::invalid_argument);
  addIPCExternalObject(config, validPath, nlohmann::json::array({ 0.0, 0.0, 0.0 }), "bad");
  EXPECT_THROW(pgo::RunIPCSim::buildShellIpcSimulation(config), std::invalid_argument);
  addIPCExternalObject(config, validPath, nlohmann::json::array({ 0.0, 0.0, 0.0 }), std::numeric_limits<double>::infinity());
  EXPECT_THROW(pgo::RunIPCSim::buildShellIpcSimulation(config), std::invalid_argument);

  config.handle()["ipc-external-objects"] = nlohmann::json::array(
    { { { "filename", validPath.generic_string() }, { "scale", 1.0 }, { "movement", { 0.0, 0.0, 0.0 } } } });
  EXPECT_NO_THROW(pgo::RunIPCSim::buildShellIpcSimulation(config));
  addIPCExternalObject(config, validPath, nlohmann::json::array({ 0.0, 0.0, 0.0 }), 1.0,
    nlohmann::json::array({ 0.0, 0.0 }));
  EXPECT_THROW(pgo::RunIPCSim::buildShellIpcSimulation(config), std::invalid_argument);
  addIPCExternalObject(config, validPath, nlohmann::json::array({ 0.0, 0.0, 0.0 }), 1.0,
    nlohmann::json::array({ 0.0, "bad", 0.0 }));
  EXPECT_THROW(pgo::RunIPCSim::buildShellIpcSimulation(config), std::invalid_argument);
  addIPCExternalObject(config, validPath, nlohmann::json::array({ 0.0, 0.0, 0.0 }), 1.0,
    nlohmann::json::array({ 0.0, std::numeric_limits<double>::infinity(), 0.0 }));
  EXPECT_THROW(pgo::RunIPCSim::buildShellIpcSimulation(config), std::invalid_argument);

  addIPCExternalObject(config, validPath, nlohmann::json::array({ 0.0, 0.1, 0.0 }));
  EXPECT_THROW(pgo::RunIPCSim::buildShellIpcSimulation(config), std::invalid_argument);
  addIPCExternalObject(config, validPath, nlohmann::json::array({ 0.0, 0.0 }));
  EXPECT_THROW(pgo::RunIPCSim::buildShellIpcSimulation(config), std::invalid_argument);
  addIPCExternalObject(config, validPath, nlohmann::json::array({ 0.0, "bad", 0.0 }));
  EXPECT_THROW(pgo::RunIPCSim::buildShellIpcSimulation(config), std::invalid_argument);
  addIPCExternalObject(config, validPath, nlohmann::json::array({ 0.0, std::numeric_limits<double>::infinity(), 0.0 }));
  EXPECT_THROW(pgo::RunIPCSim::buildShellIpcSimulation(config), std::invalid_argument);
  addIPCExternalObject(config, tempDir.path() / "missing.obj");
  EXPECT_THROW(pgo::RunIPCSim::buildShellIpcSimulation(config), std::invalid_argument);
  addIPCExternalObject(config, emptyPath);
  EXPECT_THROW(pgo::RunIPCSim::buildShellIpcSimulation(config), std::invalid_argument);
  addIPCExternalObject(config, repeatedPath);
  EXPECT_THROW(pgo::RunIPCSim::buildShellIpcSimulation(config), std::invalid_argument);
  addIPCExternalObject(config, zeroAreaPath);
  EXPECT_THROW(pgo::RunIPCSim::buildShellIpcSimulation(config), std::invalid_argument);
  addIPCExternalObject(config, nonFiniteAreaPath);
  EXPECT_THROW(pgo::RunIPCSim::buildShellIpcSimulation(config), std::invalid_argument);

  config.handle().erase("ipc-external-objects");
  config.handle()["external-objects"] = nlohmann::json::array();
  EXPECT_THROW(pgo::RunIPCSim::buildShellIpcSimulation(config), std::invalid_argument);
  config.handle().erase("external-objects");

  config.handle()["ipc-external-objects"] = "not-an-array";
  EXPECT_THROW(pgo::RunIPCSim::buildShellIpcSimulation(config), std::invalid_argument);
  config.handle()["scale"] = 0.0;
  config.handle()["ipc-external-objects"] = nlohmann::json::array();
  EXPECT_THROW(pgo::RunIPCSim::buildShellIpcSimulation(config), std::invalid_argument);
}

TEST(RunIPCSimSetupGTest, ExternalObjectDefaultsAndTypedMovementArePreserved)
{
  initializeRunIPCSimTestEnvironment();
  ScopedTempDir tempDir;
  const fs::path configPath = tempDir.path() / "shell-external-defaults.json";
  const fs::path externalPath = tempDir.path() / "external-plane.obj";
  writeTextFile(configPath, makeShellIPCConfig(tempDir.path(), 0));
  writeExternalPlane(externalPath, -0.2);

  pgo::ConfigFileJSON config;
  ASSERT_TRUE(config.open(configPath.string().c_str()));
  config.handle()["ipc-external-objects"] = nlohmann::json::array(
    { { { "filename", externalPath.generic_string() } } });

  const auto specs = pgo::RunIPCSim::parseIPCExternalObjectSpecs(config);
  ASSERT_EQ(specs.size(), 1u);
  EXPECT_EQ(specs[0].path, externalPath);
  EXPECT_DOUBLE_EQ(specs[0].scale, 1.0);
  EXPECT_TRUE(specs[0].initialTranslation.isZero(0.0));
  EXPECT_TRUE(specs[0].movement.isZero(0.0));
  EXPECT_NO_THROW(pgo::RunIPCSim::buildShellIpcSimulation(config));
}

TEST(RunIPCSimCliGTest, StaticRestartRequestIsRejected)
{
  const fs::path binary = runIPCSimBinaryPath();
  ASSERT_FALSE(binary.empty());
  ASSERT_TRUE(fs::exists(binary));
  ScopedTempDir tempDir;
  const fs::path configPath = tempDir.path() / "shell-static-restart.json";
  std::string config = makeStaticConfig(makeShellIPCConfig(tempDir.path(), 1));
  config = addBoolConfigField(std::move(config), "restart-from-u", true);
  writeTextFile(configPath, config);

  EXPECT_NE(runCommand(shellExecutable(binary) + " " + quotePath(configPath)), 0);
}

TEST(RunIPCSimCliGTest, InvalidDynamicRestartStatesAreRejected)
{
  const fs::path binary = runIPCSimBinaryPath();
  ASSERT_FALSE(binary.empty());
  ASSERT_TRUE(fs::exists(binary));
  ScopedTempDir tempDir;
  const fs::path outputDir = tempDir.path() / "shell-output";
  const fs::path configPath = tempDir.path() / "shell-malformed-restart.json";
  fs::create_directories(outputDir);
  writeTextFile(configPath, addBoolConfigField(makeShellIPCConfig(tempDir.path(), 2), "restart-from-u", true));
  const std::string command = shellExecutable(binary) + " " + quotePath(configPath);
  const fs::path restartPath = outputDir / "deform0000.u";

  writeTextFile(restartPath, "not a matrix\n");
  EXPECT_NE(runCommand(command), 0);

  ES::MXd wrongSize = ES::MXd::Zero(2, 2);
  ASSERT_EQ(ES::writeMatrix(restartPath.string().c_str(), wrongSize), 0);
  EXPECT_NE(runCommand(command), 0);

  pgo::Mesh::TriMeshGeo mesh;
  ASSERT_TRUE(mesh.load((fs::path(kShellExampleDir) / "shell.obj").string()));
  ES::MXd nonFinite = ES::MXd::Zero(mesh.numVertices() * 3, 3);
  nonFinite(0, 0) = std::numeric_limits<double>::quiet_NaN();
  ASSERT_EQ(ES::writeMatrix(restartPath.string().c_str(), nonFinite), 0);
  EXPECT_NE(runCommand(command), 0);
}
