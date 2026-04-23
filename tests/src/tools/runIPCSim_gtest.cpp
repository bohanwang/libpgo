#include <gtest/gtest.h>

#include "EigenSupport.h"
#include "barycentricCoordinates.h"
#include "configFileJSON.h"
#include "initPredicates.h"
#include "pgoLogging.h"
#include "runIPCSimSetup.h"
#include "runSimCliLogging.h"
#include "runSimVolumeMeshIO.h"
#include "triMeshGeo.h"
#include "volumetricMesh.h"

#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
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
  double ipcDhat = 0.002, double ipcKappa = 3000.0, int dumpInterval = 1)
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
  double ipcDhat = 0.002, double ipcKappa = 3000.0, int dumpInterval = 1)
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

std::string makeTetIPCConfig(const fs::path &tempDir, int numTimesteps,
  double scale = 1.0, bool includeIPCFields = true, bool ipcHeuristic = false, const std::string &material = "stable-neo",
  int dumpInterval = 1)
{
  return makeVolumeIPCConfig(tetIPCExampleDir(), "tet-mesh", tempDir / "tet-output", numTimesteps,
    scale, includeIPCFields, ipcHeuristic, material, 0.002, 3000.0, dumpInterval);
}

std::string makeCubicIPCConfig(const fs::path &tempDir, int numTimesteps,
  double scale = 1.0, bool includeIPCFields = true, bool ipcHeuristic = false, const std::string &material = "stable-neo",
  int dumpInterval = 1)
{
  return makeVolumeIPCConfig(cubicIPCExampleDir(), "cubic-mesh", tempDir / "cubic-output", numTimesteps,
    scale, includeIPCFields, ipcHeuristic, material, 0.002, 3000.0, dumpInterval);
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

TEST(RunIPCSimCliGTest, LogFlagWritesCliOutputNextToConfig)
{
  const fs::path binary = runIPCSimBinaryPath();
  ASSERT_FALSE(binary.empty());
  ASSERT_TRUE(fs::exists(binary));

  ScopedTempDir tempDir;
  const fs::path configPath = tempDir.path() / "shell-ipc.json";
  const fs::path logPath = tempDir.path() / "shell-ipc.log";

  writeTextFile(configPath, makeShellIPCConfig(tempDir.path(), 0));

  std::ostringstream command;
  command << shellExecutable(binary)
          << " --log "
          << quotePath(configPath);

  ASSERT_EQ(runCommand(command.str()), 0);
  ASSERT_TRUE(fs::exists(logPath));

  std::ifstream in(logPath);
  ASSERT_TRUE(in.is_open());
  const std::string contents((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
  EXPECT_NE(contents.find("ipc-dhat"), std::string::npos);
  EXPECT_NE(contents.find("ipc-kappa"), std::string::npos);
  EXPECT_NE(contents.find("eps_ee"), std::string::npos);
  EXPECT_NE(contents.find("slackness"), std::string::npos);
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
  const fs::path logPath = tempDir.path() / "shell-ipc-heuristic-override.log";

  writeTextFile(configPath, makeShellIPCConfig(tempDir.path(), 0, true, true, 9.0, 42.0));

  std::ostringstream command;
  command << shellExecutable(binary)
          << " --log "
          << quotePath(configPath);

  ASSERT_EQ(runCommand(command.str()), 0);
  ASSERT_TRUE(fs::exists(logPath));

  std::ifstream in(logPath);
  ASSERT_TRUE(in.is_open());
  const std::string contents((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
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
