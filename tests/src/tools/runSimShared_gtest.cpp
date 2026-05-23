#include <gtest/gtest.h>

#include "EigenSupport.h"
#include "barycentricCoordinates.h"
#include "configFileJSON.h"
#include "cubicMeshDeformationModel.h"
#include "cubicMesh.h"
#include "deformationModelAssembler.h"
#include "deformationModelEnergy.h"
#include "deformationModelManager.h"
#include "generateMassMatrix.h"
#include "pgoLogging.h"
#include "cli/cliLogging.h"
#include "setup/femSetup.h"
#include "io/volumeMeshIO.h"
#include "simulationMesh.h"
#include "tetMesh.h"
#include "legacy_penalty/triangleMeshExternalContactHandler.h"
#include "legacy_penalty/triangleMeshSelfContactHandler.h"
#include "triMeshGeo.h"
#include "volumetricMesh.h"

#include <chrono>
#include <filesystem>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>

namespace
{
namespace ES = pgo::EigenSupport;
namespace fs = std::filesystem;
using pgo::RunSim::VolumeMeshInputConfig;
using pgo::RunSim::ResolvedRunSimPaths;
using pgo::VolumetricMeshes::VolumetricMesh;
using pgo::SolidDeformationModel::CubicMeshDeformationModel;
using pgo::SolidDeformationModel::SimulationMeshType;

constexpr const char *kTetBoxVegPath = LIBPGO_TEST_TET_BOX_VEG;
constexpr const char *kTetBoxObjPath = LIBPGO_TEST_TET_BOX_OBJ;
constexpr const char *kCubicBoxVegPath = LIBPGO_TEST_CUBIC_BOX_VEG;
constexpr const char *kCubicBoxObjPath = LIBPGO_TEST_CUBIC_BOX_OBJ;
class ScopedTempDir
{
public:
  ScopedTempDir()
  {
    const auto base = fs::temp_directory_path();
    for (int attempt = 0; attempt < 32; ++attempt) {
      const auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
      path_ = base / ("libpgo-runSim-gtest-" + std::to_string(stamp) + "-" + std::to_string(std::rand()) + "-" + std::to_string(attempt));

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

fs::path tetExampleDir()
{
  return fs::path(kTetBoxVegPath).parent_path();
}

fs::path cubicExampleDir()
{
  return fs::path(kCubicBoxVegPath).parent_path();
}

std::string tetConfigPath()
{
  return (tetExampleDir() / "box.json").string();
}

std::string cubicConfigPath()
{
  return (cubicExampleDir() / "box.json").string();
}

VolumeMeshInputConfig parseConfig(const char *key, const char *filename, const std::string &configFilename)
{
  pgo::ConfigFileJSON config;
  if (!config.open(configFilename.c_str())) {
    throw std::runtime_error("Failed to open test config: " + configFilename);
  }
  config.handle().erase("tet-mesh");
  config.handle().erase("cubic-mesh");
  config.handle()[key] = filename;
  return pgo::RunSim::parseVolumeMeshInputConfig(config);
}

void expectInvalidArgumentContaining(const std::function<void()> &fn, const std::string &needle)
{
  try {
    fn();
    FAIL() << "Expected std::invalid_argument";
  }
  catch (const std::invalid_argument &err) {
    EXPECT_NE(std::string(err.what()).find(needle), std::string::npos) << err.what();
  }
}

void expectScaledBoundingBox(const VolumetricMesh &referenceMesh, const VolumetricMesh &scaledMesh, double scale)
{
  const auto referenceBB = referenceMesh.getBoundingBox();
  const auto scaledBB = scaledMesh.getBoundingBox();
  for (int axis = 0; axis < 3; axis++) {
    EXPECT_DOUBLE_EQ(scaledBB.bmin()[axis], referenceBB.bmin()[axis] * scale);
    EXPECT_DOUBLE_EQ(scaledBB.bmax()[axis], referenceBB.bmax()[axis] * scale);
  }
}

ResolvedRunSimPaths resolvePaths(const pgo::ConfigFileJSON &config, const std::string &configFilename)
{
  (void)configFilename;
  return pgo::RunSim::resolveRunSimPaths(config);
}

void expectCommonPreprocessingWorks(const std::string &configFilename, const char *meshKey, const char *meshPath,
  const char *surfacePath, VolumetricMesh::elementType expectedType)
{
  pgo::ConfigFileJSON config;
  ASSERT_TRUE(config.open(configFilename.c_str()));
  config.handle()[meshKey] = meshPath;
  config.handle()["surface-mesh"] = surfacePath;
  config.handle()["output"] = "ret-test";

  const VolumeMeshInputConfig meshConfig = pgo::RunSim::parseVolumeMeshInputConfig(config);
  const ResolvedRunSimPaths paths = resolvePaths(config, configFilename);
  std::unique_ptr<VolumetricMesh> volumetricMesh = pgo::RunSim::loadValidatedVolumeMesh(meshConfig, 1.0);
  ASSERT_NE(volumetricMesh, nullptr);
  ASSERT_EQ(volumetricMesh->getElementType(), expectedType);

  pgo::Mesh::TriMeshGeo surfaceMesh;
  ASSERT_TRUE(surfaceMesh.load(paths.surfaceMeshFilename));

  ES::VXd surfaceRestPositions(surfaceMesh.numVertices() * 3);
  for (int vi = 0; vi < surfaceMesh.numVertices(); vi++) {
    surfaceRestPositions.segment<3>(vi * 3) = surfaceMesh.pos(vi);
  }

  pgo::InterpolationCoordinates::BarycentricCoordinates bc(
    surfaceMesh.numVertices(), surfaceRestPositions.data(), volumetricMesh.get());
  EXPECT_EQ(bc.getNumLocations(), surfaceMesh.numVertices());
  EXPECT_EQ(bc.getNumElementVertices(), expectedType == VolumetricMesh::TET ? 4 : 8);

  const ES::SpMatD W = bc.generateInterpolationMatrix();
  EXPECT_EQ(W.rows(), surfaceMesh.numVertices() * 3);
  EXPECT_EQ(W.cols(), volumetricMesh->getNumVertices() * 3);
  EXPECT_GT(W.nonZeros(), 0);

  ES::SpMatD M;
  pgo::VolumetricMeshes::GenerateMassMatrix::computeMassMatrix(volumetricMesh.get(), M, true);
  EXPECT_EQ(M.rows(), volumetricMesh->getNumVertices() * 3);
  EXPECT_EQ(M.cols(), volumetricMesh->getNumVertices() * 3);
  EXPECT_GT(M.nonZeros(), 0);
}

void expectAllFinite(const ES::VXd &v)
{
  for (Eigen::Index i = 0; i < v.size(); i++) {
    EXPECT_TRUE(std::isfinite(v[i])) << "Non-finite vector entry at " << i;
  }
}

void expectAllFinite(const ES::SpMatD &m)
{
  for (Eigen::Index i = 0; i < m.nonZeros(); i++) {
    EXPECT_TRUE(std::isfinite(m.valuePtr()[i])) << "Non-finite sparse entry at " << i;
  }
}

struct ContactEmbeddingTestInput
{
  std::unique_ptr<VolumetricMesh> volumetricMesh;
  pgo::Mesh::TriMeshGeo surfaceMesh;
  std::vector<int> embeddingVertexIndices;
  std::vector<double> embeddingWeights;
  ES::SpMatD expectedEmbeddingMatrix;
  int embeddingArity = 0;
  int nDOFs = 0;
};

ContactEmbeddingTestInput makeContactEmbeddingTestInput(
  const std::string &configFilename,
  VolumetricMesh::elementType expectedType)
{
  pgo::Logging::init();

  pgo::ConfigFileJSON config;
  if (!config.open(configFilename.c_str())) {
    throw std::runtime_error("Failed to open test config: " + configFilename);
  }

  const double scale = config.getDouble("scale", 1);
  const VolumeMeshInputConfig meshConfig = pgo::RunSim::parseVolumeMeshInputConfig(config);

  ContactEmbeddingTestInput input;
  input.volumetricMesh = pgo::RunSim::loadValidatedVolumeMesh(meshConfig, scale);
  if (!input.volumetricMesh) {
    throw std::runtime_error("Failed to load volumetric mesh for contact embedding test.");
  }
  if (input.volumetricMesh->getElementType() != expectedType) {
    throw std::runtime_error("Unexpected volumetric mesh type in contact embedding test.");
  }

  const ResolvedRunSimPaths paths = resolvePaths(config, configFilename);
  if (!input.surfaceMesh.load(paths.surfaceMeshFilename)) {
    throw std::runtime_error("Failed to load surface mesh for contact embedding test.");
  }
  for (int vi = 0; vi < input.surfaceMesh.numVertices(); vi++) {
    input.surfaceMesh.pos(vi) *= scale;
  }

  ES::VXd surfaceRestPositions(input.surfaceMesh.numVertices() * 3);
  for (int vi = 0; vi < input.surfaceMesh.numVertices(); vi++) {
    surfaceRestPositions.segment<3>(vi * 3) = input.surfaceMesh.pos(vi);
  }

  pgo::InterpolationCoordinates::BarycentricCoordinates bc(
    input.surfaceMesh.numVertices(), surfaceRestPositions.data(), input.volumetricMesh.get());

  input.embeddingArity = bc.getNumElementVertices();
  input.embeddingVertexIndices = bc.getEmbeddingVertexIndices();
  input.embeddingWeights = bc.getEmbeddingWeights();
  input.expectedEmbeddingMatrix = bc.generateInterpolationMatrix();
  input.nDOFs = input.volumetricMesh->getNumVertices() * 3;

  return input;
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

TEST(RunSimVolumeMeshIOGTest, AcceptsLegacyTetMeshKey)
{
  const VolumeMeshInputConfig config = parseConfig("tet-mesh", "box.veg", tetConfigPath());
  EXPECT_EQ(config.configKey, "tet-mesh");
  EXPECT_EQ(config.meshFilename, kTetBoxVegPath);
  EXPECT_EQ(config.expectedElementType, VolumetricMesh::TET);

  pgo::VolumetricMeshes::TetMesh referenceMesh(kTetBoxVegPath);
  std::unique_ptr<VolumetricMesh> scaledMesh = pgo::RunSim::loadValidatedVolumeMesh(config, 2.0);
  ASSERT_NE(scaledMesh, nullptr);
  EXPECT_EQ(scaledMesh->getElementType(), VolumetricMesh::TET);
  expectScaledBoundingBox(referenceMesh, *scaledMesh, 2.0);
}

TEST(RunSimVolumeMeshIOGTest, AcceptsCubicMeshKey)
{
  const VolumeMeshInputConfig config = parseConfig("cubic-mesh", "box.veg", cubicConfigPath());
  EXPECT_EQ(config.configKey, "cubic-mesh");
  EXPECT_EQ(config.meshFilename, kCubicBoxVegPath);
  EXPECT_EQ(config.expectedElementType, VolumetricMesh::CUBIC);

  pgo::VolumetricMeshes::CubicMesh referenceMesh(kCubicBoxVegPath);
  std::unique_ptr<VolumetricMesh> scaledMesh = pgo::RunSim::loadValidatedVolumeMesh(config, 2.0);
  ASSERT_NE(scaledMesh, nullptr);
  EXPECT_EQ(scaledMesh->getElementType(), VolumetricMesh::CUBIC);
  expectScaledBoundingBox(referenceMesh, *scaledMesh, 2.0);
}

TEST(RunSimVolumeMeshIOGTest, RejectsMissingVolumeMeshKey)
{
  pgo::ConfigFileJSON config;
  ASSERT_TRUE(config.open(tetConfigPath().c_str()));
  config.handle().erase("tet-mesh");
  config.handle().erase("cubic-mesh");
  expectInvalidArgumentContaining(
    [&]() { (void)pgo::RunSim::parseVolumeMeshInputConfig(config); },
    "exactly one");
}

TEST(RunSimVolumeMeshIOGTest, RejectsDuplicateVolumeMeshKeys)
{
  pgo::ConfigFileJSON config;
  ASSERT_TRUE(config.open(tetConfigPath().c_str()));
  config.handle()["tet-mesh"] = kTetBoxVegPath;
  config.handle()["cubic-mesh"] = kCubicBoxVegPath;

  expectInvalidArgumentContaining(
    [&]() { (void)pgo::RunSim::parseVolumeMeshInputConfig(config); },
    "exactly one");
}

TEST(RunSimVolumeMeshIOGTest, RejectsTetFileBoundToCubicKey)
{
  const VolumeMeshInputConfig config = parseConfig("cubic-mesh", "box.veg", tetConfigPath());
  expectInvalidArgumentContaining(
    [&]() { (void)pgo::RunSim::loadValidatedVolumeMesh(config, 1.0); },
    "expects a cubic volumetric mesh");
}

TEST(RunSimVolumeMeshIOGTest, RejectsCubicFileBoundToTetKey)
{
  const VolumeMeshInputConfig config = parseConfig("tet-mesh", "box.veg", cubicConfigPath());
  expectInvalidArgumentContaining(
    [&]() { (void)pgo::RunSim::loadValidatedVolumeMesh(config, 1.0); },
    "expects a tet volumetric mesh");
}

TEST(RunSimVolumeMeshIOGTest, ResolvesCubicExamplePathsAgainstConfigDirectory)
{
  pgo::ConfigFileJSON config;
  ASSERT_TRUE(config.open(cubicConfigPath().c_str()));
  config.handle()["cubic-mesh"] = "box.veg";
  config.handle()["surface-mesh"] = "box.obj";
  config.handle()["output"] = "ret-cubic-box";
  config.handle()["fixed-vertices"] = nlohmann::json::array({ { { "filename", "fixed.txt" }, { "movement", { 0.0, 0.0, 0.0 } }, { "coeff", 1.0 } } });
  config.handle()["external-objects"] = nlohmann::json::array({ { { "filename", "../bottom.obj" }, { "movement", { 0.0, 0.0, 0.0 } } } });

  const ResolvedRunSimPaths paths = resolvePaths(config, cubicConfigPath());
  EXPECT_EQ(paths.surfaceMeshFilename, kCubicBoxObjPath);
  EXPECT_EQ(paths.outputPath, (cubicExampleDir() / "ret-cubic-box").string());
  ASSERT_EQ(paths.fixedVertexFilenames.size(), 1u);
  EXPECT_EQ(paths.fixedVertexFilenames[0], (cubicExampleDir() / "fixed.txt").string());
  ASSERT_EQ(paths.externalObjectFilenames.size(), 1u);
  EXPECT_EQ(paths.externalObjectFilenames[0], (cubicExampleDir() / "../bottom.obj").lexically_normal().string());
}

TEST(RunSimVolumeMeshIOGTest, ResolvesLegacyTetExamplePathsAgainstConfigDirectory)
{
  pgo::ConfigFileJSON config;
  ASSERT_TRUE(config.open(tetConfigPath().c_str()));
  config.handle()["tet-mesh"] = "box.veg";
  config.handle()["surface-mesh"] = "box.obj";
  config.handle()["output"] = "ret-box";
  config.handle()["external-objects"] = nlohmann::json::array({ { { "filename", "../bottom.obj" }, { "movement", { 0.0, 0.0, 0.0 } } } });

  const ResolvedRunSimPaths paths = resolvePaths(config, tetConfigPath());
  EXPECT_EQ(paths.surfaceMeshFilename, kTetBoxObjPath);
  EXPECT_EQ(paths.outputPath, (tetExampleDir() / "ret-box").string());
  ASSERT_EQ(paths.externalObjectFilenames.size(), 1u);
  EXPECT_EQ(paths.externalObjectFilenames[0], (tetExampleDir() / "../bottom.obj").lexically_normal().string());
}

TEST(RunSimCliLoggingGTest, DerivesDefaultLogPathFromConfigPath)
{
  const fs::path logPath = pgo::RunSim::deriveDefaultLogPathFromConfig(cubicExampleDir() / "bunny.json");
  EXPECT_EQ(logPath, cubicExampleDir() / "bunny.log");
}

TEST(RunSimCliLoggingGTest, RedirectsStdoutAndStderrToLogFile)
{
  const fs::path tempDir = fs::temp_directory_path() / "libpgo_runSim_cli_logging_gtest";
  fs::remove_all(tempDir);
  fs::create_directories(tempDir);

  const fs::path logPath = tempDir / "runSim.log";
  {
    pgo::RunSim::ScopedRunSimCliLogRedirect redirect(logPath.string());
    std::cout << "stdout redirect check" << std::endl;
    std::cerr << "stderr redirect check" << std::endl;
  }

  std::ifstream in(logPath);
  ASSERT_TRUE(in.is_open());
  const std::string contents((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
  EXPECT_NE(contents.find("stdout redirect check"), std::string::npos);
  EXPECT_NE(contents.find("stderr redirect check"), std::string::npos);

  std::cout << "stdout restored check" << std::endl;
  std::cerr << "stderr restored check" << std::endl;

  std::ifstream after(logPath);
  ASSERT_TRUE(after.is_open());
  const std::string afterContents((std::istreambuf_iterator<char>(after)), std::istreambuf_iterator<char>());
  EXPECT_EQ(afterContents.find("stdout restored check"), std::string::npos);
  EXPECT_EQ(afterContents.find("stderr restored check"), std::string::npos);
}

TEST(RunSimCliLoggingGTest, ResolveConfiguredLogLevelDefaultsToInfo)
{
  ScopedTempDir tempDir;
  const fs::path configPath = tempDir.path() / "loglevel-default.json";
  writeTextFile(configPath, "{\n  \"output\": \"ret\"\n}\n");

  pgo::ConfigFileJSON config;
  ASSERT_TRUE(config.open(configPath.string().c_str()));
  EXPECT_EQ(pgo::RunSim::resolveConfiguredLogLevel(config), spdlog::level::info);
}

TEST(RunSimCliLoggingGTest, ResolveConfiguredLogLevelParsesTraceAndWarn)
{
  ScopedTempDir tempDir;
  const fs::path traceConfigPath = tempDir.path() / "loglevel-trace.json";
  writeTextFile(traceConfigPath, "{\n  \"loglevel\": \"trace\",\n  \"output\": \"ret\"\n}\n");

  pgo::ConfigFileJSON traceConfig;
  ASSERT_TRUE(traceConfig.open(traceConfigPath.string().c_str()));
  EXPECT_EQ(pgo::RunSim::resolveConfiguredLogLevel(traceConfig), spdlog::level::trace);

  const fs::path warnConfigPath = tempDir.path() / "loglevel-warn.json";
  writeTextFile(warnConfigPath, "{\n  \"loglevel\": \"warn\",\n  \"output\": \"ret\"\n}\n");

  pgo::ConfigFileJSON warnConfig;
  ASSERT_TRUE(warnConfig.open(warnConfigPath.string().c_str()));
  EXPECT_EQ(pgo::RunSim::resolveConfiguredLogLevel(warnConfig), spdlog::level::warn);
}

TEST(RunSimVolumeMeshIOGTest, BuildsCommonPreprocessingForTetAndCubic)
{
  expectCommonPreprocessingWorks(tetConfigPath(), "tet-mesh", "box.veg", "box.obj", VolumetricMesh::TET);
  expectCommonPreprocessingWorks(cubicConfigPath(), "cubic-mesh", "box.veg", "box.obj", VolumetricMesh::CUBIC);
}

TEST(RunSimVolumeMeshIOGTest, InitializesCubicRuntimeMainPath)
{
  pgo::Logging::init();

  pgo::ConfigFileJSON config;
  ASSERT_TRUE(config.open(cubicConfigPath().c_str()));

  const VolumeMeshInputConfig meshConfig = pgo::RunSim::parseVolumeMeshInputConfig(config);
  std::unique_ptr<VolumetricMesh> volumetricMesh = pgo::RunSim::loadValidatedVolumeMesh(meshConfig, config.getDouble("scale", 1));
  ASSERT_NE(volumetricMesh, nullptr);
  ASSERT_EQ(volumetricMesh->getElementType(), VolumetricMesh::CUBIC);

  const auto initialized = pgo::RunSim::initializeVolumetricSimulation(
    *volumetricMesh, pgo::SolidDeformationModel::DeformationModelElasticMaterial::STABLE_NEO);

  ASSERT_NE(initialized.simMesh, nullptr);
  EXPECT_EQ(initialized.simMesh->getElementType(), SimulationMeshType::CUBIC);
  ASSERT_GT(initialized.simMesh->getNumElements(), 0);

  ASSERT_NE(initialized.dmm, nullptr);
  ASSERT_NE(initialized.assembler, nullptr);
  ASSERT_NE(initialized.elasticEnergy, nullptr);

  const auto *cubicFEM = dynamic_cast<const CubicMeshDeformationModel *>(initialized.dmm->getDeformationModel(0));
  ASSERT_NE(cubicFEM, nullptr);
  EXPECT_EQ(cubicFEM->getNumDOFs(), 24);

  EXPECT_EQ(initialized.assembler->getNumDOFs(), initialized.restPosition.size());
  EXPECT_EQ(initialized.plasticity.size(),
    initialized.simMesh->getNumElements() * initialized.dmm->getNumPlasticParameters());

  ES::VXd zero = ES::VXd::Zero(initialized.restPosition.size());
  ES::SpMatD hess;
  initialized.elasticEnergy->createHessian(hess);
  initialized.elasticEnergy->hessian(zero, hess);
  EXPECT_EQ(hess.rows(), initialized.restPosition.size());
  EXPECT_EQ(hess.cols(), initialized.restPosition.size());
  expectAllFinite(hess);

  ES::VXd grad = ES::VXd::Zero(initialized.assembler->getNumDOFs());
  initialized.assembler->computeGradient(initialized.restPosition.data(), initialized.plasticity.data(), nullptr, grad.data());
  EXPECT_EQ(grad.size(), initialized.assembler->getNumDOFs());
  expectAllFinite(grad);
}

TEST(RunSimVolumeMeshIOGTest, InitializeVolumetricSimulationCanDisableMaterialMaxStep)
{
  pgo::Logging::init();

  pgo::ConfigFileJSON config;
  ASSERT_TRUE(config.open(tetConfigPath().c_str()));

  const VolumeMeshInputConfig meshConfig = pgo::RunSim::parseVolumeMeshInputConfig(config);
  std::unique_ptr<VolumetricMesh> volumetricMesh = pgo::RunSim::loadValidatedVolumeMesh(meshConfig, config.getDouble("scale", 1));
  ASSERT_NE(volumetricMesh, nullptr);

  const auto initialized = pgo::RunSim::initializeVolumetricSimulation(
    *volumetricMesh,
    pgo::SolidDeformationModel::DeformationModelElasticMaterial::STABLE_NEO,
    pgo::SolidDeformationModel::DeformationModelPlasticMaterial::VOLUMETRIC_DOF6,
    false);

  ASSERT_NE(initialized.elasticEnergy, nullptr);
  EXPECT_FALSE(initialized.elasticEnergy->isMaterialMaxStepEnabled());
}

TEST(RunSimVolumeMeshIOGTest, TetExternalContactEmbeddingMatchesBarycentricInterpolation)
{
  const auto input = makeContactEmbeddingTestInput(tetConfigPath(), VolumetricMesh::TET);
  ASSERT_EQ(input.embeddingArity, 4);

  const std::vector<pgo::Mesh::TriMeshRef> noExternalSurfaces;
  pgo::Contact::TriangleMeshExternalContactHandler handler(
    input.surfaceMesh.positions(), input.surfaceMesh.triangles(), input.nDOFs,
    noExternalSurfaces, 1, &input.embeddingVertexIndices, &input.embeddingWeights);

  expectSparseMatrixNear(handler.getSampleEmbeddingMatrix(), input.expectedEmbeddingMatrix);
}

TEST(RunSimVolumeMeshIOGTest, TetSelfContactEmbeddingMatchesBarycentricInterpolation)
{
  const auto input = makeContactEmbeddingTestInput(tetConfigPath(), VolumetricMesh::TET);
  ASSERT_EQ(input.embeddingArity, 4);

  pgo::Contact::TriangleMeshSelfContactHandler handler(
    input.surfaceMesh.positions(), input.surfaceMesh.triangles(), input.nDOFs,
    1, &input.embeddingVertexIndices, &input.embeddingWeights);

  expectSparseMatrixNear(handler.getSampleEmbeddingMatrix(), input.expectedEmbeddingMatrix);
}

TEST(RunSimVolumeMeshIOGTest, CubicExternalContactEmbeddingMatchesBarycentricInterpolation)
{
  const auto input = makeContactEmbeddingTestInput(cubicConfigPath(), VolumetricMesh::CUBIC);
  ASSERT_EQ(input.embeddingArity, 8);

  const std::vector<pgo::Mesh::TriMeshRef> noExternalSurfaces;
  pgo::Contact::TriangleMeshExternalContactHandler handler(
    input.surfaceMesh.positions(), input.surfaceMesh.triangles(), input.nDOFs,
    noExternalSurfaces, 1, &input.embeddingVertexIndices, &input.embeddingWeights);

  expectSparseMatrixNear(handler.getSampleEmbeddingMatrix(), input.expectedEmbeddingMatrix);
}

TEST(RunSimVolumeMeshIOGTest, CubicSelfContactEmbeddingMatchesBarycentricInterpolation)
{
  const auto input = makeContactEmbeddingTestInput(cubicConfigPath(), VolumetricMesh::CUBIC);
  ASSERT_EQ(input.embeddingArity, 8);

  pgo::Contact::TriangleMeshSelfContactHandler handler(
    input.surfaceMesh.positions(), input.surfaceMesh.triangles(), input.nDOFs,
    1, &input.embeddingVertexIndices, &input.embeddingWeights);

  expectSparseMatrixNear(handler.getSampleEmbeddingMatrix(), input.expectedEmbeddingMatrix);
}

TEST(RunSimVolumeMeshIOGTest, ExternalContactRejectsMalformedEmbeddingArrays)
{
  auto input = makeContactEmbeddingTestInput(tetConfigPath(), VolumetricMesh::TET);
  ASSERT_FALSE(input.embeddingWeights.empty());
  input.embeddingWeights.pop_back();

  const std::vector<pgo::Mesh::TriMeshRef> noExternalSurfaces;
  expectInvalidArgumentContaining([&]() {
    pgo::Contact::TriangleMeshExternalContactHandler handler(
      input.surfaceMesh.positions(), input.surfaceMesh.triangles(), input.nDOFs,
      noExternalSurfaces, 1, &input.embeddingVertexIndices, &input.embeddingWeights);
  }, "embedding");
}

TEST(RunSimVolumeMeshIOGTest, SelfContactRejectsMalformedEmbeddingArrays)
{
  auto input = makeContactEmbeddingTestInput(tetConfigPath(), VolumetricMesh::TET);
  ASSERT_FALSE(input.embeddingWeights.empty());
  input.embeddingWeights.pop_back();

  expectInvalidArgumentContaining([&]() {
    pgo::Contact::TriangleMeshSelfContactHandler handler(
      input.surfaceMesh.positions(), input.surfaceMesh.triangles(), input.nDOFs,
      1, &input.embeddingVertexIndices, &input.embeddingWeights);
  }, "embedding");
}
