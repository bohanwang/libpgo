#include <gtest/gtest.h>

#include "EigenSupport.h"
#include "barycentricCoordinates.h"
#include "configFileJSON.h"
#include "cubicMesh.h"
#include "generateMassMatrix.h"
#include "runSimVolumeMeshIO.h"
#include "tetMesh.h"
#include "triMeshGeo.h"
#include "volumetricMesh.h"

#include <filesystem>
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

constexpr const char *kTetBoxVegPath = LIBPGO_TEST_TET_BOX_VEG;
constexpr const char *kTetBoxObjPath = LIBPGO_TEST_TET_BOX_OBJ;
constexpr const char *kCubicBoxVegPath = LIBPGO_TEST_CUBIC_BOX_VEG;
constexpr const char *kCubicBoxObjPath = LIBPGO_TEST_CUBIC_BOX_OBJ;

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
  config.handle()["external-objects"] = nlohmann::json::array({ { { "filename", "../../bottom.obj" }, { "movement", { 0.0, 0.0, 0.0 } } } });

  const ResolvedRunSimPaths paths = resolvePaths(config, cubicConfigPath());
  EXPECT_EQ(paths.surfaceMeshFilename, kCubicBoxObjPath);
  EXPECT_EQ(paths.outputPath, (cubicExampleDir() / "ret-cubic-box").string());
  ASSERT_EQ(paths.fixedVertexFilenames.size(), 1u);
  EXPECT_EQ(paths.fixedVertexFilenames[0], (cubicExampleDir() / "fixed.txt").string());
  ASSERT_EQ(paths.externalObjectFilenames.size(), 1u);
  EXPECT_EQ(paths.externalObjectFilenames[0], (cubicExampleDir() / "../../bottom.obj").lexically_normal().string());
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

TEST(RunSimVolumeMeshIOGTest, BuildsCommonPreprocessingForTetAndCubic)
{
  expectCommonPreprocessingWorks(tetConfigPath(), "tet-mesh", "box.veg", "box.obj", VolumetricMesh::TET);
  expectCommonPreprocessingWorks(cubicConfigPath(), "cubic-mesh", "box.veg", "box.obj", VolumetricMesh::CUBIC);
}
