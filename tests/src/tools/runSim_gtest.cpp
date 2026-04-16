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

#include <functional>
#include <memory>
#include <stdexcept>
#include <string>

namespace
{
namespace ES = pgo::EigenSupport;
using pgo::RunSim::VolumeMeshInputConfig;
using pgo::VolumetricMeshes::VolumetricMesh;

constexpr const char *kTetBoxVegPath = LIBPGO_TEST_TET_BOX_VEG;
constexpr const char *kTetBoxObjPath = LIBPGO_TEST_TET_BOX_OBJ;
constexpr const char *kCubicBoxVegPath = LIBPGO_TEST_CUBIC_BOX_VEG;
constexpr const char *kCubicBoxObjPath = LIBPGO_TEST_CUBIC_BOX_OBJ;

VolumeMeshInputConfig parseConfig(const char *key, const char *filename)
{
  pgo::ConfigFileJSON config;
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

void expectCommonPreprocessingWorks(const char *meshKey, const char *meshPath, const char *surfacePath, VolumetricMesh::elementType expectedType)
{
  const VolumeMeshInputConfig config = parseConfig(meshKey, meshPath);
  std::unique_ptr<VolumetricMesh> volumetricMesh = pgo::RunSim::loadValidatedVolumeMesh(config, 1.0);
  ASSERT_NE(volumetricMesh, nullptr);
  ASSERT_EQ(volumetricMesh->getElementType(), expectedType);

  pgo::Mesh::TriMeshGeo surfaceMesh;
  ASSERT_TRUE(surfaceMesh.load(surfacePath));

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
  const VolumeMeshInputConfig config = parseConfig("tet-mesh", kTetBoxVegPath);
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
  const VolumeMeshInputConfig config = parseConfig("cubic-mesh", kCubicBoxVegPath);
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
  expectInvalidArgumentContaining(
    [&]() { (void)pgo::RunSim::parseVolumeMeshInputConfig(config); },
    "exactly one");
}

TEST(RunSimVolumeMeshIOGTest, RejectsDuplicateVolumeMeshKeys)
{
  pgo::ConfigFileJSON config;
  config.handle()["tet-mesh"] = kTetBoxVegPath;
  config.handle()["cubic-mesh"] = kCubicBoxVegPath;

  expectInvalidArgumentContaining(
    [&]() { (void)pgo::RunSim::parseVolumeMeshInputConfig(config); },
    "exactly one");
}

TEST(RunSimVolumeMeshIOGTest, RejectsTetFileBoundToCubicKey)
{
  const VolumeMeshInputConfig config = parseConfig("cubic-mesh", kTetBoxVegPath);
  expectInvalidArgumentContaining(
    [&]() { (void)pgo::RunSim::loadValidatedVolumeMesh(config, 1.0); },
    "expects a cubic volumetric mesh");
}

TEST(RunSimVolumeMeshIOGTest, RejectsCubicFileBoundToTetKey)
{
  const VolumeMeshInputConfig config = parseConfig("tet-mesh", kCubicBoxVegPath);
  expectInvalidArgumentContaining(
    [&]() { (void)pgo::RunSim::loadValidatedVolumeMesh(config, 1.0); },
    "expects a tet volumetric mesh");
}

TEST(RunSimVolumeMeshIOGTest, BuildsCommonPreprocessingForTetAndCubic)
{
  expectCommonPreprocessingWorks("tet-mesh", kTetBoxVegPath, kTetBoxObjPath, VolumetricMesh::TET);
  expectCommonPreprocessingWorks("cubic-mesh", kCubicBoxVegPath, kCubicBoxObjPath, VolumetricMesh::CUBIC);
}
