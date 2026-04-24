#include <gtest/gtest.h>

#include "simulationMesh.h"
#include "cubicMesh.h"
#include "volumetricMeshENuMaterial.h"

#include <memory>

namespace
{
constexpr const char *kCubicBoxVegPath = LIBPGO_TEST_CUBIC_BOX_VEG;
}

TEST(SimulationMeshGTest, LoadsCubicMeshFromExampleFile)
{
  pgo::VolumetricMeshes::CubicMesh cubicMesh(kCubicBoxVegPath);
  std::unique_ptr<pgo::SolidDeformationModel::SimulationMesh> simMesh(
    pgo::SolidDeformationModel::loadCubicMesh(&cubicMesh));

  ASSERT_NE(simMesh, nullptr);
  EXPECT_EQ(simMesh->getElementType(), pgo::SolidDeformationModel::SimulationMeshType::CUBIC);
  EXPECT_EQ(simMesh->getNumVertices(), cubicMesh.getNumVertices());
  EXPECT_EQ(simMesh->getNumElements(), cubicMesh.getNumElements());
  EXPECT_EQ(simMesh->getNumElementVertices(), 8);

  for (int j = 0; j < 8; j++) {
    EXPECT_EQ(simMesh->getVertexIndex(0, j), cubicMesh.getVertexIndex(0, j));
  }

  double simPos[3];
  simMesh->getVertex(0, simPos);
  const pgo::Vec3d cubicPos = cubicMesh.getVertex(0);
  EXPECT_DOUBLE_EQ(simPos[0], cubicPos[0]);
  EXPECT_DOUBLE_EQ(simPos[1], cubicPos[1]);
  EXPECT_DOUBLE_EQ(simPos[2], cubicPos[2]);

  const auto *simMat = dynamic_cast<const pgo::SolidDeformationModel::SimulationMeshENuMaterial *>(
    simMesh->getElementMaterial(0, 0));
  ASSERT_NE(simMat, nullptr);

  const auto *cubicMat = pgo::VolumetricMeshes::downcastENuMaterial(cubicMesh.getElementMaterial(0));
  ASSERT_NE(cubicMat, nullptr);
  EXPECT_DOUBLE_EQ(simMat->getE(), cubicMat->getE());
  EXPECT_DOUBLE_EQ(simMat->getNu(), cubicMat->getNu());
}
