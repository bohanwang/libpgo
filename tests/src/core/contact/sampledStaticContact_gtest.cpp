#include <gtest/gtest.h>

#include "createTriMesh.h"
#include "pointPenetrationEnergy.h"
#include "pgoLogging.h"
#include "triangleMeshExternalContactHandler.h"

#include <memory>
#include <vector>

namespace
{
namespace ES = pgo::EigenSupport;

TEST(SampledStaticContactGTest, ContactEnergyBuiltBeforePenetrationMissesLaterContact)
{
  pgo::Logging::init();

  const std::vector<ES::V3d> surfaceVertices = {
    ES::V3d(-0.4, -0.4, 1.1),
    ES::V3d(-0.4, 0.4, 1.1),
    ES::V3d(0.4, -0.4, 1.1),
  };
  // The normal points down, opposite the outward normal of the cube's top face.
  const std::vector<ES::V3i> surfaceTriangles = { ES::V3i(0, 1, 2) };

  const pgo::Mesh::TriMeshGeo externalMesh = pgo::Mesh::createBoxMesh(
    ES::V3d(-1.0, -1.0, -1.0), ES::V3d(1.0, 1.0, 1.0));
  const std::vector<pgo::Mesh::TriMeshRef> externalMeshes = { pgo::Mesh::TriMeshRef(externalMesh) };

  pgo::Contact::TriangleMeshExternalContactHandler handler(
    surfaceVertices, surfaceTriangles, 9, externalMeshes, 1);

  ES::VXd beforeContact = ES::VXd::Zero(9);
  ES::VXd afterContact = ES::VXd::Zero(9);
  for (int vertex = 0; vertex < 3; ++vertex)
    afterContact[3 * vertex + 2] = -0.2;

  handler.execute(beforeContact.data());
  ASSERT_EQ(handler.getNumCollidingSamples(), 0);
  auto staleEnergy = handler.buildContactEnergy();
  staleEnergy->setComputePosFunction([&surfaceVertices](const ES::V3d &u, ES::V3d &p, int dofStart) {
    p = surfaceVertices[static_cast<std::size_t>(dofStart / 3)] + u;
  });
  auto *staleBuffer = staleEnergy->allocateBuffer();
  staleEnergy->setBuffer(staleBuffer);
  EXPECT_DOUBLE_EQ(staleEnergy->func(afterContact), 0.0);
  staleEnergy->freeBuffer(staleBuffer);
  staleEnergy.reset();

  handler.execute(afterContact.data());
  ASSERT_GT(handler.getNumCollidingSamples(), 0);
  auto rebuiltEnergy = handler.buildContactEnergy();
  rebuiltEnergy->setComputePosFunction([&surfaceVertices](const ES::V3d &u, ES::V3d &p, int dofStart) {
    p = surfaceVertices[static_cast<std::size_t>(dofStart / 3)] + u;
  });
  auto *rebuiltBuffer = rebuiltEnergy->allocateBuffer();
  rebuiltEnergy->setBuffer(rebuiltBuffer);
  EXPECT_GT(rebuiltEnergy->func(afterContact), 0.0);
  rebuiltEnergy->freeBuffer(rebuiltBuffer);
}
}  // namespace
