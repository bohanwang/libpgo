#include <gtest/gtest.h>

#include "ipc/assembly/ipcCollisionSurfaceAssembler.h"

#include <stdexcept>
#include <vector>

namespace
{
namespace ES = pgo::EigenSupport;
using pgo::Contact::CIPC::IPCExternalSurface;
using pgo::Contact::CIPC::assembleIPCCollisionSurface;

ES::SpMatD makeDeformableMap()
{
  std::vector<ES::TripletD> triplets = {
    { 0, 0, 1.0 },
    { 4, 1, 2.0 },
    { 8, 0, -0.5 },
  };
  ES::SpMatD map(9, 2);
  map.setFromTriplets(triplets.begin(), triplets.end());
  return map;
}

IPCExternalSurface makeExternalSurface()
{
  IPCExternalSurface surface;
  surface.vertices.resize(3, 3);
  surface.vertices << 0.0, 0.0, 1.0,
    1.0, 0.0, 1.0,
    0.0, 1.0, 1.0;
  surface.triangles.resize(1, 3);
  surface.triangles << 0, 1, 2;
  return surface;
}
}  // namespace

TEST(IPCCollisionSurfaceAssemblerGTest, ConcatenatesGeometryMapAndRoleMask)
{
  ES::MXd deformableVertices(3, 3);
  deformableVertices << 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0,
    0.0, 1.0, 0.0;
  ES::MXi deformableTriangles(1, 3);
  deformableTriangles << 0, 1, 2;
  const ES::SpMatD deformableMap = makeDeformableMap();

  const auto data = assembleIPCCollisionSurface(
    deformableVertices, deformableTriangles, deformableMap, { makeExternalSurface() });

  ASSERT_EQ(data.vertices.rows(), 6);
  ASSERT_EQ(data.triangles.rows(), 2);
  EXPECT_TRUE(data.vertices.topRows(3).isApprox(deformableVertices));
  EXPECT_EQ(data.triangles.row(1), (ES::V3i{ 3, 4, 5 }).transpose());
  EXPECT_EQ(data.displacementMap.rows(), 18);
  EXPECT_EQ(data.displacementMap.cols(), 2);
  EXPECT_TRUE(ES::MXd(data.displacementMap.topRows(9)).isApprox(ES::MXd(deformableMap)));
  EXPECT_EQ(data.displacementMap.bottomRows(9).nonZeros(), 0);
  EXPECT_EQ(data.vertexIsDeformable, (std::vector<uint8_t>{ 1, 1, 1, 0, 0, 0 }));
}

TEST(IPCCollisionSurfaceAssemblerGTest, RejectsInvalidExternalFaceIndex)
{
  ES::MXd deformableVertices = ES::MXd::Zero(3, 3);
  ES::MXi deformableTriangles(1, 3);
  deformableTriangles << 0, 1, 2;
  IPCExternalSurface external = makeExternalSurface();
  external.triangles(0, 2) = 3;

  EXPECT_THROW(
    assembleIPCCollisionSurface(
      deformableVertices, deformableTriangles, makeDeformableMap(), { external }),
    std::invalid_argument);
}
