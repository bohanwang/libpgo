#include "geometryQuery.h"
#include "triMeshGeo.h"
#include "triMeshNeighbor.h"
#include "pgoLogging.h"
#include "boundingVolumeTree.h"

#include <sstream>
#include <iomanip>
#include <numeric>

#include <fmt/format.h>

std::string printVec3(const pgo::Vec3d &v)
{
  return fmt::format("{:12.8f},{:12.8f},{:12.8f}", v[0], v[1], v[2]);
}

std::string printMat3(const pgo::Mat3d &m)
{
  return fmt::format("{:12.8f},{:12.8f},{:12.8f}\n{:12.8f},{:12.8f},{:12.8f}\n{:12.8f},{:12.8f},{:12.8f}\n",
    m(0, 0), m(0, 1), m(0, 2),
    m(1, 0), m(1, 1), m(1, 2),
    m(2, 0), m(2, 1), m(2, 2));
}

int main(int argc, char *argv[])
{
  using namespace pgo;
  using namespace Mesh;
  namespace ES = EigenSupport;

  pgo::Logging::init();

  TriMeshGeo mesh;
  if (!mesh.load(argv[1])) {
    return 1;
  }
  mesh.save("a.obj");

  TriMeshNeighbor meshNeighbor(mesh);

  double surfaceArea = 0;
  for (int ti = 0; ti < mesh.numTriangles(); ti++) {
    surfaceArea += getTriangleArea(mesh.pos(ti, 0), mesh.pos(ti, 1), mesh.pos(ti, 2));
  }
  std::cout << "Surf A" << std::endl;
  std::cout << surfaceArea << std::endl;

  std::cout << "Tri N" << std::endl;
  int ti = 0;
  std::cout << printVec3(getTriangleScaledNormal(mesh.pos(ti, 0), mesh.pos(ti, 1), mesh.pos(ti, 2))) << std::endl;
  std::cout << printVec3(getTriangleNormal(mesh.pos(ti, 0), mesh.pos(ti, 1), mesh.pos(ti, 2))) << std::endl;

  std::cout << "Tri M" << std::endl;
  std::cout << printVec3(getTriangleCenterOfMass(mesh.pos(ti, 0), mesh.pos(ti, 1), mesh.pos(ti, 2))) << std::endl;

  std::cout << "Tri Angle" << std::endl;
  std::cout << getTriangleAngle(mesh.pos(ti, 0), mesh.pos(ti, 1), mesh.pos(ti, 2)) << std::endl;
  std::cout << getTriangleAngleRobust(mesh.pos(ti, 0), mesh.pos(ti, 1), mesh.pos(ti, 2)) << std::endl;

  std::cout << "Plane Dist:" << std::endl;
  std::cout << getScaledSignedDistanceToTrianglePlane(mesh.pos(mesh.numVertices() - 1), mesh.pos(ti, 0), mesh.pos(ti, 1), mesh.pos(ti, 2)) << std::endl;
  Vec3d n = getTriangleScaledNormal(mesh.pos(ti, 0), mesh.pos(ti, 1), mesh.pos(ti, 2));
  std::cout << printVec3(getClosestPointToPlaneWithScaledNormal(mesh.pos(mesh.numVertices() - 1), n, mesh.pos(ti, 0))) << std::endl;

  n.normalize();
  std::cout << printVec3(getClosestPointToPlaneWithNormal(mesh.pos(mesh.numVertices() - 1), n, mesh.pos(ti, 0))) << std::endl;

  std::cout << "Edge Angle" << std::endl;
  int vi = 0;
  int tj = meshNeighbor.getTriangleNeighbors(ti)[vi];
  PGO_ALOG(tj >= 0);

  Vec3i tri0 = mesh.tri(ti), tri1 = mesh.tri(tj);
  PGO_ALOG(getInvertedIndex(tri1, tri0[(vi + 2) % 3]) < 0);
  PGO_ALOG(getInvertedIndex(tri1, tri0[(vi + 0) % 3]) >= 0);
  PGO_ALOG(getInvertedIndex(tri1, tri0[(vi + 1) % 3]) >= 0);
  int vj = firstNonOverlapIndex(tri1, tri0);
  PGO_ALOG(vj >= 0);
  std::cout << getTwoTriangleDihedralAngle(mesh.pos(ti, (vi + 0) % 3), mesh.pos(ti, (vi + 1) % 3), mesh.pos(ti, (vi + 2) % 3), mesh.pos(tj, vj)) << std::endl;

  std::cout << "Proj Vec" << std::endl;
  std::cout << printVec3(getProjectedVectorToPlaneWithNormal(mesh.pos(0) - mesh.pos(mesh.numVertices() - 1), n)) << std::endl;

  std::cout << "Seg Dist" << std::endl;
  std::cout << getSquaredDistanceToLineSegment(mesh.pos(mesh.numVertices() / 2), mesh.pos(0), mesh.pos(mesh.numVertices() - 1)) << std::endl;
  std::cout << printVec3(getClosestPointToLineSegment(mesh.pos(mesh.numVertices() / 2), mesh.pos(0), mesh.pos(mesh.numVertices() - 1))) << std::endl;
  Vec2d w;
  std::cout << printVec3(getClosestPointToLineSegment(mesh.pos(mesh.numVertices() / 2), mesh.pos(0), mesh.pos(mesh.numVertices() - 1), w)) << std::endl;
  std::cout << printVec3(Vec3d(w[0], w[1], 0)) << std::endl;

  std::cout << "Line Dist" << std::endl;
  std::cout << printVec3(getClosestPointToLine(mesh.pos(mesh.numVertices() / 2), mesh.pos(0), mesh.pos(mesh.numVertices() - 1))) << std::endl;

  std::cout << "Bary W" << std::endl;
  std::cout << printVec3(getBarycentricWeightProjectedOnTrianglePlane(mesh.pos(mesh.numVertices() / 2), mesh.pos(ti, 0), mesh.pos(ti, 1), mesh.pos(ti, 2))) << std::endl;

  std::cout << "Seg Seg dist" << std::endl;
  std::cout << minimalDistance2OfTwoLineSegments(mesh.pos(mesh.numVertices() / 2), mesh.pos(ti, 0), mesh.pos(ti, 1), mesh.pos(ti, 2)) << std::endl;

  std::cout << "tri dist" << std::endl;
  int f;
  Vec3d cpt;
  Vec3d w3;
  std::cout << getSquaredDistanceToTriangle(mesh.pos(mesh.numVertices() / 2), mesh.pos(ti, 0), mesh.pos(ti, 1), mesh.pos(ti, 2), f) << std::endl;
  std::cout << getSquaredDistanceToTriangle(mesh.pos(mesh.numVertices() / 2), mesh.pos(ti, 0), mesh.pos(ti, 1), mesh.pos(ti, 2), f, cpt, w3) << std::endl;
  std::cout << getSquaredDistanceToTriangle(mesh.pos(mesh.numVertices() / 2), mesh.pos(ti, 0), mesh.pos(ti, 1), mesh.pos(ti, 2)) << std::endl;
  std::cout << printVec3(cpt) << std::endl;
  std::cout << printVec3(w3) << std::endl;
  std::cout << f << std::endl;

  std::cout << printVec3(getClosestPointToTriangleWithFeature(mesh.pos(mesh.numVertices() / 2), mesh.pos(ti, 0), mesh.pos(ti, 1), mesh.pos(ti, 2), f)) << std::endl;
  std::cout << printVec3(getClosestPointToTriangleWithNormalAndFeature(mesh.pos(mesh.numVertices() / 2), mesh.pos(ti, 0), mesh.pos(ti, 1), mesh.pos(ti, 2), n, f)) << std::endl;

  std::cout << "BB" << std::endl;
  std::vector<int> idx(100);
  std::iota(idx.begin(), idx.end(), 10);
  Vec3d bmin, bmax;
  getBoundingBoxFromSelectedVertices((int)idx.size(), idx.data(), mesh.positions().data(), bmin, bmax);
  std::cout << printVec3(bmin) << std::endl;
  std::cout << printVec3(bmax) << std::endl;

  TriMeshBVTree bvTree;
  bvTree.buildByInertiaPartition(mesh);
  auto ret = bvTree.closestTriangleQuery(mesh, Vec3d(0, 0, 0));
  std::cout << "bvtree\n";
  std::cout << ret.closestPosition << ',' << ret.triID << ',' << ret.feature << ',' << printVec3(ret.triBaryWeight) << std::endl;

  Vec3d tetV[4] = {
    Vec3d(1, 0, 0),
    Vec3d(0, 0, 0),
    Vec3d(0, 0, 1),
    Vec3d(0, 1, 0),
  };

  std::cout << "Tet" << std::endl;
  std::cout << getTetDeterminant(tetV[0], tetV[1], tetV[2], tetV[3]) << std::endl;
  std::cout << getTetVolume(tetV[0], tetV[1], tetV[2], tetV[3]) << std::endl;
  std::cout << getTetSignedVolume(tetV[0], tetV[1], tetV[2], tetV[3]) << std::endl;
  std::cout << printVec3(getTetCenterOfMass(tetV[0], tetV[1], tetV[2], tetV[3])) << std::endl;

  double w4[4];
  getTetBarycentricWeights(Vec3d(0.1, 0.1, 0.1), tetV[0], tetV[1], tetV[2], tetV[3], w4);
  std::cout << w4[0] << ',' << w4[1] << ',' << w4[2] << ',' << w4[3] << std::endl;
  std::cout << getSquaredDistanceToTet(Vec3d(1.1, 1.1, 1.1), tetV[0], tetV[1], tetV[2], tetV[3]) << std::endl;

  std::cout << "Triangle Inertia" << std::endl;
  std::cout << printMat3(getPointInertiaTensor(Vec3d(1, 2, 3), 0.5)) << std::endl;

  double triTensor[6];
  getUnitMassTriangleInertiaTensorVectorAroundOrigin(mesh.pos(ti, 0), mesh.pos(ti, 1), mesh.pos(ti, 2), triTensor);
  std::cout << triTensor[0] << ',' << triTensor[1] << ',' << triTensor[2] << ',' << triTensor[3] << ',' << triTensor[4] << ',' << triTensor[5] << std::endl;
  double ar = getTriangleArea(mesh.pos(ti, 0), mesh.pos(ti, 1), mesh.pos(ti, 2));
  std::cout << "==\n";
  std::cout << printMat3(getTriangleInertiaTensorAroundOrigin(mesh.pos(ti, 0), mesh.pos(ti, 1), mesh.pos(ti, 2), ar)) << std::endl;
  std::cout << "==\n";
  std::cout << printMat3(getTriangleInertiaTensorAroundOrigin(mesh.pos(ti, 0), mesh.pos(ti, 1), mesh.pos(ti, 2))) << std::endl;

  std::cout << "==\n";
  Mat3d comI = getTriangleInertiaTensorAroundCOM(mesh.pos(ti, 0), mesh.pos(ti, 1), mesh.pos(ti, 2));
  Mat3d refI = shiftInertiaTensorAroundMassCenterToReferencePoint(comI, ar, getTriangleCenterOfMass(mesh.pos(ti, 0), mesh.pos(ti, 1), mesh.pos(ti, 2)));
  std::cout << "==\n";
  std::cout << printMat3(refI) << std::endl;

  std::cout << "Triangle Inertia" << std::endl;
  std::cout << printMat3(getTetInertiaTensorAroundCOM(tetV[0], tetV[1], tetV[2], tetV[3])) << std::endl;
  std::cout << "==\n";
  std::cout << printMat3(getTetInertiaTensorAroudOrigin(tetV[0], tetV[1], tetV[2], tetV[3], getTetDeterminant(tetV[0], tetV[1], tetV[2], tetV[3]))) << std::endl;

  Vec3d side(5, 6, 7);
  Mat3d bI = getBoxInertiaTensorAroundCOM(side, side[0] * side[1] * side[2]);
  std::cout << "==\n";
  std::cout << printMat3(bI) << std::endl;

  return 0;
}