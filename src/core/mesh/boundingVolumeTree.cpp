/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "mesh" library , Copyright (C) 2018 USC                               *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code authors: Yijing Li, Jernej Barbic                                *
 * http://www.jernejbarbic.com/vega                                      *
 *                                                                       *
 * Research: Jernej Barbic, Hongyi Xu, Yijing Li,                        *
 *           Danyong Zhao, Bohan Wang,                                   *
 *           Fun Shing Sin, Daniel Schroeder,                            *
 *           Doug L. James, Jovan Popovic                                *
 *                                                                       *
 * Funding: National Science Foundation, Link Foundation,                *
 *          Singapore-MIT GAMBIT Game Lab,                               *
 *          Zumberge Research and Innovation Fund at USC,                *
 *          Sloan Foundation, Okawa Foundation,                          *
 *          USC Annenberg Foundation                                     *
 *                                                                       *
 * This library is free software; you can redistribute it and/or         *
 * modify it under the terms of the BSD-style license that is            *
 * included with this library in the file LICENSE.txt                    *
 *                                                                       *
 * This library is distributed in the hope that it will be useful,       *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the file     *
 * LICENSE.TXT for more details.                                         *
 *                                                                       *
 *************************************************************************/

#include "predicates.h"
#include "boundingVolumeTree.h"
#include "geometryQuery.h"
#include "tribox3.h"

#include "triple.h"
#include "basicAlgorithms.h"
#include "containerHelper.h"
#include "pgoLogging.h"

#include <tbb/parallel_for.h>
#include <tbb/enumerable_thread_specific.h>

#include <iostream>
#include <vector>
#include <stack>
#include <cstring>
#include <cassert>
#include <algorithm>
#include <numeric>
#include <queue>
#include <fstream>

constexpr double eps = 1e-6;

namespace pgo
{
namespace Mesh
{

using namespace BasicAlgorithms;

template<class DistanceToBBType, class DistanceToElementType>
double BoundingVolumeTreeBase::nearestQuery(std::vector<std::tuple<double, double, int>> &nodeStack, DistanceToBBType toBB, DistanceToElementType toElement,
  double distanceHi, double distanceLo) const
{
  PGO_ALOG(nodes.empty() == false);
  PGO_ALOG(distanceLo <= distanceHi);
  // distanceHi stores the current closest distances
  // traverse the hierarchy
  double minDistance = distanceHi;

  // each entry in the stack is for one node: <higher/lower bound of distanceToNode, node ID>
  // vector<tuple<double, double, int>> nodeStack;
  auto nodeBound = toBB(nodes[0].bb);  // get the lower and higher bound of this BV node
  if (nodeBound.first == DBL_MAX)      // if there are some special situations resulting in early return
    return DBL_MAX;

  nodeStack.emplace_back(nodeBound.first, nodeBound.second, 0);

  while (nodeStack.size() > 0) {
    const auto &top = nodeStack.back();
    int nodeID = get<2>(top);
    double nodeLowerBound = get<0>(top), nodeHigherBound = get<1>(top);
    nodeStack.pop_back();
    const Node &node = nodes[nodeID];
    if (nodeLowerBound > minDistance)  // if this BV node's distance lower bound is larger than distanceHi
      continue;
    // note here we don't use if (nodeLowerBound >= minDistance), so in this way two node which share the same lower bound
    // as minDistance will be both visited

    // this node needs to be inspected further
    if (nodeHigherBound < minDistance)
      minDistance = nodeHigherBound;

    if (node.isLeaf()) {
      // leaf node, traverse all children
      for (int eleID : node.indices) {
        double dist = toElement(eleID, minDistance);
        if (dist <= minDistance) {
          minDistance = dist;
          if (minDistance <= distanceLo)  // lower bound reached
            return minDistance;
        }
      }
    }
    else  // non-leaf node
    {
      int lastEnd = nodeStack.size();
      for (int childID : node.childrenIDs) {
        if (childID < 0)
          continue;
        auto distBound = toBB(nodes[childID].bb);
        if (distBound.first > minDistance || distBound.first == DBL_MAX)  // this node will be skipped
          continue;
        nodeStack.emplace_back(distBound.first, distBound.second, childID);
      }
      sort(nodeStack.begin() + lastEnd, nodeStack.end(), std::greater<std::tuple<double, double, int>>());
    }
  }

  return minDistance;
}

template<class DistanceToBBType, class DistanceToElementType>
double BoundingVolumeTreeBase::nearestQuery(DistanceToBBType toBB, DistanceToElementType toElement,
  double distanceHi, double distanceLo) const
{
  std::vector<std::tuple<double, double, int>> nodeStack;
  return nearestQuery<DistanceToBB, DistanceToElement>(nodeStack, toBB, toElement, distanceHi, distanceLo);
}

/////////////////////////////////////////////////////////////
//                   VertexBVTree
/////////////////////////////////////////////////////////////

void VertexBVTree::clear()
{
  nodes.clear();
  numVertices = 0;
}

void VertexBVTree::buildOctree(const ArrayRef<Vec3d> vertices, int maxDepth, int maxNumVerticesPerNode)
{
  nodes.clear();
  numVertices = vertices.size();
  std::vector<int> rootIndices(vertices.size());
  std::iota(rootIndices.begin(), rootIndices.end(), 0);

  nodes.emplace_back(0, std::move(rootIndices));
  nodes[0].bb = BoundingBox(vertices);

  auto subdivide = [&](const std::vector<int> &elementList, const BoundingBox &bb,
                     std::vector<std::vector<int>> &childIDs, std::vector<BoundingBox> &subBBs) -> void {
    Vec3d center = bb.center();
    int numVertices = elementList.size();

    childIDs.resize(8);
    subBBs.resize(8);

    // for each child, generate a list of vertices that intersect the child bounding box
    for (int i = 0; i < numVertices; i++) {
      const Vec3d &v = vertices[elementList[i]];
      unsigned char code = 0;
      code += (v[0] > center[0]);
      code += ((v[1] > center[1]) << 1);
      code += ((v[2] > center[2]) << 2);
      PGO_ALOG(code < 8);
      childIDs[code].push_back(elementList[i]);
    }

    for (int i = 0; i < 8; i++) {
      if (childIDs[i].size() > 0) {
        subBBs[i] = BoundingBox(vertices.data(), childIDs[i]);
      }
    }
  };
  buildFromRoot(subdivide, maxDepth, maxNumVerticesPerNode);
}

void VertexBVTree::buildByInertiaPartition(const ArrayRef<Vec3d> vertices, int maxDepth, int maxNumVerticesPerNode)
{
  nodes.clear();
  numVertices = vertices.size();
  std::vector<int> rootIndices(vertices.size());
  std::iota(rootIndices.begin(), rootIndices.end(), 0);

  nodes.emplace_back(0, std::move(rootIndices));
  nodes[0].bb = BoundingBox(vertices);

  auto subdivide = [&](const std::vector<int> &elementList, const BoundingBox &bb,
                     std::vector<std::vector<int>> &childIDs, std::vector<BoundingBox> &subBBs) -> void {
    // PGO_ALOG(bb.verifyBox());
    elementsInertiaBinaryPartition(
      elementList, bb, childIDs, INT_MAX,
      [&](int) { return 1.0; },
      [&](int i) { return vertices[i]; },
      [&](int i, const Vec3d &reference) { return getPointInertiaTensor(vertices[i] - reference, 1.0); });

    subBBs.resize(2);
    for (int i = 0; i < sizei(childIDs); i++)
      subBBs[i] = BoundingBox(vertices.data(), childIDs[i]);
  };

  buildFromRoot(subdivide, maxDepth, maxNumVerticesPerNode);
}

void VertexBVTree::rangeQuery(const ArrayRef<Vec3d> vertices, const SimpleSphere &simpleSphere, std::vector<int> &vertexIDList) const
{
  PGO_ALOG(numVertices == vertices.size());
  auto toBB = [&](const BoundingBox &bb) -> bool {
    return simpleSphere.doesBoundingBoxIntersect(bb);
  };

  auto toEle = [&](int idx) -> bool {
    return simpleSphere.inside(vertices[idx]);
  };

  BoundingVolumeTreeBase::rangeQuery(toBB, toEle, vertexIDList);
}

void VertexBVTree::rangeQuery(const ArrayRef<Vec3d> vertices, const HalfSpace &halfSpace, std::vector<int> &vertexIDList) const
{
  PGO_ALOG(numVertices == vertices.size());
  auto toBB = [&](const BoundingBox &bb) -> bool {
    return halfSpace.intersect(bb);
  };

  auto toEle = [&](int idx) -> bool {
    return halfSpace.outside(vertices[idx]) <= 0;
  };

  BoundingVolumeTreeBase::rangeQuery(toBB, toEle, vertexIDList);
}

bool VertexBVTree::intersectAABB(const ArrayRef<Vec3d> vertices, const BoundingBox &bbQuery) const
{
  PGO_ALOG(numVertices == vertices.size());
  auto toBB = [&](const BoundingBox &bb) -> bool {
    return bbQuery.intersect(bb);
  };

  auto toEle = [&](int idx) -> bool {
    return bbQuery.checkInside(vertices[idx]);
  };

  return BoundingVolumeTreeBase::rangeQueryWithEarlyStop(toBB, toEle);
}

int VertexBVTree::getClosestVertex(const ArrayRef<Vec3d> vertices, const Vec3d &queryPosition) const
{
  PGO_ALOG(numVertices == vertices.size());
  int closestVtxID = -1;

  auto toBB = [&](const BoundingBox &bb) -> std::pair<double, double> {
    double dist2 = bb.distanceToPoint2(queryPosition) * (1 - eps);
    double furDist2 = bb.furthestDistanceToPoint2(queryPosition) * (1 + eps);
    return { dist2, furDist2 };
  };

  auto toElement = [&](int vtxID, double minDistance) -> double {
    double dist2 = (vertices[vtxID] - queryPosition).squaredNorm();
    // std::cout << vtxID << ',' << dist2 << ',' << minDistance << '\n';

    if (dist2 <= minDistance) {
      closestVtxID = vtxID;
    }
    return dist2;
  };

  nearestQuery(toBB, toElement);

  // PGO_ALOG(closestVtxID >= 0);
  return closestVtxID;
}

/////////////////////////////////////////////////////////////
//                 BoundingBoxBVTree
/////////////////////////////////////////////////////////////

void BoundingBoxBVTree::clear()
{
  numBBs = 0;
  nodes.clear();
}

void BoundingBoxBVTree::buildByInertiaPartition(const ArrayRef<BoundingBox> boundingBoxes, int maxDepth, int maxNumTrianglesPerNode)
{
  nodes.clear();
  numBBs = boundingBoxes.size();
  PGO_ALOG(numBBs > 0);

  std::vector<int> rootIndices(numBBs);
  std::iota(rootIndices.begin(), rootIndices.end(), 0);
  BoundingBox rootBB(numBBs, boundingBoxes.data());
  nodes.emplace_back(0, std::move(rootIndices));
  nodes[0].bb = rootBB;

  std::vector<double> bbVolumes(numBBs);
  std::vector<Mat3d> bbInertia(numBBs);
  for (int i = 0; i < numBBs; i++) {
    bbVolumes[i] = boundingBoxes[i].volume();
    bbInertia[i] = getBoxInertiaTensorAroundCOM(boundingBoxes[i].sides(), bbVolumes[i]);
  }

  auto subdivide = [&](const std::vector<int> &elementList, const BoundingBox &bb,
                     std::vector<std::vector<int>> &childIDs, std::vector<BoundingBox> &subBBs) -> void {
    // PGO_ALOG(bb.verifyBox());
    //    elementsInertiaPartition(elementList, bb, 4, childIDs,
    elementsInertiaBinaryPartition(
      elementList, bb, childIDs, INT_MAX,
      [&](int i) { return bbVolumes[i]; },
      [&](int i) { return boundingBoxes[i].center(); },
      [&](int i, const Vec3d &reference) { return shiftInertiaTensorAroundMassCenterToReferencePoint(bbInertia[i], bbVolumes[i],
                                             reference - boundingBoxes[i].center()); });
    subBBs.resize(2);
    for (int i = 0; i < sizei(childIDs); i++)
      subBBs[i] = BoundingBox(boundingBoxes.data(), childIDs[i]);
  };

  buildFromRoot(subdivide, maxDepth, maxNumTrianglesPerNode);
}

void BoundingBoxBVTree::getClosestBoundingBoxes(const ArrayRef<BoundingBox> boundingBoxes, const Vec3d &queryPosition,
  std::vector<int> &bbIDList) const
{
  int bbIDListLastEnd = bbIDList.size();
  PGO_ALOG(numBBs == boundingBoxes.size());

  auto toBB = [&](const BoundingBox &nodeBB) -> std::pair<double, double> {
    double dist2 = nodeBB.distanceToPoint2(queryPosition) * (1 - eps);
    double furDist2 = nodeBB.furthestDistanceToPoint2(queryPosition) * (1 + eps);
    return { dist2, furDist2 };
  };

  auto toElement = [&](int bbID, double minDistance) -> double {
    double dist2 = boundingBoxes[bbID].distanceToPoint2(queryPosition);

    if (dist2 < minDistance)             // if has smaller value
    {
      bbIDList.resize(bbIDListLastEnd);  // clear previous stored values
      bbIDList.push_back(bbID);
    }
    else if (dist2 == minDistance) {
      bbIDList.push_back(bbID);
    }
    return dist2;
  };

  // we assign the fourth parameter of nearestQuery() to be distanceLowerBound = -1.0 to avoid early termination of the function
  // once it finds minDistance == 0.0
  // In this way, we can find all the bounding boxes that contain queryPosition
  nearestQuery(toBB, toElement, DBL_MAX, -1.0);

  sortAndDeduplicate(bbIDList, bbIDListLastEnd);
}

void BoundingBoxBVTree::intersectAABB(const ArrayRef<BoundingBox> bbs, const BoundingBox &bbQuery, std::vector<int> &bbIDList) const
{
  auto toBB = [&](const BoundingBox &bb) -> bool {
    return bbQuery.intersect(bb);
  };

  auto toEle = [&](int idx) -> bool {
    return bbQuery.intersect(bbs[idx]);
  };

  BoundingVolumeTreeBase::rangeQuery(toBB, toEle, bbIDList);
}

/////////////////////////////////////////////////////////////
//                   TriMeshBVTree
/////////////////////////////////////////////////////////////

void TriMeshBVTree::clear()
{
  nodes.clear();
  triBBs.clear();
  numVertices = numTriangles = 0;
}

void TriMeshBVTree::buildExactOctree(const TriMeshRef triMesh, int maxDepth, int maxNumTrianglesPerNode)
{
  nodes.clear();
  numVertices = triMesh.numVertices();
  numTriangles = triMesh.numTriangles();
  PGO_ALOG(numTriangles > 0);

  triBBs = triMesh.getTriangleBoundingBoxes();
  std::vector<int> rootIndices(triMesh.numTriangles());
  std::iota(rootIndices.begin(), rootIndices.end(), 0);
  BoundingBox rootBB(triMesh.positions(), triMesh.triangles(), rootIndices);
  nodes.emplace_back(0, std::move(rootIndices));
  // We don't build a bounding box out of all triMesh's positions because some vertices might not
  // be used by triangles
  nodes[0].bb = rootBB;

  // cout << streamRange(rootIndices) << endl;
  // cout << nodes[0].bb << endl;
  // PGO_ALOG(nodes[0].bb.verifyBox());

  auto bbFromTris = [&](const std::vector<int> &triList) -> BoundingBox {
    // get the bounding box of all triangles inside this childNode
    auto bb = BoundingBox(triMesh.positions(), triMesh.triangles(), triList);
    // PGO_ALOG(bb.verifyBox());
    return bb;
  };
  auto triBBIntersect = [&](int triID, const BoundingBox &bb) -> bool {
    // use the exact predicate
    return intersectTriAABB(triMesh.pos(triID, 0).data(), triMesh.pos(triID, 1).data(), triMesh.pos(triID, 2).data(), bb.bmin().data(), bb.bmax().data());
  };

  auto subdivide = [&](const std::vector<int> &elementList, const BoundingBox &bb,
                     std::vector<std::vector<int>> &childIDs, std::vector<BoundingBox> &subBBs) -> void {
    // PGO_ALOG(bb.verifyBox());
    boundingBoxOctreePartition(elementList, bb, childIDs, subBBs, triBBs, bbFromTris, triBBIntersect);
  };

  buildFromRoot(subdivide, maxDepth, maxNumTrianglesPerNode);
}

void TriMeshBVTree::buildInexactOctree(const TriMeshRef triMesh, const BoundingBox &initialBoundingBox, int maxDepth, int maxNumTrianglesPerNode)
{
  nodes.clear();
  numVertices = triMesh.numVertices();
  numTriangles = triMesh.numTriangles();
  PGO_ALOG(numTriangles > 0);
  triBBs = triMesh.getTriangleBoundingBoxes();
  std::vector<int> rootIndices(triMesh.numTriangles());
  std::iota(rootIndices.begin(), rootIndices.end(), 0);
  nodes.emplace_back(0, std::move(rootIndices));
  nodes[0].bb = initialBoundingBox;

  auto bbFromTris = [&](const std::vector<int> &triList) -> BoundingBox {
    // get the bounding box of all triangles inside this childNode
    auto bb = BoundingBox(triMesh.positions(), triMesh.triangles(), triList);
    // PGO_ALOG(bb.verifyBox());
    return bb;
  };
  auto triBBIntersect = [&](int triID, const BoundingBox &bb) -> bool {
    // inexact tri/BB intersection test:
    return triBoxOverlap(bb.center().data(), bb.halfSides().data(),
      triMesh.pos(triID, 0).data(), triMesh.pos(triID, 1).data(), triMesh.pos(triID, 2).data());
  };

  auto subdivide = [&](const std::vector<int> &elementList, const BoundingBox &bb,
                     std::vector<std::vector<int>> &childIDs, std::vector<BoundingBox> &subBBs) -> void {
    // PGO_ALOG(bb.verifyBox());
    const bool shrink = false;  // we don't shrink bounding box here to prevent numerical problems rising
    boundingBoxOctreePartition(elementList, bb, childIDs, subBBs, triBBs, bbFromTris, triBBIntersect, shrink);
  };

  buildFromRoot(subdivide, maxDepth, maxNumTrianglesPerNode);
}

void TriMeshBVTree::buildByInertiaPartition(const TriMeshRef triMesh, int maxDepth, int maxNumTrianglesPerNode)
{
  nodes.clear();
  numVertices = triMesh.numVertices();
  numTriangles = triMesh.numTriangles();
  PGO_ALOG(numTriangles > 0);
  triBBs = triMesh.getTriangleBoundingBoxes();
  std::vector<int> rootIndices(triMesh.numTriangles());
  std::iota(rootIndices.begin(), rootIndices.end(), 0);
  BoundingBox rootBB(triMesh.positions(), triMesh.triangles(), rootIndices);
  nodes.emplace_back(0, std::move(rootIndices));
  // We don't build a bounding box out of all triMesh's positions because some vertices might not
  // be used by triangles
  nodes[0].bb = rootBB;

  std::vector<double> triSurfaceAreas(numTriangles);
  std::vector<Vec3d> triCenterOfMass(numTriangles);
  std::vector<Mat3d> triInertia(numTriangles);
  for (int i = 0; i < numTriangles; i++) {
    Vec3d v0 = triMesh.pos(i, 0), v1 = triMesh.pos(i, 1), v2 = triMesh.pos(i, 2);
    triSurfaceAreas[i] = getTriangleArea(v0, v1, v2);
    triCenterOfMass[i] = getTriangleCenterOfMass(v0, v1, v2);
    triInertia[i] = getTriangleInertiaTensorAroundCOM(v0, v1, v2, triCenterOfMass[i], triSurfaceAreas[i]);
  }

  auto subdivide = [&](const std::vector<int> &elementList, const BoundingBox &bb,
                     std::vector<std::vector<int>> &childIDs, std::vector<BoundingBox> &subBBs) -> void {
    // PGO_ALOG(bb.verifyBox());
    elementsInertiaBinaryPartition(
      elementList, bb, childIDs, INT_MAX,
      [&](int i) { return triSurfaceAreas[i]; },
      [&](int i) { return triCenterOfMass[i]; },
      [&](int i, const Vec3d &reference) { return shiftInertiaTensorAroundMassCenterToReferencePoint(triInertia[i], triSurfaceAreas[i],
                                             reference - triCenterOfMass[i]); });
    subBBs.resize(2);
    for (int i = 0; i < sizei(childIDs); i++)
      subBBs[i] = BoundingBox(triMesh.positions(), triMesh.triangles(), childIDs[i]);
  };

  buildFromRoot(subdivide, maxDepth, maxNumTrianglesPerNode);
}

void TriMeshBVTree::updateBoundingVolumes(const TriMeshRef triMesh)
{
  PGO_ALOG(numVertices == triMesh.numVertices());
  PGO_ALOG(numTriangles == triMesh.numTriangles());
  for (auto it = nodes.rbegin(); it != nodes.rend(); it++) {
    Node &node = *it;
    if (node.isLeaf()) {
      PGO_ALOG(node.indices.size() > 0);
      node.bb = BoundingBox(triMesh.positions(), triMesh.triangles(), node.indices);
    }
    else {
      PGO_ALOG(node.childrenIDs.size() > 0);
      node.bb = nodes[node.childrenIDs[0]].bb;
      for (size_t i = 1; i < node.childrenIDs.size(); i++)
        node.bb.expand(nodes[node.childrenIDs[i]].bb);
    }
  }

  for (int triID = 0; triID < numTriangles; triID++) {
    triBBs[triID] = BoundingBox(triMesh.positions(), triMesh.tri(triID));
  }
}

void TriMeshBVTree::selfIntersectionExact(const TriMeshRef triMesh, std::vector<std::pair<int, int>> &triangleIDList) const
{
  // cout << "selfIntersectionExact" << endl;
  PGO_ALOG(numVertices == triMesh.numVertices());
  PGO_ALOG(numTriangles == triMesh.numTriangles());
  PGO_ALOG(nodes.empty() == false);

  std::unordered_set<UEdgeKey> candidateSet;

  auto triTriProcess = [&](int triIDA, int triIDB) {
    if (overlap(triMesh.tri(triIDA), triMesh.tri(triIDB)))
      return;  // if two triangles are neighbors

    UEdgeKey ue(triIDA, triIDB);

    if (usetFind(candidateSet, ue))
      return;

    // do tri vs tri intersection
    if (triBBs[triIDA].intersect(triBBs[triIDB]) == false)
      return;

    candidateSet.insert(ue);
  };

  selfIntersectionQuery(triTriProcess);

  //  sortAndDeduplicate(candidatePairs);
  std::vector<UEdgeKey> candidatePairs(candidateSet.begin(), candidateSet.end());  // store all those candidate triangle pairs for parallel evaluations
  std::vector<char> intersected(candidatePairs.size(), 0);

  tbb::parallel_for(0, sizei(candidatePairs), [&](int i) {
    // for (size_t i = 0; i < candidatePairs.size(); ++i) {

    int triIDA = candidatePairs[i][0], triIDB = candidatePairs[i][1];
    if (intersectTriTri(triMesh.pos(triIDA, 0).data(), triMesh.pos(triIDA, 1).data(), triMesh.pos(triIDA, 2).data(),
          triMesh.pos(triIDB, 0).data(), triMesh.pos(triIDB, 1).data(), triMesh.pos(triIDB, 2).data())) {
      intersected[i] = 1;
    }
  });  // end for locations

  for (size_t i = 0; i < candidatePairs.size(); i++) {
    if (intersected[i]) {
      triangleIDList.emplace_back(candidatePairs[i][0], candidatePairs[i][1]);
    }
  }

  // cout << "selfIntersectionExact Done" << endl;
}

void TriMeshBVTree::intersectionExact(const TriMeshRef triMesh, const TriMeshRef &otherMesh,
  std::vector<std::pair<int, int>> &triangleIDList) const
{
  PGO_ALOG(numVertices == triMesh.numVertices());
  PGO_ALOG(numTriangles == triMesh.numTriangles());
  PGO_ALOG(nodes.empty() == false);

  std::vector<std::pair<int, int>> allPairList;

  struct ThreadLocalData
  {
    std::vector<int> IDlist;
    std::vector<std::pair<int, int>> pairList;
  };

  tbb::enumerable_thread_specific<ThreadLocalData> threadLocalData;
  tbb::parallel_for(tbb::blocked_range<int>(0, otherMesh.numTriangles()), [&](const tbb::blocked_range<int> &rng) {
    //  for(int oID = 0; oID < otherMesh.numTriangles(); oID++)
    auto &local = threadLocalData.local();
    auto &IDlist = local.IDlist;
    auto &pairList = local.pairList;
    for (int oID = rng.begin(); oID != rng.end(); oID++) {
      IDlist.clear();
      triangleIntersectionExact(triMesh, otherMesh.pos(oID, 0), otherMesh.pos(oID, 1), otherMesh.pos(oID, 2), IDlist);
      sortAndDeduplicate(IDlist);
      for (int triID : IDlist)
        pairList.emplace_back(triID, oID);
    }
  });
  for (const auto &local : threadLocalData)
    vectorInsertRangeBack(allPairList, local.pairList);

  sortAndDeduplicate(allPairList);
  vectorInsertRangeBack(triangleIDList, allPairList);
}

void TriMeshBVTree::intersectionExact(const TriMeshRef triMesh, const TriMeshBVTree &otherTree, const TriMeshRef &otherMesh,
  std::vector<std::pair<int, int>> &triangleIDList) const
{
  PGO_ALOG(numVertices == triMesh.numVertices());
  PGO_ALOG(numTriangles == triMesh.numTriangles());
  PGO_ALOG(nodes.empty() == false);
  PGO_ALOG(otherTree.numVertices == otherMesh.numVertices());
  PGO_ALOG(otherTree.numTriangles == otherMesh.numTriangles());
  PGO_ALOG(otherTree.nodes.empty() == false);

  std::vector<std::pair<int, int>> candidatePairs;
  auto triTriIntersect = [&](int triIDA, int triIDB) -> void {
    if (triBBs[triIDA].intersect(otherTree.triBBs[triIDB]))
      candidatePairs.emplace_back(triIDA, triIDB);
  };

  intersectionQuery(otherTree, triTriIntersect);

  // cout << "Finish finding candidate pairs: " << candidatePairs.size() << endl;
  sortAndDeduplicate(candidatePairs);

  std::vector<char> intersected(candidatePairs.size(), 0);
  tbb::parallel_for(
    tbb::blocked_range<int>(0, candidatePairs.size()), [&](const tbb::blocked_range<int> &rng) {
      for (int i = rng.begin(); i != rng.end(); ++i) {
        int triIDA = candidatePairs[i].first, triIDB = candidatePairs[i].second;
        PGO_ALOG(triIDA >= 0 && triIDA < numTriangles);
        PGO_ALOG(triIDB >= 0 && triIDB < otherMesh.numTriangles());
        if (intersectTriTri(triMesh.pos(triIDA, 0).data(), triMesh.pos(triIDA, 1).data(), triMesh.pos(triIDA, 2).data(),
              otherMesh.pos(triIDB, 0).data(), otherMesh.pos(triIDB, 1).data(), otherMesh.pos(triIDB, 2).data())) {
          intersected[i] = 1;
        }
      }
    },
    tbb::auto_partitioner());  // end for locations

  for (size_t i = 0; i < candidatePairs.size(); i++) {
    if (intersected[i])
      triangleIDList.push_back(candidatePairs[i]);
  }
}

void TriMeshBVTree::lineSegmentIntersectionExact(const TriMeshRef triMesh, Vec3d segStart, Vec3d segEnd, std::vector<int> &triangleIDList) const
{
  PGO_ALOG(numVertices == triMesh.numVertices());
  PGO_ALOG(numTriangles == triMesh.numTriangles());

  auto toBB = [&](const BoundingBox &bb) -> bool {
    return intersectSegAABB(&segStart[0], &segEnd[0], &bb.bmin()[0], &bb.bmax()[0]);
  };

  auto toEle = [&](int idx) -> bool {
    return intersectSegTri(&segStart[0], &segEnd[0], &triMesh.pos(idx, 0)[0], &triMesh.pos(idx, 1)[0], &triMesh.pos(idx, 2)[0]);
  };

  BoundingVolumeTreeBase::rangeQuery(toBB, toEle, triangleIDList);
}

// note: no tbb threading should be used in triangleIntersectionExact
// because I use tbb on intersectionExact with another TriMesh
// and this function calls triangleIntersectionExact in an parallelized loop
// nested multi-threading can be a problem when I use thread-local data
void TriMeshBVTree::triangleIntersectionExact(const TriMeshRef triMesh, Vec3d t0, Vec3d t1, Vec3d t2, std::vector<int> &triangleIDList) const
{
  PGO_ALOG(numVertices == triMesh.numVertices());
  PGO_ALOG(numTriangles == triMesh.numTriangles());

  auto toBB = [&](const BoundingBox &bb) -> bool {
    return intersectTriAABB(&t0[0], &t1[0], &t2[0], &bb.bmin()[0], &bb.bmax()[0]);
  };

  auto toEle = [&](int idx) -> bool {
    return intersectTriTri(&t0[0], &t1[0], &t2[0], &triMesh.pos(idx, 0)[0], &triMesh.pos(idx, 1)[0], &triMesh.pos(idx, 2)[0]);
  };

  BoundingVolumeTreeBase::rangeQuery(toBB, toEle, triangleIDList);
}

int TriMeshBVTree::lineSegmentFirstIntersectionPoint(const TriMeshRef triMesh, Vec3d segStart, Vec3d segEnd,
  double segWeight[2], double triWeight[3]) const
{
  PGO_ALOG(numVertices == triMesh.numVertices());
  PGO_ALOG(numTriangles == triMesh.numTriangles());

  auto toBB = [&](const BoundingBox &bb) -> bool {
    Vec3d intersectionPoint;
    return bb.lineSegmentIntersection(segStart, segEnd, &intersectionPoint);
  };

  if (segWeight)
    segWeight[0] = segWeight[1] = DBL_MAX;
  if (triWeight)
    triWeight[0] = triWeight[1] = triWeight[2] = DBL_MAX;

  int closestTriID = -1;
  double closestWeight = DBL_MAX;
  auto toEle = [&](int triID) -> void {
    TriangleBasic tb(triMesh.pos(triID, 0), triMesh.pos(triID, 1), triMesh.pos(triID, 2));

    double t = 0.0;
    double baryWeight[3];
    int collisionResult = tb.lineSegmentIntersection(segStart, segEnd, nullptr, &t, baryWeight);
    if (collisionResult > 0) {
      if (t < closestWeight) {
        closestTriID = triID;
        closestWeight = t;
        if (segWeight) {
          segWeight[0] = 1.0 - t;
          segWeight[1] = t;
        }
        if (triWeight)
          std::memcpy(triWeight, baryWeight, sizeof(double) * 3);
      }
    }
  };

  BoundingVolumeTreeBase::rangeQuery(toBB, toEle);
  return closestTriID;
}

int TriMeshBVTree::lineSegmentFirstIntersectionPointSemiExact(const TriMeshRef triMesh, Vec3d segStart, Vec3d segEnd,
  double segWeight[2], double triWeight[3], std::stack<int> *externalStack) const
{
  PGO_ALOG(numVertices == triMesh.numVertices());
  PGO_ALOG(numTriangles == triMesh.numTriangles());

  auto toBB = [&](const BoundingBox &bb) -> bool {
    return intersectSegAABB(&segStart[0], &segEnd[0], &bb.bmin()[0], &bb.bmax()[0]);
  };

  if (segWeight)
    segWeight[0] = segWeight[1] = DBL_MAX;
  if (triWeight)
    triWeight[0] = triWeight[1] = triWeight[2] = DBL_MAX;

  int closestTriID = -1;
  double closestWeight = DBL_MAX;
  auto toEle = [&](int triID) -> void {
    double alpha[5];
    if (intersectSegTri(segStart.data(), segEnd.data(), triMesh.pos(triID, 0).data(), triMesh.pos(triID, 1).data(), triMesh.pos(triID, 2).data(), alpha, alpha + 2)) {
      // intersection_point = segStart * alpha[0] + segEnd * alpha[1]
      if (alpha[1] < closestWeight) {
        closestWeight = alpha[1];
        closestTriID = triID;

        if (segWeight)
          std::memcpy(segWeight, alpha, sizeof(double) * 2);

        if (triWeight)
          std::memcpy(triWeight, alpha + 2, sizeof(double) * 3);
      }
    }
  };

  BoundingVolumeTreeBase::rangeQuery(toBB, toEle, externalStack);
  return closestTriID;
}

bool TriMeshBVTree::hasLineSegmentIntersectionExact(const TriMeshRef triMesh, Vec3d segStart, Vec3d segEnd, std::vector<int> *skippedTriIDs) const
{
  PGO_ALOG(numVertices == triMesh.numVertices());
  PGO_ALOG(numTriangles == triMesh.numTriangles());

  auto toBB = [&](const BoundingBox &bb) -> bool {
    return intersectSegAABB(&segStart[0], &segEnd[0], &bb.bmin()[0], &bb.bmax()[0]);
  };

  auto toEle = [&](int idx) -> bool {
    if (skippedTriIDs && binarySearchFound(*skippedTriIDs, idx) == true)
      return false;
    return intersectSegTri(&segStart[0], &segEnd[0], &triMesh.pos(idx, 0)[0], &triMesh.pos(idx, 1)[0], &triMesh.pos(idx, 2)[0]);
  };

  return BoundingVolumeTreeBase::rangeQueryWithEarlyStop(toBB, toEle);
}

std::vector<int> TriMeshBVTree::trianglesUnderNode(int nodeID) const
{
  std::vector<int> ret;
  std::vector<int> candidates = { nodeID };
  int candBegin = 0, candEnd = 1;
  while (candBegin < candEnd) {
    for (int i = candBegin; i < candEnd; i++) {
      int ID = candidates[i];
      if (nodes[ID].isLeaf()) {
        ret.insert(ret.end(), nodes[ID].indices.begin(), nodes[ID].indices.end());
      }
      else {
        for (int childID : nodes[ID].childrenIDs) {
          if (childID < 0)
            continue;
          candidates.push_back(childID);
        }
      }
    }
    candBegin = candEnd;
    candEnd = int(candidates.size());
  }

  sortAndDeduplicate(ret);
  return ret;
}

template<class GetTriangleDistance2>
std::pair<int, double> TriMeshBVTree::getClosestTriangleHelper(GetTriangleDistance2 getTriangleDistance2, const Vec3d &queryPosition,
  int *feature, int warmStartTriID, double dist2HigherBound, std::vector<std::tuple<double, double, int>> &nodeStack) const
{
  int closestTriID = -1, closestFeature = -1;

  auto toBB = [&](const BoundingBox &bb) -> std::pair<double, double> {
    return { bb.distanceToPoint2(queryPosition) * (1 - eps), bb.furthestDistanceToPoint2(queryPosition) * (1 + eps) };
  };

  auto toElement = [&](int triID, double minDistance) -> double {
    int feature = -1;
    double dist2 = getTriangleDistance2(queryPosition, triID, feature);

    if (dist2 <= minDistance) {
      closestTriID = triID;
      closestFeature = feature;
    }
    return dist2;
  };

  if (warmStartTriID >= 0) {
    double dist2 = getTriangleDistance2(queryPosition, warmStartTriID, closestFeature);
    if (dist2 <= dist2HigherBound) {
      closestTriID = warmStartTriID;
      dist2HigherBound = dist2;
    }
  }
  double minDist = nearestQuery(nodeStack, toBB, toElement, dist2HigherBound);

  if (feature)
    *feature = closestFeature;
  return std::make_pair(closestTriID, minDist);
}

template<class GetTriangleDistance2>
std::pair<int, double> TriMeshBVTree::getClosestTriangleHelper(GetTriangleDistance2 getTriangleDistance2, const Vec3d &queryPosition, int *feature,
  int warmStartTriID, double distance2HigherBound) const
{
  std::vector<std::tuple<double, double, int>> nodeStack;
  return getClosestTriangleHelper<GetTriangleDistance2>(getTriangleDistance2, queryPosition, feature,
    warmStartTriID, distance2HigherBound, nodeStack);
}

int TriMeshBVTree::getClosestTriangle(const TriMeshRef triMesh, const Vec3d &queryPosition, int *feature, int warmStartTriID,
  double dist2HigherBound) const
{
  PGO_ALOG(triMesh.numVertices() == numVertices && triMesh.numTriangles() == numTriangles);

  auto getTriangleDist2 = [&](const Vec3d &queryPosition, int triID, int &feature) {
    return getSquaredDistanceToTriangle(queryPosition, triMesh.pos(triID, 0), triMesh.pos(triID, 1), triMesh.pos(triID, 2), feature);
  };

  return getClosestTriangleHelper(getTriangleDist2, queryPosition, feature, warmStartTriID, dist2HigherBound).first;
}

int TriMeshBVTree::getClosestTriangle(const TriangleWithCollisionInfo *triangleWithColInfos, const Vec3d &queryPosition, int *feature,
  int warmStartTriID, double dist2HigherBound) const
{
  auto getTriangleDist2 = [&](const Vec3d &queryPosition, int triID, int &feature) {
    return triangleWithColInfos[triID].distanceToPoint2(queryPosition, &feature);
  };

  return getClosestTriangleHelper(getTriangleDist2, queryPosition, feature, warmStartTriID, dist2HigherBound).first;
}

TriMeshBVTree::ClosestTriangleQueryResult TriMeshBVTree::closestTriangleQuery(const TriMeshRef triMesh, const Vec3d &queryPosition,
  int warmStartTriID, double dist2HigherBound) const
{
  std::vector<std::tuple<double, double, int>> nodeStack;
  return closestTriangleQuery(triMesh, queryPosition, nodeStack, warmStartTriID, dist2HigherBound);
}

TriMeshBVTree::ClosestTriangleQueryResult TriMeshBVTree::closestTriangleQuery(const TriMeshRef triMesh, const Vec3d &queryPosition,
  std::vector<std::tuple<double, double, int>> &nodeStack,
  int warmStartTriID, double dist2HigherBound) const
{
  PGO_ALOG(triMesh.numVertices() == numVertices && triMesh.numTriangles() == numTriangles);

  ClosestTriangleQueryResult ret;

  auto getTriangleDist2 = [&](const Vec3d &queryPosition, int triID, int &feature) {
    return getSquaredDistanceToTriangle(queryPosition, triMesh.pos(triID, 0), triMesh.pos(triID, 1), triMesh.pos(triID, 2), feature);
  };

  auto returnedPair = getClosestTriangleHelper(getTriangleDist2, queryPosition, &(ret.feature), warmStartTriID, dist2HigherBound, nodeStack);
  ret.triID = returnedPair.first;
  ret.dist2 = returnedPair.second;

  if (ret.triID >= 0) {
    ret.closestPosition = getClosestPointToTriangleWithFeature(queryPosition,
      triMesh.pos(ret.triID, 0), triMesh.pos(ret.triID, 1), triMesh.pos(ret.triID, 2), ret.feature);

    ret.triBaryWeight = getBarycentricWeightProjectedOnTrianglePlane(ret.closestPosition,
      triMesh.pos(ret.triID, 0), triMesh.pos(ret.triID, 1), triMesh.pos(ret.triID, 2));
  }

  return ret;
}

TriMeshBVTree::ClosestTriangleQueryResult TriMeshBVTree::closestTriangleQuery(const TriangleWithCollisionInfo *triangles,
  const Vec3d &queryPosition, std::vector<std::tuple<double, double, int>> &nodeStack,
  int warmStartTriID, double dist2HigherBound) const
{
  ClosestTriangleQueryResult ret;
  auto getTriangleDist2 = [&](const Vec3d &queryPosition, int triID, int &feature) {
    return triangles[triID].distanceToPoint2(queryPosition, &feature);
  };

  auto returnedPair = getClosestTriangleHelper(getTriangleDist2, queryPosition, &(ret.feature), warmStartTriID, dist2HigherBound, nodeStack);
  ret.triID = returnedPair.first;
  ret.dist2 = returnedPair.second;
  if (ret.triID >= 0) {
    const auto &tri = triangles[ret.triID];
    ret.closestPosition = getClosestPointToTriangleWithFeature(queryPosition, tri[0], tri[1], tri[2], ret.feature);
    ret.triBaryWeight = getBarycentricWeightProjectedOnTrianglePlane(queryPosition, tri[0], tri[1], tri[2]);
  }
  return ret;
}

bool TriMeshBVTree::sanityCheck(const TriMeshRef triMesh) const
{
  PGO_ALOG(triMesh.numVertices() == numVertices && triMesh.numTriangles() == numTriangles);

  for (size_t i = 0; i < nodes.size(); i++) {
    const BoundingBox &bb = nodes[i].bb;
    for (int childID : nodes[i].childrenIDs) {
      if (childID < 0)
        continue;
      const BoundingBox &childBB = nodes[childID].bb;
      if (bb.checkInside(childBB) == false) {
        std::cout << "Error: boundingBox ID " << childID << " is not inside its parent: ID " << i << std::endl;
        bb.print();
        childBB.print();
        return false;
      }
    }
    // we don't check whether triangles are inside the box because our construction
    // allow triangles to be partially outside as long as the outside part are covered by other boxes
  }
  return true;
}

/////////////////////////////////////////////////////////////
//                   TetMeshBVTree
/////////////////////////////////////////////////////////////

void TetMeshBVTree::clear()
{
  nodes.clear();
  tetBBs.clear();
  numVertices = numTets = 0;
}

void TetMeshBVTree::buildByInertiaPartition(const TetMeshRef tetMesh, int maxDepth, int maxNumTetsPerNode)
{
  nodes.clear();
  numVertices = tetMesh.numVertices();
  numTets = tetMesh.numTets();
  PGO_ALOG(numTets > 0);
  tetBBs = tetMesh.getTetBoundingBoxes();
  std::vector<int> rootIndices;
  if (elementFlags.size()) {
    for (int i = 0; i < numTets; i++) {
      if (elementFlags[i]) {
        rootIndices.emplace_back(i);
      }
    }
  }
  else {
    rootIndices.resize(numTets);
    iota(rootIndices.begin(), rootIndices.end(), 0);
  }

  nodes.emplace_back(0, std::move(rootIndices));
  nodes[0].bb = BoundingBox(tetMesh.numVertices(), tetMesh.positions());

  std::vector<double> tetVolume(numTets);
  std::vector<Vec3d> tetCenterOfMass(numTets);
  std::vector<Mat3d> tetInertia(numTets);
  for (int i = 0; i < numTets; i++) {
    Vec3d v0 = tetMesh.pos(i, 0), v1 = tetMesh.pos(i, 1), v2 = tetMesh.pos(i, 2), v3 = tetMesh.pos(i, 3);
    tetVolume[i] = getTetVolume(v0, v1, v2, v3);
    tetCenterOfMass[i] = getTetCenterOfMass(v0, v1, v2, v3);
    tetInertia[i] = getTetInertiaTensorAroundCOM(v0, v1, v2, v3, tetCenterOfMass[i], tetVolume[i]);
  }

  auto subdivide = [&](const std::vector<int> &elementList, const BoundingBox &bb,
                     std::vector<std::vector<int>> &childIDs, std::vector<BoundingBox> &subBBs) -> void {
    // PGO_ALOG(bb.verifyBox());
    elementsInertiaBinaryPartition(
      elementList, bb, childIDs, INT_MAX,
      [&](int i) { return tetVolume[i]; },
      [&](int i) { return tetCenterOfMass[i]; },
      [&](int i, const Vec3d &reference) { return shiftInertiaTensorAroundMassCenterToReferencePoint(tetInertia[i], tetVolume[i],
                                             reference - tetCenterOfMass[i]); });
    subBBs.resize(2);
    for (int i = 0; i < sizei(childIDs); i++)
      subBBs[i] = BoundingBox(tetMesh.positions(), tetMesh.tets(), childIDs[i]);
    //    cout << "Done partitioning on " << elementList.size() << " elements, children sizes are ";
    //    for(const auto & c : childIDs)
    //      cout << c.size() << " ";
    //    cout << endl;
  };

  buildFromRoot(subdivide, maxDepth, maxNumTetsPerNode);

  elementFlags.assign(tetMesh.numTets(), 1);
}

void TetMeshBVTree::buildExactOctree(const TetMeshRef tetMesh, int maxDepth, int maxNumTetsPerNode)
{
  nodes.clear();
  numVertices = tetMesh.numVertices();
  numTets = tetMesh.numTets();
  tetBBs = tetMesh.getTetBoundingBoxes();
  std::vector<int> rootIndices(numTets);
  std::iota(rootIndices.begin(), rootIndices.end(), 0);

  nodes.emplace_back(0, std::move(rootIndices));
  nodes[0].bb = BoundingBox(tetMesh.numVertices(), tetMesh.positions());

  auto bbFromTris = [&](const std::vector<int> &tetList) -> BoundingBox {
    // get the bounding box of all triangles inside this childNode
    std::vector<int> allVtxIDsInThisNode;
    tetMesh.getVerticesInTets(tetList, allVtxIDsInThisNode);
    return BoundingBox(tetMesh.positions(), allVtxIDsInThisNode);
  };
  auto tetBBIntersect = [&](int tetID, const BoundingBox &bb) -> bool {
    return intersectTetAABB(tetMesh.pos(tetID, 0).data(), tetMesh.pos(tetID, 1).data(), tetMesh.pos(tetID, 2).data(), tetMesh.pos(tetID, 3).data(),
      bb.bmin().data(), bb.bmax().data());
  };
  auto subdivide = [&](const std::vector<int> &elementList, const BoundingBox &bb,
                     std::vector<std::vector<int>> &childIDs, std::vector<BoundingBox> &subBBs) -> void {
    boundingBoxOctreePartition(elementList, bb, childIDs, subBBs, tetBBs, bbFromTris, tetBBIntersect);
  };

  buildFromRoot(subdivide, maxDepth, maxNumTetsPerNode);
}

int TetMeshBVTree::getClosestTet(const TetMeshRef tetMesh, const Vec3d &queryPosition) const
{
  PGO_ALOG(numVertices == tetMesh.numVertices());
  PGO_ALOG(numTets == tetMesh.numTets());

  int closestTetID = -1;

  auto toBB = [&](const BoundingBox &bb) -> std::pair<double, double> {
    return { bb.distanceToPoint2(queryPosition) * (1 - eps), bb.furthestDistanceToPoint2(queryPosition) * (1 + eps) };
  };

  auto toElement = [&](int tetID, double minDistance) -> double {
    if (elementFlags[tetID] == 0) {
      return 1e300;
    }

    double dist2 = getSquaredDistanceToTet(queryPosition, tetMesh.pos(tetID, 0), tetMesh.pos(tetID, 1), tetMesh.pos(tetID, 2), tetMesh.pos(tetID, 3));
    if (dist2 <= minDistance) {
      closestTetID = tetID;
    }
    return dist2;
  };

  nearestQuery(toBB, toElement);

  PGO_ALOG(closestTetID >= 0);
  return closestTetID;
}

int TetMeshBVTree::getClosestTet(const TetMeshRef tetMesh, const Vec3d &queryPosition, std::vector<std::tuple<double, double, int>> &nodeStack) const
{
  PGO_ALOG(numVertices == tetMesh.numVertices());
  PGO_ALOG(numTets == tetMesh.numTets());

  int closestTetID = -1;

  auto toBB = [&](const BoundingBox &bb) -> std::pair<double, double> {
    return { bb.distanceToPoint2(queryPosition) * (1 - eps), bb.furthestDistanceToPoint2(queryPosition) * (1 + eps) };
  };

  auto toElement = [&](int tetID, double minDistance) -> double {
    if (elementFlags[tetID] == 0) {
      return 1e300;
    }

    double dist2 = getSquaredDistanceToTet(queryPosition, tetMesh.pos(tetID, 0), tetMesh.pos(tetID, 1), tetMesh.pos(tetID, 2), tetMesh.pos(tetID, 3));
    if (dist2 <= minDistance) {
      closestTetID = tetID;
    }
    return dist2;
  };

  nearestQuery(nodeStack, toBB, toElement);

  PGO_ALOG(closestTetID >= 0);
  return closestTetID;
}

}  // namespace Mesh
}  // namespace pgo