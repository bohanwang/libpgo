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

#pragma once

// An octree implementation using exact predicates for precise geometry queries.
// The base class is BoundingVolumeTreeBase. It has several derived classes:
// VertexBVTree, TriMeshBVTree and TetMeshBVTree.

#include "boundingVolumeTreeBase.h"
#include "boundingBox.h"
#include "simpleSphere.h"
#include "halfSpace.h"
#include "triMeshGeo.h"
#include "tetMeshGeo.h"
#include "triangle.h"

#include "arrayRef.h"

#include <cstring>
#include <vector>
#include <array>
#include <cassert>
#include <functional>
#include <cfloat>
#include <climits>
#include <stack>

namespace pgo
{
namespace Mesh
{
class TriangleWithCollisionInfo;  // a better triangle class which pre-computes necessary data for fast intersection queries
class TriangleSampler;

// bounding volume tree to query vertices

class VertexBVTree : public BoundingVolumeTreeBase
{
public:
  VertexBVTree() {}
  virtual ~VertexBVTree() { clear(); }

  // maxDepth: maximum tree depth allowed
  // build the octree
  void buildOctree(const BasicAlgorithms::ArrayRef<Vec3d> verticesList, int maxDepth = 10, int maxNumVerticesPerNode = 10);

  // Build a BV tree by partitioning the mesh along a primary axis. The primary axis is computed from doing SVD on the inertia tensor of the vertex group.
  // Each tree node will have at most two children.
  void buildByInertiaPartition(const BasicAlgorithms::ArrayRef<Vec3d> verticesList, int maxDepth = INT_MAX, int maxNumTrianglesPerNode = 5);

  // query vertices inside simpleSphere/halfSpace with double-precision
  // vertexIDList are not sorted or deduplicated
  void rangeQuery(const BasicAlgorithms::ArrayRef<Vec3d> vertices, const SimpleSphere &simpleSphere, std::vector<int> &vertexIDList) const;
  void rangeQuery(const BasicAlgorithms::ArrayRef<Vec3d> vertices, const HalfSpace &halfSpace, std::vector<int> &vertexIDList) const;
  bool intersectAABB(const BasicAlgorithms::ArrayRef<Vec3d> vertices, const BoundingBox &bb) const;

  // double precision
  int getClosestVertex(const BasicAlgorithms::ArrayRef<Vec3d> vertices, const Vec3d &queryPosition) const;

  void clear();

protected:
  int numVertices = 0;
};

// bounding volume tree for an array of bounding boxes

class BoundingBoxBVTree : public BoundingVolumeTreeBase
{
public:
  BoundingBoxBVTree() {}
  virtual ~BoundingBoxBVTree() { clear(); }

  // Build a BV tree by partitioning the boxes along a primary axis. The primary axis is computed from doing SVD on the inertia tensor of the boxes.
  // Each tree node will have at most two children.
  void buildByInertiaPartition(const BasicAlgorithms::ArrayRef<BoundingBox> boundingBoxes, int maxDepth = INT_MAX, int maxNumBBsPerNode = 5);

  // If queryPosition is inside at least one of the bounding box, return all the bounding boxes containing it,
  // otherwise, return the closest bounding box.
  // Returned bounding box IDs are pushed into bbIDList. They are sorted and has no duplicates
  // Note the point on the surface of the BoundingBox is considered as "contained"
  void getClosestBoundingBoxes(const BasicAlgorithms::ArrayRef<BoundingBox> boundingBoxes, const Vec3d &queryPosition, std::vector<int> &bbIDList) const;

  void intersectAABB(const BasicAlgorithms::ArrayRef<BoundingBox> vertices, const BoundingBox &bb, std::vector<int> &bbIDList) const;

  void clear();

protected:
  int numBBs = 0;
};

// bounding volume tree to query a triangle mesh

class TriMeshBVTree : public BoundingVolumeTreeBase
{
public:
  TriMeshBVTree() {}
  virtual ~TriMeshBVTree() { clear(); }

  // build the exact octree where the bounding box must cover all the triangles passing to it.
  // Inexact bounding volume hierarchy will result in some triangles missing from its bounding box node if they are just touching the bounding box
  // or they completely lie on the surface of the bounding box.
  // An exact bounding volume hierarchy is needed if you want exact intersection
  // An inexact bounding volume hierarchy is acceptable if you just want a rough estimation on the distance between a point and the triangle octree
  void buildExactOctree(const TriMeshRef triMesh, int maxDepth = 10, int maxNumTrianglesPerNode = 15);

  // initialBounding box should be large enough to cover triMesh
  void buildInexactOctree(const TriMeshRef triMesh, const BoundingBox &initialBoundingBox, int maxDepth, int maxNumTrianglesPerNode);

  // Build a BV tree by partitioning the mesh along a primary axis. The primary axis is computed from doing SVD on the inertia tensor of the mesh.
  // Each tree node will have at most two children.
  void buildByInertiaPartition(const TriMeshRef triMesh, int maxDepth = INT_MAX, int maxNumTrianglesPerNode = 5);

  // assuming we have built a bounding volume hierarchy and now the triMesh geometry changed,
  // change the size of each bounding volumes to accomodate deformed geometry
  void updateBoundingVolumes(const TriMeshRef triMesh);

  const std::vector<BoundingBox> &triangleBoundingBoxes() const { return triBBs; }

  int numMeshVertices() const { return numVertices; }
  int numMeshTriangles() const { return numTriangles; }

  // void rangeQuery(const SimpleSphere & simpleSphere, std::vector<int> & vertexIDList);

  // do exact self intersection
  // input triMesh should be the same mesh as the one used in buildExact()
  // output triangleID pairs are sorted and deduplicated, the triangle indices within each pair are also sorted
  void selfIntersectionExact(const TriMeshRef triMesh, std::vector<std::pair<int, int>> &triangleIDList) const;

  // do exact intersection with another triangle mesh
  // output triangleID pairs are sorted and deduplicated
  void intersectionExact(const TriMeshRef triMesh, const TriMeshBVTree &otherBVTree, const TriMeshRef &otherMesh,
    std::vector<std::pair<int, int>> &triangleIDList) const;
  void intersectionExact(const TriMeshRef triMesh, const TriMeshRef &otherMesh,
    std::vector<std::pair<int, int>> &triangleIDList) const;

  // do exact intersection with line segment / triangle
  void lineSegmentIntersectionExact(const TriMeshRef triMesh, Vec3d segStart, Vec3d segEnd, std::vector<int> &triangleIDList) const;
  void triangleIntersectionExact(const TriMeshRef triMesh, Vec3d t0, Vec3d t1, Vec3d t2, std::vector<int> &triangleIDList) const;

  // Detect whether the line segment hits triMesh exactly and return the first triangle it hits.
  // If it does, return the first triangle ID triID. Otherwise, return -1.
  // It records the weights to locate the first intersection point from segStart to segEnd by:
  // intersection_point = segStart * segWeight[0] + segEnd * segWeight[1]
  //                    = \sum_i={0,1,2} triMesh.pos(triID, i) * triWeights[i]
  // Note:
  // For SemiExact version of the function, although it is guaranteed to find a triangle it intersects exactly if
  // buildExactOctree() or buildByInertiaPartition() is called, it cannot guarantee the returned triangle is the first to be hit,
  // since the comparison on which is the first is done in double-precision.
  // For the default version, the computation is done in double-precision so it can happen that intersections are missing.
  int lineSegmentFirstIntersectionPoint(const TriMeshRef triMesh, Vec3d segStart, Vec3d segEnd,
    double segWeight[2] = nullptr, double triWeights[3] = nullptr) const;
  int lineSegmentFirstIntersectionPointSemiExact(const TriMeshRef triMesh, Vec3d segStart, Vec3d segEnd,
    double segWeight[2] = nullptr, double triWeights[3] = nullptr, std::stack<int> *externalStack = nullptr) const;

  // if skippedTriIDs is not nullptr, it must be pre-sorted
  // triangles in skippedTriIDs will be skipped during the search
  bool hasLineSegmentIntersectionExact(const TriMeshRef triMesh, Vec3d segStart, Vec3d segEnd, std::vector<int> *skippedTriIDs = nullptr) const;

  // get the closest triangle to queryPosition, double-precision
  // distance2HigherBound: search only triangles closer or eaqual to sqrt(distance2HigherBound)
  //                      if we only want to find triangles within a certain distance to queryPosition,
  //                      then this distance2HigherBound is needed
  int getClosestTriangle(const TriMeshRef triMesh, const Vec3d &queryPosition, int *feature = nullptr,
    int warmStartTriID = -1, double distance2HigherBound = DBL_MAX) const;
  int getClosestTriangle(const TriangleWithCollisionInfo *triangles, const Vec3d &queryPosition, int *feature = nullptr,
    int warmStartTriID = -1,
    double distance2HigherBound = DBL_MAX) const;

  struct ClosestTriangleQueryResult
  {
    int triID = -1;
    int feature = -1;
    double dist2 = DBL_MAX;
    Vec3d closestPosition{ 0.0, 0.0, 0.0 };
    Vec3d triBaryWeight{ 0.0, 0.0, 0.0 };  // the barycentric weight of this closest position on the triangle: triID
  };

  ClosestTriangleQueryResult closestTriangleQuery(const TriMeshRef triMesh, const Vec3d &queryPosition,
    int warmStartTriID = -1, double distance2HigherBound = DBL_MAX) const;
  ClosestTriangleQueryResult closestTriangleQuery(const TriMeshRef triMesh, const Vec3d &queryPosition,
    std::vector<std::tuple<double, double, int>> &nodeStack,  // external stack to reduce memory allocation
    int warmStartTriID = -1, double distance2HigherBound = DBL_MAX) const;
  ClosestTriangleQueryResult closestTriangleQuery(const TriangleWithCollisionInfo *triangles, const Vec3d &queryPosition,
    std::vector<std::tuple<double, double, int>> &nodeStack,  // external stack to reduce memory allocation
    int warmStartTriID = -1, double distance2HigherBound = DBL_MAX) const;

  //  void closestTriangleQueryFromTriMeshVertices(const TriangleWithCollisionInfo * triangles, const TriMeshBVTree & otherBVTree,
  //      const TriMeshRef otherTriMesh);

  // check whether bounding boxes are fine
  bool sanityCheck(const TriMeshRef triMesh) const;

  void clear();

protected:
  std::vector<int> trianglesUnderNode(int nodeID) const;

  template<class GetTriangleDistance2>
  std::pair<int, double> getClosestTriangleHelper(GetTriangleDistance2, const Vec3d &queryPosition, int *feature,
    int warmStartTriID, double distance2HigherBound, std::vector<std::tuple<double, double, int>> &nodeStack) const;

  template<class GetTriangleDistance2>
  std::pair<int, double> getClosestTriangleHelper(GetTriangleDistance2, const Vec3d &queryPosition, int *feature,
    int warmStartTriID, double distance2HigherBound) const;

  int numVertices = 0, numTriangles = 0;
  std::vector<BoundingBox> triBBs;
};

// bounding volume tree to query a tet mesh

class TetMeshBVTree : public BoundingVolumeTreeBase
{
public:
  TetMeshBVTree() {}

  void buildByInertiaPartition(const TetMeshRef tetMesh, int maxDepth = INT_MAX, int maxNumTrianglesPerNode = 5);

  void buildExactOctree(const TetMeshRef tetMesh, int maxDepth = 10, int maxNumTetsPerNode = 15);

  // get the closest tet to queryPosition, double-precision
  int getClosestTet(const TetMeshRef tetMesh, const Vec3d &queryPosition) const;
  int getClosestTet(const TetMeshRef tetMesh, const Vec3d &queryPosition, std::vector<std::tuple<double, double, int>> &nodeStack) const;

  void setElementFlags(int count, const int *flags) { elementFlags.assign(flags, flags + count); }

  // if queryPosition is inside the leaf bounding box, get the closest tet to queryPosition
  // otherwise, return -1
  //  int getClosestTetIfInsideBoundingBox(const TetMeshRef tetMesh, const Vec3d & queryPosition, int * feature = nullptr, int warmStartTriID = -1) const;

  void clear();

protected:
  int numVertices = 0, numTets = 0;
  std::vector<BoundingBox> tetBBs;
  std::vector<int> elementFlags;
};

}  // namespace Mesh
}  // namespace pgo
