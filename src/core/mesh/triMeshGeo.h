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

#include "meshLinearAlgebra.h"
#include "boundingBox.h"
#include "edgeKey.h"

#include "arrayRef.h"

#include <vector>
#include <algorithm>
#include <map>

namespace pgo
{
namespace Mesh
{

// a triangle struct to hold triangle index, its vertex indices and positions
struct IndexedTriangle
{
  int triID = -1;
  Vec3i vtxID{ -1, -1, -1 };
  Vec3d pos[3]{ { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 } };

  IndexedTriangle() {}
  IndexedTriangle(int triID, const Vec3i &vtxID, const Vec3d &p0, const Vec3d &p1, const Vec3d &p2):
    triID(triID), vtxID(vtxID), pos{ p0, p1, p2 } {}
};

// class to reference an external triangle mesh
class TriMeshRef
{
public:
  TriMeshRef() {}  // empty mesh
  TriMeshRef(int numVertices, const double *vertices, int numTriangles, const int *triangles);
  TriMeshRef(int numVertices, const Vec3d *vertices, int numTriangles, const int *triangles);
  TriMeshRef(int numVertices, const Vec3d *vertices, int numTriangles, const Vec3i *triangles);
  TriMeshRef(const std::vector<Vec3d> &vertices, const std::vector<Vec3i> &triangles);
  TriMeshRef(int numVertices, const Vec3d *vertices, const std::vector<Vec3i> &triangles);

  int numVertices() const { return numVertices_; }
  int numTriangles() const { return numTriangles_; }

  const Vec3d &pos(int vtxID) const { return positions_[vtxID]; }
  const Vec3i &tri(int triID) const { return triangles_[triID]; }

  int triVtxID(int triID, int i) const { return triangles_[triID][i]; }
  const Vec3d &pos(int triID, int i) const { return positions_[triangles_[triID][i]]; }

  IndexedTriangle getIndexedTriangle(int triID) const;

  BasicAlgorithms::ArrayRef<Vec3d> positionsRef() const { return BasicAlgorithms::makeArrayRef(numVertices_, positions_); }
  BasicAlgorithms::ArrayRef<Vec3i> trianglesRef() const { return BasicAlgorithms::makeArrayRef(numTriangles_, triangles_); }

  const Vec3d *positionsPtr() const { return positions_; }
  const Vec3i *trianglesPtr() const { return triangles_; }

  const double *positionBuffer() const { return (const double *)positions_; }
  const int *triangleBuffer() const { return (const int *)triangles_; }

  const Vec3d *positions() const { return positions_; }
  const Vec3i *triangles() const { return triangles_; }

  std::vector<Vec3d> exportPositions() const { return std::vector<Vec3d>(positions_, positions_ + numVertices_); }
  std::vector<Vec3i> exportTriangles() const { return std::vector<Vec3i>(triangles_, triangles_ + numTriangles_); }
  void exportPositions(std::vector<Vec3d> &v) const { v.insert(v.end(), positions_, positions_ + numVertices_); }
  void exportTriangles(std::vector<Vec3i> &t) const { t.insert(t.end(), triangles_, triangles_ + numTriangles_); }

  template<class TriangleIDContainer>
  void getVerticesInTriangles(const TriangleIDContainer &triangleIDs, std::vector<int> &vertexIDs) const;

  // get bounding box for every triangle
  std::vector<BoundingBox> getTriangleBoundingBoxes() const;
  // compute bounding box for all the triangles
  // return a default BoundingBox if no triangles
  BoundingBox computeTriangleBoundingBox() const;

  Vec3d computeTriangleNormal(int triID) const;

  Vec3d computeAverageVertexPosition() const;

  // compute barycenter of each triangle
  Vec3d computeTriangleCentroid(int triID) const;

  Vec3d computeAverageTriangleCentroid() const;

  // x = V - E + F
  // for closed connected orientable manifold polygon, x = 2 - 2 g, where g = genus
  // for connected orientable manifold polygon with b boundaries, x = 2 - 2 g - b
  int computeEulerCharacteristic() const;
  // it the mesh is closed edge-manifold, then E = F * 3 / 2
  int computeEulerCharacteristicAssumingClosedManifold() const;

  // compute mass, center of mass and inertia tensor of the solid interior of the mesh
  // assuming interior has unit density
  void computeSolidInertiaParameters(double &mass, Vec3d &COM, Mat3d &inertiaTensor) const;

  // vtxID: [0, TriMesh::numVertices())
  double getTriangleAngleAtVertex(int triID, int vtxID) const;
  double getTriangleAngleAtVertexRobust(int triID, int vtxID) const;
  // see geometryQuery.h for definition of feature
  Vec3d getTriangleClosestPoint(int triID, const Vec3d &queryPos, int feature) const;

  double computeSurfaceArea() const;
  void computeTriangleSurfaceAreas(double *triangleSurfaceAreas) const;
  void computeVertexSurfaceAreas(double *vertexSurfaceAreas, double *triangleSurfaceAreas = nullptr) const;

  inline Vec3d getSurfacePositionWithBarycentricWeights(int triID, const Vec3d &weight) const;
  inline Vec3d getSurfacePositionWithBarycentricWeights(const Vec3i &vtxIDs, const Vec3d &weight) const;
  inline Vec3d getSurfacePositionWithBarycentricWeights(const Vec3i &vtxIDs, const Vec3d &weight, const Vec3d *vtxDisp) const;

  // use an O(n) method to compute winding number
  double computeWindingNumber(const Vec3d &pos) const;

  // use an O(n) method to compute distance, not exact
  // return quiet nan if no triangles
  double computeSquaredDistanceToPoint(const Vec3d &pos, int *closestTriangle = nullptr, int *feature = nullptr);

  // save to obj mesh
  bool save(const std::string &filename) const;

protected:
  int numVertices_ = 0, numTriangles_ = 0;
  const Vec3d *positions_ = nullptr;
  const Vec3i *triangles_ = nullptr;
};

// class to store basic triangle mesh data
class TriMeshGeo
{
public:
  TriMeshGeo() {}  // empty mesh
  TriMeshGeo(int numVertices, const double *vertices, int numTriangles, const int *triangles);
  TriMeshGeo(int numVertices, const Vec3d *vertices, int numTriangles, const Vec3i *triangles);
  TriMeshGeo(int numVertices, const Vec3d *vertices, std::vector<Vec3i> triangles);
  TriMeshGeo(std::vector<Vec3d> vertices, std::vector<Vec3i> triangles);
  TriMeshGeo(const TriMeshRef meshRef);

  int numVertices() const { return static_cast<int>(positions_.size()); }
  int numTriangles() const { return static_cast<int>(triangles_.size()); }

  const Vec3d &pos(int vtxID) const { return positions_[vtxID]; }
  Vec3d &pos(int vtxID) { return positions_[vtxID]; }
  const Vec3i &tri(int triID) const { return triangles_[triID]; }
  Vec3i &tri(int triID) { return triangles_[triID]; }

  int triVtxID(int triID, int i) const { return triangles_[triID][i]; }
  const Vec3d &pos(int triID, int i) const { return positions_[triangles_[triID][i]]; }
  Vec3d &pos(int triID, int i) { return positions_[triangles_[triID][i]]; }

  void addPos(const Vec3d &p) { positions_.push_back(p); }
  void addTri(const Vec3i &t) { triangles_.push_back(t); }
  void addMesh(const TriMeshRef mesh);
  void reserve(int nVtx, int nTri) { positions_.reserve(nVtx), triangles_.reserve(nTri); }
  void clear()
  {
    positions_.clear();
    triangles_.clear();
  }

  const std::vector<Vec3d> &positions() const { return positions_; }
  const std::vector<Vec3i> &triangles() const { return triangles_; }
  std::vector<Vec3d> &positions() { return positions_; }
  std::vector<Vec3i> &triangles() { return triangles_; }

  TriMeshRef ref() const { return { positions_, triangles_ }; }
  // implicit conversion
  operator TriMeshRef() const { return ref(); }

  bool load(const std::string &filename);
  // save to obj mesh
  bool save(const std::string &filename) const { return ref().save(filename); }

protected:
  std::vector<Vec3d> positions_;
  std::vector<Vec3i> triangles_;
};

// =========================================================
//              Single Triangle Utilities
// =========================================================
// they will eventually go to a triangle related header

// get vertex ID opposite the edge of (e0, e1)
// assume tri is valid and contains e0 and e1
int getTriangleVertexOppositeEdge(const Vec3i &tri, int e0, int e1);
int getTriangleVertexOppositeEdge(const Vec3i &tri, const UEdgeKey &edge);

// give an uedge, get the oedge on the tri
// assume tri is valid
// return an invalid key {-1,-1} if no such uedge on tri
OEdgeKey getTriangleOEdge(const Vec3i &tri, const UEdgeKey &uedge);

// get the shared vertex IDs from two triangles, assuming input t0 and t1 are valid triangles (v[i] >= 0 && v[i] != v[j])
// if no vtx shared, return {-1, -1, -1}
// if v0 is shared, return {v0, -1, -1}
// if v0, v1 are shared, return {v0, v1, -1}
// if v0, v1, v2 are shared, return {v0, v1, v2}
Vec3i getSharedVertices(const Vec3i &t0, const Vec3i &t1);

// if there are triangles that are valid:  v[i] >= 0 && v[i] != v[j]
inline bool isTriangleValid(const Vec3i &tri);
inline bool isTriangleInvalid(const Vec3i &tri);

// =========================================================
//                  Polygon Utilities
// =========================================================
// they will eventually go to a polygon related header

// given a polygon formed by a vertexID list, triangulate it
// return an empty vector if the polygon is degenerate
std::vector<Vec3i> triangulatePolygon(const std::vector<int> &polygon);
// if polygon is degenerate (contains 1 or 2 vtx), create a degenerate triangle
// if only one vtx (v0) in polygon, return a triangle of (v0, v0, v0)
// if only two vtx (v0, v1) in polygon, return a triangle of (v0, v1, v1)
std::vector<Vec3i> triangulatePolygonRobust(const std::vector<int> &polygon);

// remove vertices of the same positions
// return new vertices
std::vector<Vec3d> removeIdenticalVertices(const std::vector<Vec3d> &vertices,
  std::vector<int> *newVtxID2OldVtxID = nullptr, std::vector<int> *oldVtxID2NewVtxID = nullptr);

// =========================================================
//                  TriMesh Utilities
// =========================================================

// if there are triangles that are invalid:  v[i] < 0 || v[i] == v[j]
bool hasInvalidTriangles(const std::vector<Vec3i> &triangles);

// return triangles that are invalid
std::vector<Vec3i> getInvalidTriangles(const std::vector<Vec3i> &triangles);

// return triangles that are valid:  v[i] >=0 && v[i] != v[j]
std::vector<Vec3i> getOnlyValidTriangles(const std::vector<Vec3i> &triangles);
// also return valid triID in input triangles
std::vector<Vec3i> getOnlyValidTriangles(const std::vector<Vec3i> &triangles, std::vector<int> &validIDs);

template<class TriangleIDContainer>
void getVerticesInSelectedTriangles(const std::vector<Vec3i> &triangles, const TriangleIDContainer &triangleIDs,
  std::vector<int> &vertexIDs);
template<class TriangleIDContainer>
void getVerticesInSelectedTriangles(const Vec3i *triangles, const TriangleIDContainer &triangleIDs,
  std::vector<int> &vertexIDs);

template<class Vec3iContainer>
void getVerticesInTriangles(const Vec3iContainer &triangles, std::vector<int> &vertexIDs);

// remove those vertices in the mesh that do not contribute to a triangle
// the returned TriMesh has different vertices, but the triangle ordering is the same as meshRef
// newVtxID2OldVtxID has size #returnedMeshVertices
// oldVtxID2NewVtxID has size meshRef.numVertices()
TriMeshGeo removeIsolatedVertices(const TriMeshRef meshRef, std::vector<int> *newVtxID2OldVtxID = nullptr,
  std::vector<int> *oldVtxID2NewVtxID = nullptr);
// for a small sub mesh, it is wasteful to construct a vector for old->new mapping
// this function uses map to realize the old->new mapping
TriMeshGeo removeIsolatedVerticesWithMap(const TriMeshRef meshRef, std::vector<int> *newVtxID2OldVtxID = nullptr,
  std::map<int, int> *oldVtxID2NewVtxID = nullptr);

// keep only one vertex if several vertices are on the same position
TriMeshGeo removeIdenticalVertices(const TriMeshRef meshRef, std::vector<int> *newVtxID2OldVtxID = nullptr,
  std::vector<int> *oldVtxID2NewVtxID = nullptr);

TriMeshGeo removeInvalidTriangles(const TriMeshRef meshRef);

// triangleIDsToRemove must be ordered
TriMeshGeo removeTriangles(const TriMeshRef meshRef, const std::vector<int> &triangleIDsToRemove);

// oldVtxID2NewVtxID is of size meshRef.numVertices(), it creates a mapping from old vtxID to the new vtxID
// it asserts the mapping of oldVtxID2NewVtxID is valid, new VtxID are continuous, ranging [0, #newVtxID)
TriMeshGeo mergeVertices(const TriMeshRef meshRef, std::vector<int> &oldVtxID2NewVtxID);

// return a sub mesh which are triangles from selected triangleIDs
// isolated vertices in the sub mesh ARE NOT removed
// call removeIsolatedVertices to remove them from the result TriMesh
template<class TriangleIDContainer>
TriMeshGeo getSubTriMeshWithSameVertices(const TriMeshRef meshRef, const TriangleIDContainer &triangleIDs);
template<class TriangleIDContainer>
void getSubMesh(const TriMeshRef meshRef, const TriangleIDContainer &triangleIDs, std::vector<Vec3i> &outputTriangles);

template<class TriangleIDContainer>
TriMeshGeo getSubTriMesh(const TriMeshRef &mesh, const TriangleIDContainer &subTriIDs,
  std::vector<int> *subVtxID2FullVtxID = nullptr,
  std::map<int, int> *fullVtxID2SubVtxID = nullptr);

// get part of the triangles to form a sub mesh
// the triangles have vertices all from sortedSubVtxIDs
TriMeshGeo getSubTriMeshOnlyOnSortedSubVertexIDs(const TriMeshRef &mesh, BasicAlgorithms::ArrayRef<int> sortedSubVtxIDs);
// return triangle IDs in the sub mesh
void getSubTriMeshOnlyOnSortedSubVertexIDs(const TriMeshRef &mesh, BasicAlgorithms::ArrayRef<int> sortedSubVtxIDs, std::vector<int> &outputTriangleIDs);

// simple merge of mesh by retaining all vertices and triangles from both meshes
TriMeshGeo mergeMesh(const TriMeshRef mesh1, const TriMeshRef mesh2);

TriMeshGeo mergeMesh(size_t numMeshes, const TriMeshRef *meshes);

// =========================================================
//                  Implementations
// =========================================================

template<class TriangleIDContainer>
void getVerticesInSelectedTriangles(const std::vector<Vec3i> &triangles, const TriangleIDContainer &triangleIDs,
  std::vector<int> &vertexIDs)
{
  return getVerticesInSelectedTriangles(triangles.data(), triangleIDs, vertexIDs);
}

template<class TriangleIDContainer>
void getVerticesInSelectedTriangles(const Vec3i *triangles, const TriangleIDContainer &triangleIDs,
  std::vector<int> &vertexIDs)
{
  int numPrevIDs = vertexIDs.size();
  for (int triID : triangleIDs) {
    for (int j = 0; j < 3; j++)
      vertexIDs.push_back(triangles[triID][j]);
  }
  // remove duplicate vertex IDs
  std::sort(vertexIDs.begin() + numPrevIDs, vertexIDs.end());
  auto newEnd = std::unique(vertexIDs.begin() + numPrevIDs, vertexIDs.end());
  vertexIDs.resize(std::distance(vertexIDs.begin(), newEnd));
}

template<class Vec3iContainer>
void getVerticesInTriangles(const Vec3iContainer &triangles, std::vector<int> &vertexIDs)
{
  int numPrevIDs = vertexIDs.size();
  for (Vec3i t : triangles) {
    for (int j = 0; j < 3; j++)
      vertexIDs.push_back(t[j]);
  }
  // remove duplicate vertex IDs
  std::sort(vertexIDs.begin() + numPrevIDs, vertexIDs.end());
  auto newEnd = std::unique(vertexIDs.begin() + numPrevIDs, vertexIDs.end());
  vertexIDs.resize(std::distance(vertexIDs.begin(), newEnd));
}

template<class TriangleIDContainer>
void TriMeshRef::getVerticesInTriangles(const TriangleIDContainer &triangleIDs, std::vector<int> &vertexIDs) const
{
  getVerticesInSelectedTriangles(triangles_, triangleIDs, vertexIDs);
}

inline Vec3d TriMeshRef::getSurfacePositionWithBarycentricWeights(int triID, const Vec3d &w) const
{
  return pos(triID, 0) * w[0] + pos(triID, 1) * w[1] + pos(triID, 2) * w[2];
}

inline Vec3d TriMeshRef::getSurfacePositionWithBarycentricWeights(const Vec3i &vtxIDs, const Vec3d &w) const
{
  Vec3d p;
  p.setZero();
  for (int i = 0; i < 3; i++) {
    if (vtxIDs[i] < 0)
      continue;
    p += pos(vtxIDs[i]) * w[i];
  }
  return p;
}

inline Vec3d TriMeshRef::getSurfacePositionWithBarycentricWeights(const Vec3i &vtxIDs, const Vec3d &w, const Vec3d *vtxDisp) const
{
  Vec3d p;
  p.setZero();
  for (int i = 0; i < 3; i++) {
    if (vtxIDs[i] < 0)
      continue;
    p += (pos(vtxIDs[i]) + vtxDisp[vtxIDs[i]]) * w[i];
  }
  return p;
}

template<class TriangleIDContainer>
TriMeshGeo getSubTriMeshWithSameVertices(const TriMeshRef meshRef, const TriangleIDContainer &triangleIDs)
{
  std::vector<Vec3i> subTris;
  for (int triID : triangleIDs) {
    subTris.push_back(meshRef.tri(triID));
  }
  return TriMeshGeo(meshRef.numVertices(), meshRef.positions(), subTris.size(), subTris.data());
}

template<class TriangleIDContainer>
void getSubMesh(const TriMeshRef meshRef, const TriangleIDContainer &triangleIDs, std::vector<Vec3i> &outputTriangles)
{
  for (int triID : triangleIDs) {
    outputTriangles.push_back(meshRef.tri(triID));
  }
}

template<class TriangleIDContainer>
TriMeshGeo getSubTriMesh(const TriMeshRef &mesh, const TriangleIDContainer &subTriIDs,
  std::vector<int> *subVtxID2FullVtxID, std::map<int, int> *fullVtxID2SubVtxID)
{
  std::vector<int> usedVtxIDs;
  for (int triID : subTriIDs) {
    for (int i = 0; i < 3; i++)
      usedVtxIDs.push_back(mesh.tri(triID)[i]);
  }
  // remove duplicated vertices
  std::sort(usedVtxIDs.begin(), usedVtxIDs.end());
  auto newEnd = std::unique(usedVtxIDs.begin(), usedVtxIDs.end());
  usedVtxIDs.resize(std::distance(usedVtxIDs.begin(), newEnd));

  std::map<int, int> oldToNewVtxMap;
  for (size_t i = 0; i < usedVtxIDs.size(); i++)
    oldToNewVtxMap[usedVtxIDs[i]] = i;

  std::vector<Vec3i> newTris(subTriIDs.size());
  int newTriID = 0;
  for (int triID : subTriIDs) {
    for (int i = 0; i < 3; i++) {
      int newID = oldToNewVtxMap[mesh.tri(triID)[i]];
      newTris[newTriID][i] = newID;
    }
    newTriID++;
  }

  std::vector<Vec3d> newPos(usedVtxIDs.size());
  for (size_t i = 0; i < usedVtxIDs.size(); i++)
    newPos[i] = mesh.pos(usedVtxIDs[i]);
  if (subVtxID2FullVtxID)
    *subVtxID2FullVtxID = std::move(usedVtxIDs);
  if (fullVtxID2SubVtxID)
    *fullVtxID2SubVtxID = std::move(oldToNewVtxMap);

  return TriMeshGeo(std::move(newPos), std::move(newTris));
}

inline bool isTriangleValid(const Vec3i &t)
{
  return t[0] >= 0 && t[1] >= 0 && t[2] >= 0 && t[0] != t[1] && t[1] != t[2] && t[2] != t[0];
}

inline bool isTriangleInvalid(const Vec3i &t)
{
  return !isTriangleValid(t);
}

}  // namespace Mesh
}  // namespace pgo
