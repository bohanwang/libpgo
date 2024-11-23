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

#include "triMeshGeo.h"
#include "triKey.h"
#include "geometryQuery.h"

#include "basicAlgorithms.h"
#include "valueIndex.h"
#include "stringHelper.h"
#include "pgoLogging.h"

#include <iostream>
#include <cstring>
#include <cassert>
#include <fstream>
#include <cmath>
#include <functional>
#include <limits>
#include <unordered_set>
#include <map>
#include <numeric>
#include <iomanip>

/////////////////////////////////////////////////////////////////
//                   TriMeshRef functions
/////////////////////////////////////////////////////////////////

namespace pgo
{
namespace Mesh
{
TriMeshRef::TriMeshRef(int numVertices, const double *vertices, int numTriangles, const int *t):
  numVertices_(numVertices), numTriangles_(numTriangles)
{
  positions_ = (const Vec3d *)vertices;
  triangles_ = (const Vec3i *)t;
}

TriMeshRef::TriMeshRef(int numVertices, const Vec3d *vertices, int numTriangles, const int *t):
  numVertices_(numVertices), numTriangles_(numTriangles)
{
  positions_ = vertices;
  triangles_ = (const Vec3i *)t;
}

TriMeshRef::TriMeshRef(int numVertices, const Vec3d *vertices, int numTriangles, const Vec3i *t):
  numVertices_(numVertices), numTriangles_(numTriangles)
{
  positions_ = vertices;
  triangles_ = t;
}

TriMeshRef::TriMeshRef(const std::vector<Vec3d> &vertices, const std::vector<Vec3i> &triangles):
  numVertices_(static_cast<int>(vertices.size())), numTriangles_(static_cast<int>(triangles.size()))
{
  positions_ = vertices.data();
  triangles_ = triangles.data();
}

TriMeshRef::TriMeshRef(int numVertices, const Vec3d *vertices, const std::vector<Vec3i> &triangles):
  numVertices_(numVertices), numTriangles_(static_cast<int>(triangles.size()))
{
  positions_ = vertices;
  triangles_ = triangles.data();
}

IndexedTriangle TriMeshRef::getIndexedTriangle(int triID) const
{
  return IndexedTriangle(triID, tri(triID), pos(triID, 0), pos(triID, 1), pos(triID, 2));
}

std::vector<BoundingBox> TriMeshRef::getTriangleBoundingBoxes() const
{
  std::vector<BoundingBox> bbs;
  bbs.reserve(numTriangles());
  for (int i = 0; i < numTriangles(); i++) {
    bbs.emplace_back(positions_, triangles_[i]);
  }
  return bbs;
}

BoundingBox TriMeshRef::computeTriangleBoundingBox() const
{
  if (numTriangles() == 0)
    return BoundingBox();

  const int *ts = (int *)triangles();
  const int *te = ts + numTriangles() * 3;
  auto r = BasicAlgorithms::makeRange(ts, te);
  return BoundingBox(positions(), r);
}

Vec3d TriMeshRef::computeTriangleNormal(int triID) const
{
  Vec3d normal = (pos(triID, 1) - pos(triID, 0)).cross(pos(triID, 2) - pos(triID, 0));
  normal.normalize();
  return normal;
}

Vec3d TriMeshRef::computeAverageVertexPosition() const
{
  Vec3d center;
  center.setZero();

  if (numVertices() > 0) {
    for (int i = 0; i < numVertices(); i++)
      center += pos(i);
    center /= numVertices();
  }
  return center;
}

Vec3d TriMeshRef::computeTriangleCentroid(int triID) const
{
  return (1.0 / 3.0) * (pos(triID, 0) + pos(triID, 1) + pos(triID, 2));
}

Vec3d TriMeshRef::computeAverageTriangleCentroid() const
{
  Vec3d centroid;
  centroid.setZero();

  if (numTriangles() == 0)
    return centroid;
  for (int i = 0; i < numTriangles(); i++) {
    centroid += computeTriangleCentroid(i);
  }
  centroid /= numTriangles();
  return centroid;
}

double TriMeshRef::getTriangleAngleAtVertex(int triID, int vtxID) const
{
  const auto &tri = triangles_[triID];
  int i = getInvertedIndex(tri, vtxID);
  PGO_ALOG(i >= 0);

  int j = (i + 1) % 3, k = (i + 2) % 3;
  return getTriangleAngle(pos(tri[i]), pos(tri[j]), pos(tri[k]));
}

double TriMeshRef::getTriangleAngleAtVertexRobust(int triID, int vtxID) const
{
  const auto &tri = triangles_[triID];
  int i = getInvertedIndex(tri, vtxID);
  PGO_ALOG(i >= 0);

  int j = (i + 1) % 3, k = (i + 2) % 3;
  return getTriangleAngleRobust(pos(tri[i]), pos(tri[j]), pos(tri[k]));
}

Vec3d TriMeshRef::getTriangleClosestPoint(int triID, const Vec3d &queryPos, int feature) const
{
  return getClosestPointToTriangleWithFeature(queryPos, pos(triID, 0), pos(triID, 1), pos(triID, 2), feature);
}

double TriMeshRef::computeSurfaceArea() const
{
  double area = 0.0;
  for (int triID = 0; triID < numTriangles(); triID++) {
    area += getTriangleArea(pos(triID, 0), pos(triID, 1), pos(triID, 2));
  }
  return area;
}

void TriMeshRef::computeTriangleSurfaceAreas(double *triangleSurfaceAreas) const
{
  for (int triID = 0; triID < numTriangles(); triID++)
    triangleSurfaceAreas[triID] = getTriangleArea(pos(triID, 0), pos(triID, 1), pos(triID, 2));
}

void TriMeshRef::computeVertexSurfaceAreas(double *vertexSurfaceAreas, double *triangleSurfaceAreas) const
{
  memset(vertexSurfaceAreas, 0, sizeof(double) * numVertices());
  for (int triID = 0; triID < numTriangles(); triID++) {
    double triangleArea = getTriangleArea(pos(triID, 0), pos(triID, 1), pos(triID, 2));
    if (triangleSurfaceAreas)
      triangleSurfaceAreas[triID] = triangleArea;

    double vtxLocalArea = triangleArea * (1 / 3.0);
    for (int i = 0; i < 3; i++)
      vertexSurfaceAreas[triVtxID(triID, i)] += vtxLocalArea;
  }
}

double TriMeshRef::computeWindingNumber(const Vec3d &queryPos) const
{
  double w = 0;
  for (int j = 0; j < numTriangles(); j++) {
    Vec3d a = pos(j, 0) - queryPos;
    Vec3d b = pos(j, 1) - queryPos;
    Vec3d c = pos(j, 2) - queryPos;
    double la = a.norm(), lb = b.norm(), lc = c.norm();
    Mat3d mat = asMat3d(a, b, c);
    double omega = 2 * std::atan2(mat.determinant(), (la * lb * lc + a.dot(b) * lc + b.dot(c) * la + c.dot(a) * lb));
    PGO_ALOG(std::isnan(omega) == false);
    w += omega;
  }

  w /= 4 * M_PI;
  return w;
}

bool saveToAscii(const TriMeshRef &mesh, const std::string &filename)
{
  // open file
  std::ofstream fout(filename.c_str());

  if (!fout) {
    std::cout << "Error: could not write to file " << filename << "." << std::endl;
    return false;
  }

  fout << "# Generated by the TriMeshRef class" << std::endl;
  fout << "# Number of vertices: " << mesh.numVertices() << std::endl;
  fout << "# Number of faces: " << mesh.numTriangles() << std::endl;

  // vertices...
  for (int i = 0; i < mesh.numVertices(); i++) {
    const Vec3d &pos = mesh.pos(i);
    fout << std::setprecision(17) << "v " << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;
  }

  // groups and faces...
  //  fout << "g default"<< endl;
  for (int i = 0; i < mesh.numTriangles(); i++) {
    const Vec3i &t = mesh.tri(i);
    fout << "f " << t[0] + 1 << " " << t[1] + 1 << " " << t[2] + 1 << std::endl;
  }

  std::cout << "Saved mesh (#v: " << mesh.numVertices() << ", #t: " << mesh.numTriangles() << ") to " << filename << "." << std::endl;
  fout.close();
  return true;
}

bool saveToBinary(const TriMeshRef &mesh, const std::string &filename)
{
  FILE *fout = fopen(filename.c_str(), "wb");
  if (fout == nullptr) {
    std::cout << "Error: could not write to file " << filename << std::endl;
    return false;
  }

  // first pass: count the total number of bytes to be written to the file
  // second pass: do the actual writing
  enum
  {
    COUNT_BYTES,
    WRITE_TO_DISK,
    NUM_PASSES
  };
  int totalPasses = NUM_PASSES;

  unsigned int totalBytes = 0;
  for (int pass = 0; pass < totalPasses; pass++) {
    size_t bytes = 0;
    unsigned int items = 0;

    // the header will be the number of bytes (including the totalbytes itself)
    if (pass == WRITE_TO_DISK && fwrite(&totalBytes, sizeof(unsigned int), 1, fout) != 1)
      return false;
    bytes += sizeof(unsigned int);

    // save the flag that determines whether to output materials or not
    int outputMaterials = 0;
    if (pass == WRITE_TO_DISK && fwrite(&outputMaterials, sizeof(int), 1, fout) != 1)
      return 1;
    bytes += sizeof(int);

    // save the number of vertices
    unsigned int numVertices = mesh.numVertices();
    if (pass == WRITE_TO_DISK && fwrite(&numVertices, sizeof(unsigned int), 1, fout) != 1)
      return 1;
    bytes += sizeof(unsigned int);

    // save vertices
    items = 3 * numVertices;
    if (pass == WRITE_TO_DISK && fwrite(mesh.positions(), sizeof(double), items, fout) != items)
      return 1;
    bytes += items * sizeof(double);

    // save the number of texture coordinates
    unsigned int numTexCoordinates = 0;
    if (pass == WRITE_TO_DISK && fwrite(&numTexCoordinates, sizeof(unsigned int), 1, fout) != 1)
      return 1;
    bytes += sizeof(unsigned int);

    // save the number of normals
    unsigned int numNormals = 0;
    if (pass == WRITE_TO_DISK && fwrite(&numNormals, sizeof(unsigned int), 1, fout) != 1)
      return 1;
    bytes += sizeof(unsigned int);

    // save the number of groups
    unsigned int numGroups = 1;
    if (pass == WRITE_TO_DISK && fwrite(&numGroups, sizeof(unsigned int), 1, fout) != 1)
      return 1;
    bytes += sizeof(unsigned int);

    // save group name
    const char *groupName = "defaultGroup";
    unsigned int strLength = static_cast<unsigned int>(strlen(groupName));
    if (pass == WRITE_TO_DISK && fwrite(&strLength, sizeof(unsigned int), 1, fout) != 1)
      return 1;
    bytes += sizeof(unsigned int);

    items = strLength;
    if (pass == WRITE_TO_DISK && fwrite(groupName, sizeof(char), items, fout) != items)
      return 1;
    bytes += items * sizeof(char);

    // save the number of faces of the current group
    unsigned int numFaces = mesh.numTriangles();
    if (pass == WRITE_TO_DISK && fwrite(&numFaces, sizeof(unsigned int), 1, fout) != 1)
      return 1;
    bytes += sizeof(unsigned int);

    // save the number of vertices of each face in current group
    const unsigned int numVtxPerFace = 3;
    if (pass == WRITE_TO_DISK) {
      for (int i = 0; i < mesh.numTriangles(); i++) {
        if (fwrite(&numVtxPerFace, sizeof(unsigned int), 1, fout) != 1)
          return 1;
      }
    }
    bytes += mesh.numTriangles() * sizeof(unsigned int);

    // save the vertices of each face
    items = mesh.numTriangles() * 3;
    std::vector<unsigned int> vtxIndices(mesh.numTriangles() * 3);
    for (int i = 0; i < mesh.numTriangles(); i++)
      for (int j = 0; j < 3; j++)
        vtxIndices[i * 3 + j] = mesh.triVtxID(i, j) + 1;  // 1-index
    if (pass == WRITE_TO_DISK && fwrite(vtxIndices.data(), sizeof(unsigned int), items, fout) != items)
      return 1;
    bytes += items * sizeof(unsigned int);

    std::vector<unsigned int> buffer(mesh.numTriangles() * 3, 0.0);

    // save texture coordinates
    if (pass == WRITE_TO_DISK && fwrite(buffer.data(), sizeof(unsigned int), items, fout) != items)
      return 1;
    bytes += items * sizeof(unsigned int);

    // save normal coordinates
    if (pass == WRITE_TO_DISK && fwrite(buffer.data(), sizeof(unsigned int), items, fout) != items)
      return 1;
    bytes += items * sizeof(unsigned int);

    if (pass == COUNT_BYTES)
      totalBytes = (unsigned int)bytes;
  }  // for pass

  fclose(fout);
  return true;
}

bool TriMeshRef::save(const std::string &filename) const
{
  if (BasicAlgorithms::iendWith(filename, ".objb"))
    return saveToBinary(*this, filename);
  if (BasicAlgorithms::iendWith(filename, ".obj"))
    return saveToAscii(*this, filename);

  printf("Unknown file extension when saving %s, try ASCII format...\n", filename.c_str());
  return saveToAscii(*this, filename);
}

/////////////////////////////////////////////////////////////////
//                   TriMesh functions
/////////////////////////////////////////////////////////////////

TriMeshGeo::TriMeshGeo(int numVertices, const double *vertices, int numTriangles, const int *t):
  positions_(numVertices), triangles_(numTriangles)
{
  memcpy(&positions_[0][0], vertices, sizeof(double) * numVertices * 3);
  memcpy(&triangles_[0][0], t, sizeof(int) * numTriangles * 3);
}

TriMeshGeo::TriMeshGeo(int numVertices, const Vec3d *vertices, int numTriangles, const Vec3i *t):
  positions_(numVertices), triangles_(numTriangles)
{
  memcpy(positions_[0].data(), vertices, sizeof(double) * numVertices * 3);
  memcpy(triangles_[0].data(), t, sizeof(int) * numTriangles * 3);
}

TriMeshGeo::TriMeshGeo(int numVertices, const Vec3d *vertices, std::vector<Vec3i> triangles):
  positions_(numVertices)
{
  memcpy(positions_[0].data(), vertices, sizeof(double) * numVertices * 3);
  triangles_ = triangles;
}

TriMeshGeo::TriMeshGeo(std::vector<Vec3d> vertices, std::vector<Vec3i> triangles):
  positions_(std::move(vertices)), triangles_(std::move(triangles))
{
}

TriMeshGeo::TriMeshGeo(const TriMeshRef meshRef):
  TriMeshGeo(meshRef.numVertices(), meshRef.positions(), meshRef.numTriangles(), meshRef.triangles())
{
}

/////////////////////////////////////////////////////////////////
//                   Other functions
/////////////////////////////////////////////////////////////////

int getTriangleVertexOppositeEdge(const Vec3i &tri, int e0, int e1)
{
  for (int i = 0; i < 3; i++) {
    if (tri[i] != e0 && tri[i] != e1)
      return tri[i];
  }
  PGO_ALOG(0);
  return -1;
}

int getTriangleVertexOppositeEdge(const Vec3i &tri, const UEdgeKey &edge)
{
  return getTriangleVertexOppositeEdge(tri, edge[0], edge[1]);
}

OEdgeKey getTriangleOEdge(const Vec3i &tri, const UEdgeKey &uedge)
{
  for (int i = 0; i < 3; i++) {
    UEdgeKey curUEdge(tri[i], tri[(i + 1) % 3]);
    if (curUEdge == uedge)
      return OEdgeKey(tri[i], tri[(i + 1) % 3]);
  }
  return OEdgeKey();
}

Vec3i getSharedVertices(const Vec3i &t0, const Vec3i &t1)
{
  Vec3i ret(-1);
  int index = 0;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      if (t0[i] == t1[j]) {
        ret[index++] = t0[i];
        break;
      }
    }
  }
  return ret;
}

std::vector<Vec3i> triangulatePolygon(const std::vector<int> &polygon)
{
  std::vector<Vec3i> ret;
  for (size_t i = 2; i < polygon.size(); i++) {
    ret.emplace_back(polygon[0], polygon[i - 1], polygon[i]);
  }
  return ret;
}

std::vector<Vec3i> triangulatePolygonRobust(const std::vector<int> &polygon)
{
  if (polygon.size() > 2)
    return triangulatePolygon(polygon);
  if (polygon.size() == 2)
    return { Vec3i(polygon[0], polygon[1], polygon[1]) };
  if (polygon.size() == 1)
    return { Vec3i(polygon[0]) };
  return {};
}

std::vector<Vec3d> removeIdenticalVertices(const std::vector<Vec3d> &vertices,
  std::vector<int> *newVtxID2OldVtxID, std::vector<int> *oldVtxID2NewVtxID)
{
  if (vertices.size() == 0) {
    if (newVtxID2OldVtxID)
      newVtxID2OldVtxID->clear();
    if (oldVtxID2NewVtxID)
      oldVtxID2NewVtxID->clear();
    return {};
  }

  std::vector<int> indices(vertices.size());  // eventually stores newVtxID->oldVtxID
  std::iota(indices.begin(), indices.end(), 0);
  std::sort(indices.begin(), indices.end(), [&](int a, int b) { return vertices[a] < vertices[b]; });

  std::vector<Vec3d> ret;
  ret.push_back(vertices[indices[0]]);
  if (newVtxID2OldVtxID)
    *newVtxID2OldVtxID = { indices[0] };
  if (oldVtxID2NewVtxID)
    oldVtxID2NewVtxID->assign(vertices.size(), -1);

  for (size_t i = 0, j = 1; j < indices.size(); j++) {
    if (vertices[indices[i]] == vertices[indices[j]]) {
      if (oldVtxID2NewVtxID)
        oldVtxID2NewVtxID->at(indices[j]) = (int)ret.size() - 1;
    }
    else {
      if (oldVtxID2NewVtxID)
        oldVtxID2NewVtxID->at(indices[j]) = (int)ret.size();
      if (newVtxID2OldVtxID)
        newVtxID2OldVtxID->push_back(indices[j]);

      ret.push_back(vertices[indices[j]]);
      i = j - 1;
    }
  }

  if (oldVtxID2NewVtxID)
    PGO_ALOG(find(oldVtxID2NewVtxID->begin(), oldVtxID2NewVtxID->end(), -1) == oldVtxID2NewVtxID->end());
  return ret;
}

bool hasInvalidTriangles(const std::vector<Vec3i> &triangles)
{
  for (const auto &t : triangles) {
    if (isTriangleInvalid(t))
      return true;
  }
  return false;
}

std::vector<Vec3i> getInvalidTriangles(const std::vector<Vec3i> &triangles)
{
  std::vector<Vec3i> ret;
  for (const auto &t : triangles) {
    if (isTriangleInvalid(t))
      ret.push_back(t);
  }
  return ret;
}

std::vector<Vec3i> getOnlyValidTriangles(const std::vector<Vec3i> &triangles)
{
  std::vector<Vec3i> ret;
  for (const auto &t : triangles) {
    if (isTriangleValid(t))
      ret.push_back(t);
  }
  return ret;
}

std::vector<Vec3i> getOnlyValidTriangles(const std::vector<Vec3i> &triangles, std::vector<int> &validIDs)
{
  validIDs.clear();
  std::vector<Vec3i> ret;
  for (size_t i = 0; i < triangles.size(); i++) {
    const auto &t = triangles[i];
    if (isTriangleValid(t)) {
      ret.push_back(t);
      validIDs.push_back(static_cast<int>(i));
    }
  }
  return ret;
}

TriMeshGeo removeIsolatedVertices(const TriMeshRef meshRef, std::vector<int> *newVtxID2OldVtxID, std::vector<int> *oldVtxID2NewVtxID)
{
  std::vector<int> usedVtxIDs;
  for (int triID = 0; triID < meshRef.numTriangles(); triID++) {
    for (int i = 0; i < 3; i++)
      usedVtxIDs.push_back(meshRef.tri(triID)[i]);
  }
  BasicAlgorithms::sortAndDeduplicate(usedVtxIDs);

  if (usedVtxIDs.size() == (size_t)meshRef.numVertices()) {
    // now usedVtxIDs contains [0, meshRef.numVertices)
    if (newVtxID2OldVtxID) {
      *newVtxID2OldVtxID = usedVtxIDs;
    }
    if (oldVtxID2NewVtxID) {
      *oldVtxID2NewVtxID = std::move(usedVtxIDs);
    }
    return TriMeshGeo(meshRef);
  }

  std::vector<int> oldToNewVtxMap(meshRef.numVertices(), -1);
  for (size_t i = 0; i < usedVtxIDs.size(); i++) {
    oldToNewVtxMap[usedVtxIDs[i]] = static_cast<int>(i);
  }

  std::vector<Vec3i> newTriangles(meshRef.numTriangles());
  for (int triID = 0; triID < meshRef.numTriangles(); triID++) {
    for (int i = 0; i < 3; i++) {
      int newID = oldToNewVtxMap[meshRef.tri(triID)[i]];
      PGO_ALOG(newID >= 0);
      newTriangles[triID][i] = newID;
    }
  }

  std::vector<Vec3d> newPos(usedVtxIDs.size());
  for (size_t i = 0; i < usedVtxIDs.size(); i++)
    newPos[i] = meshRef.pos(usedVtxIDs[i]);
  if (newVtxID2OldVtxID) {
    *newVtxID2OldVtxID = std::move(usedVtxIDs);
  }
  if (oldVtxID2NewVtxID) {
    *oldVtxID2NewVtxID = std::move(oldToNewVtxMap);
  }

  return TriMeshGeo(std::move(newPos), std::move(newTriangles));
}

TriMeshGeo removeIsolatedVerticesWithMap(const TriMeshRef meshRef, std::vector<int> *newVtxID2OldVtxID, std::map<int, int> *oldVtxID2NewVtxID)
{
  std::vector<int> usedVtxIDs;
  for (int triID = 0; triID < meshRef.numTriangles(); triID++) {
    for (int i = 0; i < 3; i++)
      usedVtxIDs.push_back(meshRef.tri(triID)[i]);
  }
  BasicAlgorithms::sortAndDeduplicate(usedVtxIDs);
  if (usedVtxIDs.size() == (size_t)meshRef.numVertices()) {
    // now usedVtxIDs contains [0, meshRef.numVertices)
    if (newVtxID2OldVtxID) {
      *newVtxID2OldVtxID = std::move(usedVtxIDs);
    }
    if (oldVtxID2NewVtxID) {
      for (int i = 0; i < meshRef.numVertices(); i++)
        (*oldVtxID2NewVtxID)[i] = i;
    }
    return TriMeshGeo(meshRef);
  }

  std::map<int, int> oldToNewVtxMap;
  for (size_t i = 0; i < usedVtxIDs.size(); i++) {
    oldToNewVtxMap[usedVtxIDs[i]] = static_cast<int>(i);
  }

  std::vector<Vec3i> newTriangles(meshRef.numTriangles());
  for (int triID = 0; triID < meshRef.numTriangles(); triID++) {
    for (int i = 0; i < 3; i++) {
      int newID = oldToNewVtxMap[meshRef.tri(triID)[i]];
      newTriangles[triID][i] = newID;
    }
  }

  std::vector<Vec3d> newPos(usedVtxIDs.size());
  for (size_t i = 0; i < usedVtxIDs.size(); i++)
    newPos[i] = meshRef.pos(usedVtxIDs[i]);
  if (newVtxID2OldVtxID) {
    *newVtxID2OldVtxID = std::move(usedVtxIDs);
  }
  if (oldVtxID2NewVtxID) {
    *oldVtxID2NewVtxID = std::move(oldToNewVtxMap);
  }

  return TriMeshGeo(std::move(newPos), std::move(newTriangles));
}

TriMeshGeo removeIdenticalVertices(const TriMeshRef meshRef, std::vector<int> *newVtxID2OldVtxID,
  std::vector<int> *oldVtxID2NewVtxID)
{
  if (meshRef.numVertices() == 0) {
    if (newVtxID2OldVtxID)
      newVtxID2OldVtxID->clear();
    if (oldVtxID2NewVtxID)
      oldVtxID2NewVtxID->clear();
    return TriMeshGeo(meshRef);
  }

  std::vector<int> indices(meshRef.numVertices());  // eventually stores newVtxID->oldVtxID
  std::iota(indices.begin(), indices.end(), 0);
  std::sort(indices.begin(), indices.end(), [&](int a, int b) { return meshRef.pos(a) < meshRef.pos(b); });

  std::vector<Vec3d> newVertices;
  newVertices.push_back(meshRef.pos(indices[0]));
  std::vector<int> old2new(meshRef.numVertices(), -1);
  if (newVtxID2OldVtxID)
    *newVtxID2OldVtxID = { indices[0] };
  old2new[indices[0]] = 0;

  for (size_t i = 0, j = 1; j < indices.size(); j++) {
    if (meshRef.pos(indices[i]) == meshRef.pos(indices[j])) {
      old2new[indices[j]] = static_cast<int>(newVertices.size()) - 1;
    }
    else {
      old2new[indices[j]] = static_cast<int>(newVertices.size());
      newVertices.push_back(meshRef.pos(indices[j]));
      if (newVtxID2OldVtxID) {
        newVtxID2OldVtxID->push_back(indices[j]);
      }
      i = j - 1;
    }
  }

  PGO_ALOG(find(old2new.begin(), old2new.end(), -1) == old2new.end());

  std::vector<Vec3i> newTri(meshRef.numTriangles());
  for (int i = 0; i < meshRef.numTriangles(); i++) {
    for (int j = 0; j < 3; j++)
      newTri[i][j] = old2new[meshRef.tri(i)[j]];
  }
  if (oldVtxID2NewVtxID)
    *oldVtxID2NewVtxID = old2new;

  return TriMeshGeo(std::move(newVertices), std::move(newTri));
}

TriMeshGeo removeInvalidTriangles(const TriMeshRef meshRef)
{
  std::vector<int> invalidTriIDs;  // invalid triangles
  for (int i = 0; i < meshRef.numTriangles(); i++) {
    if (isTriangleInvalid(meshRef.tri(i)))
      invalidTriIDs.push_back(i);
  }
  std::vector<Vec3i> triangles = meshRef.exportTriangles();
  BasicAlgorithms::removeByIndices(triangles, invalidTriIDs);
  return { meshRef.numVertices(), meshRef.positions(), std::move(triangles) };
}

TriMeshGeo removeTriangles(const TriMeshRef meshRef, const std::vector<int> &triangleIDsToRemove)
{
  std::vector<Vec3i> triangles = meshRef.exportTriangles();
  BasicAlgorithms::removeByIndices(triangles, triangleIDsToRemove);
  std::vector<Vec3d> pos = meshRef.exportPositions();
  return { std::move(pos), std::move(triangles) };
}

TriMeshGeo mergeVertices(const TriMeshRef meshRef, std::vector<int> &oldVtxID2NewVtxID)
{
  PGO_ALOG(oldVtxID2NewVtxID.size() == (size_t)meshRef.numVertices());
  std::vector<Vec3i> newTriangles = meshRef.exportTriangles();
  for (int triID = 0; triID < meshRef.numTriangles(); triID++) {
    for (int i = 0; i < 3; i++)
      newTriangles[triID][i] = oldVtxID2NewVtxID[meshRef.triVtxID(triID, i)];
  }
  std::vector<Vec3d> newVertices;
  for (int vID = 0; vID < meshRef.numVertices(); vID++) {
    int newID = oldVtxID2NewVtxID[vID];
    if (newID >= (int)newVertices.size())
      newVertices.resize(newID + 1);
    newVertices[newID] = meshRef.pos(vID);
  }
  return { std::move(newVertices), std::move(newTriangles) };
}

TriMeshGeo mergeMesh(const TriMeshRef mesh1, const TriMeshRef mesh2)
{
  std::vector<Vec3d> vertices = mesh1.exportPositions();
  std::vector<Vec3i> triangles = mesh1.exportTriangles();
  int vtxStart = static_cast<int>(vertices.size());
  int triStart = static_cast<int>(triangles.size());
  mesh2.exportPositions(vertices);
  mesh2.exportTriangles(triangles);
  for (size_t i = triStart; i < triangles.size(); i++) {
    for (int j = 0; j < 3; j++)
      triangles[i][j] += vtxStart;
  }
  return { std::move(vertices), std::move(triangles) };
}

TriMeshGeo mergeMesh(size_t numMeshes, const TriMeshRef *meshes)
{
  std::vector<Vec3d> newVtx;
  std::vector<Vec3i> newTri;

  int offset = 0;
  for (size_t i = 0; i < numMeshes; i++) {
    for (int vi = 0; vi < meshes[i].numVertices(); vi++) {
      newVtx.push_back(meshes[i].pos(vi));
    }

    for (int fi = 0; fi < meshes[i].numTriangles(); fi++) {
      Vec3i tri = meshes[i].tri(fi);
      tri[0] += offset;
      tri[1] += offset;
      tri[2] += offset;
      newTri.push_back(tri);
    }

    offset += meshes[i].numVertices();
  }

  return TriMeshGeo(std::move(newVtx), std::move(newTri));
}

double TriMeshRef::computeSquaredDistanceToPoint(const Vec3d &pt, int *closestTriangle, int *feature)
{
  static_assert(std::numeric_limits<double>::has_quiet_NaN, "No quiet NaN");
  if (numTriangles() == 0) {
    if (closestTriangle)
      *closestTriangle = -1;
    if (feature)
      *feature = -1;
    return std::numeric_limits<double>::quiet_NaN();
  }

  MinValueKey<std::pair<int, int>> vi(std::make_pair(-1, -1));
  for (int i = 0; i < numTriangles(); i++) {
    int f = -1;
    double d2 = getSquaredDistanceToTriangle(pt, pos(i, 0), pos(i, 1), pos(i, 2), f);
    vi.update(d2, std::make_pair(i, f));
  }
  PGO_ALOG(vi.key.first >= 0);
  if (closestTriangle)
    *closestTriangle = vi.key.first;
  if (feature)
    *feature = vi.key.second;
  return vi.value;
}

int TriMeshRef::computeEulerCharacteristic() const
{
  int ret = numVertices() + numTriangles();
  std::unordered_set<UEdgeKey> ues;
  for (int i = 0; i < numTriangles(); i++) {
    for (int j = 0; j < 3; j++)
      ues.emplace(triangles_[i][j], triangles_[i][(j + 1) % 3]);
  }
  ret -= static_cast<int>(ues.size());
  return ret;
}

int TriMeshRef::computeEulerCharacteristicAssumingClosedManifold() const
{
  PGO_ALOG(numTriangles() % 2 == 0);
  return numVertices() - numTriangles() / 2;
}

void TriMeshRef::computeSolidInertiaParameters(double &mass, Vec3d &COM, Mat3d &inertiaTensor) const
{
  mass = 0.0;
  COM.setZero();

  Mat3d originTensor;  // inertia tensor around origin
  originTensor.setZero();

  for (int triID = 0; triID < numTriangles_; triID++) {
    Vec3d a;
    a.setZero();

    const Vec3d &b = pos(triID, 0);
    const Vec3d &c = pos(triID, 1);
    const Vec3d &d = pos(triID, 2);
    double td = getTetDeterminant(a, b, c, d);
    mass += td;
    COM += td * getTetCenterOfMass(a, b, c, d);
    originTensor += getTetInertiaTensorAroudOrigin(a, b, c, d, td);
  }

  COM /= mass;
  mass /= 6.0;  // because tetDet = 6 * signedVolume

  inertiaTensor = shiftInertiaTensorAroundReferencePointToMassCenter(originTensor, mass, COM);
}

void TriMeshGeo::addMesh(const TriMeshRef mesh)
{
  int numPrevVtx = static_cast<int>(positions_.size());
  auto meshPos = mesh.positionsRef();
  insertRangeToEnd(positions_, meshPos);
  for (int triID = 0; triID < mesh.numTriangles(); triID++) {
    const Vec3i &tri = mesh.tri(triID);
    triangles_.emplace_back(tri[0] + numPrevVtx, tri[1] + numPrevVtx, tri[2] + numPrevVtx);
  }
}

TriMeshGeo getSubTriMeshOnlyOnSortedSubVertexIDs(const TriMeshRef &mesh, BasicAlgorithms::ArrayRef<int> sortedSubVtxIDs)
{
  std::map<int, int> oldToNewVtxMap;
  for (int i = 0; i < sortedSubVtxIDs.size(); i++)
    oldToNewVtxMap[sortedSubVtxIDs[i]] = i;

  std::vector<Vec3i> newTris;
  for (int triID = 0; triID < mesh.numTriangles(); triID++) {
    Vec3i newTri;
    int i = 0;
    for (; i < 3; i++) {
      int vtxID = mesh.triVtxID(triID, i);
      auto it = BasicAlgorithms::binarySearchFind(sortedSubVtxIDs.begin(), sortedSubVtxIDs.end(), vtxID);
      if (it != sortedSubVtxIDs.end())
        newTri[i] = static_cast<int>(std::distance(sortedSubVtxIDs.begin(), it));
      else
        break;
    }
    if (i == 3)  // all vtxIDs in the tri are in sortedSubVtxIDs
      newTris.push_back(newTri);
  }

  std::vector<Vec3d> newPos(sortedSubVtxIDs.size());
  for (int i = 0; i < sortedSubVtxIDs.size(); i++)
    newPos[i] = mesh.pos(sortedSubVtxIDs[i]);

  return TriMeshGeo(std::move(newPos), std::move(newTris));
}

void getSubTriMeshOnlyOnSortedSubVertexIDs(const TriMeshRef &mesh, BasicAlgorithms::ArrayRef<int> sortedSubVtxIDs, std::vector<int> &outputTriangleIDs)
{
  for (int triID = 0; triID < mesh.numTriangles(); triID++) {
    int i = 0;
    for (; i < 3; i++) {
      int vtxID = mesh.triVtxID(triID, i);
      if (BasicAlgorithms::binarySearchFound(sortedSubVtxIDs, vtxID) == false)
        break;
    }
    if (i == 3)  // all vtxIDs in the tri are in sortedSubVtxIDs
      outputTriangleIDs.push_back(triID);
  }
}

}  // namespace Mesh
}  // namespace pgo
