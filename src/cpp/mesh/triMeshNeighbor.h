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
#include "triMeshGeo.h"
#include "edgeKey.h"
#include "arrayRef.h"

#include <unordered_map>
#include <map>

namespace pgo
{
namespace Mesh
{

// Edge-manifold neighboring structure for vector of triangles.
// No vtx positions needed.
class TriangleNeighbor
{
public:
  TriangleNeighbor() {}  // empty
  // throw int 1 if triangles are not edge-manifold or contain invalid triangles
  TriangleNeighbor(BasicAlgorithms::ArrayRef<Vec3i> triangles);

  const Vec3i &getTriangleNeighbors(int triID) const { return triNbrs[triID]; }

  // input triangles should be the same one used to construct this neighbor class
  // return vector of <boundary triID, its boundary OEdgeKey> pairs
  std::vector<std::pair<int, OEdgeKey>> findBoundaryTriangles(BasicAlgorithms::ArrayRef<Vec3i> triangles) const;

  // input triangles should be the same one used to construct this neighbor class
  // each vector<int> in the returned data contains one loop. The ordering agrees with the boundary edge direction
  std::vector<std::vector<int>> findBoundaryLoops(BasicAlgorithms::ArrayRef<Vec3i> triangles) const;

  // get the triangles arround the boundary vertex vtxID, which starts at a triangle with OEdgeKey(prevVtxID, vtxID)
  std::vector<int> findTrianglesArroundBoundaryVertex(int prevVtxID, int vtxID,
    BasicAlgorithms::ArrayRef<Vec3i> triangles) const;

  // get the triangle with ordered edge
  // return -1 if not found
  int getTriangleAtEdge(const OEdgeKey &edge) const;

  // get the vtx opposite to the ordered edge on the triangle with the ordered edge
  // return -1 if not found
  int getVertexOppositeToEdge(const OEdgeKey &edge, BasicAlgorithms::ArrayRef<Vec3i> triangles) const;

protected:
  int numTriangles = 0;
  std::vector<Vec3i> triNbrs;                  // triID -> i[0,3) -> nbring triID at edge <tri[i], tri[(i+1)%3]>
  std::unordered_map<OEdgeKey, int> oedgeTri;  // <v0,v1> -> tri with ordered edge <v0, v1>
};

// Edge-manifold neighboring structure for triMesh
class TriMeshNeighbor
{
public:
  TriMeshNeighbor() {}  // empty
  // throw int 1 if triMesh are not edge-manifold or contain invalid triangles
  TriMeshNeighbor(TriMeshRef triMesh);

  const std::vector<int> &getVtxNearbyTriangles(int vtxID) const { return vtxNbrTri[vtxID]; }

  const Vec3i &getTriangleNeighbors(int triID) const { return triNbrs[triID]; }

  // input triangles should be the same one used to construct this neighbor class
  // return vector of <boundary triID, its boundary OEdgeKey> pairs
  std::vector<std::pair<int, OEdgeKey>> findBoundaryTriangles(BasicAlgorithms::ArrayRef<Vec3i> triangles) const;

  // input triangles should be the same one used to construct this neighbor class
  // each vector<int> in the returned data contains one loop. The ordering agrees with the boundary edge direction
  std::vector<std::vector<int>> findBoundaryLoops(BasicAlgorithms::ArrayRef<Vec3i> triangles) const;

  // for vtxID, get a mapping: nbring vtxID -> the triangleID which has the ordered edge of <vtxID, nbring vtxID>
  const std::map<int, int> &getVertexEdgeTriMap(int vtxID) const { return vtxEdgeTri[vtxID]; }

  // result are sorted
  std::vector<int> getVtxNearbyVertices(int vtxID, TriMeshRef mesh) const;
  bool areVerticesNeighbors(int vtxID0, int vtxID1) const;

  // get the triangle with ordered edge
  // return -1 if not found
  int getTriangleAtEdge(const OEdgeKey &edge) const;

  // input triangles should be the same one used to construct this neighbor class
  // get the vtx opposite to the ordered edge on the triangle with the ordered edge
  // return -1 if not found
  int getVertexOppositeToEdge(const OEdgeKey &edge, BasicAlgorithms::ArrayRef<Vec3i> triangles) const;

protected:
  int numVertices = 0, numTriangles = 0;
  std::vector<std::vector<int>> vtxNbrTri;     // vtxID -> nearby triIDs
  std::vector<Vec3i> triNbrs;                  // triID -> i[0,3) -> nbring triID at edge <tri[i], tri[(i+1)%3]>
  std::vector<std::map<int, int>> vtxEdgeTri;  // v0 -> v1 -> tri with ordered edge <v0, v1>
};

// Possibly non-manifold neighboring structure for vector of triangles
// No vtx positions needed
class NonManifoldTriangleNeighbor
{
public:
  NonManifoldTriangleNeighbor() {}  // empty
  // throw int 1 if triangles contain invalid triangles
  NonManifoldTriangleNeighbor(const std::vector<Vec3i> &triangles);

  // return -1 if not found
  std::vector<int> getTriangleAtOEdge(const OEdgeKey &oedge) const;

  // get ordered edges that are non-manifold (more than one triangles own this oedge)
  std::vector<OEdgeKey> findNonManifoldOEdges() const;

protected:
  int numTriangles = 0;
  std::unordered_map<OEdgeKey, std::vector<int>> oedgeTris;
};

// return if triangles are edge-manifold
// return false if triangles contain invalid or degenerate triangles
bool areTrianglesEdgeManifold(BasicAlgorithms::ArrayRef<Vec3i> triangles);

// return if triangles are edge-manifold and vtx-manifold
// return false if triangles contain invalid or degenerate triangles
bool areTrianglesManifold(BasicAlgorithms::ArrayRef<Vec3i> triangles);

// return vector of size #triangles, mapping: triID -> nbring triIDs sorted
// triangle connections are found by sharing edges
// work on non-manifold meshes
std::vector<std::vector<int>> getTriangleNeighborsByEdge(BasicAlgorithms::ArrayRef<Vec3i> triangles);

// return vector of size #triangles, mapping: triID -> nbring triIDs sorted
// triangle connections are found by sharing vertices
// work on non-manifold meshes
std::vector<std::vector<int>> getTriangleNeighborsByVertex(BasicAlgorithms::ArrayRef<Vec3i> triangles);

// return vector of size #vertices, mapping: vtxID -> nbring triIDs sorted
// work on non-manifold meshes
std::vector<std::vector<int>> getVertexTriangleNeighbors(BasicAlgorithms::ArrayRef<Vec3i> triangles, int numVertices);
std::map<int, std::vector<int>> getVertexTriangleNeighbors(BasicAlgorithms::ArrayRef<Vec3i> triangles);

// return the triangle indices in each connected component in triangles
// triIDs in each connected component is sorted
// the connectivity is defined by sharing (unordered) edges
// work on non-manifold meshes
std::vector<std::vector<int>> getConnectedComponentsByEdge(BasicAlgorithms::ArrayRef<Vec3i> triangles);

// return the triangle indices in each connected component in triangles
// triIDs in each connected component is sorted
// the connectivity is defined by sharing vertices
// work on non-manifold meshes
std::vector<std::vector<int>> getConnectedComponentsByVertex(BasicAlgorithms::ArrayRef<Vec3i> triangles);

// build oedgeTri: <v0,v1> -> tri with ordered edge <v0, v1>
// return false if triangles are invalid (degenerate or invalid vtx indices) or not edge-manifold
bool getOEdgeTriMap(BasicAlgorithms::ArrayRef<Vec3i> triangles, std::unordered_map<OEdgeKey, int> &oedgeTri);

// get non-manifold vertices on an edge-manifold mesh
// need input oedgeTri: <v0,v1> -> tri with ordered edge <v0, v1>
// return sorted non-manifold vertex IDs
std::vector<int> getNonManifoldVerticesOnEdgeManifoldTriangles(BasicAlgorithms::ArrayRef<Vec3i> triangle,
  const std::unordered_map<OEdgeKey, int> &oedgeTri);
// fix those non-manifold vertices by duplicating them to make the entire mesh manifold
// modify input triMesh to be manifold
// optionally output newVtx2OldVtxMap: newly generated vtxID -> original vtxID
void fixNonManifoldVerticesOnEdgeManifoldTriangles(TriMeshGeo &triMesh, const std::unordered_map<OEdgeKey, int> &oedgeTri,
  std::map<int, int> *newVtx2OldVtxMap = nullptr);

// build oedgeTris: <v0, v1> -> tris with ordered edge <v0, v1>
// work on non-manifold meshes
// return false only if triangles are invalid (degenerate or invalid vtx indices)
bool getNonManifoldOEdgeTrisMap(BasicAlgorithms::ArrayRef<Vec3i> triangles, std::unordered_map<OEdgeKey, std::vector<int>> &oedgeTris);

// On a closed manifold mesh, each edge is visited exactly twice: (i, j) and (j, i)
// Return those edges that are #(i,j) - #(j,i) != 0
// This includes boundary edges and non-manifold ones
// If (i, j) is a boundary edge, (i, j) will be returned
// If k = #(i,j) - #(j,i) > 0, then edge (i, j) will appear k times
// If k = #(i,j) - #(j,i) < 0, then edge (j, i) will appear -k times
// Useful when computing winding numbers
// interface inspired from libigl, which is subject to the terms of the Mozilla Public License, v. 2.0. You can obtain one copy at https://mozilla.org/MPL/2.0/.
std::vector<OEdgeKey> getExteriorEdges(BasicAlgorithms::ArrayRef<Vec3i> triangles);

// given a triMesh, find one vertex-connected component that is part of the outer layer of the mesh
// no manifoldness or closeness required
// return sorted triIDs for the outer component
// details: the method finds one extreme triangle (the highest triangle) and assumes this is part of the outer layer
// and then do a BFS to find the connected component with the triangle
std::vector<int> getOneOuterTriMeshConnectedComponentByVertex(TriMeshRef &triMesh);

}  // namespace Mesh
}  // namespace pgo
