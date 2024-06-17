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

#include "triMeshNeighbor.h"

#include "basicAlgorithms.h"
#include "graphSearchAlgorithms.h"
#include "containerHelper.h"
#include "valueIndex.h"
#include "pgoLogging.h"

#include <unordered_map>
#include <set>
#include <iostream>

namespace pgo
{
namespace Mesh
{

using namespace BasicAlgorithms;

std::vector<std::vector<int>> findBoundaryLoops(ArrayRef<Vec3i> triangles, ArrayRef<Vec3i> triNbrs)
{
  PGO_ALOG(triNbrs.size() == triangles.size());
  std::vector<std::vector<int>> ret;

  std::vector<bool> triVisited(triangles.size() * 3, false);
  for (int triID = 0; triID < triNbrs.size(); triID++) {
    const auto &n = triNbrs[triID];
    const auto &t = triangles[triID];
    //    cout << "visit triID " << triID << endl;
    for (int i = 0; i < 3; i++) {
      if (triVisited[triID * 3 + i])
        continue;
      triVisited[triID * 3 + i] = true;
      if (n[i] >= 0)
        continue;
      // find a boundary, let's go through the loop and find other edges in this loop
      int t1 = t[(i + 1) % 3];
      std::vector<int> loop = { t[i], t1 };

      int start = t[i];
      int end = t1;
      int curTriID = triID;
      //      cout << "start a loop at " << t[i] << " " << t1 << endl;

      int nextI = (i + 1) % 3;
      while (true) {
        // end = tri[curTriID][nextI+1]
        // check whether the edge <tri[curTriID][nextI] = end, tri[curTriID][nextI+1]> is boundary
        triVisited[curTriID * 3 + nextI] = true;
        if (triNbrs[curTriID][nextI] < 0)  // the next edge on cur triID is a boundary
        {
          end = triangles[curTriID][(nextI + 1) % 3];
          if (end == start) {
            break;
          }  // end this loop
          loop.push_back(end);
          nextI = (nextI + 1) % 3;
          //          cout << "find next vtx on triangle: " << end << endl;
          continue;
        }

        curTriID = triNbrs[curTriID][nextI];
        //        PGO_ALOG(curTriID != triID || triVisited[curTriID] == false);
        // next cutTriID has the edge <tri[old curTriID][nextI+1], tri[old curTriID][nextI] = end>
        nextI = getInvertedIndex(triangles[curTriID], end);
        PGO_ALOG(nextI >= 0);
        //        cout << "move to another triangle triID " <<  curTriID << endl;
      }

      ret.emplace_back(std::move(loop));
      break;  // because one triangle can only be in one loop
    }         // end for i
  }

  return ret;
}

std::vector<std::pair<int, OEdgeKey>> findBoundaryTriangles(ArrayRef<Vec3i> triangles, ArrayRef<Vec3i> triNbrs)
{
  PGO_ALOG(triNbrs.size() == triangles.size());
  std::vector<std::pair<int, OEdgeKey>> ret;
  for (int triID = 0; triID < triNbrs.size(); triID++) {
    const auto &n = triNbrs[triID];
    const auto &t = triangles[triID];
    for (int i = 0; i < 3; i++) {
      if (n[i] < 0)  // find a boundary
      {
        ret.push_back(std::make_pair(triID, OEdgeKey(t[i], t[(i + 1) % 3])));
      }
    }
  }
  return ret;
}

TriangleNeighbor::TriangleNeighbor(ArrayRef<Vec3i> triangles):
  numTriangles(triangles.size()),
  triNbrs(triangles.size(), Vec3i(-1,-1,-1))
{
  if (getOEdgeTriMap(triangles, oedgeTri) == false)
    throw 1;

  for (int triID = 0; triID < numTriangles; triID++) {
    for (int j = 0; j < 3; j++) {
      if (triNbrs[triID][j] < 0) {
        int v0 = triangles[triID][j];
        int v1 = triangles[triID][(j + 1) % 3];
        auto iter = oedgeTri.find(OEdgeKey(v1, v0));  // get the revserse key
        if (iter != oedgeTri.end()) {
          int triID2 = iter->second;
          triNbrs[triID][j] = triID2;
          int j2 = getInvertedIndex(triangles[triID2], v1);
          PGO_ALOG(j2 >= 0);
          triNbrs[triID2][j2] = triID;
        }
      }
    }
  }
}

int TriangleNeighbor::getTriangleAtEdge(const OEdgeKey &edge) const
{
  auto it = oedgeTri.find(edge);
  if (it == oedgeTri.end())
    return -1;
  return it->second;
}

int TriangleNeighbor::getVertexOppositeToEdge(const OEdgeKey &edge, ArrayRef<Vec3i> triangles) const
{
  PGO_ALOG(sizei(triNbrs) == triangles.size());
  auto it = oedgeTri.find(edge);
  if (it == oedgeTri.end())
    return -1;
  int triID = it->second;
  const Vec3i &tri = triangles[triID];
  for (int i = 0; i < 3; i++)
    if (tri[i] == edge[0]) {
      PGO_ALOG(tri[(i + 1) % 3] == edge[1]);
      return tri[(i + 2) % 3];
    }
  PGO_ALOG(0);  // error, the triangle does not match the values in TriangleNeighbor
  return -1;
}

std::vector<std::pair<int, OEdgeKey>> TriangleNeighbor::findBoundaryTriangles(ArrayRef<Vec3i> triangles) const
{
  return Mesh::findBoundaryTriangles(triangles, triNbrs);
}

std::vector<std::vector<int>> TriangleNeighbor::findBoundaryLoops(ArrayRef<Vec3i> triangles) const
{
  return Mesh::findBoundaryLoops(triangles, triNbrs);
}

std::vector<int> TriangleNeighbor::findTrianglesArroundBoundaryVertex(int prevVtxID, int vtxID, ArrayRef<Vec3i> triangles) const
{
  PGO_ALOG(sizei(triNbrs) == triangles.size());
  std::vector<int> ret;
  OEdgeKey prevEdge(prevVtxID, vtxID);
  auto it = oedgeTri.find(prevEdge);
  if (it == oedgeTri.end())
    return ret;
  int triID = it->second;
  do {
    ret.push_back(triID);
    int i = getInvertedIndex(triangles[triID], vtxID);
    PGO_ALOG(i >= 0);
    triID = triNbrs[triID][i];
  } while (triID >= 0);
  return ret;
}

// ==============================================================================
//                              TriMeshNeighbor
// ==============================================================================

TriMeshNeighbor::TriMeshNeighbor(TriMeshRef triMesh):
  numVertices(triMesh.numVertices()), numTriangles(triMesh.numTriangles()),
  vtxNbrTri(triMesh.numVertices()), triNbrs(triMesh.numTriangles(), Vec3i(-1, -1, -1)), vtxEdgeTri(triMesh.numVertices())
{
  for (int triID = 0; triID < numTriangles; triID++) {
    for (int j = 0; j < 3; j++) {
      int v0 = triMesh.tri(triID)[j];
      int v1 = triMesh.tri(triID)[(j + 1) % 3];
      if (v0 == v1)  // error, TriMesh is not valid!
      {
        //        cout << "Error, triangle is broken: " << triMesh.tri(triID) << endl;
        throw 1;
      }
      vtxNbrTri[v0].push_back(triID);
      if (vtxEdgeTri[v0].find(v1) != vtxEdgeTri[v0].end()) {
        //        cout << "Error: triMesh is not manifold" << endl;
        //        triMesh.save("debug.obj");
        throw 1;  // error, TriMesh is not manifold!
      }
      vtxEdgeTri[v0][v1] = triID;
    }
  }

  for (int triID = 0; triID < numTriangles; triID++) {
    for (int j = 0; j < 3; j++) {
      if (triNbrs[triID][j] < 0) {
        int v0 = triMesh.tri(triID)[j];
        int v1 = triMesh.tri(triID)[(j + 1) % 3];
        auto iter = vtxEdgeTri[v1].find(v0);
        if (iter != vtxEdgeTri[v1].end()) {
          int triID2 = iter->second;
          triNbrs[triID][j] = triID2;
          int j2 = getInvertedIndex(triMesh.tri(triID2), v1);
          PGO_ALOG(j2 >= 0);
          triNbrs[triID2][j2] = triID;
        }
      }
    }
  }
}

std::vector<std::pair<int, OEdgeKey>> TriMeshNeighbor::findBoundaryTriangles(ArrayRef<Vec3i> triangles) const
{
  return Mesh::findBoundaryTriangles(triangles, triNbrs);
}

std::vector<std::vector<int>> TriMeshNeighbor::findBoundaryLoops(ArrayRef<Vec3i> triangles) const
{
  return Mesh::findBoundaryLoops(triangles, triNbrs);
}

std::vector<int> TriMeshNeighbor::getVtxNearbyVertices(int vtxID, TriMeshRef mesh) const
{
  PGO_ALOG(numVertices == mesh.numVertices() && numTriangles == mesh.numTriangles());
  std::vector<int> ret;
  for (int triID : vtxNbrTri[vtxID]) {
    const auto &t = mesh.tri(triID);
    ret.insert(ret.end(), t.begin(), t.end());
  }
  BasicAlgorithms::sortAndDeduplicate(ret);
  return ret;
}

bool TriMeshNeighbor::areVerticesNeighbors(int vtxID0, int vtxID1) const
{
  return mapFind(vtxEdgeTri[vtxID0], vtxID1) || mapFind(vtxEdgeTri[vtxID1], vtxID0);
}

int TriMeshNeighbor::getTriangleAtEdge(const OEdgeKey &edge) const
{
  auto it = vtxEdgeTri[edge[0]].find(edge[1]);
  if (it == vtxEdgeTri[edge[0]].end())
    return -1;
  return it->second;
}

int TriMeshNeighbor::getVertexOppositeToEdge(const OEdgeKey &edge, ArrayRef<Vec3i> triangles) const
{
  PGO_ALOG(sizei(triNbrs) == triangles.size());
  auto it = vtxEdgeTri[edge[0]].find(edge[1]);
  if (it == vtxEdgeTri[edge[0]].end())
    return -1;
  int triID = it->second;
  const Vec3i &tri = triangles[triID];
  for (int i = 0; i < 3; i++)
    if (tri[i] == edge[0]) {
      PGO_ALOG(tri[(i + 1) % 3] == edge[1]);
      return tri[(i + 2) % 3];
    }
  PGO_ALOG(0);  // error, the triangle does not match the values in TriangleNeighbor
  return -1;
}

// ==============================================================================
//                          NonManifoldTriangleNeighbor
// ==============================================================================

NonManifoldTriangleNeighbor::NonManifoldTriangleNeighbor(const std::vector<Vec3i> &triangles):
  numTriangles(static_cast<int>(triangles.size()))
{
  if (getNonManifoldOEdgeTrisMap(triangles, oedgeTris) == false)
    throw 1;
}

std::vector<int> NonManifoldTriangleNeighbor::getTriangleAtOEdge(const OEdgeKey &oedge) const
{
  auto it = oedgeTris.find(oedge);
  if (it == oedgeTris.end())
    return {};
  return it->second;
}

std::vector<OEdgeKey> NonManifoldTriangleNeighbor::findNonManifoldOEdges() const
{
  std::vector<OEdgeKey> ret;
  for (const auto &p : oedgeTris) {
    if (p.second.size() > 1)
      ret.push_back(p.first);
  }
  return ret;
}
/////////////////////////////////////////////////////////////////////////
//                       Other Functions
/////////////////////////////////////////////////////////////////////////

using Neighbor = std::function<const std::vector<int> &(int nodeA)>;
std::vector<std::vector<int>> getConnectedComponents(int numNodes, Neighbor getNeighbor)
{
  std::vector<std::vector<int>> ret;
  if (numNodes == 0)
    return ret;

  std::vector<bool> visited(numNodes, false);
  std::vector<int> component;  // record traversed triangles

  for (int seedID = 0; seedID < numNodes; seedID++) {
    if (visited[seedID])
      continue;
    component.clear();
    visited[seedID] = true;
    component.push_back(seedID);
    size_t candidateBegin = 0, candidateEnd = 1;

    while (candidateBegin != candidateEnd) {
      for (size_t i = candidateBegin; i < candidateEnd; i++) {
        int node = component[i];
        const auto &nbrs = getNeighbor(node);
        for (int nbr : nbrs) {
          if (nbr < 0 || visited[nbr])
            continue;
          visited[nbr] = true;
          component.push_back(nbr);
        }
      }
      candidateBegin = candidateEnd;
      candidateEnd = component.size();
    }

    sort(component.begin(), component.end());
    PGO_ALOG(unique(component.begin(), component.end()) == component.end());
    ret.emplace_back(std::move(component));
  }
  return ret;
}

// build oedgeTri: <v0,v1> -> tri with ordered edge <v0, v1>
bool getOEdgeTriMap(ArrayRef<Vec3i> triangles, std::unordered_map<OEdgeKey, int> &oedgeTri)
{
  int numTriangles = triangles.size();
  oedgeTri.clear();

  for (int triID = 0; triID < numTriangles; triID++) {
    for (int j = 0; j < 3; j++) {
      int v0 = triangles[triID][j];
      int v1 = triangles[triID][(j + 1) % 3];
      if (v0 == v1 || v0 < 0)  // error, TriMesh is not valid!
      {
        return false;
      }
      OEdgeKey edge(v0, v1);
      if (oedgeTri.find(edge) != oedgeTri.end()) {
        return false;  // this signed edge appears twice, not manifold!
      }
      oedgeTri.insert(std::make_pair(edge, triID));
    }
  }
  return true;
}

// go through a triangle fan around vtxID on an edge-manifold mesh
// input includes vtxID, a startingTriID to start the search,
// and the localVtxID == triangles[startingTriID].getInvertedIndex(vtxID), localVtxID: [0,3)
// once a new triangle is visited, its triID and the local vtxID [0,3) of the input vtxID is passed to
// processTriangle for custom processing
void visitTriangleFanNeighboringVtx(ArrayRef<Vec3i> triangles,
  const std::unordered_map<OEdgeKey, int> &oedgeTri, int vtxID, int startingTriID, int localVtxID,
  std::function<void(int newTriID, int newLocalVtxID)> processTriangle)
{
  // visit the neighboring triangles to search for the fan shape.
  // ((lvID + edgeVtxOffset)%3, (lvID + edgeVtxOffset+1)%3) represents an OEdge on triangle triID
  // return true if we go 360 degrees around the vtx, finishing one full loop.
  // if it returns false, then we hit a boundary edge on the triangle mesh
  auto visitNbr = [&](int edgeVtxOffset) {
    //        cout << "begin search" << endl;
    // loop over triangles around this vtx, starting at the edge (vEdgeStart, vEdgeEnd)
    int vEdgeStart = triangles[startingTriID][(localVtxID + edgeVtxOffset) % 3];
    int vEdgeEnd = triangles[startingTriID][(localVtxID + edgeVtxOffset + 1) % 3];
    // edge<vEdgeStart, vEdgeEnd> is an OEdge of triID
    // since we want to find the nbring tri on this edge, we should search the reversed edge in oedgeTri
    auto it = oedgeTri.find({ vEdgeEnd, vEdgeStart });
    while (it != oedgeTri.end()) {
      int nextTriID = it->second;
      if (nextTriID == startingTriID)
        return true;  // we finished one full loop around the vtx v
      int nextlvID = getInvertedIndex(triangles[nextTriID], vtxID);
      PGO_ALOG(nextlvID >= 0);

      processTriangle(nextTriID, nextlvID);

      vEdgeStart = triangles[nextTriID][(nextlvID + edgeVtxOffset) % 3];
      vEdgeEnd = triangles[nextTriID][(nextlvID + edgeVtxOffset + 1) % 3];
      it = oedgeTri.find({ vEdgeEnd, vEdgeStart });
    }
    return false;
  };

  if (visitNbr(0) == false)  // visitNbr(0) return true only if it finishes one complete loop
    visitNbr(2);
}

// The algorithm here is that for each vtx, we visit each neighboring triangles and record this "fan shape" we visited.
// If the mesh is vtx-manifold, then there will only be one fan shape for each vtx.
// But on a non-vtx-manifold but edge-manifold mesh, one vtx can have more than one fan shape.
std::vector<int> getNonManifoldVerticesOnEdgeManifoldTriangles(ArrayRef<Vec3i> triangles,
  const std::unordered_map<OEdgeKey, int> &oedgeTri)
{
  std::vector<int> ret;
  int numTriangles = triangles.size();
  std::set<int> vtxVisited;
  std::vector<bool> triVtxVisited(numTriangles * 3, false);

  for (int triID = 0; triID < numTriangles; triID++) {
    for (int lvID = 0; lvID < 3; lvID++)  // local vtx ID
    {
      if (triVtxVisited[triID * 3 + lvID])
        continue;
      int curVtxID = triangles[triID][lvID];
      if (vtxVisited.find(curVtxID) != vtxVisited.end())  // if this v has been visited
      {
        ret.push_back(curVtxID);                          // not vtx-manifold
                                                          //        cout << "non-manifold: " << triID << " " << v << endl;
        continue;
      }
      vtxVisited.insert(curVtxID);
      triVtxVisited[triID * 3 + lvID] = true;

      //      cout << "visit triID " << triID << " v " << v << endl;
      visitTriangleFanNeighboringVtx(triangles, oedgeTri, curVtxID, triID, lvID, [&](int nextTriID, int nextlvID) {
        PGO_ALOG(triVtxVisited[nextTriID * 3 + nextlvID] == false);
        triVtxVisited[nextTriID * 3 + nextlvID] = true;
      });
    }
  }
  sortAndDeduplicate(ret);
  return ret;
}

void fixNonManifoldVerticesOnEdgeManifoldTriangles(TriMeshGeo &triMesh, const std::unordered_map<OEdgeKey, int> &oedgeTri,
  std::map<int, int> *newVtx2OldVtxMap)
{
  std::vector<int> nmVtxIDs = getNonManifoldVerticesOnEdgeManifoldTriangles(triMesh.triangles(), oedgeTri);
  if (nmVtxIDs.size() == 0)
    return;

  std::map<int, std::vector<int>> nmVtxNbringTris;
  for (int triID = 0; triID < triMesh.numTriangles(); triID++)
    for (int i = 0; i < 3; i++) {
      int vtxID = triMesh.triVtxID(triID, i);
      if (binarySearchFound(nmVtxIDs, vtxID)) {
        nmVtxNbringTris[vtxID].push_back(triID);
      }
    }

  for (int nmVtxID : nmVtxIDs) {
    const std::vector<int> &nbringTriIDs = nmVtxNbringTris[nmVtxID];
    PGO_ALOG(nbringTriIDs.size() > 1);
    std::set<int> nbringTriIDVisited;
    std::vector<std::vector<int>> triFans;  // stores the triIDs in each triangle fan around nmVtxID
    for (int i = 0; i < sizei(nbringTriIDs); i++) {
      int triID = nbringTriIDs[i];
      if (setFind(nbringTriIDVisited, triID))
        continue;
      nbringTriIDVisited.insert(triID);

      std::vector<int> newFan = { triID };
      int lvID = getInvertedIndex(triMesh.tri(triID), nmVtxID);  // local vtx ID, [0,3)
      visitTriangleFanNeighboringVtx(triMesh.triangles(), oedgeTri, nmVtxID, triID, lvID, [&](int newTriID, int) {
        PGO_ALOG(binarySearchFound(nbringTriIDs, newTriID));
        PGO_ALOG(setNotFind(nbringTriIDVisited, newTriID));
        nbringTriIDVisited.insert(newTriID);
        newFan.push_back(newTriID);
      });
      //      cout << "Found a new fan, size " << newFan.size() << endl;
      triFans.emplace_back(std::move(newFan));
    }

    // sanity check
    size_t sum = 0;
    for (const auto &fan : triFans)
      sum += fan.size();
    //    cout << "Fans:  " << triFans.size() << endl;
    //    cout << "sum " << sum << " " << nbringTriIDs.size() << endl;
    PGO_ALOG(sum == nbringTriIDs.size());
    PGO_ALOG(triFans.size() > 1);

    for (int i = 1; i < sizei(triFans); i++)  // for each addtional fan
    {
      int newVtxID = triMesh.numVertices();
      if (newVtx2OldVtxMap)
        newVtx2OldVtxMap->emplace(newVtxID, nmVtxID);
      triMesh.addPos(triMesh.pos(nmVtxID));
      for (int triID : triFans[i]) {
        int lvID = getInvertedIndex(triMesh.tri(triID), nmVtxID);
        triMesh.tri(triID)[lvID] = newVtxID;
      }
    }
  }  // for each nmVtxID
}

bool areTrianglesEdgeManifold(ArrayRef<Vec3i> triangles)
{
  std::unordered_map<OEdgeKey, int> oedgeTri;  // <v0,v1> -> tri with ordered edge <v0, v1>
  return getOEdgeTriMap(triangles, oedgeTri);
}

bool areTrianglesManifold(ArrayRef<Vec3i> triangles)
{
  //  int numTriangles = triangles.size();
  std::unordered_map<OEdgeKey, int> oedgeTri;  // <v0,v1> -> tri with ordered edge <v0, v1>

  if (getOEdgeTriMap(triangles, oedgeTri) == false)
    return false;
  // now the mesh is at least edge-manifold

  return (getNonManifoldVerticesOnEdgeManifoldTriangles(triangles, oedgeTri).size() == 0);
}

std::vector<std::vector<int>> getTriangleNeighborsByEdge(ArrayRef<Vec3i> triangles)
{
  std::vector<std::vector<int>> nbrs(triangles.size());
  std::map<UEdgeKey, std::vector<int>> uedgeMap;
  for (int triID = 0; triID < triangles.size(); triID++) {
    const Vec3i &tri = triangles[triID];
    for (int j = 0; j < 3; j++) {
      UEdgeKey edge(tri[j], tri[(j + 1) % 3]);
      uedgeMap[edge].push_back(triID);
    }
  }
  for (const auto &p : uedgeMap) {
    const auto &tris = p.second;
    for (size_t i = 0; i < tris.size(); i++)
      for (size_t j = i + 1; j < tris.size(); j++) {
        nbrs[tris[i]].push_back(tris[j]);
        nbrs[tris[j]].push_back(tris[i]);
      }
  }
  for (auto &vec : nbrs)  // sort and remove duplicates in each nbr vector
  {
    sortAndDeduplicate(vec);
  }

  return nbrs;
}

std::vector<std::vector<int>> getTriangleNeighborsByVertex(ArrayRef<Vec3i> triangles)
{
  if (triangles.size() == 0)
    return {};

  std::map<int, std::vector<int>> vtxNbringTri;
  for (int triID = 0; triID < triangles.size(); triID++) {
    for (int j = 0; j < 3; j++)
      vtxNbringTri[triangles[triID][j]].push_back(triID);
  }
  std::vector<std::vector<int>> triNbrTri(triangles.size());
  for (const auto &p : vtxNbringTri) {
    const auto &t = p.second;  // neighborhood
    for (size_t i = 0; i < t.size(); i++) {
      for (size_t j = i + 1; j < t.size(); j++) {
        triNbrTri[t[i]].push_back(t[j]);
        triNbrTri[t[j]].push_back(t[i]);
      }
    }
  }
  for (auto &t : triNbrTri) {
    sortAndDeduplicate(t);
  }

  return triNbrTri;
}

std::vector<std::vector<int>> getVertexTriangleNeighbors(ArrayRef<Vec3i> triangles, int numVertices)
{
  std::vector<std::vector<int>> vtxNbrTri(numVertices);
  for (int triID = 0; triID < triangles.size(); triID++) {
    for (int j = 0; j < 3; j++) {
      int v0 = triangles[triID][j];
      PGO_ALOG(v0 >= 0 && v0 < numVertices);
      vtxNbrTri[v0].push_back(triID);
    }
  }
  for (int vtxID = 0; vtxID < numVertices; vtxID++)
    std::sort(vtxNbrTri[vtxID].begin(), vtxNbrTri[vtxID].end());

  return vtxNbrTri;
}

std::map<int, std::vector<int>> getVertexTriangleNeighbors(ArrayRef<Vec3i> triangles)
{
  std::map<int, std::vector<int>> vtxNbrTri;
  for (int triID = 0; triID < triangles.size(); triID++) {
    for (int j = 0; j < 3; j++) {
      int v0 = triangles[triID][j];
      PGO_ALOG(v0 >= 0);
      vtxNbrTri[v0].push_back(triID);
    }
  }
  for (auto &p : vtxNbrTri)
    std::sort(p.second.begin(), p.second.end());

  return vtxNbrTri;
}

std::vector<std::vector<int>> getConnectedComponentsByEdge(ArrayRef<Vec3i> triangles)
{
  if (triangles.size() == 0)
    return {};

  auto triNbrs = getTriangleNeighborsByEdge(triangles);

  auto getNeighbor = [&](int triID) -> const std::vector<int> & {
    return triNbrs[triID];
  };

  return getConnectedComponents(triangles.size(), getNeighbor);
}

std::vector<std::vector<int>> getConnectedComponentsByVertex(ArrayRef<Vec3i> triangles)
{
  if (triangles.size() == 0)
    return {};

  auto triNbrTri = getTriangleNeighborsByVertex(triangles);

  auto getNeighbor = [&](int triID) -> const std::vector<int> & {
    return triNbrTri[triID];
  };

  return getConnectedComponents(triangles.size(), getNeighbor);
}

bool getNonManifoldOEdgeTrisMap(ArrayRef<Vec3i> triangles, std::unordered_map<OEdgeKey, std::vector<int>> &oedgeTris)
{
  int numTriangles = triangles.size();
  oedgeTris.clear();

  for (int triID = 0; triID < numTriangles; triID++) {
    for (int j = 0; j < 3; j++) {
      int v0 = triangles[triID][j];
      int v1 = triangles[triID][(j + 1) % 3];
      if (v0 == v1 || v0 < 0)  // error, TriMesh is not valid!
      {
        return false;
      }
      OEdgeKey edge(v0, v1);
      oedgeTris[edge].push_back(triID);
    }
  }
  return true;
}

std::vector<OEdgeKey> getExteriorEdges(ArrayRef<Vec3i> triangles)
{
  std::map<OEdgeKey, int> count;
  for (int triID = 0; triID < triangles.size(); triID++) {
    for (int i = 0; i < 3; i++) {
      OEdgeKey e(triangles[triID][i], triangles[triID][(i + 1) % 3]);
      auto it = count.find(e);
      if (it != count.end()) {
        it->second++;
      }
      else {
        it = count.find(e.getReversedEdgeKey());
        if (it != count.end()) {
          it->second--;
        }
        else {
          count.emplace(e, 1);
        }
      }
    }
  }
  std::vector<OEdgeKey> ret;
  for (const auto &p : count) {
    if (p.second == 0)
      continue;
    if (p.second > 0) {
      for (int i = 0; i < p.second; i++)
        ret.push_back(p.first);
    }
    else {
      OEdgeKey e = p.first.getReversedEdgeKey();
      for (int i = 0; i < (-p.second); i++) {
        ret.push_back(e);
      }
    }
  }
  return ret;
}

std::vector<int> getOneOuterTriMeshConnectedComponentByVertex(TriMeshRef &triMesh)
{
  if (triMesh.numTriangles() == 0)
    return {};

  MaxValueIndex mvi;  // get the triangle with the highest vtx, which should be an outer surface triangle
  for (int triID = 0; triID < triMesh.numTriangles(); triID++)
    for (int i = 0; i < 3; i++)
      mvi.update(triMesh.pos(triMesh.triVtxID(triID, i))[1], triID);

  // now do BFS
  auto surfaceTriNbrTri = getTriangleNeighborsByVertex(triMesh.trianglesRef());
  std::vector<bool> surfaceTriangleVisited(triMesh.numTriangles(), false);
  auto outerTriIDs = findNodeIDsFromSeedBFS(convertArrayToFunction(surfaceTriNbrTri), surfaceTriangleVisited, mvi.index);
  sortAndDeduplicate(outerTriIDs);

  return outerTriIDs;
}

}  // namespace Mesh
}  // namespace pgo
