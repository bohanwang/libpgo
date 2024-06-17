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

#include "triMeshPseudoNormal.h"
#include "triMeshNeighbor.h"

#include "basicAlgorithms.h"
#include "containerHelper.h"
#include "pgoLogging.h"

#include <tbb/parallel_for.h>

#include <limits>

namespace pgo
{
namespace Mesh
{
TriMeshPseudoNormal::TriMeshPseudoNormal(TriMeshRef triMesh, const Vec3d *triangleNormals)
{
  int ret = buildPseudoNormals(triMesh, triangleNormals);
  if (ret != 0)
    throw ret;
}

TriMeshPseudoNormal::TriMeshPseudoNormal(TriMeshRef triMesh, const TriMeshNeighbor *externalNbr, const Vec3d *extTriNormals)
{
  buildPseudoNormals(triMesh, externalNbr, extTriNormals);
}

void TriMeshPseudoNormal::updateVertexPositions(TriMeshRef triMesh, const Vec3d *extTriNormals)
{
  buildPseudoNormals(triMesh, extTriNormals);
}

void TriMeshPseudoNormal::updateVertexPositions(TriMeshRef triMesh, const TriMeshNeighbor *externalNbr, const Vec3d *extTriNormals)
{
  buildPseudoNormals(triMesh, externalNbr, extTriNormals);
}

int TriMeshPseudoNormal::buildPseudoNormals(TriMeshRef triMesh, const TriMeshNeighbor *externalNbr, const Vec3d *extTriNormals)
{
  int ret = 0;
  vtxNormals.resize(triMesh.numVertices());
  edgeNormals.resize(triMesh.numVertices());
  triNormals.resize(triMesh.numTriangles());

  const TriMeshNeighbor *nbr = externalNbr;

  // first, compute triangle normals
  if (extTriNormals)
    memcpy(triNormals[0].data(), extTriNormals, sizeof(Vec3d) * triMesh.numTriangles());
  else {
    tbb::parallel_for(0, triMesh.numTriangles(), [&](int triID) {
      triNormals[triID] = triMesh.computeTriangleNormal(triID);
      if (triNormals[triID].hasNaN()) {
        // cout << "Error: triangle has nan pseudo-normal, ID: " << triID << endl;
        triNormals[triID].setZero();
        ret = 1;
      }
    });
  }

  // compute each vertex normal and edge normal

  tbb::parallel_for(0, triMesh.numVertices(), [&](int vtxID) {
    auto &vtxNormal = vtxNormals[vtxID];
    vtxNormal.setZero();
    edgeNormals[vtxID].clear();
    for (int triID : nbr->getVtxNearbyTriangles(vtxID)) {
      double angle = triMesh.getTriangleAngleAtVertexRobust(triID, vtxID);
      vtxNormal += angle * triNormals[triID];
      for (int vtxID2 : triMesh.tri(triID)) {
        if (vtxID >= vtxID2)  // we only consider edges where vtxID < vtxID2
          continue;

        // now compute edgeNormal for edge <vtxID, vtxID2>
        auto iter = edgeNormals[vtxID].find(vtxID2);
        if (iter == edgeNormals[vtxID].end()) {
          edgeNormals[vtxID][vtxID2] = triNormals[triID];
        }
        else {
          iter->second += triNormals[triID];
          iter->second.normalize();
          if (iter->second.hasNaN()) {
            // cout << "Error: edge pseudo normal has nan, edge vtx ID " << vtxID << " " << vtxID2;
            iter->second.setZero();
            ret = 1;
          }
        }
      }
    }  // end each nbring triangle
    if (vtxNormal.squaredNorm() > std::numeric_limits<double>::epsilon()) {
      vtxNormal.normalize();
      if (vtxNormal.hasNaN()) {
        // cout << "Error: vtx pseudo normal has nan after normalization, ID " << vtxID << endl;
        vtxNormal.setZero();
        ret = 1;
      }
    }
    else if (nbr->getVtxNearbyTriangles(vtxID).size() > 0)  // this vtx has nbr triangles but its computed vtx normal has zero length
    {
      //      cout << "Error: vtx pseudo normal has nan, ID " << vtxID << endl;
      ret = 1;
    }
    // else, this vtx has zero triangle as neighbors, then we don't consider its vertex normal
  });

  return ret;
}

int TriMeshPseudoNormal::buildPseudoNormals(TriMeshRef triMesh, const Vec3d *extTriNormals)
{
  int ret = 0;
  vtxNormals.assign(triMesh.numVertices(), asVec3d(0.0));
  edgeNormals.clear();
  edgeNormals.resize(triMesh.numVertices());
  triNormals.resize(triMesh.numTriangles());

  // first, compute triangle normals
  if (extTriNormals)
    memcpy(triNormals[0].data(), extTriNormals, sizeof(Vec3d) * triMesh.numTriangles());
  else {
    tbb::parallel_for(
      0, triMesh.numTriangles(), [&](int triID) {
        triNormals[triID] = triMesh.computeTriangleNormal(triID);
      },
      tbb::static_partitioner());
  }

  // compute each vertex normal and edge normal
  for (int triID = 0; triID < triMesh.numTriangles(); triID++) {
    if (triNormals[triID].hasNaN() || triNormals[triID] == asVec3d(0.0)) {
      triNormals[triID].setZero();
      //      cout << "Error: triangle has nan pseudo-normal, ID: " << triID << endl;
      ret = 1;  // triangle normal is broken, change return value
      continue;
    }
    for (int i = 0; i < 3; i++) {
      int vtxID = triMesh.triVtxID(triID, i);
      double angle = triMesh.getTriangleAngleAtVertexRobust(triID, vtxID);
      vtxNormals[vtxID] += angle * triNormals[triID];
      int vtxID2 = triMesh.triVtxID(triID, (i + 1) % 3);
      UEdgeKey edgeKey(vtxID, vtxID2);

      BasicAlgorithms::assignOrAddToMap(edgeNormals[edgeKey[0]], edgeKey[1], triNormals[triID]);
    }
  }

  // normalize each vertex normal and edge normal
  tbb::parallel_for(0, triMesh.numVertices(), [&](int vtxID) {
    if (vtxNormals[vtxID].squaredNorm() > 0) {
      vtxNormals[vtxID].normalize();
      if (vtxNormals[vtxID].hasNaN()) {
        //        cout << "Error: vtx pseudo normal has nan after normalization, ID " << vtxID << endl;
        vtxNormals[vtxID].setZero();
        ret = 1;
      }
    }
    // else, this vtx has zero triangle as neighbors, then we don't consider its vertex normal

    for (auto &p : edgeNormals[vtxID]) {
      Vec3d &edgeNormal = p.second;
      edgeNormal.normalize();
      if (edgeNormal.hasNaN()) {
        edgeNormal.setZero();
        //        cout << "Error: edge pseudo normal has nan, edge vtx ID " << vtxID << " " << p.first;
        ret = 1;
      }
    }
  });

  return ret;
}

const Vec3d &TriMeshPseudoNormal::edgeNormal(int vtxID0, int vtxID1) const
{
  if (vtxID0 > vtxID1) {
    std::swap(vtxID0, vtxID1);
  }
  const auto &edgeMap = edgeNormals[vtxID0];
  auto iter = edgeMap.find(vtxID1);
  PGO_ALOG(iter != edgeMap.end());
  //  if (iter == edgeMap.end()) { return Vec3d(0.0); }
  return iter->second;
}

}  // namespace Mesh
}  // namespace pgo
