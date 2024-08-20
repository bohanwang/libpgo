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

#include "triMeshGeo.h"

#include <map>

namespace pgo
{
namespace Mesh
{

class TriMeshNeighbor;

// compute pseudo-normal on a triangle mesh
// assuming TriMesh is edge-manifold, required by TriMeshNeighbor
class TriMeshPseudoNormal
{
public:
  TriMeshPseudoNormal() {}  // empty
  // throw int 1 if triMesh is not valid or vtx/edge/triangle normals have NaN
  // if triMeshNbr != nullptr, use this TriMeshNeighbor to help build normals
  TriMeshPseudoNormal(TriMeshRef triMesh, const Vec3d *triangleNormals = nullptr);
  TriMeshPseudoNormal(TriMeshRef triMesh, const TriMeshNeighbor *triMeshNbr, const Vec3d *triangleNormals = nullptr);

  // initialize data with the pseudo-normals from triMesh
  // the implementation assumes the input mesh is manifold
  // return 0 if success, 1 if invalid, non-manifold or degenerate cases found
  int buildPseudoNormals(TriMeshRef triMesh, const Vec3d *triangleNormals = nullptr);
  int buildPseudoNormals(TriMeshRef triMesh, const TriMeshNeighbor *triMeshNbr, const Vec3d *triangleNormals = nullptr);

  // update the existing TriMeshPseudoNormal with new mesh positions
  void updateVertexPositions(TriMeshRef triMesh, const Vec3d *triangleNormals = nullptr);
  void updateVertexPositions(TriMeshRef triMesh, const TriMeshNeighbor *triMeshNbr, const Vec3d *triangleNormals = nullptr);

  int numVertices() const { return static_cast<int>(vtxNormals.size()); }
  int numTriangles() const { return static_cast<int>(triNormals.size()); }

  // return Vec3d(0.0) if this vtx has no nearby triangles
  const Vec3d &vtxNormal(int vtxID) const { return vtxNormals[vtxID]; }

  // assert (vtxID0, vtxID1) is a valid edge, the order of vtxID0 and vtxID1 does not matter
  const Vec3d &edgeNormal(int vtxID0, int vtxID1) const;

  const Vec3d &triNormal(int triID) const { return triNormals[triID]; }

  // closest feature:
  //  0: vertex0
  //  1: vertex1
  //  2: vertex2
  //  3: edge among 01
  //  4: edge among 12
  //  5: edge among 20
  //  6: the face itself

  // get pseudo-normal on the closest feature
  inline const Vec3d &getPseudoNormal(const Vec3i *triangles, int triID, int closestFeature) const;

  // get pseudo closet position on the feature
  // for vertices, returns the vertex itself
  // for edges, returns the midpoint of the edge (for the purposes of distance sign test, this could be any point on the edge)
  // for faces, returns the face centroid
  inline Vec3d getPseudoClosestPosition(TriMeshRef triMesh, int triID, int closestFeature) const;

protected:
  std::vector<Vec3d> vtxNormals;
  // edgeNormals: vtx0 -> vtx1 -> the peseudo normal of the edge <vtx0, vtx1>
  // note: vtx0 < vtx1
  std::vector<std::map<int, Vec3d>> edgeNormals;
  std::vector<Vec3d> triNormals;
};

inline const Vec3d &TriMeshPseudoNormal::getPseudoNormal(const Vec3i *triangles, int triID, int closestFeature) const
{
  switch (closestFeature) {
  case 0:
  case 1:
  case 2:
    return vtxNormals[triangles[triID][closestFeature]];

  case 3:
    return edgeNormal(triangles[triID][0], triangles[triID][1]);

  case 4:
    return edgeNormal(triangles[triID][1], triangles[triID][2]);

  case 5:
    return edgeNormal(triangles[triID][2], triangles[triID][0]);

  case 6:
  default:
    return triNormals[triID];
  }
}

inline Vec3d TriMeshPseudoNormal::getPseudoClosestPosition(TriMeshRef triMesh, int triID, int closestFeature) const
{
  switch (closestFeature) {
  case 0:
  case 1:
  case 2:
    return triMesh.pos(triID, closestFeature);

  case 3:
    return 0.5 * (triMesh.pos(triID, 0) + triMesh.pos(triID, 1));

  case 4:
    return 0.5 * (triMesh.pos(triID, 1) + triMesh.pos(triID, 2));

  case 5:
    return 0.5 * (triMesh.pos(triID, 2) + triMesh.pos(triID, 0));

  case 6:
  default:
    return triMesh.computeTriangleCentroid(triID);
  }
}
}  // namespace Mesh
}  // namespace pgo
