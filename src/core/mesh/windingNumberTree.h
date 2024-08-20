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

// Winding number tree acceleration implementation from the paper:
// Alec Jacobson, Ladislav Kavan, Olga Sorkine-Hornung:
// Robust Inside-Outside Segmentation using Generalized Winding Numbers, SIGGRAPH 2013
// Implementation inspired by libigl.

#include "boundingBox.h"
#include "triMeshGeo.h"

#include <map>
#include <vector>
#include <list>
#include <limits>
#include <climits>

namespace pgo
{
namespace Mesh
{
class WindingNumberTree
{
public:
  WindingNumberTree() {}

  // build the tree, vertices buffer in mesh should not be changed or deleted when this WindingNumberTree is used
  void build(const TriMeshRef mesh, int maxDepth = INT_MAX, int maxNumTrianglesPerNode = 10);

  void clear();  // clear internal data

  // compute winding number recursively
  double windingNumber(const TriMeshRef mesh, const Vec3d &p) const;

  void print() const;  // print the tree

protected:
  int numMeshVertices = 0;
  int numMeshTriangles = 0;

  struct Node
  {
    int depth = 0;
    std::vector<Vec3i> triangles;  // triangles for this node is stored
    std::vector<Vec3i> bouFaces;   // boundary faces
    BoundingBox bb;
    int childrenIDs[2] = { -1, -1 };
    Node(int depth, const TriMeshRef mesh);
    const TriMeshRef meshInNode(const TriMeshRef fullMesh) const { return TriMeshRef(fullMesh.numVertices(), fullMesh.positions(), triangles); }
    const TriMeshRef boundaryInNode(const TriMeshRef fullMesh) const { return TriMeshRef(fullMesh.numVertices(), fullMesh.positions(), bouFaces); }
  };

  int depth = 0;
  std::vector<Node> nodes;
};
}
}  // namespace pgo

