/*************************************************************************
 *                                                                       *
 *                                                                       *
 * "mesh" library , Copyright (C) 2018 USC                               *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code authors: Yijing Li, Jernej Barbic                                *
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

#include "pointInsideOutsideQuery.h"
#include "triMeshNeighbor.h"

#include <iostream>

namespace pgo::Mesh
{
PointInsideOutsideQuery::PointInsideOutsideQuery(TriMeshGeo triMesh):
  triMesh(triMesh)
{
  triBVTree.buildByInertiaPartition(triMesh);

  std::vector<std::pair<int, int>> triIDPairs;
  triBVTree.selfIntersectionExact(triMesh, triIDPairs);

  if (triIDPairs.size() == 0)  // no self-intersection
  {
    // cout << "Input surface no self-intersection" << endl;
    if (areTrianglesManifold(triMesh.triangles()) && getExteriorEdges(triMesh.triangles()).size() == 0) {
      // cout << "Input surface is manifold and without boundary" << endl;
      // cout << "Use pseuodo-normal to solve inside/outside problems" << endl;
      triNormals.buildPseudoNormals(triMesh);

      triCols.reserve(triMesh.numTriangles());
      for (int triID = 0; triID < triMesh.numTriangles(); triID++)
        triCols.emplace_back(triMesh.pos(triID, 0), triMesh.pos(triID, 1), triMesh.pos(triID, 2));
    }
  }

  if (triNormals.numVertices() == 0)  // no pseudo normals constructed, use winding number for inside/outside tests
  {
    wnTree.build(triMesh);
  }
}

bool PointInsideOutsideQuery::isPointInside(const Vec3d &pos) const
{
  if (triNormals.numVertices() > 0)  // use pseudo normals
  {
    std::vector<std::tuple<double, double, int>> nodeStack;
    auto result = triBVTree.closestTriangleQuery(triCols.data(), pos, nodeStack);
    Vec3d normal = triNormals.getPseudoNormal(triMesh.triangles().data(), result.triID, result.feature);

    return ((pos - result.closestPosition).dot(normal) < 0);  // inside
  }

  // use winding number
  return (wnTree.windingNumber(triMesh, pos) > 0.5);  // inside
}
}  // namespace pgo::Mesh
