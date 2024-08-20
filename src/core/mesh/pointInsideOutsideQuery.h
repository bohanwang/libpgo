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

#pragma once

#include "triMeshGeo.h"
#include "boundingVolumeTree.h"
#include "windingNumberTree.h"
#include "triMeshPseudoNormal.h"
#include "triangle.h"

#include <vector>

namespace pgo
{
namespace Mesh
{
// an convenient wrapper library for fast and robust inside/outside query
// initPredicate() should be called before this class is to be used
// the input mesh must have correct orientation
// if the input mesh is manifold, self-intersection-free and without boundaries, use pseudo normal for inside/outside tests
// otherwise, use winding numbers to do that

class PointInsideOutsideQuery
{
public:
  PointInsideOutsideQuery(TriMeshGeo triMesh);

  // if the point is very close or on the surface, we consider it to be outside
  bool isPointInside(const Vec3d &pos) const;

  // the class object will always build a BVTree for triMesh
  // and it can be accessed through this function:
  const TriMeshBVTree &getBVTree() const { return triBVTree; }

protected:
  TriMeshGeo triMesh;
  TriMeshBVTree triBVTree;
  WindingNumberTree wnTree;
  TriMeshPseudoNormal triNormals;
  std::vector<TriangleWithCollisionInfo> triCols;
};
}  // namespace Mesh
}  // namespace pgo