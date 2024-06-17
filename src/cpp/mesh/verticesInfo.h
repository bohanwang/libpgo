/*
  This code is based on code from the Geometric Tools library,
  which is licensed under a boost license.
  Such usage is permitted by the boost license; for details,
  please see the boost license below.
*/

// Geometric Tools, LLC
// Copyright (c) 1998-2014
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt

/*************************************************************************
 *                                                                       *
 * We release our improvements to the wildMagic code under our standard  *
 * Vega FEM license, as follows:                                         *
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "improvements to the wildMagic library" , Copyright (C) 2018 USC      *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code author: Yijing Li                                                *
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
// handles floating-point queries among vertices

#include <vector>

namespace pgo
{
namespace Mesh
{

struct VerticesInformation
{
  // The intrinsic dimension of the input set. The parameter 'epsilon'
  // to the VerticesQuery::getInformation function is used to provide a tolerance when
  // determining the dimension.
  // if dimension == 0, all points are effectively the same (or very close based on epsilon)
  // If dimension == 1, all points effectively lie on a line segment.
  // If dimension == 2, all points effectively line on a plane.
  // If dimension == 3, the points are not coplanar.
  int dimension;

  // Axis-aligned bounding box of the input set.
  double min[3], max[3];
  // mMaxRange = max(max[0]-min[0], max[1]-min[1], and max[2]-min[2].)
  double maxRange;

  // The indices that define the maximum dimensional extents.  The
  // values mExtreme[0] and mExtreme[1] are the indices for the points
  // that define the largest extent in one of the coordinate axis
  // directions.
  // If the dimension is 2, then mExtreme[2] is the index
  // for the point that generates the largest extent in the direction
  // perpendicular to the line through the points corresponding to
  // mExtreme[0] and mExtreme[1].
  // Furthermore, if the dimension is 3, then mExtreme[3] is the index
  // for the point that generates the largest extent in the direction
  // perpendicular to the triangle defined by the other extreme points.
  int extreme[4];
  // If dimenstion == 3, the tetrahedron formed by the
  // points V[extreme0], V[extreme1], V[extreme2], V[extreme3]> is
  // counterclockwise (positive) if extremeCCW == true.
  // tet is positive if its vertices satisfy:  ((v1 - v0) x (v2 - v0)) dot (v3 - v0) >= 0
  bool extremeCCW;

  // Coordinate system.
  // The origin == vertices[exterme[0].  The
  // unit-length direction vector is valid only for 0 <= i < d.
  Vec3d origin;
  Vec3d direction[3];
};

VerticesInformation getVerticesInformation(const std::vector<Vec3d> &vertices, double epsilon);
}  // namespace Mesh
}  // namespace pgo
