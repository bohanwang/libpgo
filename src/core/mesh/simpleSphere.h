/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "mesh" library , Copyright (C) 2007 CMU, 2009 MIT, 2018 USC           *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code authors: Jernej Barbic, Christopher Twigg, Daniel Schroeder      *
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

/*
  A sphere
  Author: Jernej Barbic, CMU
*/

#include "meshLinearAlgebra.h"
#include "boundingBox.h"

namespace pgo
{
namespace Mesh
{

class SimpleSphere
{
public:
  SimpleSphere(const Vec3d &center, double radius):
    center_(center), radius_(radius) {}

  SimpleSphere(double x, double y, double z, double radius):
    center_(x, y, z), radius_(radius) {}

  const Vec3d &center() const { return center_; }
  double radius() const { return radius_; }

  // does the sphere intersect the bounding box
  bool doesBoundingBoxIntersect(const BoundingBox &box) const;

  // return true if queryPoint is inside or touches the sphere
  bool inside(const Vec3d &queryPoint) const;

  double signedDistance(const Vec3d &queryPoint) const;

private:
  Vec3d center_;
  double radius_;
};
}  // namespace Mesh
}  // namespace pgo
