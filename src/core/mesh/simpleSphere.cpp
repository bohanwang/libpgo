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

#include "simpleSphere.h"

/*
  A sphere
  Author: Jernej Barbic, CMU
*/

namespace pgo::Mesh
{
// does the sphere intersect the bounding box
bool SimpleSphere::doesBoundingBoxIntersect(const BoundingBox &box) const
{
  Vec3d bmin, bmax;
  bmin = box.bmin();
  bmax = box.bmax();

  double d;

#define COORDINATE_TEST(i)  \
  d = bmin[i] - center_[i]; \
  if (d > 0)                \
    dmin += d * d;          \
  d = center_[i] - bmax[i]; \
  if (d > 0)                \
    dmin += d * d;

  double dmin = 0;
  COORDINATE_TEST(0)
  COORDINATE_TEST(1)
  COORDINATE_TEST(2)

#undef COORDINATE_TEST

  return (dmin <= radius_ * radius_);
}

bool SimpleSphere::inside(const Vec3d &queryPoint) const
{
  return (center_ - queryPoint).squaredNorm() <= radius_ * radius_;
}

double SimpleSphere::signedDistance(const Vec3d &queryPoint) const
{
  return (center_ - queryPoint).norm() - radius_;
}

}  // namespace pgo::Mesh
