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

namespace pgo
{
namespace Mesh
{

// a plane implementation for fast geometry queries
// no normalization operation is used

struct FastPlane
{
  FastPlane(const Vec3d &v0, const Vec3d &v1, const Vec3d &v2);
  // return 1 if outside according to the normal of the plane
  //       -1 if inside
  //        0 if on the plane
  int outside(const Vec3d &v) const;

  // return the scaled distance between v and the plane
  // the distance is scaled by len(dir)
  double scaledDistance(const Vec3d &v) const;

  Vec3d dir;
  double d;
};

inline FastPlane::FastPlane(const Vec3d &v0, const Vec3d &v1, const Vec3d &v2)
{
  Vec3d e1 = v1 - v0, e2 = v2 - v0;
  dir = e1.cross(e2);
  d = v0.dot(dir);
}

inline int FastPlane::outside(const Vec3d &v) const
{
  auto s = (v.dot(dir) - d);
  return (s > 0 ? 1 : (s < 0) ? -1 :
                                0);
}

inline double FastPlane::scaledDistance(const Vec3d &v) const
{
  return std::abs(v.dot(dir) - d);
}

}  // namespace Mesh
}  // namespace pgo
