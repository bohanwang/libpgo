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

namespace pgo
{
namespace Mesh
{

class HalfSpace
{
public:
  // HalfSpace satisfying: dot(normal, X) + d > 0, normal pointing outward
  // normalization of normal is not needed
  inline HalfSpace(const Vec3d &normal, double d);
  // HalfSpace on the positive side of triangle(a,b,c)
  // normal = cross(b-a, c-a), d = - dot(a, normal)
  inline HalfSpace(const Vec3d &a, const Vec3d &b, const Vec3d &c);

  // check a point is in this half space or not
  // return +1 for outside, 0 for on the boundary, -1 for inside
  inline int outside(const Vec3d &p) const;

  // return true if bb intersect or touch the half space
  inline bool intersect(const BoundingBox &bb) const;

protected:
  Vec3d normal;
  double d;
};

inline HalfSpace::HalfSpace(const Vec3d &n, double D):
  normal(n), d(D)
{
}

inline HalfSpace::HalfSpace(const Vec3d &a, const Vec3d &b, const Vec3d &c)
{
  normal = (b - a).cross(c - a);
  d = -a.dot(normal);
}

inline int HalfSpace::outside(const Vec3d &p) const
{
  double ret = normal.dot(p) + d;
  return ret > 0 ? +1 : (ret < 0 ? -1 : 0);
}

inline bool HalfSpace::intersect(const BoundingBox &bb) const
{
  const Vec3d *bound[2] = { &bb.bmin(), &bb.bmax() };
  for (int i = 0; i < 8; i++) {
    Vec3d v((*bound[(i & 4) >> 2])[0], (*bound[(i & 2) >> 1])[1], (*bound[i & 1])[3]);
    if (outside(v) <= 0)
      return true;  // if v is inside the half space
  }
  return false;
}

}  // namespace Mesh
}  // namespace pgo
