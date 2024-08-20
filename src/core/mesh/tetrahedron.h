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

#include <array>

namespace pgo
{
namespace Mesh
{
class Tetrahedron
{
public:
  Tetrahedron(const Vec3d &v0, const Vec3d &v1, const Vec3d &v2, const Vec3d &v3):
    v{ { v0, v1, v2, v3 } } {}

  const Vec3d &operator[](int i) const { return v[i]; }
  Vec3d &operator[](int i) { return v[i]; }

  const Vec3d *begin() const { return &v[0]; }
  Vec3d *begin() { return &v[0]; }

  const Vec3d *end() const { return &v[0] + 4; }
  Vec3d *end() { return &v[0] + 4; }

protected:
  std::array<Vec3d, 4> v;
};

}  // namespace Mesh
}  // namespace pgo
