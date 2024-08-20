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

#include "tetMeshGeo.h"

#include <functional>
#include <vector>

namespace pgo
{
namespace Mesh
{
// use flood-fill to remove tets which considered outside when embedding a closed surface.
// isTetBoundary returns true if tetID embeds the surface. It will get called at most once per tet
// isTetOuter returns true if tetID is outside the surface. One way to detect this fast is to use winding number. It will get called at most once per tet
// return a vector of size numTets, each bool represents whether tetID is outer
// note: isTetBoundary and isTetOuter should be thread-safe
std::vector<bool> labelOuterTets(const TetMeshRef &tetMesh, const TetNeighbor &tetNeighbor,
  std::function<bool(int tetID)> isTetBoundary, std::function<bool(int tetID)> isTetOuter);
}  // namespace Mesh
}  // namespace pgo
