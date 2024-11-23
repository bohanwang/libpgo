/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "mesh" library , Copyright (C) 2007 CMU, 2009 MIT, 2018 USC           *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code authors: Somya Sharma, Jernej Barbic                             *
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

/*
  Tests if a 3D triangle overlaps with a 3D box.

  INPUT: center of box, box half sizes (in each of the three dimensions), and the three triangle vertices v0, v1, v2.
  OUTPUT: whether triangle intersects the box or not
  Note: all entries in boxHalfSize must be >= 0.
  Note: lower-left-front corner is boxCenter - boxHalfSize, upper-right-back corner is boxCenter + boxHalfSiz
e.
*/

#pragma once

namespace pgo
{
namespace Mesh
{
bool triBoxOverlap(const double boxcenter[3], const double boxhalfsize[3], const double trivert0[3], const double trivert1[3], const double trivert2[3]);
}
}  // namespace pgo
