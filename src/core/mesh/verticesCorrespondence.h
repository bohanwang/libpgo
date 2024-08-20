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
#include "meshLinearAlgebra.h"
#include "arrayRef.h"

namespace pgo
{
namespace Mesh
{
// given two groups of vertices, one of them is the subset of the other
// this function finds the correspondence between them based on their positions
// return the sub2full vtx ID map, which should be pre-allocated as a buffer of size subGroup.size()
// function also returns the error on the correspondences:
// pair<int, double> where int is the subGroup ID that has the largest distance towards corresponding fullGroup vtx
// and double is this distance
std::pair<int, double> findVertexCorrespondences(BasicAlgorithms::ArrayRef<Vec3d> fullGroup, BasicAlgorithms::ArrayRef<Vec3d> subGroup, int *sub2fullMap);

}  // namespace Mesh
}  // namespace pgo
