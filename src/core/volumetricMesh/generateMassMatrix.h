/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "volumetricMesh" library , Copyright (C) 2007 CMU, 2009 MIT, 2018 USC *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code author: Jernej Barbic                                            *
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
  Computes the mass matrix for the given volumetric mesh.
  See also volumetricMesh.h .
*/

#pragma once

#include "volumetricMesh.h"
#include "EigenDef.h"

namespace pgo
{
namespace VolumetricMeshes
{
namespace GenerateMassMatrix
{
// If inflate3Dim flag is on, each matrix element z will be augmented
// to a 3x3 z*I matrix (causing mtx dimensions to grow by a factor of 3).
// In order to get the true 3D mass matrix, set inflate3Dim to true (output
// matrix will be 3*numVertices x 3*numVertices).
// In order to save some space, set it to false (output matrix will be
// numVertices x numVertices).
void computeMassMatrix(const VolumetricMesh *volumetricMesh, EigenSupport::SpMatD &massMatrix, bool inflate3Dim = false, const double *elementWeight = nullptr);
// computes the mass belonging to each vertex, by lumping the mass matrix
// masses has size of volumetricMesh->getNumVertices() if (inflat3Dim == false) or 3 * n if (inflate3Dim == true)
void computeVertexMasses(const VolumetricMesh *volumetricMesh, double *masses, bool inflate3Dim = false);
// This function divides the mass of one element equally to assign to each vertex of the element,
// therefore it is faster than computeVertexMasses which computes all element mass matrices.
// At least for TetMesh and CubicMeshes, this function behaves the same as computeVertexMasses.
void computeVertexMassesByAveragingNeighboringElements(const VolumetricMesh *volumetricMesh, double *masses, bool inflate3Dim = false);
};  // namespace GenerateMassMatrix

}  // namespace VolumetricMeshes
}  // namespace pgo
