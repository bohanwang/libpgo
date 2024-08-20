/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "volumetricMesh" library , Copyright (C) 2007 CMU, 2009 MIT, 2018 USC *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code authors: Jernej Barbic, Yijing Li                                *
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
  Creates the surface mesh of the given volumetric mesh.
  Note: interior volumetric mesh vertices are kept in the surface mesh (as isolated vertices).
        So, the vertex set of the volumetric mesh is identical to the surface mesh vertex set,
        with the same order.
*/

#pragma once

#include "volumetricMesh.h"
#include "EigenDef.h"

#include <vector>

namespace pgo
{
namespace VolumetricMeshes
{
namespace GenerateSurfaceMesh
{

// The output surface mesh is a triangle mesh for tet meshes and can also be a quad mesh for cubic meshes.
// All vertices of the volumetric mesh are included in the output surface mesh, even if they don't touch any surface triangle
//   (e.g., vertices in the interior of the volumetric mesh).
// triangulate: specifies whether output mesh should be quads or triangles (in the case of cubic volumetric meshes).
// allElementFaces = true: build ALL faces for ALL elements in the mesh
//                 = false: build only surface faces of the mesh (default)
void computeMesh(const VolumetricMesh *volumetricMesh, std::vector<EigenSupport::V3d> &vertices, std::vector<std::vector<int>> &faces,
  bool triangulate = false, bool allElementFaces = false);
};  // namespace GenerateSurfaceMesh
}  // namespace VolumetricMeshes
}  // namespace pgo