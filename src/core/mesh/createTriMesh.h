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

#include "triMeshGeo.h"

namespace pgo
{
namespace Mesh
{
// =========================================================
//              Create simple TriMeshes
// =========================================================

// create an AABB mesh
TriMeshGeo createBoxMesh(const Vec3d &bmin, const Vec3d &bmax);
TriMeshGeo createCubeMesh(const Vec3d &center, double size);

// the ordering of v0-v3 affects the normals of the output mesh
// if v0-v3 have positive tet orientation, the faces in the output mesh have outward normals
TriMeshGeo createTetSurfaceMesh(const Vec3d &v0, const Vec3d &v1, const Vec3d &v2, const Vec3d &v3);

// create a mesh consisting of only one triangle
TriMeshGeo createSingleTriangleMesh(const TriMeshRef mesh, int triID);
TriMeshGeo createSingleTriangleMesh(const Vec3d &v0, const Vec3d &v1, const Vec3d &v2);

// create a cylinder without cap, centered at origin and axis at y-direction
// subdivisionAxis: #cut on the circle
// subdivisionHeight: #subdivision on the height direction
TriMeshGeo createCylinderWallMesh(double radius, double height, int subdivisionAxis, int subdivisionHeight);
TriMeshGeo createCylinderWallMesh(double radius1, double radius2, double height, int subdivisionAxis, int subdivisionHeight);

// create a cylinder mesh centered at origin and axis at y-direction
TriMeshGeo createCylinderMesh(double radius, double height, int subdivisionAxis, int subdivisionHeight);

TriMeshGeo createTorus(int radialResolution, int tubularResolution, double radius, double thickness);

TriMeshGeo createSphereMesh(double radius, int axisSubdivisions, int heightSubdivisions);
}  // namespace Mesh
}  // namespace pgo
