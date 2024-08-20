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

// Predicates for robust geometric queries.
// Interface to Shewchuk's exact predicates code.

#include "initPredicates.h"

namespace pgo
{
namespace Mesh
{

// !!!!!!!!!! IMPORTANT !!!!!!!!!!
// You must first call "initPredicates" defined in initPredicates.h to initialize variables used for exact predicates.
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!

// Return a positive value if the points pa, pb, and pc occur
// in counterclockwise order; a negative value if they occur
// in clockwise order; and zero if they are collinear.  The
// result is also a rough approximation of twice the signed
// area of the triangle defined by the three points.
double orient2d(const double pa[2], const double pb[2], const double pc[2]);
// double orient2dfast(const double pa[2], const double pb[2], const double pc[2]);

// Return a positive value if the point pd lies below the
// plane passing through pa, pb, and pc; "below" is defined so
// that pa, pb, and pc appear in counterclockwise order when
// viewed from above the plane.  Returns a negative value if
// pd lies above the plane.  Returns zero if the points are
// coplanar.  The result is also a rough approximation of six
// times the signed volume of the tetrahedron defined by the
// four points.
double orient3d(const double pa[3], const double pb[3], const double pc[3], const double pd[3]);
// double orient3dfast(const double pa[3], const double pb[3], const double pc[3], const double pd[3]);

// Return a positive value if the point pd lies inside the
// circle passing through pa, pb, and pc; a negative value if
// it lies outside; and zero if the four points are cocircular.
// The points pa, pb, and pc must be in counterclockwise
// order, or the sign of the result will be reversed.
double incircle(const double pa[2], const double pb[2], const double pc[2], const double pd[2]);
// double incirclefast(const double pa[2], const double pb[2], const double pc[2], const double pd[2]);

// Return a positive value if the point pe lies inside the
// sphere passing through pa, pb, pc, and pd; a negative value
// if it lies outside; and zero if the five points are
// cospherical.  The points pa, pb, pc, and pd must be ordered
// so that they have a positive orientation (as defined by
// orient3d()), or the sign of the result will be reversed.
double insphere(const double pa[3], const double pb[3], const double pc[3], const double pd[3], const double pe[3]);
// double inspherefast(const double pa[3], const double pb[3], const double pc[3], const double pd[3], const double pe[3]);

// assuming pa, pb, pc and pd are coplanar
// return +1 if pd lies inside the circle passing through pa, pb, pc on the plane
// return -1 if pd lies outside the circle passing through pa, pb, pc on the plane
// return 0 if they are cocircular
int inCircumsphereOnPlane(const double pa[3], const double pb[3], const double pc[3], const double pd[3]);

bool pointInTet(const double point[3], const double teta[3], const double tetb[3], const double tetc[3], const double tetd[3]);

bool intersectSegTri(const double sega[3], const double segb[3], const double tria[3], const double trib[3], const double tric[3]);
bool intersectSegTri(const double sega[3], const double segb[3], const double tria[3], const double trib[3], const double tric[3],
  double segWeight[2], double triangleWeight[3]);

bool intersectTriTri(const double pa[3], const double pb[3], const double pc[3], const double qa[3], const double qb[3], const double qc[3]);

bool intersectTriTet(const double tria[3], const double trib[3], const double tric[3],
  const double teta[3], const double tetb[3], const double tetc[3], const double tetd[3]);

bool intersectTetTet(const double tet1a[3], const double tet1b[3], const double tet1c[3], const double tet1d[3],
  const double tet2a[3], const double tet2b[3], const double tet2c[3], const double tet2d[3]);

// line segment completely inside AABB also considered as intersecting
bool intersectSegAABB(const double sa[3], const double sb[3], const double bmin[3], const double bmax[3]);

// triangle completely inside AABB also considered as intersecting
bool intersectTriAABB(const double tria[3], const double trib[3], const double tric[3], const double bmin[3], const double bmax[3]);

// tet completely inside AABB also considered as intersecting
bool intersectTetAABB(const double teta[3], const double tetb[3], const double tetc[3], const double tetd[3], const double bmin[3], const double bmax[3]);

bool isTriangleDegenerate(const double ta[3], const double tb[3], const double tc[3]);

bool intersectSegSeg2d(const double sa[2], const double sb[2], const double ta[2], const double tb[2]);
// sw[2] and tw[2] are weight of the intersection
bool intersectSegSeg2d(const double sa[2], const double sb[2], const double ta[2], const double tb[2], double sw[2], double tw[2]);
}  // namespace Mesh
}  // namespace pgo
