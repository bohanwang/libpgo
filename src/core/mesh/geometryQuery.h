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

#include <cmath>

namespace pgo
{
namespace Mesh
{

inline double rad2deg(double x)
{
  return x * (180.0 / M_PI);
}
inline double deg2rad(double x)
{
  return x * (M_PI / 180.0);
}

// get angle at v0: /_v1v0v2
double getTriangleAngle(const Vec3d &v0, const Vec3d &v1, const Vec3d &v2);

// get angle between vec0 and vec1, ranging in [0, pi]
double getVectorAngle(const Vec3d &vec0, const Vec3d vec1);

// get the scaled normal of triangle <v0, v1, v2>: cross(v1-v0, v2-v0) = 2 * tri_area * n
// its length is twice the triangle area
inline Vec3d getTriangleScaledNormal(const Vec3d &v0, const Vec3d &v1, const Vec3d &v2)
{
  return (v1 - v0).cross(v2 - v0);
}

// get the normal of triangle <v0, v1, v2>
inline Vec3d getTriangleNormal(const Vec3d &v0, const Vec3d &v1, const Vec3d &v2)
{
  return getTriangleScaledNormal(v0, v1, v2).normalized();
}

inline double getTriangleArea(const Vec3d &v0, const Vec3d &v1, const Vec3d &v2)
{
  return 0.5 * getTriangleScaledNormal(v0, v1, v2).norm();
}
inline Vec3d getTriangleCenterOfMass(const Vec3d &v0, const Vec3d &v1, const Vec3d &v2)
{
  return (v0 + v1 + v2) * (1.0 / 3.0);
}

// return dot (queryPoint - v0, cross(v1-v0, v2-v0) )
// the returned value is positive if queryPoint is above triangle, negative if under and zero if on the plane where the triangle lies
// the returned value is 2 * tri_area * signed_distance
inline double getScaledSignedDistanceToTrianglePlane(const Vec3d &queryPoint, const Vec3d &v0, const Vec3d &v1, const Vec3d &v2)
{
  return (queryPoint - v0).dot(getTriangleScaledNormal(v0, v1, v2));
}

// robust computation of the angle at v0 in a triangle v0, v1, v2
// if any two points have the same position, it treats the triangle as a degenerate one with angle 0, 90 and 90
// if all three points have the same position, it treats the triangle angles as 60, 60, 60
double getTriangleAngleRobust(const Vec3d &v0, const Vec3d &v1, const Vec3d &v2);

// get the dihedral angle between triangle (e0, e1, t0) and (e0, e1, t1), triangle ordering does not matter
// assuming the largest dihedral angle <= PI, return value of [0, PI]
double getTwoTriangleDihedralAngle(const Vec3d &e0, const Vec3d &e1, const Vec3d &t0, const Vec3d &t1);

// assert len(scaledNormal) > 0.0
Vec3d getClosestPointToPlaneWithScaledNormal(const Vec3d &queryPoint, const Vec3d &scaledNormal, const Vec3d &planeStart);
// asssert len(normal) == 1.0
Vec3d getClosestPointToPlaneWithNormal(const Vec3d &queryPoint, const Vec3d &normal, const Vec3d &planeStart);

// project vector vec onto plane with normal "planeNormal"
Vec3d getProjectedVectorToPlaneWithNormal(const Vec3d &vec, const Vec3d &planeNormal);

Vec3d getClosestPointToLineSegment(const Vec3d &queryPoint, const Vec3d &segStart, const Vec3d &segEnd);
Vec3d getClosestPointToLineSegment(const Vec3d &queryPoint, const Vec3d &segStart, const Vec3d &segEnd, Vec2d &barycentricWeight);

// Note here: the points on the line is: p = lineStart + k * lineVector
// lineVector does not need to be normalized
Vec3d getClosestPointToLine(const Vec3d &queryPoint, const Vec3d &lineStart, const Vec3d &lineVector);

// assume triangle (t0, t1, t2) is not degenerate
// queryPoint does not need to be on the triangle, the computed weight is the weight for the projected version of queryPoint on the triangle plane
Vec3d getBarycentricWeightProjectedOnTrianglePlane(const Vec3d &queryPoint, const Vec3d &t0, const Vec3d &t1, const Vec3d &t2);

double minimalDistance2OfTwoLineSegments(const Vec3d &p0, const Vec3d &p1, const Vec3d &q0, const Vec3d &q1);

// also input the closest feature to the query point:
//  0: vertex0
//  1: vertex1
//  2: vertex2
//  3: edge among 01
//  4: edge among 12
//  5: edge among 20
//  6: the face itself
Vec3d getClosestPointToTriangleWithFeature(const Vec3d &queryPoint, const Vec3d &t0, const Vec3d &t1, const Vec3d &t2, int feature);
Vec3d getClosestPointToTriangleWithNormalAndFeature(const Vec3d &queryPoint, const Vec3d &t0, const Vec3d &t1, const Vec3d &t2, const Vec3d &normal, int feature);

double getSquaredDistanceToTriangle(const Vec3d &queryPoint, const Vec3d &t0, const Vec3d &t1, const Vec3d &t2, int &feature);
double getSquaredDistanceToTriangle(const Vec3d &queryPoint, const Vec3d &t0, const Vec3d &t1, const Vec3d &t2, int &feature,
  Vec3d &closestPointOnTriangle, Vec3d &barycentricWeight);
inline double getSquaredDistanceToTriangle(const Vec3d &queryPoint, const Vec3d &t0, const Vec3d &t1, const Vec3d &t2)
{
  int f = 0;
  return getSquaredDistanceToTriangle(queryPoint, t0, t1, t2, f);
}

double getSquaredDistanceToLineSegment(const Vec3d &queryPoint, const Vec3d &segStart, const Vec3d &segEnd);

double whetherTriangleIntersectBoundingBox(const Vec3d &t0, const Vec3d &t1, const Vec3d &t2, const Vec3d &bmin, const Vec3d &bmax);

// triVtxIDs[i] can be < 0, in this case we skip this index
double interpolateValueOnTriangle(const Vec3i &triVtxIDs, const Vec3d &triVtxWeights, const double *vtxDataVector);

// ==============================================================================
//                                  Tetrahedra
// ==============================================================================

// computes det(A), for the 4x4 matrix A
//     [ 1 a ]
// A = [ 1 b ]
//     [ 1 c ]
//     [ 1 d ]
// It can also be computed as det(A) = dot(d - a, cross(b - a, c - a))
// When det(A) > 0, the tet has positive orientation.
// When det(A) = 0, the tet is degenerate.
// When det(A) < 0, the tet has negative orientation.
// The orientation can also be determined as:
// if a is under the plane of the triangle formed by <b, c, d>, then it has positive orientation
inline double getTetDeterminant(const Vec3d &a, const Vec3d &b, const Vec3d &c, const Vec3d &d)
{
  return (d - a).dot((b - a).cross(c - a));
}

// get the volume of the tet
inline double getTetVolume(const Vec3d &a, const Vec3d &b, const Vec3d &c, const Vec3d &d)
{
  return (1.0 / 6 * fabs(getTetDeterminant(a, b, c, d)));
}

// get the volume of the tet combined with the sign from tet determinant
inline double getTetSignedVolume(const Vec3d &a, const Vec3d &b, const Vec3d &c, const Vec3d &d)
{
  return 1.0 / 6 * getTetDeterminant(a, b, c, d);
}

inline Vec3d getTetCenterOfMass(const Vec3d &a, const Vec3d &b, const Vec3d &c, const Vec3d &d)
{
  return (a + b + c + d) * 0.25;
}

// compute barycentric weights of tet <a,b,c,d> for queryPoint
void getTetBarycentricWeights(const Vec3d &queryPoint, const Vec3d &a, const Vec3d &b, const Vec3d &c, const Vec3d &d, double weight[4]);

double getSquaredDistanceToTet(const Vec3d &queryPoint, const Vec3d &a, const Vec3d &b, const Vec3d &c, const Vec3d &d);

// assuming numVtx > 0
void getBoundingBoxFromSelectedVertices(int numVtx, const int *vtxIndices, const Vec3d *vertices, Vec3d &bmin, Vec3d &bmax);

// ==============================================================================
//                              Inertia Tensors
// ==============================================================================

// note: for point at origin, its inertia tensor w.r.t origin is a zero matrix
// vectorBetweenMassPointAndReferencePoint can be computed from mass-point to the reference point, or in reverse.
// It does not matter.
Mat3d getPointInertiaTensor(const Vec3d &vectorBetweenMassPointAndReferencePoint, double mass);

// use parallel axis theorem to compute a new inertia tensor from a reference point given the inertia tensor around the center of mass
// note vectorBetweenMassCenterAndReferencePoint can be computed either by from center of mass to reference point, or by
// from reference point to center of mass. It does not matter in parallel axis theorem.
inline Mat3d shiftInertiaTensorAroundMassCenterToReferencePoint(const Mat3d &inertiaTensorAroundMassCenter, double mass,
  const Vec3d &vectorBetweenMassCenterAndReferencePoint)
{
  return inertiaTensorAroundMassCenter + getPointInertiaTensor(vectorBetweenMassCenterAndReferencePoint, mass);
}

inline Mat3d shiftInertiaTensorAroundReferencePointToMassCenter(const Mat3d &inertiaTensorAroundReferencePoint, double mass,
  const Vec3d &vectorBetweenMassCenterAndReferencePoint)
{
  return inertiaTensorAroundReferencePoint - getPointInertiaTensor(vectorBetweenMassCenterAndReferencePoint, mass);
}

// ----------------------------- Triangle Tensor --------------------------------

// assume triangle has mass of 1
// tensorVec is the upper-triangular part of the symmetric inertia tensor
void getUnitMassTriangleInertiaTensorVectorAroundOrigin(const Vec3d &t0, const Vec3d &t1, const Vec3d &t2, double tensorVec[6]);

// assume triangle has uniform density of 1
// compute its inertia tensor around the origin
Mat3d getTriangleInertiaTensorAroundOrigin(const Vec3d &t0, const Vec3d &t1, const Vec3d &t2, double triangleSurfaceArea);

inline Mat3d getTriangleInertiaTensorAroundOrigin(const Vec3d &t0, const Vec3d &t1, const Vec3d &t2)
{
  return getTriangleInertiaTensorAroundOrigin(t0, t1, t2, getTriangleArea(t0, t1, t2));
}
// AroundCOM: compute triangle inertia tensor around its center of mass
inline Mat3d getTriangleInertiaTensorAroundCOM(const Vec3d &t0, const Vec3d &t1, const Vec3d &t2, const Vec3d &triangleCenterOfMass,
  double triangleSurfaceArea)
{
  return getTriangleInertiaTensorAroundOrigin(t0 - triangleCenterOfMass, t1 - triangleCenterOfMass, t2 - triangleCenterOfMass,
    triangleSurfaceArea);
}

inline Mat3d getTriangleInertiaTensorAroundCOM(const Vec3d &t0, const Vec3d &t1, const Vec3d &t2)
{
  return getTriangleInertiaTensorAroundCOM(t0, t1, t2, getTriangleCenterOfMass(t0, t1, t2), getTriangleArea(t0, t1, t2));
}

// assume triangle has uniform density of 1
inline void getTriangleInertiaTensorVectorAroundOrigin(const Vec3d &t0, const Vec3d &t1, const Vec3d &t2, double triangleArea,
  double tensorVec[6])
{
  getUnitMassTriangleInertiaTensorVectorAroundOrigin(t0, t1, t2, tensorVec);
  for (int i = 0; i < 6; i++)
    tensorVec[i] *= triangleArea;
}

inline void getTriangleInertiaTensorVectorAroundOrigin(const Vec3d &t0, const Vec3d &t1, const Vec3d &t2, double tensorVec[6])
{
  return getTriangleInertiaTensorVectorAroundOrigin(t0, t1, t2, getTriangleArea(t0, t1, t2), tensorVec);
}

// ----------------------------------- Tet Tensor ----------------------------------

// assume tet has uniform density of 1
// AroundCOM: compute tet inertia tensor around its center of mass
// Note: don't know how the orientation of the tet affect this inertia tensor
Mat3d getTetInertiaTensorAroundCOM(const Vec3d &a, const Vec3d &b, const Vec3d &c, const Vec3d &d, const Vec3d &tetCenterOfMass,
  double tetDeterminant);
inline Mat3d getTetInertiaTensorAroundCOM(const Vec3d &a, const Vec3d &b, const Vec3d &c, const Vec3d &d)
{
  return getTetInertiaTensorAroundCOM(a, b, c, d, getTetCenterOfMass(a, b, c, d), getTetDeterminant(a, b, c, d));
}

Mat3d getTetInertiaTensorAroudOrigin(const Vec3d &a, const Vec3d &b, const Vec3d &c, const Vec3d &d, double tetDeterminant);

// --------------------------------- Box Tensor ------------------------------------

// assuming box has has uniform density of 1, side = bmax - bmin
Mat3d getBoxInertiaTensorAroundCOM(const Vec3d &side, double volume);

}  // namespace Mesh
}  // namespace pgo
