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

#include "geometryQuery.h"
#include "tribox3.h"
#include "triangle.h"

#include "basicAlgorithms.h"
#include "pgoLogging.h"

#include <algorithm>
#include <limits>

namespace pgo
{
namespace Mesh
{

double getTriangleAngle(const Vec3d &v0, const Vec3d &v1, const Vec3d &v2)
{
  Vec3d e1 = v1 - v0, e2 = v2 - v0;
  double cosAngle = e1.dot(e2) / (e1.norm() * e2.norm());
  cosAngle = std::clamp(cosAngle, -1.0, 1.0);
  return std::acos(cosAngle);
}

double getVectorAngle(const Vec3d &vec1, const Vec3d vec2)
{
  double cosAngle = vec1.dot(vec2) / std::sqrt(vec1.squaredNorm() * vec2.squaredNorm());
  cosAngle = std::clamp(cosAngle, -1.0, 1.0);
  return std::acos(cosAngle);
}

double getTriangleAngleRobust(const Vec3d &v0, const Vec3d &v1, const Vec3d &v2)
{
  Vec3d e1 = v1 - v0;
  Vec3d e2 = v2 - v0;
  double l2e1 = e1.squaredNorm(), l2e2 = e2.squaredNorm();
  double alpha = 0.0;
  if (l2e1 > 0 && l2e2 > 0) {
    double cosAlpha = e1.dot(e2) / std::sqrt(l2e1 * l2e2);
    cosAlpha = std::clamp(cosAlpha, -1.0, 1.0);
    alpha = std::acos(cosAlpha);
  }
  else if (l2e1 == 0 && l2e2 == 0)
    alpha = M_PI / 3;
  else
    alpha = M_PI / 2;
  return alpha;
}

double getTwoTriangleDihedralAngle(const Vec3d &e0, const Vec3d &e1, const Vec3d &t0, const Vec3d &t1)
{
  Vec3d lineVector = e1 - e0;  // the line where the edge (e0, e1) is on
  double lineLen2 = lineVector.squaredNorm();
  double invLineLen2 = 1.0 / lineLen2;

  // compute the projection of t0 onto the line: p0
  Vec3d d0 = t0 - e0;
  Vec3d p0 = e0 + (lineVector.dot(d0) * invLineLen2) * lineVector;  // the equation is gotten from the function getClosestPointToLine()
  // now compute the projection of t1 onto the line: p1
  Vec3d d1 = t1 - e0;
  Vec3d p1 = e0 + (lineVector.dot(d1) * invLineLen2) * lineVector;
  return getVectorAngle(t0 - p0, t1 - p1);
}

// double getTwoTriangleInnerDihedralAngle(const IndexedTriangle & it0, const IndexedTriangle & it1)
//{
//   Vec3d n0 = getTriangleScaledNormal(it0.pos[0], it0.pos[1], it0.pos[2]);
//   Vec3d n1 = getTriangleScaledNormal(it1.pos[0], it1.pos[1], it1.pos[2]);
//   Vec3i e = getSharedVertices(it0.vtxID, it1.vtxID);
//   PGO_ALOG(e[0] >= 0 && e[1] >= 0 && e[2] < 0); // PGO_ALOG shared is an edge
//   Vec3d t0, t1, edgePos;
//   for(int i = 0; i < 3; i++)
//     if (it0.vtxID[i] != e[0] && it0.vtxID[i] != e[1]) {
//       t0 = it0.pos[i];
//       edgePos = it0.pos[(i+1)%3]; // one pos on the joint edge
//       break;
//     }
//
//   for(int i = 0; i < 3; i++)
//     if (it1.vtxID[i] != e[0] && it1.vtxID[i] != e[1]) {
//       t1 = it1.pos[i];
//       break;
//     }
//
//   double normalAngle = getVectorAngle(n0, n1);
//   double out0 = dot(t0 - edgePos, n1); // whether t0 is outside n1
//   double out1 = dot(t1 - edgePos, n0); // whether t1 is outside n0
//   if (out0 > 0 && out1 > 0) { return M_PI + normalAngle; }
//   else if (out0 < 0 && out1 < 0) { return M_PI - normalAngle; }
// }

// PGO_ALOG len2(normal) > 0
Vec3d getClosestPointToPlaneWithScaledNormal(const Vec3d &queryPoint, const Vec3d &scaledNormal, const Vec3d &start)
{
  Vec3d diff = queryPoint - start;
  double d = scaledNormal.dot(diff);
  double normalL2 = scaledNormal.squaredNorm();
  PGO_ALOG(normalL2 > 0);
  return queryPoint - d * scaledNormal / normalL2;
}

Vec3d getClosestPointToPlaneWithNormal(const Vec3d &queryPoint, const Vec3d &normal, const Vec3d &start)
{
  Vec3d diff = queryPoint - start;
  double d = normal.dot(diff);
  return queryPoint - d * normal;
}

Vec3d getProjectedVectorToPlaneWithNormal(const Vec3d &vec, const Vec3d &planeNormal)
{
  double d = planeNormal.dot(vec);
  return vec - d * planeNormal;
}

Vec3d getClosestPointToLine(const Vec3d &queryPoint, const Vec3d &lineStart, const Vec3d &lineVector)
{
  double d = lineVector.dot(queryPoint - lineStart);
  double lineLen2 = lineVector.squaredNorm();
  return lineStart + (d / lineLen2) * lineVector;
}

Vec3d getClosestPointToLineSegment(const Vec3d &queryPoint, const Vec3d &lineStart, const Vec3d &lineEnd)
{
  Vec3d lineVec = lineEnd - lineStart;
  double d = lineVec.dot(queryPoint - lineStart);

  // the closest point is lineStart
  if (d <= 0)
    return lineStart;

  double lineLen2 = lineVec.squaredNorm();
  // the closest point is lineEnd
  if (d >= lineLen2)
    return lineEnd;

  return lineStart + (d / lineLen2) * lineVec;
}

Vec3d getClosestPointToLineSegment(const Vec3d &queryPoint, const Vec3d &lineStart, const Vec3d &lineEnd, Vec2d &barycentricWeight)
{
  Vec3d lineVec = lineEnd - lineStart;

  double d = lineVec.dot(queryPoint - lineStart);

  // the closest point is lineStart
  if (d <= 0) {
    barycentricWeight = Vec2d(1.0, 0.0);
    return lineStart;
  }

  double lineLen2 = lineVec.squaredNorm();
  // the closest point is lineEnd
  if (d >= lineLen2) {
    barycentricWeight = Vec2d(0.0, 1.0);
    return lineEnd;
  }

  barycentricWeight[1] = (d / lineLen2);
  barycentricWeight[0] = 1.0 - barycentricWeight[1];
  return lineStart + barycentricWeight[1] * lineVec;
}

Vec3d getBarycentricWeightProjectedOnTrianglePlane(const Vec3d &queryPoint, const Vec3d &t0, const Vec3d &t1, const Vec3d &t2)
{
  // w[0] * t0 + w[1] * t1 + w[2] * t2 = q', q' is the projected point of queryPoint on the triangle
  // since w[0] = 1.0 - w[2] - w[2],
  // w[1] ( t1 - t0 ) + w[2] ( t2 - t0 ) = q' - t0
  Vec3d vec1 = t1 - t0, vec2 = t2 - t0, diff = queryPoint - t0;
  // now: w[1] vec1 + w[2] vec2 = diff', again diff' here is the projected version
  // to find a way for unprojected version, we note that the equation above is true if the equation is projected onto vec1 or vec2:
  // dot(w[1] vec1 + w[2] vec2, vec1) = dot(diff', vec1) = dot(diff, vec1)
  // dot(w[1] vec1 + w[2] vec2, vec2) = dot(diff', vec2) = dot(diff, vec2)
  // so:
  // w[1] dot(vec1, vec1) + w[2] dot(vec2, vec1) = dot(diff, vec1)
  // w[1] dot(vec1, vec2) + w[2] dot(vec2, vec2) = dot(diff, vec2)
  double d11 = vec1.dot(vec1);
  double d12 = vec1.dot(vec2);
  double d22 = vec2.dot(vec2);
  double d31 = diff.dot(vec1);
  double d32 = diff.dot(vec2);

  // now:
  // w[1] d11 + w[2] d12 = d31
  // w[1] d12 + w[2] d22 = d32
  // to solve this 2x2 system:
  double invdenom = 1.0 / (d11 * d22 - d12 * d12);
  Vec3d w;
  w[1] = (d22 * d31 - d12 * d32) * invdenom;
  w[2] = (d11 * d32 - d12 * d31) * invdenom;
  w[0] = 1.0 - w[1] - w[2];
  return w;
}

Vec3d getClosestPointToTriangleWithFeature(const Vec3d &queryPoint, const Vec3d &t0, const Vec3d &t1, const Vec3d &t2, int feature)
{
  if (feature == 0)
    return t0;
  else if (feature == 1)
    return t1;
  else if (feature == 2)
    return t2;
  else if (feature == 3)  // edge 01
  {
    return getClosestPointToLineSegment(queryPoint, t0, t1);
  }
  else if (feature == 4)  // edge 12
  {
    return getClosestPointToLineSegment(queryPoint, t1, t2);
  }
  else if (feature == 5)  // edge 20
  {
    return getClosestPointToLineSegment(queryPoint, t2, t0);
  }
  // else, feature == 6, on triangle
  Vec3d scaledNormal = getTriangleScaledNormal(t0, t1, t2);
  return getClosestPointToPlaneWithScaledNormal(queryPoint, scaledNormal, t0);
}

Vec3d getClosestPointToTriangleWithNormalAndFeature(const Vec3d &queryPoint, const Vec3d &t0, const Vec3d &t1, const Vec3d &t2, const Vec3d &normal, int feature)
{
  if (feature == 0)
    return t0;
  else if (feature == 1)
    return t1;
  else if (feature == 2)
    return t2;
  else if (feature == 3)  // edge 01
  {
    return getClosestPointToLineSegment(queryPoint, t0, t1);
  }
  else if (feature == 4)  // edge 12
  {
    return getClosestPointToLineSegment(queryPoint, t1, t2);
  }
  else if (feature == 5)  // edge 20
  {
    return getClosestPointToLineSegment(queryPoint, t2, t0);
  }
  // else, feature == 6, on triangle
  return getClosestPointToPlaneWithNormal(queryPoint, normal, t0);
}

// scaledTriangleNormal must not be zero-length, but need not be unit length
// is queryPoint to the left of the edge (edgeStart -> edgeEnd)
bool isToLeftOfTriangleEdge(const Vec3d &queryPoint, const Vec3d &scaledTriangleNormal, const Vec3d &edgeStart, const Vec3d &edgeEnd)
{
  double d = (edgeEnd - edgeStart).cross(queryPoint - edgeStart).dot(scaledTriangleNormal);
  return d > 0;
}

// delta: vector from an arbitrary point on the plane to the query point
double getSquaredDistanceToPlaneWithScaledNormalAndDelta(const Vec3d &scaledPlaneNormal, const Vec3d &delta)
{
  double d = scaledPlaneNormal.dot(delta);
  double normalLen2 = scaledPlaneNormal.squaredNorm();
  return (d * d) / normalLen2;
}

double getSquaredDistanceToLineSegment(const Vec3d &queryPoint, const Vec3d &lineStart, const Vec3d &lineEnd)
{
  Vec3d lineVec = lineEnd - lineStart;

  Vec3d segStart2Query = queryPoint - lineStart;
  double d = lineVec.dot(segStart2Query);

  // return distance from segStart to query point
  if (d <= 0)
    return segStart2Query.squaredNorm();

  double lineLen2 = lineVec.squaredNorm();
  // return distance from segEnd to query point
  if (d > lineLen2)
    return (queryPoint - lineEnd).squaredNorm();

  // return distance from the infinite line to query point
  return lineVec.cross(segStart2Query).squaredNorm() / lineLen2;
}

// also returns the closest feature to the query point:
//  0: vertex0
//  1: vertex1
//  2: vertex2
//  3: edge among 01
//  4: edge among 12
//  5: edge among 20
//  6: the face itself
double getSquaredDistanceToTriangle(const Vec3d &queryPoint, const Vec3d &vertex0, const Vec3d &vertex1, const Vec3d &vertex2, int &feature)
{
  // TriangleWithCollisionInfo tc(vertex0, vertex1, vertex2);
  // return tc.distanceToPoint2(queryPoint, &feature);

  Vec3d scaledNormal = getTriangleScaledNormal(vertex0, vertex1, vertex2);

  if (scaledNormal.cwiseAbs().sum() > std::numeric_limits<double>::epsilon() &&
    isToLeftOfTriangleEdge(queryPoint, scaledNormal, vertex0, vertex1) &&
    isToLeftOfTriangleEdge(queryPoint, scaledNormal, vertex1, vertex2) &&
    isToLeftOfTriangleEdge(queryPoint, scaledNormal, vertex2, vertex0)) {
    // the closest point on triangle to queryPoint is the same closet point on the plane where the triangle lies to queryPoint
    feature = 6;
    return getSquaredDistanceToPlaneWithScaledNormalAndDelta(scaledNormal, queryPoint - vertex0);
  }
  else  // the projection of the queryPoint onto the triangle plane is outside the triangle, or the triangle is degenerate
  {     // then we query the closest distance from the query point to all the three edges
    double d0 = getSquaredDistanceToLineSegment(queryPoint, vertex1, vertex0);
    double d1 = getSquaredDistanceToLineSegment(queryPoint, vertex2, vertex1);
    double d2 = getSquaredDistanceToLineSegment(queryPoint, vertex0, vertex2);

    std::pair<double, int> sortBuffer[3] = { { d0, 0 }, { d1, 1 }, { d2, 2 } };
    std::sort(sortBuffer, sortBuffer + 3);
    if (sortBuffer[0].first == sortBuffer[1].first) {
      // closest feature is a vertex
      const int edgeIDToFeatureMap[3][3] = { { -1, 1, 0 }, { 1, -1, 2 }, { 0, 2, -1 } };
      feature = edgeIDToFeatureMap[sortBuffer[0].second][sortBuffer[1].second];
      PGO_ALOG(feature >= 0);
    }
    else {  // closest feature is an edge
      feature = sortBuffer[0].second + 3;
    }
    return sortBuffer[0].first;
  }
}

double getSquaredDistanceToTriangle(const Vec3d &queryPoint, const Vec3d &vertex0, const Vec3d &vertex1, const Vec3d &vertex2, int &feature,
  Vec3d &closestPointOnTriangle, Vec3d &barycentricWeight)
{
  Vec3d scaledNormal = getTriangleScaledNormal(vertex0, vertex1, vertex2);

  if (scaledNormal.cwiseAbs().sum() > std::numeric_limits<double>::epsilon() &&
    isToLeftOfTriangleEdge(queryPoint, scaledNormal, vertex0, vertex1) &&
    isToLeftOfTriangleEdge(queryPoint, scaledNormal, vertex1, vertex2) &&
    isToLeftOfTriangleEdge(queryPoint, scaledNormal, vertex2, vertex0)) {
    // the closest point on triangle to queryPoint is the same closet point on the plane where the triangle lies to queryPoint
    feature = 6;
    barycentricWeight = getBarycentricWeightProjectedOnTrianglePlane(queryPoint, vertex0, vertex1, vertex2);
    closestPointOnTriangle = barycentricWeight[0] * vertex0 + barycentricWeight[1] * vertex1 + barycentricWeight[2] * vertex2;
    return (closestPointOnTriangle - queryPoint).squaredNorm();
  }
  else  // the projection of the queryPoint onto the triangle plane is outside the triangle, or the triangle is degenerate
  {     // then we query the closest distance from the query point to all the three edges
    double d0 = getSquaredDistanceToLineSegment(queryPoint, vertex1, vertex0);
    double d1 = getSquaredDistanceToLineSegment(queryPoint, vertex2, vertex1);
    double d2 = getSquaredDistanceToLineSegment(queryPoint, vertex0, vertex2);

    std::pair<double, int> sortBuffer[3] = { { d0, 0 }, { d1, 1 }, { d2, 2 } };
    std::sort(sortBuffer, sortBuffer + 3);
    if (sortBuffer[0].first == sortBuffer[1].first) {
      // closest feature is a vertex
      const int edgeIDToFeatureMap[3][3] = { { -1, 1, 0 }, { 1, -1, 2 }, { 0, 2, -1 } };
      feature = edgeIDToFeatureMap[sortBuffer[0].second][sortBuffer[1].second];
      PGO_ALOG(feature >= 0);
      barycentricWeight.setZero();
      barycentricWeight[feature] = 1.0;
      const Vec3d *v[3] = { &vertex0, &vertex1, &vertex2 };
      closestPointOnTriangle = *(v[feature]);
    }
    else {       // closest feature is an edge
      feature = sortBuffer[0].second + 3;
      Vec2d ew;  // weight on edge
      if (feature == 3) {
        closestPointOnTriangle = getClosestPointToLineSegment(queryPoint, vertex0, vertex1, ew);
        barycentricWeight[0] = ew[0];
        barycentricWeight[1] = ew[1];
        barycentricWeight[2] = 0.0;
      }
      else if (feature == 4) {
        closestPointOnTriangle = getClosestPointToLineSegment(queryPoint, vertex1, vertex2, ew);
        barycentricWeight[0] = 0.0;
        barycentricWeight[1] = ew[0];
        barycentricWeight[2] = ew[1];
      }
      else  // feature == 5
      {
        closestPointOnTriangle = getClosestPointToLineSegment(queryPoint, vertex2, vertex0, ew);
        barycentricWeight[0] = ew[1];
        barycentricWeight[1] = 0.0;
        barycentricWeight[2] = ew[0];
      }
    }
    return sortBuffer[0].first;
  }
}

double whetherTriangleIntersectBoundingBox(const Vec3d &t0, const Vec3d &t1, const Vec3d &t2, const Vec3d &bmin, const Vec3d &bmax)
{
  const Vec3d center = 0.5 * (bmin + bmax);
  const Vec3d halfSides = 0.5 * (bmax - bmin);

  return triBoxOverlap(center.data(), halfSides.data(), t0.data(), t1.data(), t2.data());
}

double interpolateValueOnTriangle(const Vec3i &triVtxIDs, const Vec3d &triVtxWeights, const double *vtxDataVector)
{
  double v = 0.0;
  for (int d = 0; d < 3; d++)
    if (triVtxIDs[d] >= 0)
      v += vtxDataVector[triVtxIDs[d]] * triVtxWeights[d];
  return v;
}

Mat3d getPointInertiaTensor(const Vec3d &position, double mass)
{
  //  double correction[6] =
  //       { b*b + c*c, -a*b, -a*c,
  //               a*a + c*c, -b*c,
  //                     a*a + b*b };
  double a = position[0];
  double b = position[1];
  double c = position[2];
  return mass * asMat3d(b * b + c * c, -a * b, -a * c, -a * b, a * a + c * c, -b * c, -a * c, -b * c, a * a + b * b);
}

void getUnitMassTriangleInertiaTensorVectorAroundOrigin(const Vec3d &v0, const Vec3d &v1, const Vec3d &v2, double t[6])
{
  t[0] = (v0[1] * v0[1] + v0[2] * v0[2] + v1[1] * v1[1] + v1[2] * v1[2] +
           v1[1] * v2[1] + v2[1] * v2[1] + v0[1] * (v1[1] + v2[1]) +
           v1[2] * v2[2] + v2[2] * v2[2] + v0[2] * (v1[2] + v2[2])) /
    6;

  t[1] = (-2 * v1[0] * v1[1] - v1[1] * v2[0] - v0[1] * (v1[0] + v2[0]) - v1[0] * v2[1] - 2 * v2[0] * v2[1] -
           v0[0] * (2 * v0[1] + v1[1] + v2[1])) /
    12;

  t[2] = (-2 * v1[0] * v1[2] - v1[2] * v2[0] - v0[2] * (v1[0] + v2[0]) -
           v1[0] * v2[2] - 2 * v2[0] * v2[2] -
           v0[0] * (2 * v0[2] + v1[2] + v2[2])) /
    12;

  t[3] = (v0[0] * v0[0] + v0[2] * v0[2] + v1[0] * v1[0] + v1[2] * v1[2] +
           v1[0] * v2[0] + v2[0] * v2[0] + v0[0] * (v1[0] + v2[0]) +
           v1[2] * v2[2] + v2[2] * v2[2] + v0[2] * (v1[2] + v2[2])) /
    6;

  t[4] = (-2 * v1[1] * v1[2] - v1[2] * v2[1] -
           v0[2] * (v1[1] + v2[1]) - v1[1] * v2[2] - 2 * v2[1] * v2[2] -
           v0[1] * (2 * v0[2] + v1[2] + v2[2])) /
    12;

  t[5] = (v0[0] * v0[0] + v0[1] * v0[1] + v1[0] * v1[0] + v1[1] * v1[1] +
           v1[0] * v2[0] + v2[0] * v2[0] + v0[0] * (v1[0] + v2[0]) +
           v1[1] * v2[1] + v2[1] * v2[1] + v0[1] * (v1[1] + v2[1])) /
    6;
}

Mat3d getTriangleInertiaTensorAroundOrigin(const Vec3d &v0, const Vec3d &v1, const Vec3d &v2, double surfaceArea)
{
  double t[6];
  getTriangleInertiaTensorVectorAroundOrigin(v0, v1, v2, surfaceArea, t);

  Mat3d IT;
  IT(0, 0) = t[0];
  IT(0, 1) = t[1];
  IT(0, 2) = t[2];

  IT(1, 0) = t[1];
  IT(1, 1) = t[3];
  IT(1, 2) = t[4];

  IT(2, 0) = t[2];
  IT(2, 1) = t[4];
  IT(2, 2) = t[5];

  return IT;
}

void getTetBarycentricWeights(const Vec3d &queryPoint, const Vec3d &a, const Vec3d &b, const Vec3d &c, const Vec3d &d, double weight[4])
{
  //       |x1 y1 z1 1|         |x  y  z  1|        |x1 y1 z1 1|        |x1 y1 z1 1|        |x1 y1 z1 1|
  //  D0 = |x2 y2 z2 1|   D1 =  |x2 y2 z2 1|   D2 = |x  y  z  1|   D3 = |x2 y2 z2 1|   D4 = |x2 y2 z2 1|
  //       |x3 y3 z3 1|         |x3 y3 z3 1|        |x3 y3 z3 1|        |x  y  z  1|        |x3 y3 z3 1|
  //       |x4 y4 z4 1|         |x4 y4 z4 1|        |x4 y4 z4 1|        |x4 y4 z4 1|        |x  y  z  1|
  //  wi = Di / D0

  double tetDet = getTetDeterminant(a, b, c, d);

  for (int i = 0; i < 4; i++) {
    // compute D[i+1]
    Vec3d buf[4] = { a, b, c, d };
    buf[i] = queryPoint;
    double D = getTetDeterminant(buf[0], buf[1], buf[2], buf[3]);
    weight[i] = D / tetDet;
  }
}

double getSquaredDistanceToTet(const Vec3d &queryPoint, const Vec3d &a, const Vec3d &b, const Vec3d &c, const Vec3d &d)
{
  //  const int OTetKey::tetFaceIndex[4][3] =  { { 1, 2, 3 }, { 0, 3, 2 }, { 0, 1, 3 }, { 0, 2, 1 } };
  Vec3d diff = queryPoint - a;
  Vec3d scaledTriNormal = getTriangleScaledNormal(a, d, c);  // the tet face opposite vertex b
  double scaledSignedDistance = diff.dot(scaledTriNormal);   // the scaled signed distance from queryPoint to the face
  if (scaledSignedDistance > 0)                              // it the vertex is outside this tet face,
                                                             // then the distance to the tet is the distance to this face
    return getSquaredDistanceToTriangle(queryPoint, a, d, c);

  scaledTriNormal = getTriangleScaledNormal(a, b, d);  // tet face opposite vertex c
  scaledSignedDistance = diff.dot(scaledTriNormal);
  if (scaledSignedDistance > 0)
    return getSquaredDistanceToTriangle(queryPoint, a, b, d);

  scaledTriNormal = getTriangleScaledNormal(a, c, b);  // tet face opposite vertex d
  scaledSignedDistance = diff.dot(scaledTriNormal);
  if (scaledSignedDistance > 0)
    return getSquaredDistanceToTriangle(queryPoint, a, c, b);

  diff = queryPoint - b;
  scaledTriNormal = getTriangleScaledNormal(b, c, d);  // tet face opposite vertex a
  scaledSignedDistance = diff.dot(scaledTriNormal);
  if (scaledSignedDistance > 0)
    return getSquaredDistanceToTriangle(queryPoint, b, c, d);

  return 0.0;  // otherwise, queryPoint is inside the tet
}

Mat3d getTetInertiaTensorAroudOrigin(const Vec3d &a, const Vec3d &b, const Vec3d &c, const Vec3d &d, double tetDeterminant)
{
  double absdetJ = tetDeterminant;

  double x1 = a[0], x2 = b[0], x3 = c[0], x4 = d[0];
  double y1 = a[1], y2 = b[1], y3 = c[1], y4 = d[1];
  double z1 = a[2], z2 = b[2], z3 = c[2], z4 = d[2];

  double A = absdetJ * (y1 * y1 + y1 * y2 + y2 * y2 + y1 * y3 + y2 * y3 + y3 * y3 + y1 * y4 + y2 * y4 + y3 * y4 + y4 * y4 + z1 * z1 + z1 * z2 + z2 * z2 + z1 * z3 + z2 * z3 + z3 * z3 + z1 * z4 + z2 * z4 + z3 * z4 + z4 * z4) / 60.0;

  double B = absdetJ * (x1 * x1 + x1 * x2 + x2 * x2 + x1 * x3 + x2 * x3 + x3 * x3 + x1 * x4 + x2 * x4 + x3 * x4 + x4 * x4 + z1 * z1 + z1 * z2 + z2 * z2 + z1 * z3 + z2 * z3 + z3 * z3 + z1 * z4 + z2 * z4 + z3 * z4 + z4 * z4) / 60.0;

  double C = absdetJ * (x1 * x1 + x1 * x2 + x2 * x2 + x1 * x3 + x2 * x3 + x3 * x3 + x1 * x4 + x2 * x4 + x3 * x4 + x4 * x4 + y1 * y1 + y1 * y2 + y2 * y2 + y1 * y3 + y2 * y3 + y3 * y3 + y1 * y4 + y2 * y4 + y3 * y4 + y4 * y4) / 60.0;

  double Ap = absdetJ * (2 * y1 * z1 + y2 * z1 + y3 * z1 + y4 * z1 + y1 * z2 + 2 * y2 * z2 + y3 * z2 + y4 * z2 + y1 * z3 + y2 * z3 + 2 * y3 * z3 + y4 * z3 + y1 * z4 + y2 * z4 + y3 * z4 + 2 * y4 * z4) / 120.0;

  double Bp = absdetJ * (2 * x1 * z1 + x2 * z1 + x3 * z1 + x4 * z1 + x1 * z2 + 2 * x2 * z2 + x3 * z2 + x4 * z2 + x1 * z3 + x2 * z3 + 2 * x3 * z3 + x4 * z3 + x1 * z4 + x2 * z4 + x3 * z4 + 2 * x4 * z4) / 120.0;

  double Cp = absdetJ * (2 * x1 * y1 + x2 * y1 + x3 * y1 + x4 * y1 + x1 * y2 + 2 * x2 * y2 + x3 * y2 + x4 * y2 + x1 * y3 + x2 * y3 + 2 * x3 * y3 + x4 * y3 + x1 * y4 + x2 * y4 + x3 * y4 + 2 * x4 * y4) / 120.0;

  return asMat3d(A, -Bp, -Cp, -Bp, B, -Ap, -Cp, -Ap, C);
}

Mat3d getTetInertiaTensorAroundCOM(const Vec3d &t0, const Vec3d &t1, const Vec3d &t2, const Vec3d &t3, const Vec3d &tetCenterOfMass,
  double tetDeterminant)
{
  Vec3d a = t0 - tetCenterOfMass;
  Vec3d b = t1 - tetCenterOfMass;
  Vec3d c = t2 - tetCenterOfMass;
  Vec3d d = t3 - tetCenterOfMass;

  return getTetInertiaTensorAroudOrigin(a, b, c, d, tetDeterminant);
}

Mat3d getBoxInertiaTensorAroundCOM(const Vec3d &side, double volume)
{
  // For a solid cuboid of width w, height h, depth d, and mass m
  //      [ h^2 + d^2                     ]
  // m/12 [           w^2 + d^2           ]
  //      [                     w^2 + h^2 ]
  double mover12 = volume / 12.0;
  return asMat3d(mover12 * side[1] * side[1] + side[2] * side[2], 0.0, 0.0,
    0.0, mover12 * side[0] * side[0] + side[2] * side[2], 0.0,
    0.0, 0.0, mover12 * side[0] * side[0] + side[1] * side[1]);
}

void getBoundingBoxFromSelectedVertices(int numVtx, const int *vtxIndices, const Vec3d *vertices, Vec3d &bmin, Vec3d &bmax)
{
  bmin = bmax = vertices[vtxIndices[0]];
  for (int j = 1; j < numVtx; j++) {
    const Vec3d &v = vertices[vtxIndices[j]];
    for (int k = 0; k < 3; k++) {
      if (v[k] < bmin[k])
        bmin[k] = v[k];
      else if (v[k] > bmax[k])
        bmax[k] = v[k];
    }
  }
}

}  // namespace Mesh
}  // namespace pgo