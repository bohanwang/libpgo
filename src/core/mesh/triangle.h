/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "mesh" library , Copyright (C) 2007 CMU, 2009 MIT, 2018 USC           *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code authors: Jernej Barbic, Christopher Twigg, Daniel Schroeder      *
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
/*
  A triangle
  Jernej Barbic, CMU
*/

#include "meshLinearAlgebra.h"
#include "boundingBox.h"

#include <vector>
#include <iostream>
#include <cstring>
#include <cmath>
#include <type_traits>

namespace pgo
{
namespace Mesh
{
class TriangleBasic
{
public:
  TriangleBasic(const Vec3d &first, const Vec3d &second, const Vec3d &third)
  {
    vertex[0] = first;
    vertex[1] = second;
    vertex[2] = third;
  }

  // an invalid triangle
  TriangleBasic()
  {
    vertex[0].setZero();
    vertex[1].setZero();
    vertex[2].setZero();
  }

  // accessors
  inline const Vec3d &first() const { return vertex[0]; }
  inline const Vec3d &second() const { return vertex[1]; }
  inline const Vec3d &third() const { return vertex[2]; }
  inline const Vec3d &operator[](size_t i) const { return vertex[i]; }
  inline int index() const { return index_; }

  // index can be used to keep track of what triangle this is
  inline void setIndex(int index_) { this->index_ = index_; }

  // squared 3d distance to a point
  // slower than the implementation in TriangleWithCollisionInfo, but TriangleBasic does not store pre-computed data
  double distanceToPoint2(const Vec3d &point) const;
  inline double distanceToPoint(const Vec3d &point) const { return std::sqrt(TriangleBasic::distanceToPoint2(point)); }

  // return the triangle normal
  Vec3d triangleNormal() const { return (vertex[1] - vertex[0]).cross(vertex[2] - vertex[0]).normalized(); }

  // whether a point is inside/outside/on the triangle
  // return: 1 = point is outside the triangle (assuming normal points outward)
  //        -1 = point is inside
  //         0 = point is on the triangle
  int pointSide(const Vec3d &point) const;

  bool doesIntersectBox(const BoundingBox &bbox) const;

  // int lineSegmentIntersection(Vec3d segmentStart, Vec3d segmentEnd, Vec3d * intersectionPoint);
  //    Output: intersection point (when it exists)
  //    Return: -1 = triangle is degenerate (a segment or point)
  //             0 = disjoint (no intersect)
  //             1 = intersect in unique point I1
  //             2 = are in the same plane
  int lineSegmentIntersection(const Vec3d &segmentStart, const Vec3d &segmentEnd, Vec3d *intersectionPoint = nullptr, double *t = nullptr, double barycentricWeight[3] = nullptr) const;

  Vec3d getBarycentricLocation(double alpha, double beta, double gamma) const { return alpha * vertex[0] + beta * vertex[1] + gamma * vertex[2]; }

protected:
  Vec3d vertex[3];
  int index_ = -1;
};

class TriangleWithCollisionInfo : public TriangleBasic
{
public:
  TriangleWithCollisionInfo(const Vec3d &first_g, const Vec3d &second_g, const Vec3d &third_g):
    TriangleBasic(first_g, second_g, third_g) { ComputeCollisionInfo(); }

  // an invalid triangle
  inline TriangleWithCollisionInfo();

  // squared 3d distance to a point
  // also returns the closest feature to the query point:
  //  0: vertex0
  //  1: vertex1
  //  2: vertex2
  //  3: edge among 01
  //  4: edge among 12
  //  5: edge among 20
  //  6: the face itself
  double distanceToPoint2(const Vec3d &point, int *closestFeature, double *alpha = nullptr, double *beta = nullptr, double *gamma = nullptr) const;  // also returns the barycentric coordinates of the closest point
  inline double distanceToPoint2(const Vec3d &point) const
  {
    int f = -1;
    return distanceToPoint2(point, &f);
  }
  inline double distanceToPoint(const Vec3d &point) const { return std::sqrt(TriangleWithCollisionInfo::distanceToPoint2(point)); }

protected:
  // note: the following collision detection parameters are pre-computed with respect to a permuted set of triangle indices (a cycle)
  Mat3d Q = Mat3d::Zero();
  Vec3d x0 = Vec3d::Zero();
  double sidea = 0, sideb = 0, sidec = 0, area = 0;
  Vec3d S1 = Vec3d::Zero(), S2 = Vec3d::Zero(), N11 = Vec3d::Zero(), N12 = Vec3d::Zero(), N21 = Vec3d::Zero(), N22 = Vec3d::Zero();

  int permutation[6];  // used internally so that first triangle edge is always the longest
  // transformation(x) equals Q * x + x0;

  void ComputeCollisionInfo();
};

class TriangleWithCollisionInfoAndPseudoNormals : public TriangleWithCollisionInfo
{
public:
  // pseudoNormals is a caller-provided list of 7 normals, indexed the same as closest features below
  TriangleWithCollisionInfoAndPseudoNormals(const Vec3d &first_g, const Vec3d &second_g, const Vec3d &third_g, const Vec3d pseudoNormals[7]);

  // an invalid triangle
  inline TriangleWithCollisionInfoAndPseudoNormals();

  inline const Vec3d &pseudoNormal(int closestFeature) const { return pseudoNormal_[closestFeature]; }
  inline Vec3d interpolatedVertexPseudoNormal(double alpha, double beta, double gamma) const { return (alpha * pseudoNormal_[0] + beta * pseudoNormal_[1] + gamma * pseudoNormal_[2]).normalized(); }

  // for vertices, returns the vertex itself
  // for edges, returns the midpoint of the edge (for the purposes of distance sign test, this could be any point on the edge)
  // for faces, returns the face centroid
  inline const Vec3d &pseudoClosestPosition(int closestFeature) const { return pseudoClosestPosition_[closestFeature]; }

protected:
  Vec3d pseudoNormal_[7];
  Vec3d pseudoClosestPosition_[7];
};

inline TriangleWithCollisionInfo::TriangleWithCollisionInfo()
{
  memset(permutation, 0, sizeof(int) * 6);
}

inline TriangleWithCollisionInfoAndPseudoNormals::TriangleWithCollisionInfoAndPseudoNormals()
{
  for (int i = 0; i < 7; i++) {
    pseudoNormal_[i].setZero();
    pseudoClosestPosition_[i].setZero();
  }
}

template<typename T>
concept IsTriangle = is_decay_same<T, TriangleBasic> || is_decay_same<T, TriangleWithCollisionInfo> || is_decay_same<T, TriangleWithCollisionInfoAndPseudoNormals>;

// makes the triangle list unique, using the "index_" field of the <TriangleClass>
template<class TriangleClass>
  requires(IsTriangle<TriangleClass>)
void makeUniqueList(const std::vector<TriangleClass *> &triangleList, std::vector<TriangleClass *> &uniqueList);

template<class TriangleClass>
  requires(IsTriangle<TriangleClass>)
void makeUniqueList(std::vector<TriangleClass *> &triangleList);  // overwrites triangleList with the unique list

}  // namespace Mesh
}  // namespace pgo
