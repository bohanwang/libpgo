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
//  Bounding Box
//  Author: Jernej Barbic, CMU

#include "meshLinearAlgebra.h"

#include <vector>
#include <cfloat>
#include <cmath>
#include <iostream>
#include <algorithm>

namespace pgo
{
namespace Mesh
{

class TriangleBasic;
class TriangleWithCollisionInfo;
class TriangleWithCollisionInfoAndPseudoNormals;

// axis-aligned bounding box
// the default BoundingBox considers touching as intersection
// it also stores center and halfSides for faster queries
class BoundingBox
{
public:
  BoundingBox():
    bmin_(0.0, 0.0, 0.0), bmax_(1.0, 1.0, 1.0), center_(0.5, 0.5, 0.5), halfSides_(0.5, 0.5, 0.5) {}
  BoundingBox(Vec3d bmin_g, Vec3d bmax_g):
    bmin_(bmin_g), bmax_(bmax_g) { updateData(); }

  // an array of vertices
  BoundingBox(int numVertices, const Vec3d *vertices);
  // an array of bounding boxes
  BoundingBox(int numBBs, const BoundingBox *bbs);

  // a container holding vertices
  template<class Vec3dContainer>
  explicit BoundingBox(const Vec3dContainer &vertices);

  // an array of all the vertices and an int container holding selected indices of the vertices
  template<class IntContainer>
  explicit BoundingBox(const Vec3d *allVertices, const IntContainer &vertexIDs);

  // an array of all the vertices, an array of all the triangles and an int container holding selected indices of the triangles
  template<class IntContainer>
  explicit BoundingBox(const Vec3d *allVertices, const Vec3i *allTriangles, const IntContainer &triIDs);

  // an array of all the vertices, an array of all the tets and an int container holding selected indices of the tets
  template<class IntContainer>
  explicit BoundingBox(const Vec3d *allVertices, const Vec4i *allTets, const IntContainer &tetIDs);

  // an array of all the bounding boxes and an int container holding selected indces of the bounding boxes
  template<class IntContainer>
  explicit BoundingBox(const BoundingBox *allBBs, const IntContainer &bbIDs);

  explicit BoundingBox(const std::vector<BoundingBox> &bbs):
    BoundingBox(static_cast<int>(bbs.size()), bbs.data()) {}
  explicit BoundingBox(const std::vector<TriangleBasic> &bbs);
  explicit BoundingBox(const std::vector<TriangleWithCollisionInfo> &bbs);
  explicit BoundingBox(const std::vector<TriangleWithCollisionInfoAndPseudoNormals> &bbs);

  // accessors
  const Vec3d &bmin() const { return bmin_; }
  const Vec3d &bmax() const { return bmax_; }

  const Vec3d &center() const { return center_; }
  const Vec3d &halfSides() const { return halfSides_; }

  double diameter() const { return 2.0 * halfSides_.norm(); }
  Vec3d sides() const { return bmax_ - bmin_; }
  double volume() const { return halfSides_[0] * halfSides_[1] * halfSides_[2] * 8; }

  // mutators
  void setbmin(const Vec3d &bmin_g)
  {
    bmin_ = bmin_g;
    updateData();
  }

  void setbmin(double x, double y, double z)
  {
    bmin_ = Vec3d(x, y, z);
    updateData();
  }

  void setbmax(const Vec3d &bmax_g)
  {
    bmax_ = bmax_g;
    updateData();
  }

  void setbmax(double x, double y, double z)
  {
    bmax_ = Vec3d(x, y, z);
    updateData();
  }

  void setbminmax(const Vec3d &bmin, const Vec3d &bmax)
  {
    bmin_ = bmin;
    bmax_ = bmax;
    updateData();
  }

  double distanceToPoint(const Vec3d &point) const { return std::sqrt(distanceToPoint2(point)); }
  double furthestDistanceToPoint(const Vec3d &point) const { return std::sqrt(furthestDistanceToPoint2(point)); }

  // get squared distance
  // computing squared distance is much faster than computing distance because the later requires std::sqrt()
  double distanceToPoint2(const Vec3d &point) const;
  double furthestDistanceToPoint2(const Vec3d &point) const;

  // return true if point is inside or touching the bounding box
  bool checkInside(const Vec3d &point) const;
  // return true if bb is completely inside or touching from inside the bounding box
  bool checkInside(const BoundingBox &bb) const;

  // sanity check bmin <= bmax
  bool verifyBox() const;

  // expands from the center
  // factor of 1.0 indicates no expansion
  void expand(double expansionFactor);
  // expand the bounding box to include this point
  void expand(const Vec3d &point);
  // expand the bounding box to include another bounding box
  void expand(const BoundingBox &bb);

  void regularize();  // expand the box into one with all sides equal

  bool lineSegmentIntersection(const Vec3d &segmentStart, const Vec3d &segmentEnd, Vec3d *intersection) const;
  bool intersect(const BoundingBox &bb) const;

  // return the index of the longest side (bmax[i] - bmin[i]) and its value
  std::pair<int, double> longestSide() const;

  friend inline std::ostream &operator<<(std::ostream &o, const BoundingBox &bb)
  {
    return o << "bmin:" << bb.bmin_.transpose() << "\nbmax:" << bb.bmax_.transpose();
  }

  void print() const { std::cout << (*this) << std::endl; }

  // The ordering of children is as follows:
  // child ID -> < b2_b1_b0 > (binary notation), each bi is 0 or 1
  // bi = 0: in dimension i, this child is on the negative side
  // bi - 1: in dimension i, this child is on the positive side
  void createChildBoundingBoxes(BoundingBox children[8]) const;

  // get intersection of two bounding boxes
  // will return an invalid bounding box if this and bb does not intersect
  BoundingBox getIntersection(const BoundingBox &bb) const;

protected:
  template<class Triangle>
  void buildFromTriangles(const std::vector<Triangle> &tripool);

  inline void updateOnePoint(const Vec3d &p);               // update only bmin_ bmax_ from one point
  inline void updateOneBoundingBox(const BoundingBox &bb);  // update only bmin_ bmax_ from one bounding box
  inline void updateData();                                 // updates center and half-sides
  Vec3d bmin_, bmax_;
  Vec3d center_, halfSides_;
};

// axis-aligned bounding box for scenerios which require little memory usage and computation
// the LightBoundingBox does not consider touching as intersection
// it only stores bmin and bmax

class LightBoundingBox
{
public:
  LightBoundingBox():
    bmin_(DBL_MAX, DBL_MAX, DBL_MAX), bmax_(-DBL_MAX, -DBL_MAX, -DBL_MAX) {}
  LightBoundingBox(Vec3d bmin_g, Vec3d bmax_g):
    bmin_(bmin_g), bmax_(bmax_g) {}

  // an array of vertices
  LightBoundingBox(int numVertices, const Vec3d *vertices);

  // a container holding vertices
  template<class Vec3dContainer>
  explicit LightBoundingBox(const Vec3dContainer &vertices);

  // an array of all the vertices and an int container holding selected indices of the vertices
  template<class IntContainer>
  explicit LightBoundingBox(const Vec3d *allVertices, const IntContainer &vertexIDs);

  // an array of all the vertices, an array of all the triangles and an int container holding selected indices of the triangles
  template<class IntContainer>
  explicit LightBoundingBox(const Vec3d *allVertices, const Vec3i *allTriangles, const IntContainer &triIDs);

  // an array of all the vertices, an array of all the tets and an int container holding selected indices of the tets
  template<class IntContainer>
  explicit LightBoundingBox(const Vec3d *allVertices, const Vec4i *allTets, const IntContainer &tetIDs);

  // accessors
  const Vec3d &bmin() const { return bmin_; }
  const Vec3d &bmax() const { return bmax_; }

  Vec3d center() const { return (bmin_ + bmax_) / 2.0; }
  Vec3d halfSides() const { return (bmax_ - bmin_) / 2.0; }
  Vec3d sides() const { return bmax_ - bmin_; }

  double diameter() const { return (bmax_ - bmin_).norm(); }
  double volume() const
  {
    Vec3d s = sides();
    return s[0] * s[1] * s[2];
  }

  // mutators
  void setbmin(const Vec3d &bmin_g) { bmin_ = bmin_g; }
  void setbmin(double x, double y, double z) { bmin_ = Vec3d(x, y, z); }
  void setbmax(const Vec3d &bmax_g) { bmax_ = bmax_g; }
  void setbmax(double x, double y, double z) { bmax_ = Vec3d(x, y, z); }
  void setbminmax(const Vec3d &bmin, const Vec3d &bmax)
  {
    bmin_ = bmin;
    bmax_ = bmax;
  }

  double distanceToPoint(const Vec3d &point) const { return std::sqrt(distanceToPoint2(point)); }
  double furthestDistanceToPoint(const Vec3d &point) const { return std::sqrt(furthestDistanceToPoint2(point)); }

  // get squared distance
  // computing squared distance is much faster than computing distance because the later requires std::sqrt()
  double distanceToPoint2(const Vec3d &point) const;
  double furthestDistanceToPoint2(const Vec3d &point) const;

  // return true if point is inside or touching the bounding box
  inline bool checkInside(const Vec3d &point) const;
  // return true if bb is completely inside from inside the bounding box
  inline bool checkInside(const LightBoundingBox &bb) const;
  // return true if bb is completely inside or touching from inside the bounding box
  inline bool checkInsideOrTouching(const LightBoundingBox &bb) const;

  // sanity check bmin <= bmax
  bool verifyBox() const;

  // expands from the center
  // factor of 1.0 indicates no expansion
  inline void expand(double expansionFactor);
  // expand the bounding box to include this point
  inline void expand(const Vec3d &point);
  // expand the bounding box to include another bounding box
  inline void expand(const LightBoundingBox &bb);

  void regularize();  // expand the box into one with all sides equal

  bool lineSegmentIntersection(const Vec3d &segmentStart, const Vec3d &segmentEnd, Vec3d *intersection) const;
  inline bool intersect(const LightBoundingBox &bb) const;

  // return the index of the longest side (bmax[i] - bmin[i]) and its value
  std::pair<int, double> longestSide() const;

  friend inline std::ostream &operator<<(std::ostream &o, const LightBoundingBox &bb)
  {
    return o << "bmin:" << bb.bmin_.transpose() << "\nbmax:" << bb.bmax_.transpose();
  }

  inline void print() const { std::cout << (*this) << std::endl; }

  // get intersection of two bounding boxes
  // will return an invalid bounding box if this and bb does not intersect
  inline LightBoundingBox getIntersection(const LightBoundingBox &bb) const;

  Vec3d bmin_, bmax_;
};

///////////////////////////////////////////////////////////////////////
//                  Below are implementations
///////////////////////////////////////////////////////////////////////

// ================================================
// Implementation for BoundingBox

template<class Vec3dContainer>
BoundingBox::BoundingBox(const Vec3dContainer &vertices)
{
  // set bmin_, bmax_
  bmin_.setConstant(+DBL_MAX);
  bmax_.setConstant(-DBL_MAX);
  for (const Vec3d &p : vertices) {
    updateOnePoint(p);
  }
  updateData();
}

template<class IntContainer>
BoundingBox::BoundingBox(const Vec3d *allVertices, const IntContainer &vertexIDs)
{
  bmin_.setConstant(+DBL_MAX);
  bmax_.setConstant(-DBL_MAX);
  for (int ID : vertexIDs) {
    updateOnePoint(allVertices[ID]);
  }
  updateData();
}

template<class IntContainer>
BoundingBox::BoundingBox(const Vec3d *allVertices, const Vec3i *allTriangles, const IntContainer &triIDs)
{
  bmin_.setConstant(+DBL_MAX);
  bmax_.setConstant(-DBL_MAX);
  for (int triID : triIDs)
    for (int vtxID : allTriangles[triID]) {
      updateOnePoint(allVertices[vtxID]);
    }
  updateData();
}

template<class IntContainer>
BoundingBox::BoundingBox(const Vec3d *allVertices, const Vec4i *allTets, const IntContainer &tetIDs)
{
  bmin_.setConstant(+DBL_MAX);
  bmax_.setConstant(-DBL_MAX);
  for (int tetID : tetIDs)
    for (int vtxID : allTets[tetID]) {
      updateOnePoint(allVertices[vtxID]);
    }
  updateData();
}

template<class IntContainer>
BoundingBox::BoundingBox(const BoundingBox *allBBs, const IntContainer &bbIDs)
{
  bmin_.setConstant(+DBL_MAX);
  bmax_.setConstant(-DBL_MAX);
  for (int ID : bbIDs) {
    updateOneBoundingBox(allBBs[ID]);
  }
  updateData();
}

inline void BoundingBox::updateOnePoint(const Vec3d &p)
{
  for (int i = 0; i < 3; i++) {
    bmin_[i] = (std::min)(bmin_[i], p[i]);
    bmax_[i] = (std::max)(bmax_[i], p[i]);
  }
}

inline void BoundingBox::updateOneBoundingBox(const BoundingBox &bb)
{
  for (int j = 0; j < 3; j++) {
    bmin_[j] = (std::min)(bb.bmin()[j], bmin_[j]);
    bmax_[j] = (std::max)(bb.bmax()[j], bmax_[j]);
  }
}

inline void BoundingBox::updateData()
{
  center_ = 0.5 * (bmin_ + bmax_);
  halfSides_ = 0.5 * (bmax_ - bmin_);
}

// ================================================
// Implementation for LightBoundingBox

template<class Vec3dContainer>
LightBoundingBox::LightBoundingBox(const Vec3dContainer &vertices):
  LightBoundingBox()
{
  for (const Vec3d &p : vertices)
    expand(p);
}

template<class IntContainer>
LightBoundingBox::LightBoundingBox(const Vec3d *allVertices, const IntContainer &vertexIDs):
  LightBoundingBox()
{
  for (int ID : vertexIDs)
    expand(allVertices[ID]);
}

template<class IntContainer>
LightBoundingBox::LightBoundingBox(const Vec3d *allVertices, const Vec3i *allTriangles, const IntContainer &triIDs):
  LightBoundingBox()
{
  for (int triID : triIDs)
    for (int vtxID : allTriangles[triID])
      expand(allVertices[vtxID]);
}

template<class IntContainer>
LightBoundingBox::LightBoundingBox(const Vec3d *allVertices, const Vec4i *allTets, const IntContainer &tetIDs):
  LightBoundingBox()
{
  for (int tetID : tetIDs)
    for (int vtxID : allTets[tetID])
      expand(allVertices[vtxID]);
}

inline bool LightBoundingBox::checkInside(const Vec3d &p) const
{
  return bmin_[0] < p[0] && p[0] < bmax_[0] && bmin_[1] < p[1] && p[1] < bmax_[1] && bmin_[2] < p[2] && p[2] < bmax_[2];
}

inline bool LightBoundingBox::checkInside(const LightBoundingBox &b) const
{
  return bmin_[0] < b.bmin_[0] && bmin_[1] < b.bmin_[1] && bmin_[2] < b.bmin_[2] &&
    bmax_[0] > b.bmax_[0] && bmax_[1] > b.bmax_[1] && bmax_[2] > b.bmax_[2];
}

inline bool LightBoundingBox::checkInsideOrTouching(const LightBoundingBox &b) const
{
  return bmin_[0] <= b.bmin_[0] && bmin_[1] <= b.bmin_[1] && bmin_[2] <= b.bmin_[2] &&
    bmax_[0] >= b.bmax_[0] && bmax_[1] >= b.bmax_[1] && bmax_[2] >= b.bmax_[2];
}

// should this be turned into a self-modifying function?
inline void LightBoundingBox::expand(double expansionFactor)
{
  Vec3d center_ = center();
  Vec3d halfSides_ = halfSides();
  bmin_ = center_ - expansionFactor * halfSides_;
  bmax_ = center_ + expansionFactor * halfSides_;
}

inline void LightBoundingBox::expand(const LightBoundingBox &bb)
{
  for (int j = 0; j < 3; j++) {
    bmin_[j] = (std::min)(bb.bmin()[j], bmin_[j]);
    bmax_[j] = (std::max)(bb.bmax()[j], bmax_[j]);
  }
}

inline bool LightBoundingBox::intersect(const LightBoundingBox &bb) const
{
  return (bmin_[0] < bb.bmax_[0]) && (bmax_[0] > bb.bmin_[0]) &&
    (bmin_[1] < bb.bmax_[1]) && (bmax_[1] > bb.bmin_[1]) &&
    (bmin_[2] < bb.bmax_[2]) && (bmax_[2] > bb.bmin_[2]);
}

inline LightBoundingBox LightBoundingBox::getIntersection(const LightBoundingBox &bb) const
{
  Vec3d newBmin = bmin_, newBmax = bmax_;
  for (int i = 0; i < 3; i++) {
    if (bmin_[i] < bb.bmin_[i])
      newBmin[i] = bb.bmin_[i];
    if (bmax_[i] > bb.bmax_[i])
      newBmax[i] = bb.bmax_[i];
  }
  return LightBoundingBox(newBmin, newBmax);
}

inline void LightBoundingBox::expand(const Vec3d &p)
{
  for (int i = 0; i < 3; i++) {
    bmin_[i] = (std::min)(bmin_[i], p[i]);
    bmax_[i] = (std::max)(bmax_[i], p[i]);
  }
}

}  // namespace Mesh
}  // namespace pgo