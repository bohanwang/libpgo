/*
  This code is based on code from the Geometric Tools library,
  which is licensed under a boost license.
  Such usage is permitted by the boost license; for details,
  please see the boost license below.
*/

// Geometric Tools, LLC
// Copyright (c) 1998-2014
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt

/*************************************************************************
 *                                                                       *
 * We release our improvements to the wildMagic code under our standard  *
 * Vega FEM license, as follows:                                         *
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "improvements to the wildMagic library" , Copyright (C) 2018 USC      *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code author: Yijing Li                                                *
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

#include "verticesInfo.h"

#include "pgoLogging.h"

#include <cstring>
#include <iostream>
#include <cmath>

pgo::Mesh::VerticesInformation pgo::Mesh::getVerticesInformation(const std::vector<Vec3d> &vertices, double epsilon)
{
  VerticesInformation info;
  assert(epsilon >= 0 && vertices.size() > 0);
  info.extremeCCW = false;

  // Compute the axis-aligned bounding box for the vertices.  Keep track
  // of the indices in the vertices for the current min and max.
  int indexMin[3] = { 0, 0, 0 }, indexMax[3] = { 0, 0, 0 };

  (EigenSupport::Mp<Vec3d>(info.min)) = vertices[0];
  (EigenSupport::Mp<Vec3d>(info.max)) = vertices[0];

  for (size_t i = 1; i < vertices.size(); ++i) {
    for (int j = 0; j < 3; ++j) {
      if (vertices[i][j] < info.min[j]) {
        info.min[j] = vertices[i][j];
        indexMin[j] = (int)i;
      }
      else if (vertices[i][j] > info.max[j]) {
        info.max[j] = vertices[i][j];
        indexMax[j] = (int)i;
      }
    }
  }

  // Determine the maximum range for the bounding box.
  info.maxRange = info.max[0] - info.min[0];
  info.extreme[0] = indexMin[0];
  info.extreme[1] = indexMax[0];
  double range = info.max[1] - info.min[1];
  if (range > info.maxRange) {
    info.maxRange = range;
    info.extreme[0] = indexMin[1];
    info.extreme[1] = indexMax[1];
  }
  range = info.max[2] - info.min[2];
  if (range > info.maxRange) {
    info.maxRange = range;
    info.extreme[0] = indexMin[2];
    info.extreme[1] = indexMax[2];
  }

  // The origin is either the point of minimum x-value, point of
  // minimum y-value, or point of minimum z-value.
  info.origin = vertices[info.extreme[0]];

  // Test whether the point set is (nearly) a point.
  if (info.maxRange < epsilon) {
    info.dimension = 0;
    for (int j = 0; j < 3; ++j) {
      info.extreme[j + 1] = info.extreme[0];
      info.direction[j].setZero();
    }
    return info;
  }

  // Test whether the point set is (nearly) a line segment.
  info.direction[0] = vertices[info.extreme[1]] - info.origin;
  info.direction[0].normalize();
  double maxDistance = 0.;
  info.extreme[2] = info.extreme[0];
  for (size_t i = 0; i < vertices.size(); ++i) {
    Vec3d diff = vertices[i] - info.origin;
    double dotProduct = info.direction[0].dot(diff);
    Vec3d proj = diff - dotProduct * info.direction[0];
    double distance = proj.norm();
    if (distance > maxDistance) {
      maxDistance = distance;
      info.extreme[2] = (int)i;
    }
  }

  if (maxDistance < epsilon * info.maxRange) {
    info.dimension = 1;
    info.extreme[2] = info.extreme[1];
    info.extreme[3] = info.extreme[1];
    return info;
  }

  // Test whether the point set is (nearly) a planar polygon.
  info.direction[1] = vertices[info.extreme[2]] - info.origin;
  info.direction[1] -= info.direction[0].dot(info.direction[1]) * info.direction[0];
  info.direction[1].normalize();
  info.direction[2] = info.direction[0].cross(info.direction[1]);
  maxDistance = 0;
  double maxSign = 0;
  info.extreme[3] = info.extreme[0];
  for (size_t i = 0; i < vertices.size(); ++i) {
    Vec3d diff = vertices[i] - info.origin;
    double distance = info.direction[2].dot(diff);
    int sign = distance >= 0 ? 1 : -1;
    distance = std::abs(distance);
    if (distance > maxDistance) {
      maxDistance = distance;
      maxSign = sign;
      info.extreme[3] = (int)i;
    }
  }

  if (maxDistance < epsilon * info.maxRange) {
    info.dimension = 2;
    info.extreme[3] = info.extreme[2];
    return info;
  }

  info.dimension = 3;
  assert(epsilon == 0 || maxSign != 0);
  info.extremeCCW = (maxSign >= 0);

  return info;
}
