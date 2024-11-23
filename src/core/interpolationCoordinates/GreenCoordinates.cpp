/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "interpolationCoordinates" library , Copyright (C) 2018 USC           *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code authors: Yijing Li, Jernej Barbic                                *
 * http://www.jernejbarbic.com/vega                                      *
 *                                                                       *
 * Funding: National Science Foundation                                  *
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

#include "GreenCoordinates.h"
#include "triMeshGeo.h"
#include "geometryQuery.h"
#include "triMeshPseudoNormal.h"
#include "meshLinearAlgebra.h"

#include <cfloat>
#include <climits>
#include <cstring>
#include <memory>

#include <tbb/parallel_for.h>

using namespace pgo;
using namespace pgo::InterpolationCoordinates;

namespace ES = pgo::EigenSupport;

inline int sign(double v)
{
  return (v >= 0) ? 1 : -1;
}

GreenCoordinates::GreenCoordinates(int numLocations, const double *locations, const Mesh::TriMeshRef &cage)
{
  this->numLocations = numLocations;

  // cage->buildFaceNormals();
  numVertices = cage.numVertices();
  numTriangles = cage.numTriangles();

  vtxWeights.resize(numLocations * numVertices);
  nmlWeights.resize(numLocations * numTriangles);
  Vec3d zero = asVec3d(0.);

  for (int fi = 0; fi < cage.numTriangles(); fi++) {
    triangleArea.emplace_back(Mesh::getTriangleArea(cage.pos(fi, 0), cage.pos(fi, 1), cage.pos(fi, 2)));
    triangles.emplace_back(cage.triVtxID(fi, 0));
    triangles.emplace_back(cage.triVtxID(fi, 1));
    triangles.emplace_back(cage.triVtxID(fi, 2));
  }

  vertices.resize(3 * numVertices);
  for (int i = 0; i < cage.numVertices(); i++) {
    ES::Mp<ES::V3d>(vertices.data() + i * 3) = cage.pos(i);
  }

  Mesh::TriMeshPseudoNormal meshNormal;
  meshNormal.buildPseudoNormals(cage);

  tbb::parallel_for(tbb::blocked_range<int>(0, numLocations), [&](const tbb::blocked_range<int> &rng) {
    for (int i = rng.begin(); i != rng.end(); ++i) {
      int globalFaceCount = 0;
      Vec3d n(locations + 3 * i);

      for (int fi = 0; fi < cage.numTriangles(); fi++) {
        Vec3d normal = meshNormal.triNormal(fi);
        Vec3d v[3];
          int s[3] = { 1, 1, 1 };
          double I[3], II[3];
          Vec3d q[3];
          Vec3d N[3];
          
          for (int k = 0; k < 3; k++) {
            v[k] = cage.pos(fi, k);
            v[k] = v[k] - n;
          }
          Vec3d p = v[0].dot(normal) * normal;

          for (int k = 0; k < 3; k++) {
            int next = ((k == 2) ? 0 : (k + 1));
            
            s[k] = sign((v[k] - p).cross(v[next] - p).dot(normal));
            I[k] = GCTriInt(p, v[k], v[next], zero);
            II[k] = GCTriInt(zero, v[next], v[k], zero);
            q[k] = v[next].cross(v[k]);
            N[k] = ((q[k].squaredNorm() < std::numeric_limits<double>::epsilon()) ? zero : q[k].normalized());
          }

          double Iall = (-1) * std::abs(s[0] * I[0] + s[1] * I[1] + s[2] * I[2]);
          nmlWeights[i * numTriangles + globalFaceCount] = (-1) * Iall;
          Vec3d w = normal * Iall + N[0] * II[0] + N[1] * II[1] + N[2] * II[2];

          if ((w).squaredNorm() > std::numeric_limits<double>::epsilon()) {
            for (int k = 0; k < 3; k++) {
              int next = ((k == 2) ? 0 : (k + 1));

              int pidx = cage.triVtxID(fi, k);
              vtxWeights[i * numVertices + pidx] += ((N[next].dot( v[k]) == 0) ? 0 : N[next].dot( w) / N[next].dot( v[k]));
            }
          }

          globalFaceCount++;
      }
    } }, tbb::static_partitioner());
}

void GreenCoordinates::deform(const double *verticesDisp, double *locationDisp) const
{
  tbb::parallel_for(tbb::blocked_range<int>(0, numLocations), [&](const tbb::blocked_range<int> &rng) {
    for (int i = rng.begin(); i != rng.end(); ++i) {
      Vec3d newDisp = asVec3d(0.);

      // add vtx weights
      for (int j = 0; j < numVertices; j++) {
        Vec3d v(verticesDisp + 3 * j);
        newDisp += v * vtxWeights[i * numVertices + j];
      }

      // add normal weights
      for (int j = 0; j < numTriangles; j++) {
        if (triangleArea[j] == 0)
          continue;
        Vec3d v[3], u[3];
        for (unsigned int k = 0; k < 3; k++) {
          int idx = triangles[3 * j + k];
          u[k] = asVec3d(&vertices[3 * idx]);
          v[k] = u[k] + Vec3d(verticesDisp + 3 * idx);
        }
        Vec3d oldVec[2] = { u[1] - u[0], u[2] - u[0] };
        Vec3d newVec[2] = { v[1] - v[0], v[2] - v[0] };

        Vec3d newNormal = newVec[0].cross(newVec[1]);

        if (newNormal.squaredNorm() > 0)
          newNormal.normalize();

        double s = sqrt(newVec[0].squaredNorm() * oldVec[1].squaredNorm() - 2 * newVec[0].dot(newVec[1]) *oldVec[0].dot( oldVec[1]) + newVec[1].squaredNorm() * oldVec[0].squaredNorm()) / (sqrt(8) * triangleArea[j]);

        newDisp += newNormal * s * nmlWeights[i * numTriangles + j];
      }
      memcpy(locationDisp + 3 * i, &newDisp[0], sizeof(double) * 3);
    } }, tbb::static_partitioner());  // end for locations
}

double GreenCoordinates::GCTriInt(Vec3d &p, Vec3d &v1, Vec3d &v2, Vec3d &n)
{
  double alpha = Mesh::getTriangleAngleRobust(v1, v2, p);
  double beta = Mesh::getTriangleAngleRobust(p, v1, v2);
  double sqrtlambda = std::abs((p - v1).norm() * sin(alpha));
  double lambda = sqrtlambda * sqrtlambda;
  double sqrtc = (p - n).norm();
  double c = sqrtc * sqrtc;
  double t[2] = { M_PI - alpha, M_PI - alpha - beta };
  double I[2];
  for (int i = 0; i < 2; i++) {
    double S = sin(t[i]);
    double C = cos(t[i]);
    I[i] = (-0.5) * sign(S) * (2 * sqrtc * atan(sqrtc * C / sqrt(lambda + S * S * c)) + sqrtlambda * log(2 * sqrtlambda * S * S * (1 - 2 * c * C / (c * (1 + C) + lambda + sqrt(lambda * lambda + lambda * c * S * S))) / ((1 - C) * (1 - C))));
  }
  return std::abs(I[0] - I[1] - sqrtc * beta) * (-1) / (4 * M_PI);
}
