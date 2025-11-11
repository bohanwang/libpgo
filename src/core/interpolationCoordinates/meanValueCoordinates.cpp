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

#include "meanValueCoordinates.h"
#include "basicAlgorithms.h"
#include "EigenSupport.h"
#include "meshLinearAlgebra.h"

#include <cfloat>
#include <climits>
#include <memory>
#include <cstring>

#include <tbb/parallel_for.h>

using namespace pgo;
using namespace pgo::InterpolationCoordinates;

namespace ES = pgo::EigenSupport;

MeanValueCoordinates::MeanValueCoordinates(int numLocations_, const double *locations_, int numCageVertices_, const double *cageVertices_,
  int numCageTriangles_, const int *cageTriangles_)
{
  this->numLocations = numLocations_;
  this->numCageVertices = numCageVertices_;
  this->numCageTriangles = numCageTriangles_;

  // compute weights
  double epsilon = 1e-10;
  weights.resize(numCageVertices_ * numLocations_);

  tbb::parallel_for(tbb::blocked_range<int>(0, numLocations), [&](const tbb::blocked_range<int> &rng) {
    std::vector<ES::V3d> u(numCageVertices);
    std::vector<double> d(numCageVertices);
    for (int i = rng.begin(); i != rng.end(); ++i) {
      bool lieInTriangle = false;
      ES::V3d x = asVec3d(locations_ + 3 * i);
      bool closeToVertex = false;
      for (int j = 0; j < numCageVertices_; j++) {
        ES::V3d p = asVec3d(cageVertices_ + 3 * j);
        d[j] = (p - x).norm();
        if (d[j] < epsilon) {
          weights[i * numCageVertices + j] = 1.0;
          closeToVertex = true;
          break;
        }
        // cout << "p = " << p << endl;
        u[j] = (p - x) / d[j];
      }
      if (closeToVertex)
        continue;  // simple case, done

      double W = 0;
      for (int j = 0; j < numCageTriangles; j++) {
        double l[3] = { 0, 0, 0 };
        const int *pi = cageTriangles_ + 3 * j;
        double theta[3];
        for (int k = 0; k < 3; k++) {
          // int
          l[k] = (u[pi[(k + 1) % 3]] - u[pi[(k + 2) % 3]]).norm();
          theta[k] = 2 * std::asin(std::clamp(l[k] / 2, -1.0, 1.0));
        }
        // cout << theta[0] << " " << theta[1] << " " << theta[2] << endl;

        double h = (theta[0] + theta[1] + theta[2]) / 2.;
        if ((M_PI - h) < epsilon) {
          // x lies on t, use 2D barycentric coord.
          double w[3];
          W = 0;
          for (int k = 0; k < 3; k++) {
            w[k] = sin(theta[k]) * d[pi[(k + 2) % 3]] * d[pi[(k + 1) % 3]];
            W += w[k];
          }
          lieInTriangle = true;
          memset(&weights[i * numCageVertices], 0, sizeof(double) * numCageVertices);
          for (int k = 0; k < 3; k++) {
            weights[i * numCageVertices + pi[k]] = w[k] / W;
          }
          break;
        }
        double c[3];
        for (int k = 0; k < 3; k++) {
          c[k] = (2 * sin(h) * sin(h - theta[k])) / (sin(theta[(k + 1) % 3]) * sin(theta[(k + 2) % 3])) - 1;
        }
        // cout << c[0] << " " << c[1] << " " << c[2] << endl;

        ES::M3d m = asMat3d(u[pi[0]], u[pi[1]], u[pi[2]]);
        int sdm = BasicAlgorithms::signAsInt(m.determinant());
        // PRINT(dm);
        double s[3];
        bool onSamePlane = false;
        for (int k = 0; k < 3; k++) {
          s[k] = sdm * BasicAlgorithms::sqrtSafe(1 - c[k] * c[k]);
          if (fabs(s[k]) <= epsilon) {
            // x lies outside t on the same plane, ignore t
            onSamePlane = true;
            break;
          }
        }
        if (onSamePlane)
          continue;

        for (int k = 0; k < 3; k++) {
          double w = (theta[k] - c[(k + 1) % 3] * theta[(k + 2) % 3] - c[(k + 2) % 3] * theta[(k + 1) % 3]) / (d[pi[k]] * sin(theta[(k + 1) % 3]) * s[(k + 2) % 3]);
          weights[i * numCageVertices + pi[k]] += w;
          W += w;
        }
      }  // end for triangle

      if (lieInTriangle)
        continue;

      for (int j = 0; j < numCageVertices; j++) {
        weights[i * numCageVertices + j] /= W;
      }
    } }, tbb::static_partitioner());  // end for locations
}

void MeanValueCoordinates::deform(const double *cageDisp, double *locationDisp) const
{
  memset(locationDisp, 0, sizeof(double) * 3 * numLocations);

  tbb::parallel_for(tbb::blocked_range<int>(0, numLocations), [&](const tbb::blocked_range<int> &rng) {
    for (int i = rng.begin(); i != rng.end(); ++i) {
      double *l = locationDisp + 3 * i;
      for (int j = 0; j < numCageVertices; j++) {
        double w = weights[i * numCageVertices + j];
        for (int k = 0; k < 3; k++)
          l[k] += w * cageDisp[3 * j + k];
      }
    } }, tbb::static_partitioner());
}
