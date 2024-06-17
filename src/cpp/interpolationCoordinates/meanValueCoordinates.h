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

#pragma once

#include "interpolationCoordinatesBase.h"
#include <vector>

namespace pgo
{
namespace InterpolationCoordinates
{

// implemnted mean-value coordinates according to
// Tao Ju, 2005, Mean Value Coordinates for Closed Triangular Meshes
class MeanValueCoordinates : public InterpolationCoordinatesBase
{
public:
  MeanValueCoordinates() {}
  MeanValueCoordinates(int numLocations, const double *locations, int numCageVertices, const double *cageVertices,
    int numCageTriangles, const int *cageTriangles);
  virtual ~MeanValueCoordinates() {}

  virtual void deform(const double *verticesDisp, double *locationDisp) const override;

  const double *getWeights(int embeddedVtx) const { return &weights[embeddedVtx * numCageVertices]; }

protected:
  // std::vector<double> locations, vertices;
  // std::vector<int> triangles;
  std::vector<double> weights;
  int numLocations = 0, numCageVertices = 0, numCageTriangles = 0;
};

}  // namespace InterpolationCoordinates
}  // namespace pgo
