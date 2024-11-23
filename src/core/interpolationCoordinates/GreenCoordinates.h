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
#include "triMeshGeo.h"

#include <vector>

namespace pgo
{
namespace InterpolationCoordinates
{

// Green Coordinates, implemented according to
// Yaron Lipman, 2008, Green Coordinates

class GreenCoordinates : public InterpolationCoordinatesBase
{
public:
  // NOTICE: cage must be triangular mesh and have face normals
  GreenCoordinates(int numLocations, const double *locations, const Mesh::TriMeshRef &mesh);

  virtual ~GreenCoordinates() {}

  virtual void deform(const double *verticesDisp, double *locationDisp) const override;

protected:
  double GCTriInt(Vec3d &p, Vec3d &v1, Vec3d &v2, Vec3d &n);

  int numLocations = 0, numVertices = 0, numTriangles = 0;
  std::vector<double> vtxWeights, nmlWeights;
  std::vector<double> triangleArea;
  std::vector<double> vertices;
  std::vector<int> triangles;
};

}
}  // namespace pgo
