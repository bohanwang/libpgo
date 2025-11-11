/*
author: Bohan Wang
copyright to USC,MIT,NUS
*/

#pragma once

#include "meshLinearAlgebra.h"

namespace pgo
{
namespace SolidDeformationModel
{
void computeTetMeshOccupation(int numTetVertices, const Vec3d *tetVertices, int numTets, const Vec4i *tets,
  int numSurfaceVertices, const Vec3d *surfaceVertices, int numTriangles, const Vec3i *triangles,
  double minThreshold, int sampleCount, double *weights);

void computeTetMeshOccupation(int numTetVertices, const Vec3d *tetVertices, int numTets, const Vec4i *tets,
  const Vec3d bb[2], double minThreshold, int sampleCount, double *weights);
}  // namespace SolidDeformationModel
}  // namespace pgo