/*
author: Bohan Wang
copyright to USC,MIT
*/

#pragma once

#include "vec3d.h"
#include "vec3i.h"

class TetMesh;

namespace pgo
{
namespace SolidDeformationModel
{
class ElementLocalDirectionSolverData;

class ElementLocalDirectionSolver
{
public:
  ElementLocalDirectionSolver();
  ~ElementLocalDirectionSolver();

  void setTetMesh(const TetMesh *tetmesh);
  void setSurfaceMesh(size_t numVertices, const Vec3d *vertices, size_t numTriangles, const Vec3i *triangles, const Vec3d *normals = nullptr);
  void setDirichletBoundaryCoundition(size_t numVertices, const int *vertexIDs, const double *vertexValues);
  void setNeumannBoundaryCoundition(size_t numVertices, const int *vertexIDs, const double *values);

  enum SolverMode
  {
    SM_PURE_PENALITY,
    SM_HARD_LAPLACIAN,
    SM_HARD_CONSTRAINTS
  };
  void solve(SolverMode sm, double WL = 100.0, double WDC = 100.0, double WNC = 1.0);
  const double *getScalarField() const;
  const double *getGradientField() const;

protected:
  ElementLocalDirectionSolverData *data;
};
}  // namespace SolidDeformationModel
}  // namespace pgo
