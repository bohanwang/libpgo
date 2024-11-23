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
class ElementLocalDirectionData;

class ElementLocalDirection
{
public:
  ElementLocalDirection();
  ~ElementLocalDirection();

  void setWorkingDir(const char *filename);
  void addMuscle(size_t numVertices, const Vec3d *vertices, size_t numTriangles,
    const Vec3i *triangles, const Vec3d *normals = nullptr);
  void addWrapperTetMesh(const TetMesh *tetMesh);
  void setOutputScale(double s) { outputScale = s; }

  void compute();
  void computeNonEmbedding(size_t numZeroVertices, const int *zeroVertices,
    size_t numOneVertices, const int *oneVertices, const double weights[3] = nullptr);

  Vec3d getFiberDirection(int elementID) const;
  Vec3d getVertexFiberDirection(int vertexID) const;
  size_t getMuscleClosedMeshNumVertices(int mi) const;
  const Vec3d *getMuscleClosedMeshVertices(int mi) const;
  size_t getMuscleClosedMeshNumTriangles(int mi) const;
  const Vec3i *getMuscleClosedMeshTriangles(int mi) const;

  void save(const char *filenamePrefix) const;
  void saveSurfaceVertices(const char *filenamePrefix) const;
  int load(const char *filenamePrefix) const;

protected:
  ElementLocalDirectionData *data;
  double outputScale = 1000.0;
};
}  // namespace SolidDeformationModel
}  // namespace pgo
