/*
author: Bohan Wang
copyright to USC,MIT
*/

#pragma once

namespace pgo
{

namespace Mesh
{
class TriMeshGeo;
class TetMeshGeo;
}  // namespace Mesh

namespace VolumetricMeshes
{
class TetMesh;
}

namespace SolidDeformationModel
{
class SimulationMeshMaterial
{
public:
  SimulationMeshMaterial() {}
  virtual ~SimulationMeshMaterial() {}

  virtual SimulationMeshMaterial *clone() const = 0;
};

class SimulationMeshENuMaterial : public SimulationMeshMaterial
{
public:
  SimulationMeshENuMaterial() {}
  SimulationMeshENuMaterial(double E_, double nu_):
    E(E_), nu(nu_) {}
  virtual ~SimulationMeshENuMaterial() {}

  void setE(double E_) { E = E_; }
  void setNu(double nu_) { nu = nu_; }

  double getMuLame() const { return E / (2 * (1 + nu)); }
  double getLambdaLame() const { return (nu * E) / ((1 + nu) * (1 - 2 * nu)); }

  double getE() const { return E; }
  double getNu() const { return nu; }

  virtual SimulationMeshMaterial *clone() const override { return new SimulationMeshENuMaterial(E, nu); }

protected:
  double E = 6e3, nu = 0.4;
};

class SimulationMeshENuhMaterial : public SimulationMeshENuMaterial
{
public:
  SimulationMeshENuhMaterial() {}
  SimulationMeshENuhMaterial(double E_, double nu_, double h_):
    SimulationMeshENuMaterial(E_, nu_), h(h_) {}

  virtual ~SimulationMeshENuhMaterial() {}

  void seth(double h_) { h = h_; }
  double geth() const { return h; }

  virtual SimulationMeshMaterial *clone() const override { return new SimulationMeshENuhMaterial(E, nu, h); }

protected:
  double h = 1e-4;
};

enum class SimulationMeshType
{
  TET,
  CUBIC,
  TRIANGLE,
  EDGE_QUAD,
  SHELL,
};

class SimulationMeshImpl;

class SimulationMesh
{
public:
  SimulationMesh(int numVertices, const double *vertexPositions,
    int numElements, int numElementVertices, const int *elementVertexIndices,
    const int *elementMaterialIndices, int numMaterials, const SimulationMeshMaterial *const *materials,
    SimulationMeshType meshType);
  ~SimulationMesh();

  int getNumElements() const;
  int getNumVertices() const;
  int getNumElementVertices() const;
  int getVertexIndex(int ele, int j) const;
  const int *getVertexIndices(int ele) const;

  void getVertex(int vi, double pos[3]) const;
  void getVertex(int ele, int j, double pos[3]) const;

  void assignElementUVs(const double *uvs);
  bool hasElementUV() const;
  void getElementUV(int ele, int j, double uv[2]) const;

  SimulationMeshType getElementType() const;

  const SimulationMeshMaterial *getElementMaterial(int ele) const;

protected:
  SimulationMeshImpl *impl;
};

SimulationMesh *loadTetMesh(const VolumetricMeshes::TetMesh *tetmesh);
SimulationMesh *loadTriMesh(const Mesh::TriMeshGeo &triMeshGeo, const SimulationMeshMaterial *mat, int toTriangle);
SimulationMesh *loadTriMesh(const Mesh::TriMeshGeo &triMeshGeo, int numMaterials, const SimulationMeshMaterial *const *const mat, const int *materialIndices, int toTriangle);
SimulationMesh *loadShellMesh(const Mesh::TriMeshGeo &triMeshGeo, const SimulationMeshMaterial *mat);

void computeTriangleUV(SimulationMesh *mesh, double scaleFactor);
}  // namespace SolidDeformationModel
}  // namespace pgo
