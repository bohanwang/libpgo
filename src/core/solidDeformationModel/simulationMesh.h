/*
author: Bohan Wang
copyright to USC,MIT,NUS
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
  SimulationMeshENuMaterial(double E_, double nu_, double J_ = 10000):
    E(E_), nu(nu_), J(J_) {}
  virtual ~SimulationMeshENuMaterial() {}

  void setE(double E_) { E = E_; }
  void setNu(double nu_) { nu = nu_; }
  void setJ(double J_) { J = J_; }

  double getMuLame() const { return E / (2 * (1 + nu)); }
  double getLambdaLame() const { return (nu * E) / ((1 + nu) * (1 - 2 * nu)); }

  double getE() const { return E; }
  double getNu() const { return nu; }
  double getCompressionRatio() const { return J; }

  virtual SimulationMeshMaterial *clone() const override { return new SimulationMeshENuMaterial(E, nu, J); }

protected:
  double E = 6e3, nu = 0.4;
  double J = 10000;
};

class SimulationMeshENuhMaterial : public SimulationMeshENuMaterial
{
public:
  SimulationMeshENuhMaterial() {}
  SimulationMeshENuhMaterial(double E_, double nu_, double h_, double J_ = 10000):
    SimulationMeshENuMaterial(E_, nu_, J_), h(h_) {}

  virtual ~SimulationMeshENuhMaterial() {}

  void seth(double h_) { h = h_; }
  double geth() const { return h; }

  virtual SimulationMeshMaterial *clone() const override { return new SimulationMeshENuhMaterial(E, nu, h, J); }

protected:
  double h = 1e-4;
};

class SimulationMeshHillMaterial : public SimulationMeshMaterial
{
public:
  SimulationMeshHillMaterial() {}
  SimulationMeshHillMaterial(double E_act_, double gamma_, double lo_):
    E_act(E_act_), gamma(gamma_), lo(lo_) {}
  virtual ~SimulationMeshHillMaterial() {}

  void setEact(double E_) { E_act = E_; }
  void setGamma(double gamma_) { gamma = gamma_; }
  void setLo(double lo_) { lo = lo_; }

  double getEact() const { return E_act; }
  double getGamma() const { return gamma; }
  double getLo() const { return lo; }

  virtual SimulationMeshMaterial *clone() const override { return new SimulationMeshHillMaterial(E_act, gamma, lo); }

protected:
  double E_act = 0.1e6, gamma = 1, lo = 0.6;
};

class SimulationMeshMooneyRivlinMaterial : public SimulationMeshMaterial
{
public:
  SimulationMeshMooneyRivlinMaterial(int N, int M, const double *Cpq, const double *D_)
  {
    C = new double[(N + 1) * (N + 1)];
    D = new double[M];

    this->N = N;
    this->M = M;

    for (int p = 0; p <= N; p++) {
      for (int q = 0; q <= N; q++) {
        getC(p, q) = Cpq[q * (N + 1) + p];
      }
    }

    for (int i = 0; i < M; i++) {
      getD(i) = D_[i];
    }
  }

  SimulationMeshMooneyRivlinMaterial(int N, int M, double E, double nu)
  {
    C = new double[(N + 1) * (N + 1)];
    D = new double[M];

    this->N = N;
    this->M = M;

    for (int p = 0; p <= N; p++) {
      for (int q = 0; q <= N; q++) {
        getC(p, q) = 0;
      }
    }

    for (int i = 0; i < M; i++) {
      getD(i) = 0;
    }

    double bulkModulus = E / (3 * (1 - 2 * nu));
    double shearModulus = E / (2 * (1 + nu));

    getC(1, 0) = shearModulus * 0.5;
    getD(0) = 2.0 / bulkModulus;
  }

  virtual ~SimulationMeshMooneyRivlinMaterial()
  {
    delete[] C;
    delete[] D;
  }

  void setC(int p, int q, double val)
  {
    getC(p, q) = val;
  }

  void setD(int i, double val)
  {
    getD(i) = val;
  }

  double &getC(int p, int q)
  {
    return C[q * (N + 1) + p];
  }

  const double &getC(int p, int q) const
  {
    return C[q * (N + 1) + p];
  }

  double &getD(int i)
  {
    return D[i];
  }

  const double &getD(int i) const
  {
    return D[i];
  }

  const double *getC() const { return C; }
  const double *getD() const { return D; }

  int getM() const { return M; }
  int getN() const { return N; }

  virtual SimulationMeshMaterial *clone() const override
  {
    return new SimulationMeshMooneyRivlinMaterial(N, M, C, D);
  }

protected:
  int M, N;
  double *C, *D;
};

class SimulationMeshMooneyRivlinhMaterial : public SimulationMeshMooneyRivlinMaterial
{
public:
  SimulationMeshMooneyRivlinhMaterial(int N_, int M_, const double *Cpq_, const double *D_, double h_):
    SimulationMeshMooneyRivlinMaterial(N_, M_, Cpq_, D_), h(h_) {}
  SimulationMeshMooneyRivlinhMaterial(int N_, int M_, double E_, double nu_, double h_):
    SimulationMeshMooneyRivlinMaterial(N_, M_, E_, nu_), h(h_) {}
  virtual ~SimulationMeshMooneyRivlinhMaterial() {}

  void seth(double h_) { h = h_; }
  double geth() const { return h; }

  virtual SimulationMeshMaterial *clone() const override { return new SimulationMeshMooneyRivlinhMaterial(N, M, C, D, h); }

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

  const SimulationMeshMaterial *getElementMaterial(int ele, int j) const;
  SimulationMeshMaterial *getElementMaterial(int ele, int j);
  int getElementNumMaterials(int ele) const;

  void setMaterial(int matID, const SimulationMeshMaterial *mat);

protected:
  SimulationMeshImpl *impl;
};

SimulationMesh *loadTetMesh(const VolumetricMeshes::TetMesh *tetmesh);

SimulationMesh *loadTriMesh(const Mesh::TriMeshGeo &triMeshGeo, const SimulationMeshMaterial *mat, int toTriangle);
SimulationMesh *loadTriMesh(const Mesh::TriMeshGeo &triMeshGeo, int numMaterials, const SimulationMeshMaterial *const *const mat, const int *materialIndices, int toTriangle);

SimulationMesh *loadShellMesh(const Mesh::TriMeshGeo &triMeshGeo, const SimulationMeshMaterial *mat);
SimulationMesh *loadShellMesh(const Mesh::TriMeshGeo &triMeshGeo, const int *elementMaterialIndices, const SimulationMeshMaterial *const *mat);

void computeTriangleUV(SimulationMesh *mesh, double scaleFactor);
}  // namespace SolidDeformationModel
}  // namespace pgo
