#pragma once

#include "pgo_c.h"

struct TetMeshGeo
{
  TetMeshGeo():
    handle(nullptr) {}

  TetMeshGeo(pgoTetMeshGeoStructHandle h):
    handle(h) {}

  pgoTetMeshGeoStructHandle handle;
};

struct TriMeshGeo
{
  TriMeshGeo():
    handle(nullptr) {}

  TriMeshGeo(pgoTriMeshGeoStructHandle h):
    handle(h) {}

  pgoTriMeshGeoStructHandle handle;
};

struct TetMesh
{
  TetMesh():
    handle(nullptr) {}

  TetMesh(pgoTetMeshStructHandle h):
    handle(h) {}

  pgoTetMeshStructHandle handle;
};

struct SmoothRSEnergy
{
  SmoothRSEnergy():
    handle(nullptr) {}

  SmoothRSEnergy(pgoSmoothRSEnergyStructHandle h):
    handle(h) {}

  pgoSmoothRSEnergyStructHandle handle;
};

struct SparseMatrix
{
  SparseMatrix():
    handle(nullptr) {}

  SparseMatrix(pgoSparseMatrixStructHandle h):
    handle(h) {}

  pgoSparseMatrixStructHandle handle;
};

class pypgoInit
{
public:
  pypgoInit()
  {
    pgo_init();
    // code initialization
  }
  ~pypgoInit()
  {
    // finalize
  }
};