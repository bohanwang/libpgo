/*
author: Bohan Wang
copyright to USC, MIT
*/

#pragma once

#include "potentialEnergy.h"
#include "EigenSupport.h"
#include "tetMeshGeo.h"

#include <memory>

class TetMesh;

namespace pgo
{
namespace PredefinedPotentialEnergies
{
struct SmoothRSEnergyBuffer;

class SmoothRSEnergy : public NonlinearOptimization::PotentialEnergy
{
public:
  SmoothRSEnergy(const Mesh::TetMeshGeo &tetMesh, const EigenSupport::SpMatD &G, const EigenSupport::SpMatD &L, const double coeffs[2], int isPosition = 1, int dofOffsets = 0);
  virtual ~SmoothRSEnergy() {}

  virtual double func(EigenSupport::ConstRefVecXd x) const;
  virtual void gradient(EigenSupport::ConstRefVecXd x, EigenSupport::RefVecXd grad) const;
  virtual void hessian(EigenSupport::ConstRefVecXd x, EigenSupport::SpMatD &hess) const;

  virtual void createHessian(EigenSupport::SpMatD &hess) const { hess = hessTemplate; }

  virtual void getDOFs(std::vector<int> &dofs) const { dofs = allDOFs; }
  virtual int getNumDOFs() const { return (int)allDOFs.size(); }

  static void convertRowMajorFG(const EigenSupport::SpMatD &GRowMajor, EigenSupport::SpMatD &GColMajor);

protected:
  void getu(const EigenSupport::ConstRefVecXd x, EigenSupport::VXd &u) const;

  EigenSupport::M3d getMat(const EigenSupport::VXd &matVec, int ei) const;
  EigenSupport::V9d packMat(const EigenSupport::M3d &m) const;

  void computeFSR(const EigenSupport::VXd &u, EigenSupport::VXd &F, EigenSupport::VXd &S, EigenSupport::VXd &R) const;

  template<int blockDim, typename MatrixType>
  void computeBlockSparsePattern(const EigenSupport::SpMatD &LTL, std::unordered_map<std::pair<int, int>, MatrixType, EigenSupport::IntPairHash, EigenSupport::IntPairEqual> &blocks) const;

  const Mesh::TetMeshGeo &tetMesh;
  EigenSupport::SpMatD LTL9;
  EigenSupport::SpMatD G;
  EigenSupport::VXd restPositions;

  double coeffR, coeffS;
  int inputIsPosition;
  int dofOffsets;

  EigenSupport::SpMatD hessTemplate, dF;
  std::vector<EigenSupport::M9i> dFblocks;
  EigenSupport::EntryMap hessMap;
  std::vector<int> allDOFs;

  EigenSupport::SymbolicMmData *dFGMapping, *LTLdFGMapping, *h1Mapping, *h2Mapping;
  EigenSupport::SpMatI h1ToHess, h2ToHess;

  std::shared_ptr<SmoothRSEnergyBuffer> evaluationBuf;
};
}  // namespace PredefinedPotentialEnergies
}  // namespace pgo
