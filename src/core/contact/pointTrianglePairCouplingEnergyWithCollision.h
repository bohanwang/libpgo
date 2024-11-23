#pragma once

#include "contactEnergyUtilities.h"

#include "potentialEnergy.h"
#include "triMeshGeo.h"
#include "triMeshNeighbor.h"

#include <array>
#include <vector>
#include <functional>
#include <unordered_map>
#include <tuple>

namespace pgo
{
namespace Contact
{
class PointTrianglePairCouplingEnergyWithCollisionBuffer;

class PointTrianglePairCouplingEnergyWithCollision : public NonlinearOptimization::PotentialEnergy
{
public:
  PointTrianglePairCouplingEnergyWithCollision(int numPairs, int numObjects,
    const std::array<int, 4> *const objectIDs,
    const std::array<int, 4> *const pointTrianglePairs,
    const int *const triangleIDs,
    const int *const objectDOFOffsets,
    const Mesh::TriMeshRef *const surfaceMeshes,
    const Mesh::TriMeshNeighbor *const surfaceMeshNeighbors,
    const std::vector<int> *const vertexIDToSampleIDs,
    const EigenSupport::SpMatD *const sampleEmbeddingWeights,
    const std::vector<double> *const sampleWeights);

  virtual ~PointTrianglePairCouplingEnergyWithCollision() {}

  PointTrianglePairCouplingEnergyWithCollisionBuffer *allocateBuffer() const;
  void freeBuffer(PointTrianglePairCouplingEnergyWithCollisionBuffer *buf) const;
  void setBuffer(PointTrianglePairCouplingEnergyWithCollisionBuffer *b) { buf = b; }

  void setToPosFunction(PosFunction f) { toPosFunc = f; }
  void setToLastPosFunction(PosFunction f) { toLastPosFunc = f; }

  void setCoeff(double c) { coeff = c; }
  void setFrictionCoeff(double c) { frictionCoeff = c; }
  void setTimestep(double t) { timestep = t; }
  void setVelEps(double eps_) { eps = eps_; }

  virtual double func(EigenSupport::ConstRefVecXd x) const override;
  virtual void gradient(EigenSupport::ConstRefVecXd x, EigenSupport::RefVecXd grad) const override;
  virtual void hessian(EigenSupport::ConstRefVecXd x, EigenSupport::SpMatD &hess) const override;

  virtual void createHessian(EigenSupport::SpMatD &hess_) const override { hess_ = hessTemplate; }

  virtual void getDOFs(std::vector<int> &dofs) const override { dofs = allDOFs; }
  virtual int getNumDOFs() const override { return (int)allDOFs.size(); }

  virtual int isQuadratic() const { return 0; }
  virtual int hasDenseHessian() const { return 0; }
  virtual int hasHessianVector() const { return 0; }

  void computeClosestPosition(const double *const x);

  const EigenSupport::VXd &getFrictionForce() const { return ff; };
  const EigenSupport::VXd &getNonFrictionForce() const { return fn; };

  const double &getFrictionEnergy() const { return eng_ff; }
  const double &getNonFrictionEnergy() const { return eng_fn; }

protected:
  std::tuple<bool, int> checkContact(EigenSupport::ConstRefVecXd x, int pi, EigenSupport::V12d &xlocal, int checkNeighboringTriangles) const;

  EigenSupport::V3d computePosition(EigenSupport::ConstRefVecXd x, int objID, int vid) const;
  EigenSupport::V3d computeLastPosition(EigenSupport::ConstRefVecXd x, int objID, int vid) const;

  int getTriangleObjectID(int pi) const { return objectIDs[pi][1]; }
  int getPointObjectID(int pi) const { return objectIDs[pi][0]; }

  double frictionPotential(const EigenSupport::V12d &x, const EigenSupport::V12d &xlast, const EigenSupport::V3d &n, const EigenSupport::V3d &w, double eps, double timestep) const;
  void frictionPotentialGradient(const EigenSupport::V12d &x, const EigenSupport::V12d &xlast, const EigenSupport::V3d &n, const EigenSupport::V3d &w, double eps, double timestep, EigenSupport::V12d &grad) const;
  void frictionPotentialHessian(const EigenSupport::V12d &x, const EigenSupport::V12d &xlast, const EigenSupport::V3d &n, const EigenSupport::V3d &w, double eps, double timestep, EigenSupport::M12d &hess) const;

  int numPairs;
  int numObjects;

  const std::array<int, 4> *const objectIDs;
  const std::array<int, 4> *const pointTrianglePairs;
  const int *const triangleIDs;

  const int *const objectDOFOffsets;
  const Mesh::TriMeshRef *const surfaceMeshesRest;
  const std::vector<int> *const vertexIDToSampleIDs;
  const EigenSupport::SpMatD *const sampleEmbeddingWeights;
  const std::vector<double> *const sampleWeights;

  EigenSupport::EigenArray<EigenSupport::V3d> barycentricWeights, normals, closestPositions;
  std::vector<int> contactStatus;
  std::vector<int> contactDistances;
  std::vector<std::vector<int>> neighboringTriangles;
  std::vector<double> contactForceMags;
  mutable EigenSupport::VXd fn, ff;
  mutable double eng_fn = 0, eng_ff = 0;

  EigenSupport::EigenArray<EigenSupport::M12d> hessianBlocks;
  EigenSupport::EigenArray<EigenSupport::V12d> gradientBlocks;  // any quadratic energy can be written as 0.5 xT * hess * x + grad(0)^T x
  EigenSupport::SpMatD hessTemplate;
  EigenSupport::EntryMap entryMap;

  PointTrianglePairCouplingEnergyWithCollisionBuffer *buf = nullptr;
  PosFunction toPosFunc = [](const EigenSupport::V3d &x, EigenSupport::V3d &p, int) {
    p = x;
  };

  PosFunction toLastPosFunc = [](const EigenSupport::V3d &x, EigenSupport::V3d &p, int) {
    p = x;
  };

  std::vector<int> allDOFs;
  double coeff = 1.0;
  double frictionCoeff = 0;
  double frictionContactForceMag = 1.0;
  double timestep = 0;
  double eps = 0;
};
}  // namespace Contact
}  // namespace pgo