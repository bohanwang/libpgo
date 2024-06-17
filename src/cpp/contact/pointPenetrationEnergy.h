#pragma once

#include "contactEnergyUtilities.h"

#include "potentialEnergy.h"
#include "triMeshGeo.h"
#include "boundingVolumeTree.h"
#include "triMeshPseudoNormal.h"

namespace pgo
{
namespace Contact
{
struct PointPenetrationEnergyBuffer;

class PointPenetrationEnergy : public NonlinearOptimization::PotentialEnergy
{
public:
  PointPenetrationEnergy(int numPoints, int nAll, const double *constraintCoeffs, const double *constraintTargetPositions, const double *constraintNormals,
    const std::vector<std::vector<int>> &barycentricIdx, const std::vector<std::vector<double>> &barycentricWeights, double bcCoeffs,
    int useN, int checkPenetration, const std::vector<Mesh::TriMeshGeo> *externalMeshes = nullptr, const std::vector<Mesh::TriMeshBVTree> *externalMeshBVTrees = nullptr,
    const std::vector<Mesh::TriMeshPseudoNormal> *externalMeshNormals = nullptr);

  PointPenetrationEnergyBuffer *allocateBuffer() const;
  void freeBuffer(PointPenetrationEnergyBuffer *buf) const;
  void setBuffer(PointPenetrationEnergyBuffer *b) { buf = b; }

  void setComputePosFunction(PosFunction f) { posFunc = f; }
  void setComputeLastPosFunction(PosFunction f) { lastPosFunc = f; }

  virtual double func(EigenSupport::ConstRefVecXd u) const override;
  virtual void gradient(EigenSupport::ConstRefVecXd u, EigenSupport::RefVecXd grad) const override;
  virtual void hessian(EigenSupport::ConstRefVecXd u, EigenSupport::SpMatD &hess) const override;

  void setCoeff(double v) { coeffAll = v; }
  void setFrictionCoeff(double c) { frictionCoeff = c; }
  void setTimestep(double t) { timestep = t; }
  void setVelEps(double eps_) { eps = eps_; }

  void computeHessian();

  virtual void createHessian(EigenSupport::SpMatD &h) const override { h = hessianConstant; }
  virtual void getDOFs(std::vector<int> &dofs) const override { dofs = allDOFs; }
  virtual int getNumDOFs() const override { return (int)hessianConstant.rows(); }

protected:
  bool isInside(const EigenSupport::V3d &p, const EigenSupport::V3d p0, const EigenSupport::V3d n) const;

  double frictionPotential(const EigenSupport::V3d &x, const EigenSupport::V3d &xlast, const EigenSupport::V3d &n, double eps, double timestep) const;
  void frictionPotentialGradient(const EigenSupport::V3d &x, const EigenSupport::V3d &xlast, const EigenSupport::V3d &n, double eps, double timestep, EigenSupport::V3d &grad) const;
  void frictionPotentialHessian(const EigenSupport::V3d &x, const EigenSupport::V3d &xlast, const EigenSupport::V3d &n, double eps, double timestep, EigenSupport::M3d &hess) const;

  int nAll;
  int numPoints;
  const double *constraintCoeffs;
  const double *constraintTargetPositions;
  const double *constraintNormals;
  const std::vector<std::vector<int>> &barycentricIdx;
  const std::vector<std::vector<double>> &barycentricWeights;

  double coeffAll = 1.0;
  int useNormal = 1;
  int checkPenetration = 1;
  PosFunction posFunc, lastPosFunc;

  double frictionCoeff = 0;
  double frictionContactForceMag = 1.0;
  double timestep = 0;
  double eps = 0;  
  int saveDebugInfo = 0;

  const std::vector<Mesh::TriMeshGeo> *externalMeshes = nullptr;
  const std::vector<Mesh::TriMeshBVTree> *externalMeshBVTrees = nullptr;
  const std::vector<Mesh::TriMeshPseudoNormal> *externalMeshNormals = nullptr;

  std::vector<int> allDOFs;
  EigenSupport::SpMatD hessianConstant;
  EigenSupport::EntryMap hessianEntryMap;

  PointPenetrationEnergyBuffer *buf;
};

}  // namespace Contact
}  // namespace pgo