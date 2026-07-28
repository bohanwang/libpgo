/*
copyright to Bohan Wang
*/

#pragma once

#include "EigenDef.h"
#include "ipc/geometry/ipcBarrier.h"
#include "ipc/geometry/ipcCCD.h"
#include "ipc/geometry/ipcDistancePrimitives.h"
#include "ipc/geometry/ipcHessianProjection.h"
#include "ipc/core/surfaceIPCPairs.h"
#include "ipc/topology/surfaceIPCTopology.h"
#include "potentialEnergy.h"

#include <vector>
#include <array>
#include <cmath>
#include <algorithm>

namespace pgo
{
namespace Contact
{
namespace CIPC
{

using namespace pgo::EigenSupport;

// =========================================================================
//  Surface-space IPC core
// =========================================================================
class SurfaceIPCCore
{
public:
  struct Parameters
  {
    double dhat = 1e-1;
    double kappa = 0.1;
    double eps_ee = 0.0;
    double slackness = 1.0;
    // Newton solves use a blockwise PSD approximation by default. Disable
    // this only when the unprojected, energy-consistent barrier Hessian is required.
    bool projectHessianToPSD = true;
  };

  SurfaceIPCCore() = default;
  explicit SurfaceIPCCore(const Parameters &params) { setParameters(params); }
  SurfaceIPCCore(const SurfaceIPCCore &other);
  SurfaceIPCCore &operator=(const SurfaceIPCCore &other);

  void setParameters(const Parameters &params);
  Parameters getParameters() const;

  void setMesh(const MXd &V, const MXi &F);

  double computeEnergy(EigenSupport::ConstRefVecXd x_surf) const;
  void computeGradient(EigenSupport::ConstRefVecXd x_surf, EigenSupport::RefVecXd g_surf) const;
  void computeHessian(EigenSupport::ConstRefVecXd x_surf, EigenSupport::SpMatD &H_surf) const;
  void computeAll(EigenSupport::ConstRefVecXd x_surf, double &energy, VXd &g_surf, SpMatD &H_surf) const;
  void prepareForSurfacePositions(EigenSupport::ConstRefVecXd x_surf) const;
  bool isPreparedFor(EigenSupport::ConstRefVecXd x_surf) const;
  void invalidatePreparedState() const;
  double computeEnergyWithPreparedPairs() const;
  void computeGradientWithPreparedPairs(EigenSupport::RefVecXd g_surf) const;
  void computeHessianWithPreparedPairs(EigenSupport::SpMatD &H_surf) const;
  void computeAllWithPreparedPairs(double &energy, VXd &g_surf, SpMatD &H_surf) const;
  double computeMaxStepSize(EigenSupport::ConstRefVecXd x_surf, EigenSupport::ConstRefVecXd dx_surf) const;

  const std::vector<PTPair> &getPTPairs() const { return ptPairs_; }
  const std::vector<EEPair> &getEEPairs() const { return eePairs_; }

  int getNumSurfaceVertices() const { return topology_.numVerts; }
  int getNumSurfaceDOFs() const { return topology_.numSurfaceDOFs(); }

private:
  void findCollisionPairs(const VXd &positions) const;
  void requirePreparedState() const;
  void validateSurfaceState(EigenSupport::ConstRefVecXd x_surf, const char *argumentName) const;
  void validateSurfaceGradient(EigenSupport::RefVecXd g_surf) const;

  static V3d vtx(const VXd &x, int i)
  {
    return x.segment<3>(3 * i);
  }

  double dhat = 1e-1;
  double kappa = 0.1;
  double eps_ee = 0.0;
  double slackness = 1.0;
  bool projectHessianToPSD = true;
  bool hasMesh_ = false;
  SurfaceIPCTopology topology_;
  mutable std::vector<PTPair> ptPairs_;
  mutable std::vector<EEPair> eePairs_;
  mutable bool hasPreparedState_ = false;
  mutable VXd preparedPositions_;
};

}  // namespace CIPC
}  // namespace Contact
}  // namespace pgo
