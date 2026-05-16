/*
copyright to Bohan Wang
*/

#pragma once

#include "EigenDef.h"
#include "ipc/core/surfaceIPCPairs.h"
#include "ipc/core/surfaceIPCPreparedState.h"
#include "ipc/external/obstacleSurface.h"
#include "ipc/topology/surfaceIPCTopology.h"
#include "solveDiagnostics.h"

#include <cstdint>
#include <memory>
#include <vector>

namespace pgo
{
namespace Contact
{
namespace CIPC
{

// =========================================================================
//  Surface-space IPC core
// =========================================================================
class SurfaceIPCCore
{
public:
  struct Parameters
  {
    double dhat          = 1e-1;
    double dhat_external = 1e-1;  // external pair activation distance; defaults to dhat for self-only numerical invariance
    double kappa         = 0.1;
    double eps_ee        = 0.0;
    double slackness     = 1.0;
  };

  SurfaceIPCCore() = default;
  explicit SurfaceIPCCore(const Parameters &params) { setParameters(params); }
  SurfaceIPCCore(const SurfaceIPCCore &other);
  SurfaceIPCCore &operator=(const SurfaceIPCCore &other);

  void setParameters(const Parameters &params);
  Parameters getParameters() const;

  void setMesh(const EigenSupport::MXd &V, const EigenSupport::MXi &F);
  void prepareForSurfacePositions(EigenSupport::ConstRefVecXd x_surf) const;

  double computeEnergy(EigenSupport::ConstRefVecXd x_surf) const;
  void computeGradient(EigenSupport::ConstRefVecXd x_surf, EigenSupport::RefVecXd g_surf) const;
  void computeHessian(EigenSupport::ConstRefVecXd x_surf, EigenSupport::SpMatD &H_surf) const;
  void computeAll(EigenSupport::ConstRefVecXd x_surf, double &energy, EigenSupport::VXd &g_surf, EigenSupport::SpMatD &H_surf) const;
  double computeEnergyWithPreparedPairs() const;
  void computeGradientWithPreparedPairs(EigenSupport::RefVecXd g_surf) const;
  void computeHessianWithPreparedPairs(EigenSupport::SpMatD &H_surf) const;
  void computeAllWithPreparedPairs(double &energy, EigenSupport::VXd &g_surf, EigenSupport::SpMatD &H_surf) const;
  NonlinearOptimization::MaxStepResult computeMaxStepLimit(EigenSupport::ConstRefVecXd x_surf, EigenSupport::ConstRefVecXd dx_surf) const;

  SurfaceIPCPreparedState& preparedState() const { return preparedState_; }
  const SurfaceIPCTopology& topology() const { return topology_; }

  // Obstacle (external) registration
  int32_t addObstacleSurface(std::shared_ptr<ObstacleSurface> obs);
  void    clearObstacleSurfaces();
  void    updateObstacleStage(double tStart, double tEnd);

private:
  void findCollisionPairs(const EigenSupport::VXd &positions) const;

  double dhat = 1e-1;
  double dhat_external = 1e-1;
  double kappa = 0.1;
  double eps_ee = 0.0;
  double slackness = 1.0;
  SurfaceIPCTopology topology_;
  mutable SurfaceIPCPreparedState preparedState_;
  std::vector<std::shared_ptr<ObstacleSurface>> obstacles_;
};

}  // namespace CIPC
}  // namespace Contact
}  // namespace pgo
