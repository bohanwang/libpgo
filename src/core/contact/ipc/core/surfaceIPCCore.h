/*
copyright to Bohan Wang
*/

#pragma once

#include "EigenDef.h"
#include "ipc/core/surfaceIPCActiveSet.h"
#include "ipc/core/surfaceIPCPairs.h"
#include "ipc/external/obstacleSurface.h"
#include "ipc/topology/surfaceIPCTopology.h"
#include "solveDiagnostics.h"

#include <cstdint>
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
  SurfaceIPCCore(const Parameters &params, std::vector<ObstacleSurface> obstacles)
  {
    setParameters(params);
    setObstacles(std::move(obstacles));
  }
  SurfaceIPCCore(const SurfaceIPCCore &other);
  SurfaceIPCCore &operator=(const SurfaceIPCCore &other);

  void setParameters(const Parameters &params);
  Parameters getParameters() const;

  void setMesh(const EigenSupport::MXd &V, const EigenSupport::MXi &F);
  SurfaceIPCActiveSet buildActiveSet(EigenSupport::ConstRefVecXd x_surf) const;

  double computeEnergy(EigenSupport::ConstRefVecXd x_surf) const;
  void computeGradient(EigenSupport::ConstRefVecXd x_surf, EigenSupport::RefVecXd g_surf) const;
  void computeHessian(EigenSupport::ConstRefVecXd x_surf, EigenSupport::SpMatD &H_surf) const;
  void computeAll(EigenSupport::ConstRefVecXd x_surf, double &energy, EigenSupport::VXd &g_surf, EigenSupport::SpMatD &H_surf) const;
  double computeEnergy(const SurfaceIPCActiveSet &activeSet) const;
  void computeGradient(const SurfaceIPCActiveSet &activeSet, EigenSupport::RefVecXd g_surf) const;
  void computeHessian(const SurfaceIPCActiveSet &activeSet, EigenSupport::SpMatD &H_surf) const;
  void computeAll(const SurfaceIPCActiveSet &activeSet, double &energy, EigenSupport::VXd &g_surf, EigenSupport::SpMatD &H_surf) const;
  NonlinearOptimization::MaxStepResult computeMaxStepLimit(EigenSupport::ConstRefVecXd x_surf, EigenSupport::ConstRefVecXd dx_surf) const;

  const SurfaceIPCTopology& topology() const { return topology_; }

  void updateObstacleStage(double tStart, double tEnd);

private:
  void setObstacles(std::vector<ObstacleSurface> obstacles);

  double dhat = 1e-1;
  double dhat_external = 1e-1;
  double kappa = 0.1;
  double eps_ee = 0.0;
  double slackness = 1.0;
  SurfaceIPCTopology topology_;
  std::vector<ObstacleSurface> obstacles_;
};

}  // namespace CIPC
}  // namespace Contact
}  // namespace pgo
