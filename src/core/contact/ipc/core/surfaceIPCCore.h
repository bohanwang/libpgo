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
    // CCD minimum-separation thickness (xi). 0 = classic contact-at-zero CCD.
    // If > 0, max-step finds the first time distance drops to ccd_thickness
    // and broad-phase swept AABBs are inflated by ccd_thickness on both sides
    // to stay sound under this contact definition.
    double ccd_thickness = 0.0;
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

  // Sample all registered obstacles at absolute time t. The obstacle pose is
  // treated as fixed during the subsequent solve / line-search; intra-frame
  // swept obstacle CCD is not performed.
  // Obstacles marked static via markObstacleStatic() are skipped — their
  // cache was built once during the markObstacleStatic call.
  void setObstacleTime(double t);

  // Mark the obstacle at the given slot as having a static (time-invariant)
  // pose, and immediately call update(0.0) to populate its cache exactly once.
  // Subsequent setObstacleTime() calls skip this obstacle.
  void markObstacleStatic(int32_t objectId);

private:
  void setObstacles(std::vector<ObstacleSurface> obstacles);

  double dhat = 1e-1;
  double dhat_external = 1e-1;
  double kappa = 0.1;
  double eps_ee = 0.0;
  double slackness = 1.0;
  double ccd_thickness = 0.0;
  SurfaceIPCTopology topology_;
  std::vector<ObstacleSurface> obstacles_;
  std::vector<bool> staticObstacles_;
};

}  // namespace CIPC
}  // namespace Contact
}  // namespace pgo
