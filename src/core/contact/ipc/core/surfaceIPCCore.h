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
#include <atomic>
#include <cstdint>

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
  double computeMaxStepSize(EigenSupport::ConstRefVecXd x_surf, EigenSupport::ConstRefVecXd dx_surf) const;
  std::int64_t getContactClampCount() const { return contactClampCount_.load(std::memory_order_relaxed); }
  double getMinContactFeasibleAlphaThisSolve() const { return minContactFeasibleAlphaThisSolve_.load(std::memory_order_relaxed); }
  void resetContactMaxStepStats() const;

  const std::vector<PTPair> &getPTPairs() const { return ptPairs_; }
  const std::vector<EEPair> &getEEPairs() const { return eePairs_; }

  int getNumSurfaceVertices() const { return topology_.numVerts; }
  int getNumSurfaceDOFs() const { return topology_.numSurfaceDOFs(); }

private:
  void findCollisionPairs(const VXd &positions) const;

  static V3d vtx(const VXd &x, int i)
  {
    return x.segment<3>(3 * i);
  }

  double dhat = 1e-1;
  double kappa = 0.1;
  double eps_ee = 0.0;
  double slackness = 1.0;
  SurfaceIPCTopology topology_;
  mutable std::atomic<std::int64_t> contactClampCount_{0};
  mutable std::atomic<double> minContactFeasibleAlphaThisSolve_{1.0};
  mutable std::vector<PTPair> ptPairs_;
  mutable std::vector<EEPair> eePairs_;
};

}  // namespace CIPC
}  // namespace Contact
}  // namespace pgo
