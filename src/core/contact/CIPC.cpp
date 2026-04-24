/*
copyright to Bohan Wang
*/

#include "CIPC.h"

#include "scopedProfileSection.h"
#include "ipc/profiling/surfaceIPCProfiling.h"

#include "pgoLogging.h"

#include <numeric>
#include <stdexcept>
#include <vector>

namespace pgo
{
namespace Contact
{
namespace CIPC
{

void CIPCPotentialEnergy::setMesh(const MXd &V, const MXi &F)
{
  allDOFs_.resize(static_cast<std::size_t>(V.rows()) * 3);
  std::iota(allDOFs_.begin(), allDOFs_.end(), 0);

  restPosition.resize(V.rows() * 3);
  for (int i = 0; i < V.rows(); ++i)
    restPosition.segment<3>(3 * i) = V.row(i).transpose();

  core.setMesh(V, F);
}

void CIPCPotentialEnergy::syncCoreParametersFromWrapper() const
{
  Profiling::ScopedProfileSection scopedProfile(SurfaceIPCProfileSections::kWrapperSync);

  SurfaceIPCCore::Parameters params;
  params.dhat = dhat;
  params.kappa = kappa;
  params.eps_ee = eps_ee;
  params.slackness = slackness;
  core.setParameters(params);
}

VXd CIPCPotentialEnergy::toSurfacePositions(EigenSupport::ConstRefVecXd x) const
{
  return isInputDisp ? VXd(restPosition + x) : VXd(x);
}

VXd CIPCPotentialEnergy::toSurfaceDisplacements(EigenSupport::ConstRefVecXd dx) const
{
  return VXd(dx);
}

double CIPCPotentialEnergy::computeFloorEnergy(const VXd &x_surf) const
{
  if (!useFloor)
    return 0.0;

  Profiling::ScopedProfileSection scopedProfile(SurfaceIPCProfileSections::kFloorPostPass);

  double energy = 0.0;
  for (int vi = 0; vi < x_surf.size() / 3; ++vi) {
    const double dz = x_surf[3 * vi + 2] - floorHeight;
    if (dz < 0.0)
      energy += 0.5 * floorKappa * dz * dz;
  }
  return energy;
}

void CIPCPotentialEnergy::addFloorGradient(const VXd &x_surf, EigenSupport::RefVecXd grad) const
{
  if (!useFloor)
    return;

  Profiling::ScopedProfileSection scopedProfile(SurfaceIPCProfileSections::kFloorPostPass);

  if (grad.size() != x_surf.size())
    throw std::runtime_error("Gradient vector has wrong size");

  for (int vi = 0; vi < x_surf.size() / 3; ++vi) {
    const double dz = x_surf[3 * vi + 2] - floorHeight;
    if (dz < 0.0)
      grad[3 * vi + 2] += floorKappa * dz;
  }
}

void CIPCPotentialEnergy::addFloorHessian(const VXd &x_surf, EigenSupport::SpMatD &hess) const
{
  if (!useFloor)
    return;

  Profiling::ScopedProfileSection scopedProfile(SurfaceIPCProfileSections::kFloorPostPass);

  std::vector<TripletD> floorTriplets;
  floorTriplets.reserve(static_cast<std::size_t>(x_surf.size() / 3));
  for (int vi = 0; vi < x_surf.size() / 3; ++vi) {
    const double dz = x_surf[3 * vi + 2] - floorHeight;
    if (dz < 0.0) {
      const int row = 3 * vi + 2;
      floorTriplets.emplace_back(row, row, floorKappa);
    }
  }

  if (floorTriplets.empty())
    return;

  SpMatD floorHess(x_surf.size(), x_surf.size());
  floorHess.setFromTriplets(floorTriplets.begin(), floorTriplets.end());
  hess += floorHess;
}

double CIPCPotentialEnergy::func(EigenSupport::ConstRefVecXd x) const
{
  syncCoreParametersFromWrapper();
  const VXd x_surf = toSurfacePositions(x);
  return core.computeEnergy(x_surf) + computeFloorEnergy(x_surf);
}

void CIPCPotentialEnergy::gradient(EigenSupport::ConstRefVecXd x, EigenSupport::RefVecXd grad) const
{
  syncCoreParametersFromWrapper();
  const VXd x_surf = toSurfacePositions(x);
  core.computeGradient(x_surf, grad);
  addFloorGradient(x_surf, grad);
}

void CIPCPotentialEnergy::hessian(EigenSupport::ConstRefVecXd, EigenSupport::SpMatD &) const
{
  throw std::runtime_error("CIPCPotentialEnergy::hessian() should not be called directly. Use hessianDirect() instead.");
}

void CIPCPotentialEnergy::createHessian(EigenSupport::SpMatD &) const
{
  throw std::runtime_error("CIPCPotentialEnergy::createHessian() should not be called directly. Use hessianDirect() instead.");
}

double CIPCPotentialEnergy::computeMaxStepSize(EigenSupport::ConstRefVecXd x, EigenSupport::ConstRefVecXd dx) const
{
  syncCoreParametersFromWrapper();
  const VXd x_surf = toSurfacePositions(x);
  const VXd dx_surf = toSurfaceDisplacements(dx);
  return core.computeMaxStepSize(x_surf, dx_surf);
}

void CIPCPotentialEnergy::hessianDirect(EigenSupport::ConstRefVecXd x, EigenSupport::SpMatD &hess) const
{
  syncCoreParametersFromWrapper();
  const VXd x_surf = toSurfacePositions(x);
  core.computeHessian(x_surf, hess);
  if (auto logger = Logging::lgr(); logger)
    SPDLOG_LOGGER_INFO(logger, "Computing Hessian with {} PT pairs and {} EE pairs", core.getPTPairs().size(), core.getEEPairs().size());
  addFloorHessian(x_surf, hess);
}

}  // namespace CIPC
}  // namespace Contact
}  // namespace pgo
