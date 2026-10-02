/*
copyright to Bohan Wang
*/

#pragma once

#include "ipc/core/surfaceIPCCore.h"
#include "cudaPotentialEnergy.h"

#include <memory>

namespace pgo::Contact::CIPC
{

// Separate lagged Coulomb potential. updateFriction freezes contact geometry
// before a Newton solve. Keep the original physical-step reference unchanged
// across lagging iterations. Prescribed vertices use their current placement
// in both reference and lagged states, matching EmbeddedSurfaceIPCPotentialEnergy.
class CudaSurfaceIPCFrictionPotentialEnergy final : public NonlinearOptimization::CudaPotentialEnergy
{
public:
  CudaSurfaceIPCFrictionPotentialEnergy(const EigenSupport::MXd &surfaceRestVertices,
    const EigenSupport::MXi &surfaceTriangles, const EigenSupport::SpMatD &surfaceFromSimulationDispMap,
    const std::vector<uint8_t> &vertexIsDeformableMask, const SurfaceIPCCore::Parameters &parameters);
  ~CudaSurfaceIPCFrictionPotentialEnergy() override;

  void updateFriction(const double *deviceReference, const double *deviceLagged, double timestep);
  void updateExternalSurface(int firstVertex, const EigenSupport::MXd &positions);
  void clear();
  std::size_t getNumFrictionPairs() const;

  int getNumDOFs() const override;
  double func(const double *deviceX) override;
  EnergyChange energyChange(const double *deviceBase, const double *deviceTrial) override;
  void addGradient(const double *deviceX, double *deviceGradient) override;
  void prepareHessian(const double *deviceX) override;
  void addHessianVector(const double *deviceV, double *deviceOut) const override;
  void addDiagonal(double *deviceDiagonal) const override;
  // This convex tangential potential has no additional feasible-domain bound.
  // The accompanying barrier energy supplies collision certification.
  double computeMaxStepSize(const double *, const double *) const override { return 1.0; }

private:
  class Impl;
  std::unique_ptr<Impl> impl_;
};

}  // namespace pgo::Contact::CIPC
