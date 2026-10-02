/*
copyright to Bohan Wang
*/

#pragma once

#include "ipc/core/surfaceIPCCore.h"
#include "cudaPotentialEnergy.h"

#include <cstddef>
#include <memory>

namespace pgo::Contact::CIPC
{

// GPU counterpart of EmbeddedSurfaceIPCPotentialEnergy. The sparse map acts on
// simulation displacements; external vertices have zero map rows. Topology and
// rest measures are uploaded once. This initial implementation is frictionless.
class CudaEmbeddedSurfaceIPCPotentialEnergy final : public NonlinearOptimization::CudaPotentialEnergy
{
public:
  CudaEmbeddedSurfaceIPCPotentialEnergy(
    const EigenSupport::MXd &surfaceRestVertices,
    const EigenSupport::MXi &surfaceTriangles,
    const EigenSupport::SpMatD &surfaceFromSimulationDispMap,
    const SurfaceIPCCore::Parameters &ipcParams = {});
  CudaEmbeddedSurfaceIPCPotentialEnergy(
    const EigenSupport::MXd &surfaceRestVertices,
    const EigenSupport::MXi &surfaceTriangles,
    const EigenSupport::SpMatD &surfaceFromSimulationDispMap,
    const std::vector<uint8_t> &vertexIsDeformableMask,
    const SurfaceIPCCore::Parameters &ipcParams = {});
  ~CudaEmbeddedSurfaceIPCPotentialEnergy() override;

  int getNumDOFs() const override;
  double func(const double *deviceX) override;
  EnergyChange energyChange(const double *deviceBase, const double *deviceTrial) override;
  void addGradient(const double *deviceX, double *deviceGradient) override;
  void prepareHessian(const double *deviceX) override;
  void prepareExactHessian(const double *deviceX) override;
  void addHessianVector(const double *deviceV, double *deviceOut) const override;
  void addDiagonal(double *deviceDiagonal) const override;
  double computeMaxStepSize(const double *deviceX, const double *deviceDx) const override;

  std::size_t getNumActivePairs() const;

private:
  void prepareHessianImpl(const double *deviceX, bool projectPSD);
  class Impl;
  std::unique_ptr<Impl> impl_;
};

}  // namespace pgo::Contact::CIPC
