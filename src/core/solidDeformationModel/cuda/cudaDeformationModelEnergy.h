#pragma once

#include "cudaPotentialEnergy.h"

#include <memory>

namespace pgo::SolidDeformationModel
{
class SimulationMesh;
enum class DeformationModelElasticMaterial;

// Stable Neo-Hookean or StVK-vol tetrahedral/trilinear cubic FEM with identity
// plastic deformation. Cubic elements use the CPU's 2x2x2 Gauss quadrature and
// project the 9x9 material tangent independently at every quadrature point in
// prepareHessian(); StVK-vol adds the unprojected volume tangent afterward,
// matching CPU combined-material semantics. prepareExactHessian() retains both
// unmodified tangents.
// Enabled material feasibility defaults to the CPU's det(F)>1e-8 checks at
// each Gauss point along the displacement chord. WholeCell additionally
// certifies this bound throughout each trilinear cell; uncertain certificates
// fail closed. Neither mode changes the quadrature energy or its derivatives.
// Zero-direction queries check the initial state. Disabling material max step
// disables both motion guards; WholeCell still requires certified rest cells.
// x is displacement from the mesh's rest positions. Mesh and material data are
// copied at construction; the mesh need not outlive this energy.
class CudaDeformationModelEnergy final : public NonlinearOptimization::CudaPotentialEnergy
{
public:
  enum class CubicValidity
  {
    Gauss,     // CPU-compatible sampling; corner/interior inversions may remain.
    WholeCell  // Additional certificate over the entire trilinear cell.
  };

  explicit CudaDeformationModelEnergy(const SimulationMesh &mesh, bool enableMaterialMaxStep = true,
    CubicValidity cubicValidity = CubicValidity::Gauss);
  CudaDeformationModelEnergy(const SimulationMesh &mesh, DeformationModelElasticMaterial material,
    bool enableMaterialMaxStep = true, CubicValidity cubicValidity = CubicValidity::Gauss);
  ~CudaDeformationModelEnergy() override;
  CudaDeformationModelEnergy(const CudaDeformationModelEnergy &) = delete;
  CudaDeformationModelEnergy &operator=(const CudaDeformationModelEnergy &) = delete;

  int getNumDOFs() const override;
  double func(const double *deviceX) override;
  EnergyChange energyChange(const double *deviceBase, const double *deviceTrial) override;
  void addGradient(const double *deviceX, double *deviceGradient) override;
  void prepareHessian(const double *deviceX) override;
  void prepareExactHessian(const double *deviceX) override;
  void addHessianVector(const double *deviceV, double *deviceOut) const override;
  void addDiagonal(double *deviceDiagonal) const override;
  double computeMaxStepSize(const double *deviceX, const double *deviceDx) const override;

private:
  void prepareHessian(const double *deviceX, bool projectMaterialTangent);
  class Impl;
  std::unique_ptr<Impl> data;
};
}  // namespace pgo::SolidDeformationModel
