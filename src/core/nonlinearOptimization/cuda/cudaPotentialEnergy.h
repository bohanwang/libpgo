#pragma once

#include <cmath>

namespace pgo::NonlinearOptimization
{

// All vectors reside on the active CUDA device; operations use its default
// stream. Gradients, Hessian-vector products and diagonals ADD into the output.
// prepareHessian() caches the operator at x until the next preparation. Energy
// and max-step evaluations may return host scalars, but do not move DOF arrays
// or Hessians to the host. Coordinates follow the owning simulation's convention.
class CudaPotentialEnergy
{
public:
  struct EnergyChange
  {
    double value;
    // Absolute arithmetic scale used to estimate roundoff in value. Analytic
    // differences should use their local terms, excluding constant offsets.
    double roundoffScale;
  };

  virtual ~CudaPotentialEnergy() = default;
  virtual int getNumDOFs() const = 0;
  virtual double func(const double *deviceX) = 0;
  // Change between the actual represented endpoints. Implementations may form
  // trial-base on the GPU to avoid subtracting large, nearly equal energies.
  // The fallback preserves compatibility for energies without a local formula.
  virtual EnergyChange energyChange(const double *deviceBase, const double *deviceTrial)
  {
    const double base = func(deviceBase);
    const double trial = func(deviceTrial);
    return { trial - base, std::abs(base) + std::abs(trial) };
  }
  virtual void addGradient(const double *deviceX, double *deviceGradient) = 0;
  virtual void prepareHessian(const double *deviceX) = 0;
  // Optional unprojected tangent. Newton may try it after slow convergence,
  // then restore prepareHessian() if its linear solve detects bad curvature.
  // Energies whose normal operator is exact need no separate implementation.
  virtual void prepareExactHessian(const double *deviceX) { prepareHessian(deviceX); }
  virtual void addHessianVector(const double *deviceV, double *deviceOut) const = 0;
  virtual void addDiagonal(double *deviceDiagonal) const = 0;
  // Conservative fraction in [0, 1]. One certifies the complete original
  // chord, including any mapping/interpolation arithmetic required by the
  // energy. A smaller result is only a certified prefix. Zero is a valid
  // bound and must not be raised. A zero direction certifies the initial
  // state: return one only if all enabled model guards accept it. Disabled
  // guards must return one. Newton queries this before initial convergence.
  virtual double computeMaxStepSize(const double *deviceX, const double *deviceDx) const = 0;
};

}  // namespace pgo::NonlinearOptimization
