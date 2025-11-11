/*
author: Bohan Wang
copyright to USC,MIT,NUS
*/

#include "invariantBasedMaterialStVK.h"
#include "simulationMesh.h"

#include <cmath>

using namespace pgo::SolidDeformationModel;

InvariantBasedMaterialStVK::InvariantBasedMaterialStVK(double E, double nu, double compressionRatio)
{
  SimulationMeshENuMaterial mat(E, nu);
  lambda = mat.getLambdaLame();
  mu = mat.getMuLame();

  if (compressionRatio > 0) {
    coeffJ = compressionRatio * E / (1.0 - 2.0 * nu);
  }
  else {
    coeffJ = 0;
  }
}

double InvariantBasedMaterialStVK::compute_psi(const double invariants[3]) const
{
  double IC = invariants[0];
  double IIC = invariants[1];
  //double IIIC = invariants[2]; // not needed for StVK

  double energy = 0.125 * lambda * (IC - 3.0) * (IC - 3.0) + 0.25 * mu * (IIC - 2.0 * IC + 3.0);

  if (coeffJ > 0) {
    double IIIC = invariants[2];
    double J = std::sqrt(IIIC);

    if (J < 1) {
      energy += -coeffJ * (J - 1.0) * (J - 1.0) * (J - 1.0) / 2592.0;
      // (J-1)^3 = (I^(1/2) - 1)^3
      // 3 * (I^(1/2) - 1)^2 1/2 I^(-1/2)
      // 6 * (I^(1/2 - 1)) (1/2) I^(-1/2) (1/2) I^(-1/2) + 3 * (I^(1/2) - 1)^2 (1/2) (-1/2) I^(-3/2)
      // 6/4 * (J - 1) / J + 3 (J - 1)^2 (-1/4) / J^3
    }
  }

  return energy;
}

void InvariantBasedMaterialStVK::compute_dpsi_dI(const double invariants[3], double gradient[3]) const
{
  //printf("Entered StVKIsotropicMaterial::ComputeEnergyGradient\n");

  double IC = invariants[0];
  gradient[0] = 0.25 * lambda * (IC - 3.0) - 0.5 * mu;
  gradient[1] = 0.25 * mu;
  gradient[2] = 0.0;

  if (coeffJ > 0) {
    double IIIC = invariants[2];
    double J = sqrt(IIIC);

    if (J < 1) {
      gradient[2] += -coeffJ * (J - 1.0) * (J - 1.0) / (1728.0 * J);
    }
  }
}
void InvariantBasedMaterialStVK::compute_d2psi_dI2(const double invariants[3], double hessian[6]) const
{
  // 11
  hessian[0] = 0.25 * lambda;
  // 12
  hessian[1] = 0.0;
  // 13
  hessian[2] = 0.0;
  // 22
  hessian[3] = 0.0;
  // 23
  hessian[4] = 0.0;
  // 33
  hessian[5] = 0.0;

  if (coeffJ > 0) {
    double IIIC = invariants[2];
    double J = sqrt(IIIC);

    if (J < 1.0) {
      hessian[5] += coeffJ * (1.0 - J) * (1.0 + J) / (3456.0 * J * J * J);
    }
  }
}
