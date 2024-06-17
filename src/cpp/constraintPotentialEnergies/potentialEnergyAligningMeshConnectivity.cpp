/*
author: Bohan Wang
copyright to USC, MIT
*/

#include "potentialEnergyAligningMeshConnectivity.h"

#include "pgoLogging.h"

void pgo::ConstraintPotentialEnergies::PotentialEnergyAligningMeshConnectivity::setDOFs(const std::vector<int> &dofs)
{
  PGO_ALOG((int)dofs.size() == (int)allDOFs.size());
  allDOFs = dofs;
}