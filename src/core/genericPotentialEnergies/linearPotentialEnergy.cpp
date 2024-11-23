/*
author: Bohan Wang
copyright to USC, MIT
*/

#include "linearPotentialEnergy.h"
#include "pgoLogging.h"

#include <numeric>

using namespace pgo::PredefinedPotentialEnergies;

LinearPotentialEnergy::LinearPotentialEnergy(const EigenSupport::VXd &b_):
  b(b_)
{
  allDOFs.assign(b.size(), 0);
  std::iota(allDOFs.begin(), allDOFs.end(), 0);
}

void LinearPotentialEnergy::setDOFs(const std::vector<int> &dofs)
{
  PGO_ALOG((int)dofs.size() == (int)b.size());
  allDOFs = dofs;
}