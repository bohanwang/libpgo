/*
author: Bohan Wang
copyright to USC, MIT
*/

#include "laplacianProblem.h"

#include "potentialEnergies.h"
#include "quadraticPotentialEnergy.h"
#include "minimizeEnergy.h"

#include <numeric>
#include <iostream>

using namespace pgo;
using namespace pgo::PredefinedPotentialEnergies;

namespace ES = pgo::EigenSupport;

LaplacianProblem::LaplacianProblem(const EigenSupport::SpMatD &L, int biLaplacian, const std::vector<int> &fixedDOFs,
  const EigenSupport::VXd &fixedValues, double variableLow, double variableHi)
{
  if (biLaplacian) {
    // sys = L.transpose() * L;
    ES::mm(L, L, sys, 1);
  }
  else
    sys = L;

  xlow.resize(sys.rows());
  xhi.resize(sys.rows());

  setVariableConstraints(fixedDOFs, fixedValues, variableLow, variableHi);
}

void LaplacianProblem::setVariableConstraints(const std::vector<int> &fixedDOFs, const EigenSupport::VXd &fixedValues, double variableLow, double variableHi)
{
  for (Eigen::Index i = 0; i < sys.rows(); i++) {
    xlow[i] = variableLow;
    xhi[i] = variableHi;
  }

  for (int dof : fixedDOFs) {
    xlow[dof] = fixedValues[dof];
    xhi[dof] = fixedValues[dof];
  }

  xinit = fixedValues;
}

void LaplacianProblem::solve(EigenSupport::RefVecXd xfinal, int numIter, double eps, int verbose)
{
  std::shared_ptr<QuadraticPotentialEnergy> energy = std::make_shared<QuadraticPotentialEnergy>(sys);

  std::shared_ptr<NonlinearOptimization::PotentialEnergies> energyAll = std::make_shared<NonlinearOptimization::PotentialEnergies>(sys.rows());
  energyAll->addPotentialEnergy(energy);
  energyAll->init();

  ES::VXd lambda, g, clow, chi;

  xfinal = xinit;
  int ret = NonlinearOptimization::EnergyOptimizer::minimizeUsingKnitro(xfinal, energyAll, xlow, xhi, lambda, g,
    nullptr, clow, chi, numIter, eps, verbose, "");
  std::cout << "Solver ret: " << ret << std::endl;
}