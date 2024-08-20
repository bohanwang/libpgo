/*
author: Bohan Wang
copyright to USC, MIT
*/

#pragma once

#include "EigenSupport.h"

namespace pgo
{
namespace PredefinedPotentialEnergies
{
class LaplacianProblem
{
public:
  LaplacianProblem(const EigenSupport::SpMatD &L, int biLaplacian, const std::vector<int> &fixedDOFs,
    const EigenSupport::VXd &fixedValues, double variableLow, double variableHi);

  void setVariableConstraints(const std::vector<int> &fixedDOFs, const EigenSupport::VXd &fixedValues, double variableLow, double variableHi);

  void solve(EigenSupport::RefVecXd xfinal, int numIter, double eps, int verbose);

protected:
  EigenSupport::SpMatD sys;
  EigenSupport::VXd xlow, xhi, xinit;
};

}  // namespace NonlinearOptimization
}  // namespace pgo