#pragma once

#include "EigenDef.h"

namespace pgo
{
namespace EigenSupport
{
int solveGenEigShInv(const SpMatD &A, const SpMatD &B, int numEigenvalues, VXd &eigenvalues, MXd &eigenvectors, double sigma, int verbose);
}
}  // namespace pgo