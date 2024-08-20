/*
author: Bohan Wang
copyright to USC, MIT
*/

#include "quadraticPotentialEnergy.h"
#include "pgoLogging.h"
#include "EigenSupport.h"

#include <tbb/enumerable_thread_specific.h>

#include <numeric>
#include <vector>
#include <iostream>

using namespace pgo;
using namespace pgo::PredefinedPotentialEnergies;

namespace ES = pgo::EigenSupport;

namespace pgo::PredefinedPotentialEnergies
{
class QuadraticPotentialEnergyCache
{
public:
  QuadraticPotentialEnergyCache(int n)
  {
    temp.setZero(n);
  }

  ES::VXd temp;
};

}  // namespace pgo::NonlinearOptimization

QuadraticPotentialEnergy::QuadraticPotentialEnergy(const ES::SpMatD &A_):
  A(A_), b(nullptr), isInParentheses(0)
{
  cache = std::make_shared<QuadraticPotentialEnergyCache>(A.rows());

  allDOFs.resize(A.rows());
  std::iota(allDOFs.begin(), allDOFs.end(), 0);
}

QuadraticPotentialEnergy::QuadraticPotentialEnergy(const ES::SpMatD &A_, int):
  A(ATA), b(nullptr), isInParentheses(1)
{
  ES::mm(A_, A_, ATA, 1);

  cache = std::make_shared<QuadraticPotentialEnergyCache>(A.rows());

  allDOFs.resize(A.rows());
  std::iota(allDOFs.begin(), allDOFs.end(), 0);
}

QuadraticPotentialEnergy::QuadraticPotentialEnergy(const ES::SpMatD &A_, const double *W, int):
  A(ATA), b(nullptr), isInParentheses(1)
{
  ES::aba(A_, W, ATA, 1);

  cache = std::make_shared<QuadraticPotentialEnergyCache>(A.rows());

  allDOFs.resize(A.rows());
  std::iota(allDOFs.begin(), allDOFs.end(), 0);
}

QuadraticPotentialEnergy::QuadraticPotentialEnergy(const ES::SpMatD &A_, const ES::VXd &b_):
  A(A_), b(&b_), isInParentheses(0)
{
  cache = std::make_shared<QuadraticPotentialEnergyCache>(A.rows());

  allDOFs.resize(A.rows());
  std::iota(allDOFs.begin(), allDOFs.end(), 0);
}

QuadraticPotentialEnergy::QuadraticPotentialEnergy(const ES::SpMatD &A_, const ES::VXd &b_, int):
  A(ATA), b(&bTA), isInParentheses(1)
{
  ES::mm(A_, A_, ATA, 1);

  bTA.resize(A.rows());
  ES::mv(A_, b_, bTA, 1);

  c = b_.dot(b_) * 0.5;

  cache = std::make_shared<QuadraticPotentialEnergyCache>(A.rows());

  allDOFs.resize(A.rows());
  std::iota(allDOFs.begin(), allDOFs.end(), 0);
}

QuadraticPotentialEnergy::QuadraticPotentialEnergy(const ES::SpMatD &A_, const ES::VXd &b_, const double *W, int):
  A(ATA), b(&bTA), isInParentheses(1)
{
  // = 1/2 x^T (A^T W A) x + b^T W Ax + 1/2 b^T W b
  ES::aba(A_, W, ATA, 1);

  ES::VXd Wb = b_;
  for (ES::IDX i = 0; i < A_.rows(); i++) {
    Wb[i] = b_[i] * W[i];
  }

  bTA.resize(A.rows());
  ES::mv(A_, Wb, bTA, 1);

  c = b_.dot(Wb) * 0.5;

  cache = std::make_shared<QuadraticPotentialEnergyCache>(A.rows());

  allDOFs.resize(A.rows());
  std::iota(allDOFs.begin(), allDOFs.end(), 0);
}

void QuadraticPotentialEnergy::setDOFs(const std::vector<int> &dofs)
{
  PGO_ALOG((int)dofs.size() == (int)A.rows());
  allDOFs = dofs;
}

double QuadraticPotentialEnergy::func(ES::ConstRefVecXd x) const
{
  ES::VXd &temp = cache->temp;
  double energy = ES::vTMv(A, x, temp) * 0.5;
  // std::cout << energy << ',';

  energy += c;
  // std::cout << energy << ',';
  if (b) {
    energy += (*b).dot(x);
  }
  // std::cout << energy << ',';
  return energy;
}

void QuadraticPotentialEnergy::gradient(ES::ConstRefVecXd x, ES::RefVecXd grad) const
{
  ES::mv(A, x, grad);

  if (b) {
    grad += (*b);
  }
}

void QuadraticPotentialEnergy::hessian(ES::ConstRefVecXd, ES::SpMatD &hess) const
{
  memcpy(hess.valuePtr(), A.valuePtr(), sizeof(double) * A.nonZeros());
}

void QuadraticPotentialEnergy::hessianVector(EigenSupport::ConstRefVecXd, EigenSupport::ConstRefVecXd vec, EigenSupport::RefVecXd hVec) const
{
  ES::mv(A, vec, hVec);
}

void QuadraticPotentialEnergy::gradientComponent(ES::SpMatD *A_, ES::VXd *b_) const
{
  if (A_)
    *A_ = A;

  if (b_) {
    if (b) {
      *b_ = *b;
    }
    else {
      b_->setZero(getNumDOFs());
    }
  }
}