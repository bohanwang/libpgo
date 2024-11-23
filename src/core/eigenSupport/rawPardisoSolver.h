#pragma once

#include <Eigen/Sparse>

#include <array>

class PardisoSolverEigen
{
public:
  PardisoSolverEigen();
  ~PardisoSolverEigen();

  typedef Eigen::SparseMatrix<double, Eigen::RowMajor> SpMatD;
  typedef Eigen::Index IDX;

  int stage1(const SpMatD &A);
  int stage2(const SpMatD &A);
  int stage3(int nrhs, const double *rhs, double *x);

protected:
  std::ptrdiff_t findEntryOffset(const SpMatD &A, int row, int col);

  std::array<int, 64> iparm;
  std::array<double, 64> dparm;
  std::array<void *, 64> pt;

  int maxfct, mnum, phase, error, msglvl, solver;
  int mtype = -2; /* Real symmetric matrix */
  int n;

  std::vector<int> ia, ja;
  std::vector<double> entries;
  std::vector<std::ptrdiff_t> entryIndex;

  int facted = 0;
};