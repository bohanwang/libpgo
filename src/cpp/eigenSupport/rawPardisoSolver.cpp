#include "rawPardisoSolver.h"

#include <thread>
#include <iostream>
#include <algorithm>
#include <cstdio>

/* PARDISO prototype. */
void pardisoinit(void *, int *, int *, int *, double *, int *);
void pardiso(void *, int *, int *, int *, int *, int *,
  double *, int *, int *, int *, int *, int *,
  int *, double *, double *, int *, double *);

PardisoSolverEigen::PardisoSolverEigen()
{
}

PardisoSolverEigen::~PardisoSolverEigen()
{
  /* -------------------------------------------------------------------- */
  /* ..  Termination and release of memory.                               */
  /* -------------------------------------------------------------------- */
  phase = -1; /* Release internal memory. */

  int idum;
  double ddum;
  int error = 0;
  int nrhs = 0;

  pardiso(pt.data(), &maxfct, &mnum, &mtype, &phase,
    &n, &ddum, ia.data(), ja.data(), &idum, &nrhs,
    iparm.data(), &msglvl, &ddum, &ddum, &error, dparm.data());
}

int PardisoSolverEigen::stage1(const SpMatD &A)
{
  facted = 1;

  int n = (int)A.rows();
  int solver = 0; /* use sparse direct solver */
  int error = 0;
  pardisoinit(pt.data(), &mtype, &solver, iparm.data(), dparm.data(), &error);

  if (error != 0) {
    if (error == -10)
      printf("No license file found \n");
    if (error == -11)
      printf("License is expired \n");
    if (error == -12)
      printf("Wrong username or hostname \n");
    return 1;
  } else
    printf("[PARDISO]: License check was successful ... \n");

  iparm[2] = (int)std::thread::hardware_concurrency();

  maxfct = 1; /* Maximum number of numerical factorizations.  */
  mnum = 1;   /* Which factorization to use. */

  msglvl = 0; /* Print statistical information  */
  error = 0;  /* Initialize error flag */

  /* -------------------------------------------------------------------- */
  /* ..  Convert matrix from 0-based C-notation to Fortran 1-based        */
  /*     notation.                                                        */
  /* -------------------------------------------------------------------- */
  ia.reserve(A.rows() + 1);
  ja.reserve(A.nonZeros());
  entries.reserve(A.nonZeros());
  entryIndex.reserve(A.nonZeros());

  for (IDX i = 0; i < A.rows(); i++) {
    ia.push_back((int)ja.size() + 1);

    for (SpMatD::InnerIterator it(A, i); it; ++it) {
      if (it.col() < it.row()) {
        continue;
      }

      std::ptrdiff_t offset = findEntryOffset(A, (int)it.row(), (int)it.col());
      if (offset < 0)
        abort();

      entryIndex.push_back(offset);
      ja.push_back((int)it.col() + 1);
      entries.push_back(it.value());
    }
  }
  ia.push_back((int)ja.size() + 1);

  /* -------------------------------------------------------------------- */
  /* ..  Reordering and Symbolic Factorization.  This step also allocates */
  /*     all memory that is necessary for the factorization.              */
  /* -------------------------------------------------------------------- */
  int idum;
  double ddum;
  int nrhs = 0;

  phase = 11;
  pardiso(pt.data(), &maxfct, &mnum, &mtype, &phase,
    &n, entries.data(), ia.data(), ja.data(), &idum, &nrhs,
    iparm.data(), &msglvl, &ddum, &ddum, &error, dparm.data());

  if (error != 0) {
    printf("\nERROR during symbolic factorization: %d", error);
    return 1;
  }
  printf("\nReordering completed ... ");
  printf("\nNumber of nonzeros in factors  = %d", iparm[17]);
  printf("\nNumber of factorization MFLOPS = %d", iparm[18]);

  return 0;
}

int PardisoSolverEigen::stage2(const SpMatD &A)
{
  /* -------------------------------------------------------------------- */
  /* ..  Numerical factorization.                                         */
  /* -------------------------------------------------------------------- */
  phase = 22;

  int idum;
  double ddum;
  int nrhs = 0;
  int error = 0;

  pardiso(pt.data(), &maxfct, &mnum, &mtype, &phase,
    &n, entries.data(), ia.data(), ja.data(), &idum, &nrhs,
    iparm.data(), &msglvl, &ddum, &ddum, &error, dparm.data());

  if (error != 0) {
    printf("\nERROR during numerical factorization: %d", error);
    exit(2);
  }
  printf("\nFactorization completed ...\n ");

  return 0;
}

int PardisoSolverEigen::stage3(int nrhs, const double *rhs, double *x)
{
  /* -------------------------------------------------------------------- */
  /* ..  Back substitution and iterative refinement.                      */
  /* -------------------------------------------------------------------- */
  phase = 33;

  iparm[7] = 1; /* Max numbers of iterative refinement steps. */

  int error = 0;
  int idum;

  pardiso(pt.data(), &maxfct, &mnum, &mtype, &phase,
    &n, entries.data(), ia.data(), ja.data(), &idum, &nrhs,
    iparm.data(), &msglvl, const_cast<double *>(rhs), x, &error, dparm.data());

  if (error != 0) {
    printf("\nERROR during solution: %d", error);
    return 1;
  }

  printf("\nSolve completed ... ");

  return 0;
}

std::ptrdiff_t PardisoSolverEigen::findEntryOffset(const SpMatD &A, int row, int col)
{
  const SpMatD::StorageIndex *rowBegin = A.innerIndexPtr() + A.outerIndexPtr()[row];
  const SpMatD::StorageIndex *rowEnd = A.innerIndexPtr() + A.outerIndexPtr()[row + 1];

  auto it = std::lower_bound(rowBegin, rowEnd, col);
  if (it == rowEnd || *it != col) {
    std::cerr << "Cannot find the corresponding index!";
    abort();
  }

  return it - A.innerIndexPtr();
}