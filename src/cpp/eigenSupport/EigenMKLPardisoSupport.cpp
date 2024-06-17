#include "EigenMKLPardisoSupport.h"

#include <iostream>
#include <chrono>

using namespace pgo;
using namespace pgo::EigenSupport;

using hclock = std::chrono::high_resolution_clock;
using hclockPt = hclock::time_point;

inline double dura(const hclockPt &t1, const hclockPt &t2)
{
  return static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count()) / 1e6;
}

EigenMKLPardisoSupport::EigenMKLPardisoSupport(const SpMatD &Ain, MatrixType mt, ReorderingType rt, int di, int ml, int mr,
  int tm, int sm, int inputMatrixIsUpper):
  mtype(mt),
  rtype(rt), directIterative(di), msgLvl(ml), maxNumRefinementSteps(mr), transposeMatrix(tm),
  solverMode(sm)
{
  static_assert(sizeof(MKL_INT) == sizeof(int));

  // initialize Pardiso parameters
  maxfct = 1;  // Maximum number of numerical factorizations.
  mnum = 1;    // Which factorization to use.

  if (inputMatrixIsUpper == 0 &&
    mtype != EigenMKLPardisoSupport::MatrixType::REAL_STRUCTURAL_SYM &&
    mtype != EigenMKLPardisoSupport::MatrixType::REAL_UNSYM) {
    mapAMatrix(Ain);
  }

  setParam();

  memset(pointers.data(), 0, sizeof(pointers));
  perm.assign(Ain.rows(), 0);

  n = static_cast<MKL_INT>(Ain.rows());
}

EigenMKLPardisoSupport::~EigenMKLPardisoSupport()
{
  MKL_INT phase = -1;
  MKL_INT mi_mtype = static_cast<MKL_INT>(mtype);
  MKL_INT nrhs = 0;
  MKL_INT error = 0;

  pardiso(pointers.data(), &maxfct, &mnum, &mi_mtype, &phase,
    &n, nullptr, nullptr, nullptr, perm.data(), &nrhs,
    iparm.data(), &msgLvl, nullptr, nullptr, &error);

  if (msgLvl)
    std::cout << "Released pardiso. " << error << std::endl;
}

void EigenMKLPardisoSupport::mapAMatrix(const SpMatD &inA)
{
  std::vector<TripletD> entries;

  for (IDX rowi = 0; rowi < inA.rows(); rowi++) {
    for (SpMatD::InnerIterator it(inA, rowi); it; ++it) {
      if (it.col() >= it.row()) {
        entries.emplace_back((int)it.row(), (int)it.col(), it.value());
      }
    }
  }

  A.resize(inA.rows(), inA.cols());
  A.setFromTriplets(entries.begin(), entries.end());

  big2Small(inA, A, 0, 0, AMapping, 1);
}

void EigenMKLPardisoSupport::setParam()
{
  for (auto &v : iparm)
    v = 0;

  iparm[0] = 1;                            // Do not use the solver default values (use custom values, provided below)
  iparm[1] = static_cast<MKL_INT>(rtype);  // matrix re-ordering algorithm
  iparm[2] = 0;                            // unused

  // use iterative-direct algorithm if requested
  if (directIterative) {
    if (mtype == MatrixType::REAL_SPD)
      iparm[3] = 62;  // matrix is symmetric positive-definite; use CGS iteration for symmetric positive-definite matrices
    else
      iparm[3] = 61;  // use CGS iteration
  }
  else
    iparm[3] = 0;

  iparm[4] = 0;                      // No user fill-in reducing permutation
  iparm[5] = 0;                      // Write solution into x
  iparm[6] = 0;                      // Output: number of iterative refinement steps performed
  iparm[7] = maxNumRefinementSteps;  // Max numbers of iterative refinement steps (used during the solving stage). Value of 0 (default) means: The solver automatically performs two steps of iterative refinement when perturbed pivots are obtained during the numerical factorization.
  iparm[8] = 0;                      // Reserved. Must set to 0.

  // Pivoting perturbation; the values below are Pardiso's default values
  // Pivoting only applies to REAL_UNSYM and REAL_SYM_INDEFINITE
  if (mtype == MatrixType::REAL_UNSYM)
    iparm[9] = 13;  // For non-symmetric matrices, perturb the pivot elements with 1E-13
  else
    iparm[9] = 13;   // Use 1.0E-8 for symmetric indefinite matrices

  // Scaling and matching. The following below are the Pardiso defaults.
  if (mtype == MatrixType::REAL_UNSYM) {
    iparm[10] = 1;  // enable scaling
    iparm[12] = 1;  // enable matching
  }
  else {
    iparm[10] = 1;  // disable scaling
    iparm[12] = 1;  // disable matching
  }
  if (transposeMatrix)
    iparm[11] = 2;  // Solve with transposed or conjugate transposed matrix A. Not in use here.
  else
    iparm[11] = 0;

  iparm[13] = 0;                  // Output: Number of perturbed pivots
  iparm[14] = 0;                  // Output: Peak memory on symbolic factorization (in KB)
  iparm[15] = 0;                  // Output: Permanent memory on symbolic factorization (in KB)
  iparm[16] = 0;                  // Output: Size of factors/Peak memory on numerical factorization and solution (in KB)

  iparm[17] = -1;                 // Output: Report the number of non-zero elements in the factors.
  iparm[18] = 0;                  // Report number of floating point operations (in 10^6 floating point operations) that are necessary to factor the matrix A. Disabled.
  iparm[19] = 0;                  // Output: Report CG/CGS diagnostics.
  iparm[20] = 1;                  // Pivoting for symmetric indefinite matrices: Apply 1x1 and 2x2 Bunch-Kaufman pivoting during the factorization process.
  iparm[21] = 0;                  // Output: Inertia: number of positive eigenvalues.
  iparm[22] = 0;                  // Output: Inertia: number of negative eigenvalues.

  iparm[23] = 0;                 // Parallel factorization control. Use default. 10: Intel oneAPI Math Kernel Library PARDISO uses an improved two-level factorization algorithm for nonsymmetric matrices.
  iparm[24] = 0;                  // Parallel forward/backward solve control. Intel MKL PARDISO uses a parallel algorithm for the solve step.

  iparm[25] = 0;                  // unused
  iparm[26] = 0;                  // matrix checker; 0: not check;
  iparm[27] = 0;                  // 0: double; 1: float;

  iparm[28] = 0;                  // unused;
  iparm[29] = 0;                  // output;
  iparm[30] = 0;                  // partial solve;

  iparm[31] = 0;                  // unused;
  iparm[32] = 0;                  // unused;

  iparm[33] = 0;                  // Optimal number of OpenMP threads for conditional numerical reproducibility (CNR) mode.
  iparm[34] = 1;                  // 0: one-based indexing; 1: zero-based indexing;
  iparm[35] = 0;                  // 0: Do not compute Schur complement.
  iparm[36] = 0;                  // 0: CSR format

  iparm[37] = 0;                  // unsued;
  iparm[38] = 0;                  // Do not use low rank update functionality.;

  iparm[39] = 0;                  // unsued;
  iparm[40] = 0;                  // unsued;
  iparm[41] = 0;                  // unsued;

  iparm[42] = 0;                  // Do not compute the diagonal of inverse matrix.

  for (int i = 43; i <= 54; i++)  // unsed;
    iparm[i] = 0;

  iparm[55] = 0;           // Internal function used to work with pivot and calculation of diagonal arrays turned off.

  iparm[56] = 0;           // unsued;
  iparm[57] = 0;           // unsued;
  iparm[58] = 0;           // unsued;

  iparm[59] = solverMode;  // solver mode: 0: in-core; 2: out-of-core;

  iparm[60] = 0;           // unsued;
  iparm[61] = 0;           // unsued;

  iparm[62] = 0;           // output: Size of the minimum OOC memory for numerical factorization and solution.
  iparm[63] = 0;           // unsued;
}

int EigenMKLPardisoSupport::analyze(const SpMatD &Ain)
{
  MKL_INT phase = 11;
  MKL_INT mi_mtype = static_cast<MKL_INT>(mtype);
  MKL_INT nrhs = 0;
  MKL_INT error = 0;

  hclockPt t1 = hclock::now();

  if (AMapping.nonZeros()) {
    transferBigToSmall(Ain, A, AMapping, 1);
  }

  const SpMatD *Aptr = nullptr;
  if (AMapping.nonZeros()) {
    Aptr = &A;
  }
  else {
    Aptr = &Ain;
  }

  if (msgLvl) {
    std::cout << "Analyzing matrix (phase=11)..." << std::endl;
    std::cout << "Matrix size: " << n << std::endl;
  }

  pardiso(pointers.data(), &maxfct, &mnum, &mi_mtype, &phase,
    &n, Aptr->valuePtr(), Aptr->outerIndexPtr(), Aptr->innerIndexPtr(), perm.data(), &nrhs,
    iparm.data(), &msgLvl, nullptr, nullptr, &error);

  hclockPt t2 = hclock::now();

  if (msgLvl) {
    std::cout << "Time cost: " << dura(t1, t2) << "s" << std::endl;
  }

  if (error != 0) {
    if (msgLvl) {
      std::cerr << getErrorMessage(error) << std::endl;
    }

    return error;
  }

  return 0;
}

int EigenMKLPardisoSupport::factorize(const SpMatD &Ain)
{
  MKL_INT phase = 22;
  MKL_INT mi_mtype = static_cast<MKL_INT>(mtype);
  MKL_INT nrhs = 0;
  MKL_INT error = 0;

  if (directIterative) {
    return 0;
  }

  // compute the factorization
  if (msgLvl) {
    std::cout << "Factorizing matrix (phase=22)..." << std::endl;
    std::cout << "Matrix size: " << n << std::endl;
  }

  hclockPt t1 = hclock::now();

  if (AMapping.nonZeros()) {
    transferBigToSmall(Ain, A, AMapping, 1);
  }

  const SpMatD *Aptr = nullptr;
  if (AMapping.nonZeros()) {
    Aptr = &A;
  }
  else {
    Aptr = &Ain;
  }

  // factorize
  pardiso(pointers.data(), &maxfct, &mnum, &mi_mtype, &phase,
    &n, Aptr->valuePtr(), Aptr->outerIndexPtr(), Aptr->innerIndexPtr(), perm.data(), &nrhs,
    iparm.data(), &msgLvl, nullptr, nullptr, &error);

  hclockPt t2 = hclock::now();

  if (msgLvl) {
    std::cout << "Time cost: " << dura(t1, t2) << "s" << std::endl;
  }

  if (error != 0) {
    if (msgLvl) {
      std::cerr << getErrorMessage(error) << std::endl;
    }

    return error;
  }

  return 0;
}

int EigenMKLPardisoSupport::solve(const SpMatD &Ain, double *x, double *rhs, int nrhs_)
{
  MKL_INT phase = 33;
  MKL_INT mi_mtype = static_cast<MKL_INT>(mtype);
  MKL_INT nrhs = static_cast<MKL_INT>(nrhs_);
  MKL_INT error = 0;

  if (directIterative != 0) {
    phase = 23;
  }

  // compute the factorization
  if (msgLvl) {
    std::cout << "Solving matrix (phase=" << phase << ")..." << std::endl;
    std::cout << "# rhs: " << nrhs << std::endl;
  }

  hclockPt t1 = hclock::now();

  if (AMapping.nonZeros()) {
    transferBigToSmall(Ain, A, AMapping, 1);
  }

  const SpMatD *Aptr = nullptr;
  if (AMapping.nonZeros()) {
    Aptr = &A;
  }
  else {
    Aptr = &Ain;
  }

  pardiso(pointers.data(), &maxfct, &mnum, &mi_mtype, &phase,
    &n, Aptr->valuePtr(), Aptr->outerIndexPtr(), Aptr->innerIndexPtr(), perm.data(), &nrhs,
    iparm.data(), &msgLvl, rhs, x, &error);

  hclockPt t2 = hclock::now();

  if (msgLvl) {
    std::cout << "Time cost: " << dura(t1, t2) << "s" << std::endl;
  }

  if (error != 0) {
    if (msgLvl) {
      std::cerr << getErrorMessage(error) << std::endl;
    }

    return error;
  }

  return 0;
}

int EigenMKLPardisoSupport::solve(double *x, double *rhs, int nrhs_)
{
  MKL_INT phase = 33;
  MKL_INT mi_mtype = static_cast<MKL_INT>(mtype);
  MKL_INT nrhs = static_cast<MKL_INT>(nrhs_);
  MKL_INT error = 0;

  if (directIterative != 0) {
    phase = 23;
  }

  // compute the factorization
  if (msgLvl) {
    std::cout << "Solving matrix (phase=" << phase << ")..." << std::endl;
    std::cout << "# rhs: " << nrhs << std::endl;
  }

  hclockPt t1 = hclock::now();

  const SpMatD *Aptr = &A;
  pardiso(pointers.data(), &maxfct, &mnum, &mi_mtype, &phase,
    &n, Aptr->valuePtr(), Aptr->outerIndexPtr(), Aptr->innerIndexPtr(), perm.data(), &nrhs,
    iparm.data(), &msgLvl, rhs, x, &error);

  hclockPt t2 = hclock::now();

  if (msgLvl) {
    std::cout << "Time cost: " << dura(t1, t2) << "s" << std::endl;
  }

  if (error != 0) {
    if (msgLvl) {
      std::cerr << getErrorMessage(error) << std::endl;
    }

    return error;
  }

  return 0;
}

int EigenMKLPardisoSupport::forward(const SpMatD &Ain, double *x, double *rhs, int nrhs_)
{
  if (directIterative != 0) {
    std::cerr << "Iterative solver. directly exit" << std::endl;
    return 1;
  }

  MKL_INT phase = 331;
  MKL_INT mi_mtype = static_cast<MKL_INT>(mtype);
  MKL_INT nrhs = static_cast<MKL_INT>(nrhs_);
  MKL_INT error = 0;

  // compute the factorization
  if (msgLvl) {
    std::cout << "Solving matrix (phase=" << phase << ")..." << std::endl;
    std::cout << "# rhs: " << nrhs << std::endl;
  }

  hclockPt t1 = hclock::now();

  if (AMapping.nonZeros()) {
    transferBigToSmall(Ain, A, AMapping, 1);
  }

  const SpMatD *Aptr = nullptr;
  if (AMapping.nonZeros()) {
    Aptr = &A;
  }
  else {
    Aptr = &Ain;
  }

  pardiso(pointers.data(), &maxfct, &mnum, &mi_mtype, &phase,
    &n, Aptr->valuePtr(), Aptr->outerIndexPtr(), Aptr->innerIndexPtr(), perm.data(), &nrhs,
    iparm.data(), &msgLvl, rhs, x, &error);

  hclockPt t2 = hclock::now();

  if (msgLvl) {
    std::cout << "Time cost: " << dura(t1, t2) << "s" << std::endl;
  }

  if (error != 0) {
    if (msgLvl) {
      std::cerr << getErrorMessage(error) << std::endl;
    }

    return error;
  }

  return 0;
}

int EigenMKLPardisoSupport::backward(const SpMatD &Ain, double *x, double *rhs, int nrhs_)
{
  if (directIterative != 0) {
    std::cerr << "Iterative solver. directly exit" << std::endl;
    return 1;
  }

  MKL_INT phase = 333;
  MKL_INT mi_mtype = static_cast<MKL_INT>(mtype);
  MKL_INT nrhs = static_cast<MKL_INT>(nrhs_);
  MKL_INT error = 0;

  // compute the factorization
  if (msgLvl) {
    std::cout << "Solving matrix (phase=" << phase << ")..." << std::endl;
    std::cout << "# rhs: " << nrhs << std::endl;
  }

  hclockPt t1 = hclock::now();

  if (AMapping.nonZeros()) {
    transferBigToSmall(Ain, A, AMapping, 1);
  }

  const SpMatD *Aptr = nullptr;
  if (AMapping.nonZeros()) {
    Aptr = &A;
  }
  else {
    Aptr = &Ain;
  }

  pardiso(pointers.data(), &maxfct, &mnum, &mi_mtype, &phase,
    &n, Aptr->valuePtr(), Aptr->outerIndexPtr(), Aptr->innerIndexPtr(), perm.data(), &nrhs,
    iparm.data(), &msgLvl, rhs, x, &error);

  hclockPt t2 = hclock::now();

  if (msgLvl) {
    std::cout << "Time cost: " << dura(t1, t2) << "s" << std::endl;
  }

  if (error != 0) {
    if (msgLvl) {
      std::cerr << getErrorMessage(error) << std::endl;
    }

    return error;
  }

  return 0;
}

std::string EigenMKLPardisoSupport::getErrorMessage(int errorCode) const
{
  auto iter = errorMessages.find(errorCode);
  if (iter != errorMessages.end())
    return iter->second;
  else
    return "Unknown error.";
}

std::map<int, std::string> EigenMKLPardisoSupport::errorMessages = {
  { 0, "Success" },
  { -1, "Input inconsistent" },
  { -2, "Not enough memory" },
  { -3, "Reordering problem" },
  { -4, "Zero pivot, numerical factorization or iterative refinement problem.If the error appears during the solution phase, try to change the pivoting perturbation(iparm[9]) and also increase the number of iterative refinement steps.If it does not help, consider changing the scaling, matchingand pivoting options(iparm[10], iparm[12], iparm[20])" },
  { -5, "Unclassified(internal) error" },
  { -6, "Reordering failed(matrix types 11 and 13 only)" },
  { -7, "Diagonal matrix is singular" },
  { -8, "32-bit integer overflow problem" },
  { -9, "Not enough memory for OOC" },
  { -10, "Error opening OOC files" },
  { -11, "Read/write error with OOC files" },
  { -12, "(pardiso_64 only)pardiso_64 called from 32-bit library" },
  { -13, "Interrupted by the(user-defined) mkl_progress function" },
  { -15, "Internal error which can appear for iparm[23] = 10 and iparm[12] = 1. Try switch matching off(set iparm[12] = 0 and rerun.)" },
};