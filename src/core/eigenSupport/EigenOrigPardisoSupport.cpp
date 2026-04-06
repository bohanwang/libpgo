#if defined(EIGEN_USE_MKL_ALL)
#  undef EIGEN_USE_MKL_ALL
#endif

#include "EigenOrigPardisoSupport.h"

#include <iostream>
#include <chrono>
#include <thread>

using namespace pgo;
using namespace pgo::EigenSupport;

using hclock = std::chrono::high_resolution_clock;
using hclockPt = hclock::time_point;

inline double dura(const hclockPt &t1, const hclockPt &t2)
{
  return static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count()) / 1e6;
}

/* Original PARDISO prototypes (Panua Technologies) */
/* PARDISO prototype. */
extern "C" void pardisoinit(void *, int *, int *, int *, double *, int *);
extern "C" void pardiso(void *, int *, int *, int *, int *, int *,
  double *, int *, int *, int *, int *, int *,
  int *, double *, double *, int *, double *);
extern "C" void pardiso_chkmatrix(int *, int *, double *, int *, int *, int *);
extern "C" void pardiso_chkvec(int *, int *, double *, int *);
extern "C" void pardiso_printstats(int *, int *, double *, int *, int *, int *,
  double *, int *);

EigenOrigPardisoSupport::EigenOrigPardisoSupport(const SpMatD &Ain, MatrixType mt, ReorderingType rt, int di, int ml, int mr,
  int tm, int sm, int inputMatrixIsUpper):
  mtype(mt),
  rtype(rt), directIterative(di), msgLvl(ml), maxNumRefinementSteps(mr), transposeMatrix(tm),
  solverMode(sm)
{
  pt.fill(nullptr);
  iparm.fill(0);
  dparm.fill(0.0);

  maxfct = 1;
  mnum = 1;
  n = static_cast<int>(Ain.rows());

  // For symmetric matrices stored as lower triangular input, extract upper triangular
  if (inputMatrixIsUpper == 0 &&
    mtype != EigenOrigPardisoSupport::MatrixType::REAL_STRUCTURAL_SYM &&
    mtype != EigenOrigPardisoSupport::MatrixType::REAL_UNSYM) {
    mapAMatrix(Ain);
  }

  // Initialize PARDISO handle and set default parameters
  int solver = directIterative ? 1 : 0;
  int error = 0;
  int mi_mtype = static_cast<int>(mtype);

  pardisoinit(pt.data(), &mi_mtype, &solver, iparm.data(), dparm.data(), &error);

  if (error != 0) {
    if (error == -10)
      std::cerr << "[PARDISO] No license file found.\n";
    else if (error == -11)
      std::cerr << "[PARDISO] License is expired.\n";
    else if (error == -12)
      std::cerr << "[PARDISO] Wrong username or hostname.\n";
    else
      std::cerr << "[PARDISO] pardisoinit error: " << error << "\n";
  }
  else {
    if (msgLvl)
      std::cout << "[PARDISO] License check was successful.\n";
  }

  setParam();

  // Build initial 1-based CSR from whichever matrix we will use
  const SpMatD &Aref = AMapping.nonZeros() ? A : Ain;
  buildCSR(Aref);
}

EigenOrigPardisoSupport::~EigenOrigPardisoSupport()
{
  int phase = -1;
  int mi_mtype = static_cast<int>(mtype);
  int nrhs = 0;
  int error = 0;
  int idum = 0;
  double ddum = 0.0;

  pardiso(pt.data(), &maxfct, &mnum, &mi_mtype, &phase,
    &n, &ddum, ia.data(), ja.data(), &idum, &nrhs,
    iparm.data(), &msgLvl, &ddum, &ddum, &error, dparm.data());

  if (msgLvl)
    std::cout << "Released pardiso. " << error << std::endl;
}

void EigenOrigPardisoSupport::mapAMatrix(const SpMatD &inA)
{
  std::vector<TripletD> triplets;

  for (IDX rowi = 0; rowi < inA.rows(); rowi++) {
    for (SpMatD::InnerIterator it(inA, rowi); it; ++it) {
      if (it.col() >= it.row()) {
        triplets.emplace_back((int)it.row(), (int)it.col(), it.value());
      }
    }
  }

  A.resize(inA.rows(), inA.cols());
  A.setFromTriplets(triplets.begin(), triplets.end());

  big2Small(inA, A, 0, 0, AMapping, 1);
}

void EigenOrigPardisoSupport::buildCSR(const SpMatD &inA)
{
  // Build 1-based CSR arrays from an upper-triangular row-major sparse matrix
  ia.clear();
  ja.clear();
  entries.clear();

  ia.reserve(inA.rows() + 1);
  ja.reserve(inA.nonZeros());
  entries.reserve(inA.nonZeros());

  for (IDX row = 0; row < inA.rows(); row++) {
    ia.push_back((int)ja.size() + 1);  // 1-based row pointer
    for (SpMatD::InnerIterator it(inA, row); it; ++it) {
      ja.push_back((int)it.col() + 1);  // 1-based column index
      entries.push_back(it.value());
    }
  }
  ia.push_back((int)ja.size() + 1);
}

void EigenOrigPardisoSupport::updateCSRValues(const SpMatD &inA)
{
  // Update only the values, assuming sparsity pattern is unchanged
  int idx = 0;
  for (IDX row = 0; row < inA.rows(); row++) {
    for (SpMatD::InnerIterator it(inA, row); it; ++it) {
      entries[idx++] = it.value();
    }
  }
}

void EigenOrigPardisoSupport::setParam()
{
  // iparm[0]=1 means we provide our own parameters (already set by pardisoinit default, overriding below)
  iparm[0] = 1;
  iparm[1] = static_cast<int>(rtype);  // matrix re-ordering algorithm
  iparm[2] = 64;                       // number of threads

  iparm[3] = 0;                      // default value for CG
  iparm[4] = 0;                      // No user fill-in permutation
  iparm[5] = 0;                      // Write solution into x
  iparm[6] = 0;                      // Output: number of iterative refinement steps performed
  iparm[7] = maxNumRefinementSteps;  // Max iterative refinement steps
  iparm[8] = 0;                      // Tolerance level for the relative residual.

  if (mtype == MatrixType::REAL_UNSYM)
    iparm[9] = 13;  // Pivot perturbation: 1e-13 for unsymmetric
  else
    iparm[9] = 8;  // Pivot perturbation: 1e-8 for symmetric indefinite

  // 10: Scaling vectors.
  // 12: Improved accuracy using (non-)symmetric weighted matchings.
  if (mtype == MatrixType::REAL_UNSYM) {
    iparm[10] = 1;  // Enable scaling
    iparm[12] = 1;  // Enable matching
  }
  else {
    iparm[10] = 1;  // Scaling
    iparm[12] = 1;  // Matching 1, 2
  }

  if (transposeMatrix)
    iparm[11] = 2;  // Solve with transposed matrix
  else
    iparm[11] = 0;

  iparm[13] = 0;   // Output: number of perturbed pivots
  iparm[14] = 0;   // Output: peak memory on symbolic factorization (KB)
  iparm[15] = 0;   // Output: permanent memory on symbolic factorization (KB)
  iparm[16] = 0;   // Output: peak memory on numerical factorization and solution (KB)
  iparm[17] = -1;  // Output: number of non-zeros in factors (-1 = report)
  iparm[18] = 0;   // Output: number of MFLOPS for factorization (disabled)
  iparm[19] = 0;   // Output: CG/CGS diagnostics
  iparm[20] = 3;   // Pivoting: 1x1 and 2x2 Bunch-Kaufman pivoting for symmetric indefinite
  iparm[21] = 0;   // Output: inertia - number of positive eigenvalues
  iparm[22] = 0;   // Output: inertia - number of negative eigenvalues
  iparm[23] = 1;   // Parallel factorization control
  iparm[24] = 1;   // Parallel forward/backward solve control

  iparm[25] = 0;   // Splitting of Forward/Backward Solve.
  iparm[26] = 1;   // Recomputation of Scaling.
  iparm[27] = 0;   //  Parallel Reordering for METIS
  iparm[28] = 0;   // Switch between 32-bit and 64-bit factorization.
  iparm[29] = 80;  // Control the size of the supernodes.
  iparm[30] = 0;   //  Partial solve for sparse right-hand sides and sparse solution.
  iparm[31] = 0;   // Use the multi-recursive iterative linear solver.
  iparm[32] = 0;   // Determinant of a matrix.
  iparm[33] = 1;   //  Identical solution independent on the number of processors.
  iparm[34] = 0;   // unused
  iparm[35] = 0;   //  Selected inversion for A−1
  iparm[36] = 0;   // selected inversion for A;
  iparm[37] = 0;   // Schur-complement computation.
  iparm[38] = 0;   //  Nonzeros in Schur-complement
  iparm[39] = 0;   // Permutation
  iparm[40] = 0;   // Incremental Update.

  for (int i = 41; i < 49; i++) {
    iparm[i] = 0;
  }

  iparm[49] = 0;  // out of core.
  iparm[50] = 0;  // Use parallel distributed-memory solver.
  iparm[51] = 1;  // Number of compute nodes for the distributed-memory parallel solver.
  iparm[52] = 0;  // Block size
  iparm[53] = 0;  // Dense columns
  for (int i = 54; i < 64; i++) {
    iparm[i] = 0;
  }
}

int EigenOrigPardisoSupport::analyze(const SpMatD &Ain)
{
  hclockPt t1 = hclock::now();

  if (AMapping.nonZeros()) {
    transferBigToSmall(Ain, A, AMapping, 1);
    updateCSRValues(A);
  }
  else {
    updateCSRValues(Ain);
  }

  if (msgLvl) {
    std::cout << "Analyzing matrix (phase=11)..." << std::endl;
    std::cout << "Matrix size: " << n << std::endl;
  }

  int phase = 11;
  int mi_mtype = static_cast<int>(mtype);
  int nrhs = 0;
  int error = 0;
  int idum = 0;
  double ddum = 0.0;

  pardiso(pt.data(), &maxfct, &mnum, &mi_mtype, &phase,
    &n, entries.data(), ia.data(), ja.data(), &idum, &nrhs,
    iparm.data(), &msgLvl, &ddum, &ddum, &error, dparm.data());

  hclockPt t2 = hclock::now();

  if (msgLvl) {
    std::cout << "Time cost: " << dura(t1, t2) << "s" << std::endl;
  }

  if (error != 0) {
    if (msgLvl)
      std::cerr << getErrorMessage(error) << std::endl;
    return error;
  }

  return 0;
}

int EigenOrigPardisoSupport::factorize(const SpMatD &Ain)
{
  if (directIterative) {
    return 0;
  }

  if (msgLvl) {
    std::cout << "Factorizing matrix (phase=22)..." << std::endl;
    std::cout << "Matrix size: " << n << std::endl;
  }

  hclockPt t1 = hclock::now();

  if (AMapping.nonZeros()) {
    transferBigToSmall(Ain, A, AMapping, 1);
    updateCSRValues(A);
  }
  else {
    updateCSRValues(Ain);
  }

  int phase = 22;
  int mi_mtype = static_cast<int>(mtype);
  int nrhs = 0;
  int error = 0;
  int idum = 0;
  double ddum = 0.0;

  pardiso(pt.data(), &maxfct, &mnum, &mi_mtype, &phase,
    &n, entries.data(), ia.data(), ja.data(), &idum, &nrhs,
    iparm.data(), &msgLvl, &ddum, &ddum, &error, dparm.data());

  hclockPt t2 = hclock::now();

  if (msgLvl) {
    std::cout << "Time cost: " << dura(t1, t2) << "s" << std::endl;
  }

  if (error != 0) {
    if (msgLvl)
      std::cerr << getErrorMessage(error) << std::endl;
    return error;
  }

  return 0;
}

int EigenOrigPardisoSupport::solve(const SpMatD &Ain, double *x, double *rhs, int nrhs_)
{
  int phase = directIterative ? 23 : 33;
  int mi_mtype = static_cast<int>(mtype);
  int nrhs = nrhs_;
  int error = 0;
  int idum = 0;

  if (msgLvl) {
    std::cout << "Solving matrix (phase=" << phase << ")..." << std::endl;
    std::cout << "# rhs: " << nrhs << std::endl;
  }

  hclockPt t1 = hclock::now();

  if (AMapping.nonZeros()) {
    transferBigToSmall(Ain, A, AMapping, 1);
    updateCSRValues(A);
  }
  else {
    updateCSRValues(Ain);
  }

  pardiso(pt.data(), &maxfct, &mnum, &mi_mtype, &phase,
    &n, entries.data(), ia.data(), ja.data(), &idum, &nrhs,
    iparm.data(), &msgLvl, rhs, x, &error, dparm.data());

  hclockPt t2 = hclock::now();

  if (msgLvl)
    std::cout << "Time cost: " << dura(t1, t2) << "s" << std::endl;

  if (error != 0) {
    if (msgLvl)
      std::cerr << getErrorMessage(error) << std::endl;
    return error;
  }

  return 0;
}

int EigenOrigPardisoSupport::solve(double *x, double *rhs, int nrhs_)
{
  int phase = directIterative ? 23 : 33;
  int mi_mtype = static_cast<int>(mtype);
  int nrhs = nrhs_;
  int error = 0;
  int idum = 0;

  if (msgLvl) {
    std::cout << "Solving matrix (phase=" << phase << ")..." << std::endl;
    std::cout << "# rhs: " << nrhs << std::endl;
  }

  hclockPt t1 = hclock::now();

  pardiso(pt.data(), &maxfct, &mnum, &mi_mtype, &phase,
    &n, entries.data(), ia.data(), ja.data(), &idum, &nrhs,
    iparm.data(), &msgLvl, rhs, x, &error, dparm.data());

  hclockPt t2 = hclock::now();

  if (msgLvl)
    std::cout << "Time cost: " << dura(t1, t2) << "s" << std::endl;

  if (error != 0) {
    if (msgLvl)
      std::cerr << getErrorMessage(error) << std::endl;
    return error;
  }

  return 0;
}

int EigenOrigPardisoSupport::forward(const SpMatD &Ain, double *x, double *rhs, int nrhs_)
{
  if (directIterative) {
    std::cerr << "Iterative solver. directly exit" << std::endl;
    return 1;
  }

  int phase = 331;
  int mi_mtype = static_cast<int>(mtype);
  int nrhs = nrhs_;
  int error = 0;
  int idum = 0;

  if (msgLvl) {
    std::cout << "Solving matrix (phase=" << phase << ")..." << std::endl;
    std::cout << "# rhs: " << nrhs << std::endl;
  }

  hclockPt t1 = hclock::now();

  if (AMapping.nonZeros()) {
    transferBigToSmall(Ain, A, AMapping, 1);
    updateCSRValues(A);
  }
  else {
    updateCSRValues(Ain);
  }

  pardiso(pt.data(), &maxfct, &mnum, &mi_mtype, &phase,
    &n, entries.data(), ia.data(), ja.data(), &idum, &nrhs,
    iparm.data(), &msgLvl, rhs, x, &error, dparm.data());

  hclockPt t2 = hclock::now();

  if (msgLvl)
    std::cout << "Time cost: " << dura(t1, t2) << "s" << std::endl;

  if (error != 0) {
    if (msgLvl)
      std::cerr << getErrorMessage(error) << std::endl;
    return error;
  }

  return 0;
}

int EigenOrigPardisoSupport::backward(const SpMatD &Ain, double *x, double *rhs, int nrhs_)
{
  if (directIterative) {
    std::cerr << "Iterative solver. directly exit" << std::endl;
    return 1;
  }

  int phase = 333;
  int mi_mtype = static_cast<int>(mtype);
  int nrhs = nrhs_;
  int error = 0;
  int idum = 0;

  if (msgLvl) {
    std::cout << "Solving matrix (phase=" << phase << ")..." << std::endl;
    std::cout << "# rhs: " << nrhs << std::endl;
  }

  hclockPt t1 = hclock::now();

  if (AMapping.nonZeros()) {
    transferBigToSmall(Ain, A, AMapping, 1);
    updateCSRValues(A);
  }
  else {
    updateCSRValues(Ain);
  }

  pardiso(pt.data(), &maxfct, &mnum, &mi_mtype, &phase,
    &n, entries.data(), ia.data(), ja.data(), &idum, &nrhs,
    iparm.data(), &msgLvl, rhs, x, &error, dparm.data());

  hclockPt t2 = hclock::now();

  if (msgLvl)
    std::cout << "Time cost: " << dura(t1, t2) << "s" << std::endl;

  if (error != 0) {
    if (msgLvl)
      std::cerr << getErrorMessage(error) << std::endl;
    return error;
  }

  return 0;
}

std::string EigenOrigPardisoSupport::getErrorMessage(int errorCode) const
{
  auto iter = errorMessages.find(errorCode);
  if (iter != errorMessages.end())
    return iter->second;
  else
    return "Unknown error.";
}

std::map<int, std::string> EigenOrigPardisoSupport::errorMessages = {
  { 0, "Success" },
  { -1, "Input inconsistent" },
  { -2, "Not enough memory" },
  { -3, "Reordering problem" },
  { -4, "Zero pivot, numerical factorization or iterative refinement problem" },
  { -5, "Unclassified (internal) error" },
  { -6, "Reordering failed (matrix types 11 and 13 only)" },
  { -7, "Diagonal matrix is singular" },
  { -8, "32-bit integer overflow problem" },
  { -9, "Not enough memory for OOC" },
  { -10, "Error opening OOC files" },
  { -11, "Read/write error with OOC files" },
  { -12, "pardiso_64 called from 32-bit library" },
  { -13, "Interrupted by the user-defined progress function" },
};
