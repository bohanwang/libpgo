#pragma once

#include "EigenSupport.h"

#include <array>
#include <map>
#include <string>

namespace pgo
{
namespace EigenSupport
{
class EigenOrigPardisoSupport
{
public:
  // The constructor initializes the PARDISO handle and performs symbolic factorization.
  // Only the topology of A matters for the constructor. A is not modified.
  // Note: after calling the constructor, you must call "factorize" to perform numerical factorization.
  //  "mtype" gives the matrix type:
  //  = 1   structurally symmetric matrix
  //  = 2   symmetric positive-definite matrix
  //  = -2  symmetric indefinite matrix
  //  = 11  unsymmetric matrix
  enum class MatrixType : int
  {
    REAL_STRUCTURAL_SYM = 1,
    REAL_SPD = 2,
    REAL_SYM_INDEFINITE = -2,
    REAL_UNSYM = 11
  };

  // Matrix re-ordering is specified as follows:
  // = 0   minimum degree ordering
  // = 2   nested dissection algorithm from the METIS package
  // = 3   parallel (OpenMP) version of nested dissection
  enum class ReorderingType : int
  {
    MINIMUM_DEGREE_ORDERING = 0,
    NESTED_DISSECTION_4 = 2,
    NESTED_DISSECTION_5 = 3,
    AMD = 4
  };

  // must have: numThreads >= 1
  // "directIterative" specifies whether the multi-recursive iterative solver is used (solver=1)
  EigenOrigPardisoSupport(const SpMatD &A, MatrixType mtype = MatrixType::REAL_SYM_INDEFINITE,
    ReorderingType rtype = ReorderingType::NESTED_DISSECTION_4, 
    int directIterative = 0, int msgLevel = 0, int maxNumRefinement = 0,
    int transposeMatrix = 0, int solverMode = 0, int inputMatrixIsUpper = 0);
  EigenOrigPardisoSupport(const EigenOrigPardisoSupport &other) = delete;
  EigenOrigPardisoSupport(EigenOrigPardisoSupport &&other) = delete;

  ~EigenOrigPardisoSupport();

  void setMessageLevel(int lvl) { msgLvl = lvl; }
  std::array<int, 64> &getiparam() { return iparm; }

  // Reproducible solves. PARDISO's user guide (IPARM(34)) states that the
  // solver computes the identical solution independent of the number of
  // processors when IPARM(34) = 1, the METIS reordering (IPARM(2) = 2) is
  // used together with the two-level parallel factorization (IPARM(24) = 1)
  // and the parallel solve (IPARM(25) = 1), and IPARM(3) holds the upper
  // limit of processors for which the results must coincide; the runtime
  // may then use any OMP_NUM_THREADS up to that limit. Sequential METIS
  // reordering (IPARM(28) = 0) keeps the ordering reproducible as well.
  // `threadLimit` is that upper limit; it is raised to the OpenMP team size
  // when the runtime would otherwise use more threads than the limit, since
  // IPARM(3) below the team size corrupts the factorization. Off: IPARM(3)
  // follows the OpenMP team size and IPARM(34) = 0. The default is on with a
  // limit of 64. Call before analyze().
  void setDeterministic(bool on, int threadLimit = 64);
  bool isDeterministic() const { return deterministic; }
  int getThreadLimit() const { return threadLimit; }

  int analyze(const SpMatD &A);
  int factorize(const SpMatD &A);

  // solve: A * x = rhs, using the previously computed matrix factorization
  // rhs is not modified
  int solve(double *x, double *rhs, int nrhs);
  int solve(const SpMatD &A, double *x, double *rhs, int nrhs);

  int forward(const SpMatD &A, double *x, double *rhs, int nrhs);
  int backward(const SpMatD &A, double *x, double *rhs, int nrhs);

protected:
  void setParam();
  void applyParallelParams();
  void mapAMatrix(const SpMatD &inA);
  void buildCSR(const SpMatD &inA);
  void updateCSRValues(const SpMatD &inA);
  std::string getErrorMessage(int errorCode) const;

  std::array<int, 64> iparm;
  std::array<double, 64> dparm;
  std::array<void *, 64> pt;

  MatrixType mtype;
  ReorderingType rtype;
  int directIterative;
  int msgLvl;
  int maxNumRefinementSteps;
  int transposeMatrix;
  int solverMode;
  bool deterministic = true;
  int threadLimit = 64;

  int maxfct = 1;
  int mnum = 1;
  int n;

  // Internal 1-based CSR storage (original PARDISO uses Fortran 1-based indexing)
  std::vector<int> ia, ja;
  std::vector<double> entries;

  // Mapping from full matrix to upper-triangular stored matrix
  SpMatD A;
  SpMatI AMapping;

  static std::map<int, std::string> errorMessages;
};
}  // namespace EigenSupport
}  // namespace pgo
