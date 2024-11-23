#pragma once

#include "EigenSupport.h"

#include <mkl_pardiso.h>

#include <array>
#include <map>
#include <string>

namespace pgo
{
namespace EigenSupport
{
class EigenMKLPardisoSupport
{
public:
  // The constructor computes the permutation to re-order A, and performs symbolic factorization.
  // Only the topology of A matters for the constructor. A is not modified.
  // Note: after calling the constructor, you must call "FactorMatrix" to perform numerical factorization.
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
  // = 3   parallel (OpenMP) version of nested dissection; it can decrease the computation time on multi-core computers, especially when the constructor takes a long time
  enum class ReorderingType : int
  {
    MINIMUM_DEGREE_ORDERING = 0,
    NESTED_DISSECTION = 2,
    PARALLEL_NESTED_DISSECTION = 3
  };

  // must have: numThreads >= 1
  // "directIterative" specifies whether a direct-iterative procedure is used (see Intel MKL's documentation)
  EigenMKLPardisoSupport(const SpMatD &A, MatrixType mtype = MatrixType::REAL_SYM_INDEFINITE,
    ReorderingType rtype = ReorderingType::NESTED_DISSECTION, int directIterative = 0, int msgLevel = 0, int maxNumRefinement = 0,
    int transposeMatrix = 0, int solverMode = 0, int inputMatrixIsUpper = 0);
  EigenMKLPardisoSupport(const EigenMKLPardisoSupport &other) = delete;
  EigenMKLPardisoSupport(EigenMKLPardisoSupport &&other) = delete;

  ~EigenMKLPardisoSupport();

  void setMessageLevel(int lvl) { msgLvl = static_cast<MKL_INT>(lvl); }
  std::array<MKL_INT, 64> &getiparam() { return iparm; }

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
  void mapAMatrix(const SpMatD &inA);
  std::string getErrorMessage(int errorCode) const;

  std::array<MKL_INT, 64> iparm;
  std::array<void *, 64> pointers;

  std::vector<MKL_INT> perm;

  MatrixType mtype;
  ReorderingType rtype;
  MKL_INT directIterative;
  MKL_INT msgLvl;
  MKL_INT maxNumRefinementSteps;
  MKL_INT transposeMatrix;
  MKL_INT solverMode;

  MKL_INT maxfct = 1;
  MKL_INT mnum = 1;
  MKL_INT n;

  SpMatD A;
  SpMatI AMapping;

  static std::map<int, std::string> errorMessages;
};
}  // namespace EigenSupport
}  // namespace pgo