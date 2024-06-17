/*
author: Bohan Wang
copyright to USC, MIT
*/

#pragma once

#include "EigenDef.h"

#include <initializer_list>
#include <vector>
#include <array>
#include <stdexcept>
#include <utility>
#include <fstream>
#include <iostream>

#if defined(PGO_HAS_MKL)
extern "C"
{
  struct sparse_matrix;
  typedef struct sparse_matrix *sparse_matrix_t;
}
#endif

namespace pgo
{
namespace EigenSupport
{

// ======================================================================
// for dense matrix
void tensorProduct(M3d &r, const V3d &v1, const V3d &v2);
M3d tensorProduct(const V3d &v1, const V3d &v2);
template<typename DerivedA, typename DerivedB, typename DerivedC>
void tensorProduct(Eigen::MatrixBase<DerivedC> &r, const Eigen::MatrixBase<DerivedA> &v1, const Eigen::MatrixBase<DerivedB> &v2);
M3d skewMatrix(const V3d &v);

inline int ELT(int rows, int i, int j)
{
  return j * rows + i;
}

// ======================================================================
// for sparse matrix
#if defined(PGO_HAS_MKL)
struct SparseMatrixMKL
{
  sparse_matrix_t mklHandle;
};
#endif

SpMatD createSparseMatrix(ConstRefMatXd A);
void createWeightMatrix(int numRows, int numCols, int numVertices, int numWeightsPerVertex,
  const int *vertexIndices, const int *interpolationVertexIndices, const double *interpolationWeights, int inflate, SpMatD &A);
SpMatD createWeightMatrix(int numRows, int numCols, int numVertices, int numWeightsPerVertex,
  const int *vertexIndices, const int *interpolationVertexIndices, const double *interpolationWeights, int inflate);

#if defined(PGO_HAS_MKL)
SparseMatrixMKL createSparseMatrixMKL(SpMatD &A);
void destroySparseMatrixMKL(SparseMatrixMKL &m);
#endif

void expand3(const SpMatD &A, SpMatD &A3d);
void expandN(const SpMatD &A, SpMatD &A3d, int nd);

template<typename... T>
SpMatD concatCol(T &&...Ms);

template<typename... T>
SpMatD concatRow(T &&...Ms);

template<typename Iterator>
SpMatD concatCol(Iterator beg, Iterator end);

template<typename Iterator>
SpMatD concatRow(Iterator beg, Iterator end);

template<typename... T>
void mergeSparseMatrix(SpMatD &Adst, T &&...Ms);

template<typename MatrixType>
void zero(MatrixType &m);

// IO
int saveSparseMatrix(const char *filename, const SpMatD &A);
int loadSparseMatrix(const char *filename, SpMatD &A);

// utility
std::ptrdiff_t findEntryOffset(const SpMatD &A, int row, int col);
void buildEntryMap(const SpMatD &A, EntryMap &entryMap);

// transpose
void transposeMapping(const SpMatD &A, const SpMatD &AT, SpMatI &mapping);
void transposeTransfer(const SpMatD &A, SpMatD &AT, const SpMatI &mapping, int parallel = 0);

// mapping
void small2Big(const SpMatD &Asmall, const SpMatD &Abig, int rowStart, int colStart, SpMatI &mapping);
void small2Big(const SpMatD &Asmall, const SpMatD &Abig, const std::vector<int> &dofs, SpMatI &mapping);
void small2Big(const SpMatD &Asmall, const SpMatD &Abig, const std::vector<int> &dofsRow, const std::vector<int> &dofsCol, SpMatI &mapping);

void big2Small(const SpMatD &Abig, const SpMatD &Asmall, int rowStart, int colStart, SpMatI &mapping, int allowDifferentTopo = 0);
void big2Small(const SpMatD &Abig, const SpMatD &Asmall, const std::vector<int> &dofsRow, const std::vector<int> &dofsCol, SpMatI &mapping, int allowDifferentTopo = 0);

void removeRowsCols(const SpMatD &Abig, const SpMatD &Asmall, const std::vector<int> &removedDofs, SpMatI &mapping);
void removeRowsCols(const SpMatD &Abig, const std::vector<int> &removedDofs, SpMatD &mat);
void removeRowsCols(const SpMatD &Abig, const std::vector<int> &removedRowDofs, const std::vector<int> &removedColDofs, SpMatD &mat);

void selectRowsCols(const SpMatD &Abig, const std::vector<int> &selectedDofs, SpMatD &mat);
void selectRowsCols(const SpMatD &Abig, const std::vector<int> &selectedRowDofs, const std::vector<int> &selectedColDofs, SpMatD &mat);

void removeRows(const SpMatD &Abig, const SpMatD &Asmall, const std::vector<int> &removedDofs, SpMatI &mapping);
void removeRows(int dofAll, const std::vector<int> &removedDofs, std::vector<int> &big2smallMapping, std::vector<int> &small2bigMapping);

void addSmallToBig(double coeff, const SpMatD &Asmall, SpMatD &Abig, double beta, const SpMatI &mapping, int parallel = 0);
void addBigToSmall(double coeff, const SpMatD &Abig, SpMatD &Asmall, double beta, const SpMatI &mapping, int parallel = 0);

void transferSmallToBig(const SpMatD &Asmall, SpMatD &Abig, const SpMatI &mapping, int parallel = 0);
void transferSmallToBig(ConstRefVecXd Asmall, RefVecXd Abig, const std::vector<int> &mapping, int parallel = 0);
void transferBigToSmall(const SpMatD &Abig, SpMatD &Asmall, const SpMatI &mapping, int parallel = 0);
void transferBigToSmall(ConstRefVecXd Abig, RefVecXd Asmall, const std::vector<int> &mapping, int parallel = 0);

// sparse operation
void mv(const SpMatD &A, ConstRefVecXd v, RefVecXd b, int transpose = 0);
void mv(const SpMatD &A, ConstRefVecXd v, RefVecXd b, double alpha, double beta, int transpose = 0);

#if defined(PGO_HAS_MKL)
void mv(SparseMatrixMKL A, ConstRefVecXd v, RefVecXd b, int transpose = 0);
void mv(SparseMatrixMKL A, ConstRefVecXd v, RefVecXd b, double alpha, double beta, int transpose = 0);
#endif

void mm(const SpMatD &A, ConstRefMatXd B, RefMatXd C, int transpose = 0);
void mm(const SpMatD &A, ConstRefMatXd B, RefMatXd C, double alpha, double beta, int transpose = 0);

#if defined(PGO_HAS_MKL)
void mm(SparseMatrixMKL A, ConstRefMatXd B, RefMatXd C, int transpose = 0);
void mm(SparseMatrixMKL A, ConstRefMatXd B, RefMatXd C, double alpha, double beta, int transpose = 0);
void mm(SparseMatrixMKL A, SparseMatrixMKL B, SparseMatrixMKL C, int transpose = 0);
#endif

void mm(const SpMatD &A, const SpMatD &B, SpMatD &C, int transpose = 0);

void symbolicMm(const SpMatD &A, const SpMatD &B, SpMatD &C, SymbolicMmData **dat, int transpose = 0);
void mm(const SpMatD &A, const SpMatD &B, const SymbolicMmData *dat, SpMatD &C, int transpose = 0);
void destroySymbolicMmData(SymbolicMmData *dat);

// In the following three routines,
// when transpose=1 (default), it means doing A^T M A,
// when transpose=0, it means doing A M A^T,
// The routine may not be fast compared to Mm if B is too big
void aTa(const SpMatD &A, SpMatD &C, int transpose = 1);
// B must be symmetric, and can optionally be represented by upper triangle part
void aba(const SpMatD &A, const SpMatD &B, SpMatD &C, int transpose = 1);
void aba(const SpMatD &A, const double *W, SpMatD &C, int transpose = 1);

double vTMv(const SpMatD &A, ConstRefVecXd v, RefVecXd b, int transpose = 0);

void testSparseRoutine(int dim1, int dim2);

// dense utilities
M3Xd toM3Xd(const double *data, size_t size);
MX3d toMX3d(const double *data, size_t size);
VXd toVXd(const double *data, size_t size);

template<typename Derived>
V3d rotationToVec3(const Eigen::MatrixBase<Derived> &R);

void rotationOnlySVD(const M3d &A, M3d *UOut, M3d *VOut, V3d *SOut);

void rigidFit(const M3Xd &p0, const M3Xd &p1, M3d &R, V3d &t, const double *weights = nullptr);

template<typename Derived>
int readMatrix(const char *filename, Eigen::PlainObjectBase<Derived> &mat, int overwrite = 1);
template<typename Derived>
int writeMatrix(const char *filename, const Eigen::PlainObjectBase<Derived> &mat);
}  // namespace EigenSupport

template<typename Derived>
EigenSupport::V3d EigenSupport::rotationToVec3(const Eigen::MatrixBase<Derived> &R)
{
  Eigen::Quaterniond q;
  q = R;

  double cosThetaHalf = q.w();
  double angle = 2 * acos(cosThetaHalf);

  V3d axis(q.x(), q.y(), q.z());
  axis.normalize();

  return axis * angle;
}

template<typename... T>
void EigenSupport::mergeSparseMatrix(SpMatD &Adst, T &&...Ms)
{
  static_assert(
    std::is_same<
      std::integer_sequence<bool, true, std::is_same<std::decay_t<T>, SpMatD>::value...>,
      std::integer_sequence<bool, std::is_same<std::decay_t<T>, SpMatD>::value..., true>>::value,
    "T should be SpMatD");

  std::array<const SpMatD *, sizeof...(Ms)> matrices = { &Ms... };

  std::vector<TripletD> entries;
  Eigen::Index numRows = 0, numCols = 0;
  for (const SpMatD *mptr : matrices) {
    const SpMatD &m = *mptr;
    if (numRows == 0 && numCols == 0) {
      numRows = m.rows();
      numCols = m.cols();
    }

    if (numRows != m.rows() || numCols != m.cols())
      throw std::domain_error("Different sparse matrix size");

    for (Eigen::Index outeri = 0; outeri < m.outerSize(); outeri++) {
      for (SpMatD::InnerIterator it(m, outeri); it; ++it) {
        entries.emplace_back((int)it.row(), (int)it.col(), 1.0);
      }
    }
  }

  Adst.resize(numRows, numCols);
  Adst.setFromTriplets(entries.begin(), entries.end());
}

template<typename... T>
EigenSupport::SpMatD EigenSupport::concatCol(T &&...Ms)
{
  static_assert(
    std::is_same<
      std::integer_sequence<bool, true, std::is_same<std::decay_t<T>, SpMatD>::value...>,
      std::integer_sequence<bool, std::is_same<std::decay_t<T>, SpMatD>::value..., true>>::value,
    "T should be SpMatD");

  std::array<const SpMatD *, sizeof...(Ms)> matrices = { &Ms... };

  std::vector<TripletD> entries;
  Eigen::Index numRows = 0, numCols = 0;
  for (const SpMatD *mptr : matrices) {
    const SpMatD &m = *mptr;
    if (numRows == 0) {
      numRows = m.rows();
    }

    if (numRows != m.rows())
      throw std::domain_error("Different sparse matrix rows");

    for (Eigen::Index outeri = 0; outeri < m.outerSize(); outeri++) {
      for (SpMatD::InnerIterator it(m, outeri); it; ++it) {
        entries.emplace_back((int)it.row(), (int)it.col() + (int)numCols, 1.0);
      }
    }

    numCols += (int)m.cols();
  }

  SpMatD A(numRows, numCols);
  A.setFromTriplets(entries.begin(), entries.end());

  return A;
}

template<typename Iterator>
EigenSupport::SpMatD EigenSupport::concatCol(Iterator beg, Iterator end)
{
  std::vector<TripletD> entries;
  Eigen::Index numRows = 0, numCols = 0;
  for (Iterator it = beg; it != end; ++it) {
    const SpMatD &m = *it;
    if (numRows == 0) {
      numRows = m.rows();
    }

    if (numRows != m.rows())
      throw std::invalid_argument("Different sparse matrix rows");

    for (Eigen::Index outeri = 0; outeri < m.outerSize(); outeri++) {
      for (SpMatD::InnerIterator it(m, outeri); it; ++it) {
        entries.emplace_back((int)it.row(), (int)it.col() + (int)numCols, 1.0);
      }
    }

    numCols += (int)m.cols();
  }

  SpMatD A(numRows, numCols);
  A.setFromTriplets(entries.begin(), entries.end());

  return A;
}

template<typename... T>
EigenSupport::SpMatD EigenSupport::concatRow(T &&...Ms)
{
  static_assert(
    std::is_same<
      std::integer_sequence<bool, true, std::is_same<std::decay_t<T>, SpMatD>::value...>,
      std::integer_sequence<bool, std::is_same<std::decay_t<T>, SpMatD>::value..., true>>::value,
    "T should be SpMatD");

  std::array<const SpMatD *, sizeof...(Ms)> matrices = { &Ms... };

  std::vector<TripletD> entries;
  Eigen::Index numRows = 0, numCols = 0;
  for (const SpMatD *mptr : matrices) {
    const SpMatD &m = *mptr;
    if (numCols == 0) {
      numCols = m.cols();
    }

    if (numCols != m.cols())
      throw std::domain_error("Different sparse matrix cols");

    for (Eigen::Index outeri = 0; outeri < m.outerSize(); outeri++) {
      for (SpMatD::InnerIterator it(m, outeri); it; ++it) {
        entries.emplace_back((int)it.row() + (int)numRows, (int)it.col(), 1.0);
      }
    }

    numRows += (int)m.rows();
  }

  SpMatD A(numRows, numCols);
  A.setFromTriplets(entries.begin(), entries.end());

  return A;
}

template<typename Iterator>
EigenSupport::SpMatD EigenSupport::concatRow(Iterator beg, Iterator end)
{
  std::vector<TripletD> entries;
  Eigen::Index numRows = 0, numCols = 0;
  for (Iterator it = beg; it != end; ++it) {
    const SpMatD &m = *it;
    if (numCols == 0) {
      numCols = m.cols();
    }

    if (numCols != m.cols())
      throw std::domain_error("Different sparse matrix cols");

    for (Eigen::Index outeri = 0; outeri < m.outerSize(); outeri++) {
      for (SpMatD::InnerIterator it(m, outeri); it; ++it) {
        entries.emplace_back((int)it.row() + (int)numRows, (int)it.col(), 1.0);
      }
    }

    numRows += (int)m.rows();
  }

  SpMatD A(numRows, numCols);
  A.setFromTriplets(entries.begin(), entries.end());

  return A;
}

template<>
inline void EigenSupport::zero(EigenSupport::SpMatD &m)
{
  Eigen::Index nnz = m.nonZeros();
  memset(m.valuePtr(), 0, sizeof(double) * nnz);
}

template<typename MatrixType>
void EigenSupport::zero(MatrixType &m)
{
  m.setZero();
  abort();
}

template<typename DerivedA, typename DerivedB, typename DerivedC>
void EigenSupport::tensorProduct(Eigen::MatrixBase<DerivedC> &ret, const Eigen::MatrixBase<DerivedA> &v1, const Eigen::MatrixBase<DerivedB> &v2)
{
  if (v1.rows() != v2.rows() || v1.cols() != 1 || v2.cols() != 1 ||
    ret.rows() != v1.rows() || ret.cols() != ret.rows())
    throw std::domain_error("Wrong size of the arguments");

  for (Eigen::Index c = 0; c < ret.cols(); c++) {
    for (Eigen::Index r = 0; r < ret.rows(); r++) {
      ret(r, c) = v1(r, 0) * v2(c, 0);
    }
  }
}

template<typename Derived>
int EigenSupport::readMatrix(const char *filename, Eigen::PlainObjectBase<Derived> &mat, int overwrite)
{
  std::ifstream infile(filename, std::ios_base::binary);
  if (!infile)
    return 1;

  int nRows, nCols, entrySize;
  if (!infile.read((char *)&nRows, sizeof(int))) {
    std::cerr << "Cannot read #rows." << std::endl;
    return 1;
  }

  if (!infile.read((char *)&nCols, sizeof(int))) {
    std::cerr << "Cannot read #cols." << std::endl;
    return 1;
  }

  if (!infile.read((char *)&entrySize, sizeof(int))) {
    std::cerr << "Cannot read entry size." << std::endl;
    return 1;
  }

  if (entrySize != int(sizeof(typename Eigen::PlainObjectBase<Derived>::Scalar))) {
    std::cerr << "entry size mis matched." << std::endl;
  }

  if (overwrite) {
    mat.resize(nRows, nCols);
    if (!infile.read((char *)mat.data(), sizeof(typename Eigen::PlainObjectBase<Derived>::Scalar) * mat.size())) {
      std::cerr << "Cannot read data" << std::endl;
      return 1;
    }

    return 0;
  }
  else {
    if (nRows != int(mat.rows()) || nCols != int(mat.cols())) {
      std::cerr << "#rows or #cols mismatched. File is (" << nRows << ',' << nCols << ") while input is (" << mat.rows() << ',' << mat.cols() << ")" << std::endl;
      return 1;
    }

    if (!infile.read((char *)mat.data(), sizeof(typename Eigen::PlainObjectBase<Derived>::Scalar) * mat.size())) {
      return 1;
    }

    return 0;
  }
}

template<typename Derived>
int EigenSupport::writeMatrix(const char *filename, const Eigen::PlainObjectBase<Derived> &mat)
{
  std::ofstream outfile(filename, std::ios_base::binary);
  if (!outfile)
    return 1;

  int nRows = int(mat.rows()), nCols = int(mat.cols());
  int entrySize = int(sizeof(typename Eigen::PlainObjectBase<Derived>::Scalar));
  if (!outfile.write((const char *)&nRows, sizeof(int))) {
    std::cerr << "Cannot write #rows." << std::endl;
    return 1;
  }

  if (!outfile.write((const char *)&nCols, sizeof(int))) {
    std::cerr << "Cannot write #cols." << std::endl;
    return 1;
  }

  if (!outfile.write((const char *)&entrySize, sizeof(int))) {
    std::cerr << "Cannot write #cols." << std::endl;
    return 1;
  }

  if (!outfile.write((char *)mat.data(), sizeof(typename Eigen::PlainObjectBase<Derived>::Scalar) * mat.size())) {
    return 1;
  }

  return 0;
}

}  // namespace pgo