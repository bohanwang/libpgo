/*
author: Bohan Wang
copyright to USC, MIT
*/

#ifdef __GNUC__
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wint-in-bool-context"
#endif

#include "EigenSupport.h"

#if defined(PGO_HAS_MKL)
#  include <mkl.h>
#  define MKL_INSPECTOR_EXECUTOR
#endif

#include <tbb/parallel_for.h>
#include <tbb/enumerable_thread_specific.h>

#include <vector>
#include <iostream>
#include <numeric>

#ifndef EIGEN_SUPPORT_INLINE
#  define EIGEN_SUPPORT_INLINE
#endif

EIGEN_SUPPORT_INLINE pgo::EigenSupport::SpMatD pgo::EigenSupport::createSparseMatrix(ConstRefMatXd A)
{
  std::vector<TripletD> entries;
  entries.reserve(A.size());

  for (IDX coli = 0; coli < A.cols(); coli++) {
    for (IDX rowi = 0; rowi < A.rows(); rowi++) {
      entries.emplace_back(rowi, coli, A(rowi, coli));
    }
  }

  SpMatD ret(A.rows(), A.cols());
  ret.setFromTriplets(entries.begin(), entries.end());

  return ret;
}

EIGEN_SUPPORT_INLINE void pgo::EigenSupport::createWeightMatrix(int numRows, int numCols, int numVertices, int numWeightsPerVertex,
  const int *vtxIdx, const int *interpolationVertexIndices, const double *interpolationWeights, int inflate, SpMatD &A)
{
  std::vector<int> vertexIndices;
  if (vtxIdx == nullptr) {
    vertexIndices.resize(numVertices);
    std::iota(vertexIndices.begin(), vertexIndices.end(), 0);
  }
  else {
    vertexIndices.assign(vtxIdx, vtxIdx + numVertices);
  }

  std::vector<TripletD> entries;
  if (inflate) {
    for (int vi = 0; vi < numVertices; vi++) {
      for (int vj = 0; vj < numWeightsPerVertex; vj++) {
        entries.emplace_back(vertexIndices[vi] * 3, interpolationVertexIndices[vi * numWeightsPerVertex + vj] * 3, interpolationWeights[vi * numWeightsPerVertex + vj]);
        entries.emplace_back(vertexIndices[vi] * 3 + 1, interpolationVertexIndices[vi * numWeightsPerVertex + vj] * 3 + 1, interpolationWeights[vi * numWeightsPerVertex + vj]);
        entries.emplace_back(vertexIndices[vi] * 3 + 2, interpolationVertexIndices[vi * numWeightsPerVertex + vj] * 3 + 2, interpolationWeights[vi * numWeightsPerVertex + vj]);
      }
    }
  }
  else {
    for (int vi = 0; vi < numVertices; vi++) {
      for (int vj = 0; vj < numWeightsPerVertex; vj++) {
        entries.emplace_back(vertexIndices[vi], interpolationVertexIndices[vi * numWeightsPerVertex + vj], interpolationWeights[vi * numWeightsPerVertex + vj]);
        // std::cout << entries.size() << ',' << vertexIndices[vi] << ',' << interpolationVertexIndices[vi * numWeightsPerVertex + vj] << ',' << interpolationWeights[vi * numWeightsPerVertex + vj] << '\n';
      }
    }
  }
  // std::cout << entries.size() << std::endl;

  A.resize(numRows, numCols);
  A.setFromTriplets(entries.begin(), entries.end());
}

EIGEN_SUPPORT_INLINE pgo::EigenSupport::SpMatD pgo::EigenSupport::createWeightMatrix(int numRows, int numCols, int numVertices, int numWeightsPerVertex,
  const int *vtxIdx, const int *interpolationVertexIndices, const double *interpolationWeights, int inflate)
{
  SpMatD W;
  pgo::EigenSupport::createWeightMatrix(numRows, numCols, numVertices, numWeightsPerVertex, vtxIdx, interpolationVertexIndices, interpolationWeights, inflate, W);
  std::cout << W.nonZeros() << ',' << W.rows() << ',' << W.cols() << std::endl;

  return W;
}

#if defined(PGO_HAS_MKL)
EIGEN_SUPPORT_INLINE pgo::EigenSupport::SparseMatrixMKL pgo::EigenSupport::createSparseMatrixMKL(SpMatD &A)
{
  SparseMatrixMKL sm;

  mkl_sparse_d_create_csr(&sm.mklHandle, SPARSE_INDEX_BASE_ZERO, (int)A.rows(), (int)A.cols(),
    (int *)A.outerIndexPtr(), (int *)A.outerIndexPtr() + 1,
    (int *)A.innerIndexPtr(), (double *)A.valuePtr());

  return sm;
}

EIGEN_SUPPORT_INLINE void pgo::EigenSupport::destroySparseMatrixMKL(SparseMatrixMKL &m)
{
  mkl_sparse_destroy(m.mklHandle);
  m.mklHandle = nullptr;
}
#endif

EIGEN_SUPPORT_INLINE void pgo::EigenSupport::expand3(const SpMatD &A, SpMatD &A3d)
{
  std::vector<TripletD> entries;
  entries.reserve(A.nonZeros() * 3);
  for (IDX i = 0; i < A.outerSize(); i++) {
    for (SpMatD::InnerIterator it(A, i); it; ++it) {
      entries.emplace_back((int)it.row() * 3, (int)it.col() * 3, it.value());
      entries.emplace_back((int)it.row() * 3 + 1, (int)it.col() * 3 + 1, it.value());
      entries.emplace_back((int)it.row() * 3 + 2, (int)it.col() * 3 + 2, it.value());
    }
  }

  A3d.resize(A.rows() * 3, A.cols() * 3);
  A3d.setFromTriplets(entries.begin(), entries.end());
}

EIGEN_SUPPORT_INLINE void pgo::EigenSupport::expandN(const SpMatD &A, SpMatD &A3d, int nd)
{
  std::vector<TripletD> entries;
  entries.reserve(A.nonZeros() * nd);
  for (IDX i = 0; i < A.outerSize(); i++) {
    for (SpMatD::InnerIterator it(A, i); it; ++it) {
      for (int dof = 0; dof < nd; dof++) {
        entries.emplace_back((int)it.row() * nd + dof, (int)it.col() * nd + dof, it.value());
      }
    }
  }

  A3d.resize(A.rows() * nd, A.cols() * nd);
  A3d.setFromTriplets(entries.begin(), entries.end());
}

EIGEN_SUPPORT_INLINE void pgo::EigenSupport::tensorProduct(M3d &r, const V3d &v1, const V3d &v2)
{
  r(0, 0) = v1[0] * v2[0];
  r(0, 1) = v1[0] * v2[1];
  r(0, 2) = v1[0] * v2[2];

  r(1, 0) = v1[1] * v2[0];
  r(1, 1) = v1[1] * v2[1];
  r(1, 2) = v1[1] * v2[2];

  r(2, 0) = v1[2] * v2[0];
  r(2, 1) = v1[2] * v2[1];
  r(2, 2) = v1[2] * v2[2];
}

EIGEN_SUPPORT_INLINE pgo::EigenSupport::M3d pgo::EigenSupport::tensorProduct(const V3d &v1, const V3d &v2)
{
  M3d r;

  r(0, 0) = v1[0] * v2[0];
  r(0, 1) = v1[0] * v2[1];
  r(0, 2) = v1[0] * v2[2];

  r(1, 0) = v1[1] * v2[0];
  r(1, 1) = v1[1] * v2[1];
  r(1, 2) = v1[1] * v2[2];

  r(2, 0) = v1[2] * v2[0];
  r(2, 1) = v1[2] * v2[1];
  r(2, 2) = v1[2] * v2[2];

  return r;
}

EIGEN_SUPPORT_INLINE pgo::EigenSupport::M3d pgo::EigenSupport::skewMatrix(const V3d &v)
{
  M3d skew;
  skew << 0, -v[2], v[1],
    v[2], 0, -v[0],
    -v[1], v[0], 0;

  return skew;
}

namespace pgo::EigenSupport
{
#if defined(PGO_HAS_MKL)

EIGEN_SUPPORT_INLINE sparse_matrix_t fromEigenMatrix(const SpMatD &A)
{
  // The input arrays provided are left unchanged except for the call to mkl_sparse_order,
  sparse_matrix_t AM;
  sparse_status_t ret;
  ret = mkl_sparse_d_create_csr(&AM, SPARSE_INDEX_BASE_ZERO, (int)A.rows(), (int)A.cols(),
    (int *)A.outerIndexPtr(), (int *)A.outerIndexPtr() + 1,
    (int *)A.innerIndexPtr(), (double *)A.valuePtr());
  if (ret != SPARSE_STATUS_SUCCESS) {
    std::cout << "Encounter sparse matrix creation Error: " << ret << std::endl;
    abort();
  }

  return AM;
}

EIGEN_SUPPORT_INLINE void toEigenMatrix(sparse_matrix_t CM, SpMatD &C, int onlyUpper = 0)
{
  int rows;
  int cols;
  int *rowStart;
  int *rowEnd;
  int *colIndics;
  double *values;
  sparse_index_base_t indexing;

  sparse_status_t ret = mkl_sparse_d_export_csr(CM, &indexing, &rows, &cols, &rowStart, &rowEnd, &colIndics, &values);
  if (ret != SPARSE_STATUS_SUCCESS) {
    abort();
  }

  std::vector<TripletD> entries;
  entries.reserve(rowEnd[rows - 1]);

  int offset = 0;
  if (indexing == SPARSE_INDEX_BASE_ONE)
    offset = 1;

  int maxColIndex = 0;
  for (int r = 0; r < rows; r++) {
    int row = r;
    for (int c = rowStart[r]; c < rowEnd[r]; c++) {
      int col = colIndics[c] - offset;
      double val = values[c];

      entries.emplace_back(row, col, val);
      maxColIndex = std::max(maxColIndex, col);

      if (onlyUpper && row != col) {
        entries.emplace_back(col, row, val);
      }
    }
  }
  // std::cout << "Max col index: " << maxColIndex << "; return # cols:" << cols << std::endl;
  // if (maxColIndex + 1 < cols)
  //   cols = maxColIndex + 1;
  C.resize(rows, cols);
  C.setFromTriplets(entries.begin(), entries.end());
}
#endif
}  // namespace pgo::EigenSupport

EIGEN_SUPPORT_INLINE void pgo::EigenSupport::mm(const SpMatD &A, const SpMatD &B, SpMatD &C, int transpose)
{
#if defined(MKL_INSPECTOR_EXECUTOR)
  sparse_operation_t opt = SPARSE_OPERATION_NON_TRANSPOSE;
  if (transpose)
    opt = SPARSE_OPERATION_TRANSPOSE;

  sparse_matrix_t AM = fromEigenMatrix(A);
  sparse_matrix_t BM = fromEigenMatrix(B);

  sparse_matrix_t CM;
  sparse_status_t ret = mkl_sparse_spmm(opt, AM, BM, &CM);
  if (ret != SPARSE_STATUS_SUCCESS) {
    std::cout << "Encounter sparse matrix multiplication Error: " << ret << std::endl;
    abort();
  }

  toEigenMatrix(CM, C);

  mkl_sparse_destroy(CM);
  mkl_sparse_destroy(BM);
  mkl_sparse_destroy(AM);
#else
  if (transpose)
    C = A.transpose() * B;
  else
    C = A * B;
#endif
}

EIGEN_SUPPORT_INLINE void pgo::EigenSupport::mv(const SpMatD &A, ConstRefVecXd v, RefVecXd b, int transpose)
{
#if defined(MKL_INSPECTOR_EXECUTOR)
  sparse_matrix_t AM = fromEigenMatrix(A);
  sparse_operation_t opt = SPARSE_OPERATION_NON_TRANSPOSE;
  if (transpose)
    opt = SPARSE_OPERATION_TRANSPOSE;

  matrix_descr desc;
  desc.type = SPARSE_MATRIX_TYPE_GENERAL;
  desc.mode = SPARSE_FILL_MODE_LOWER;
  desc.diag = SPARSE_DIAG_NON_UNIT;

  mkl_sparse_d_mv(opt, 1.0, AM, desc, v.data(), 0.0, b.data());

  mkl_sparse_destroy(AM);
#else
  if (transpose) {
    b.noalias() = A.transpose() * v;
  }
  else {
    b.noalias() = A * v;
  }
#endif
}

EIGEN_SUPPORT_INLINE void pgo::EigenSupport::mv(const SpMatD &A, ConstRefVecXd v, RefVecXd b, double alpha, double beta, int transpose)
{
#if defined(MKL_INSPECTOR_EXECUTOR)
  sparse_matrix_t AM = fromEigenMatrix(A);

  sparse_operation_t opt = SPARSE_OPERATION_NON_TRANSPOSE;
  if (transpose)
    opt = SPARSE_OPERATION_TRANSPOSE;

  matrix_descr desc;
  desc.type = SPARSE_MATRIX_TYPE_GENERAL;
  desc.mode = SPARSE_FILL_MODE_LOWER;
  desc.diag = SPARSE_DIAG_NON_UNIT;

  mkl_sparse_d_mv(opt, alpha, AM, desc, v.data(), beta, b.data());

  mkl_sparse_destroy(AM);
#else
  if (transpose) {
    b *= beta;
    b += A.transpose() * v * alpha;
  }
  else {
    b *= beta;
    b += A * v * alpha;
  }
#endif
}

#if defined(PGO_HAS_MKL)
EIGEN_SUPPORT_INLINE void pgo::EigenSupport::mv(SparseMatrixMKL A, ConstRefVecXd v, RefVecXd b, int transpose)
{
  sparse_operation_t opt = SPARSE_OPERATION_NON_TRANSPOSE;
  if (transpose)
    opt = SPARSE_OPERATION_TRANSPOSE;

  matrix_descr desc;
  desc.type = SPARSE_MATRIX_TYPE_GENERAL;
  desc.mode = SPARSE_FILL_MODE_LOWER;
  desc.diag = SPARSE_DIAG_NON_UNIT;

  mkl_sparse_d_mv(opt, 1.0, A.mklHandle, desc, v.data(), 0.0, b.data());
}

EIGEN_SUPPORT_INLINE void pgo::EigenSupport::mv(SparseMatrixMKL A, ConstRefVecXd v, RefVecXd b, double alpha, double beta, int transpose)
{
  sparse_operation_t opt = SPARSE_OPERATION_NON_TRANSPOSE;
  if (transpose)
    opt = SPARSE_OPERATION_TRANSPOSE;

  matrix_descr desc;
  desc.type = SPARSE_MATRIX_TYPE_GENERAL;
  desc.mode = SPARSE_FILL_MODE_LOWER;
  desc.diag = SPARSE_DIAG_NON_UNIT;

  mkl_sparse_d_mv(opt, alpha, A.mklHandle, desc, v.data(), beta, b.data());
}
#endif

EIGEN_SUPPORT_INLINE void pgo::EigenSupport::mm(const SpMatD &A, ConstRefMatXd B, RefMatXd C, int transpose)
{
#if defined(MKL_INSPECTOR_EXECUTOR)
  // std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

  sparse_matrix_t AM = fromEigenMatrix(A);

  sparse_operation_t opt = SPARSE_OPERATION_NON_TRANSPOSE;
  if (transpose)
    opt = SPARSE_OPERATION_TRANSPOSE;

  matrix_descr desc;
  desc.type = SPARSE_MATRIX_TYPE_GENERAL;
  desc.mode = SPARSE_FILL_MODE_LOWER;
  desc.diag = SPARSE_DIAG_NON_UNIT;

  mkl_sparse_d_mm(opt, 1.0, AM, desc, SPARSE_LAYOUT_COLUMN_MAJOR, B.data(), (int)B.cols(), (int)B.rows(), 0.0, C.data(), (int)C.rows());

  mkl_sparse_destroy(AM);

  // std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

#  if 0
  if (transpose) {
    MXd C1 = C;

    std::chrono::high_resolution_clock::time_point t3 = std::chrono::high_resolution_clock::now();

    std::vector<const int *> colEndPtr(A.rows());
    for (IDX row = 0; row < A.rows(); row++) {
      int rowEnd = A.outerIndexPtr()[row + 1];
      colEndPtr[row] = A.innerIndexPtr() + rowEnd;
    }

    using Buf = std::vector<const int*, tbb::cache_aligned_allocator<const int*>>;
    tbb::enumerable_thread_specific<Buf> tls(Buf(A.rows()));

    tbb::parallel_for((IDX)0, B.cols(),
      [&](int coli) {
        auto &colStartPtr = tls.local();

        for (IDX row = 0; row < A.rows(); row++) {
          int rowStart = A.outerIndexPtr()[row];
          colStartPtr[row] = A.innerIndexPtr() + rowStart;
        }

        for (IDX colj = 0; colj < A.cols(); colj++) {
          double sum = 0;

          for (IDX row = 0; row < A.rows(); row++) {
            while (colStartPtr[row] < colEndPtr[row] && *colStartPtr[row] < colj)
              colStartPtr[row]++;

            if (colStartPtr[row] < colEndPtr[row] && *colStartPtr[row] == colj) {
              std::ptrdiff_t offset = colStartPtr[row] - A.innerIndexPtr();
              sum += B(row, coli) * A.valuePtr()[offset];
            }
          }
          C1(colj, coli) = sum;
        }
      });

    std::chrono::high_resolution_clock::time_point t4 = std::chrono::high_resolution_clock::now();

    std::cout << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count()
              << "," << std::chrono::duration_cast<std::chrono::microseconds>(t4 - t3).count() << std::endl;

    std::cout << "###########################" << A.cols() << ',' << B.cols() << std::endl;
    std::cin.get();

    std::cout << (C - C1).norm() / C.norm() << std::endl;
  }
#  endif

#else
  if (transpose) {
    C.noalias() = A.transpose() * B;
  }
  else {
    C.noalias() = A * B;
  }
#endif
}

EIGEN_SUPPORT_INLINE void pgo::EigenSupport::mm(const SpMatD &A, ConstRefMatXd B, RefMatXd C, double alpha, double beta, int transpose)
{
#if defined(MKL_INSPECTOR_EXECUTOR)
  sparse_matrix_t AM = fromEigenMatrix(A);

  sparse_operation_t opt = SPARSE_OPERATION_NON_TRANSPOSE;
  if (transpose)
    opt = SPARSE_OPERATION_TRANSPOSE;

  matrix_descr desc;
  desc.type = SPARSE_MATRIX_TYPE_GENERAL;
  desc.mode = SPARSE_FILL_MODE_LOWER;
  desc.diag = SPARSE_DIAG_NON_UNIT;

  mkl_sparse_d_mm(opt, alpha, AM, desc, SPARSE_LAYOUT_COLUMN_MAJOR, B.data(), (int)B.cols(), (int)B.rows(), beta, C.data(), (int)C.rows());

  mkl_sparse_destroy(AM);
#else
  if (transpose) {
    C *= beta;
    C += A.transpose() * B * alpha;
  }
  else {
    C *= beta;
    C += A * B * alpha;
  }
#endif
}

#if defined(PGO_HAS_MKL)
EIGEN_SUPPORT_INLINE void pgo::EigenSupport::mm(SparseMatrixMKL A, ConstRefMatXd B, RefMatXd C, int transpose)
{
  sparse_operation_t opt = SPARSE_OPERATION_NON_TRANSPOSE;
  if (transpose)
    opt = SPARSE_OPERATION_TRANSPOSE;

  matrix_descr desc;
  desc.type = SPARSE_MATRIX_TYPE_GENERAL;
  desc.mode = SPARSE_FILL_MODE_LOWER;
  desc.diag = SPARSE_DIAG_NON_UNIT;

  mkl_sparse_d_mm(opt, 1.0, A.mklHandle, desc, SPARSE_LAYOUT_COLUMN_MAJOR, B.data(), (int)B.cols(), (int)B.rows(), 0.0, C.data(), (int)C.rows());
}

EIGEN_SUPPORT_INLINE void pgo::EigenSupport::mm(SparseMatrixMKL A, ConstRefMatXd B, RefMatXd C, double alpha, double beta, int transpose)
{
  sparse_operation_t opt = SPARSE_OPERATION_NON_TRANSPOSE;
  if (transpose)
    opt = SPARSE_OPERATION_TRANSPOSE;

  matrix_descr desc;
  desc.type = SPARSE_MATRIX_TYPE_GENERAL;
  desc.mode = SPARSE_FILL_MODE_LOWER;
  desc.diag = SPARSE_DIAG_NON_UNIT;

  mkl_sparse_d_mm(opt, alpha, A.mklHandle, desc, SPARSE_LAYOUT_COLUMN_MAJOR, B.data(), (int)B.cols(), (int)B.rows(), beta, C.data(), (int)C.rows());
}
#endif

EIGEN_SUPPORT_INLINE double pgo::EigenSupport::vTMv(const SpMatD &A, ConstRefVecXd v, RefVecXd b, int transpose)
{
#if defined(MKL_INSPECTOR_EXECUTOR)
  sparse_matrix_t AM = fromEigenMatrix(A);

  sparse_operation_t opt = SPARSE_OPERATION_NON_TRANSPOSE;
  if (transpose)
    opt = SPARSE_OPERATION_TRANSPOSE;

  matrix_descr desc;
  desc.type = SPARSE_MATRIX_TYPE_GENERAL;
  desc.mode = SPARSE_FILL_MODE_LOWER;
  desc.diag = SPARSE_DIAG_NON_UNIT;

  double ret;
  mkl_sparse_d_dotmv(opt, 1.0, AM, desc, v.data(), 0.0, b.data(), &ret);

  mkl_sparse_destroy(AM);

  return ret;
#else
  if (transpose) {
    b.noalias() = A.transpose() * v;
  }
  else {
    b.noalias() = A * v;
  }

  return v.dot(b);
#endif
}

#if defined(PGO_HAS_MKL)
EIGEN_SUPPORT_INLINE void pgo::EigenSupport::mm(SparseMatrixMKL A, SparseMatrixMKL B, SparseMatrixMKL C, int transpose)
{
  sparse_operation_t opt = SPARSE_OPERATION_NON_TRANSPOSE;
  if (transpose)
    opt = SPARSE_OPERATION_TRANSPOSE;

  sparse_status_t ret = mkl_sparse_spmm(opt, A.mklHandle, B.mklHandle, &C.mklHandle);
  if (ret != SPARSE_STATUS_SUCCESS) {
    abort();
  }
}
#endif

EIGEN_SUPPORT_INLINE void pgo::EigenSupport::aba(const SpMatD &A, const SpMatD &B, SpMatD &C, int transpose)
{
#if defined(MKL_INSPECTOR_EXECUTOR)
  sparse_operation_t opt = SPARSE_OPERATION_NON_TRANSPOSE;
  if (transpose)
    opt = SPARSE_OPERATION_TRANSPOSE;

  sparse_status_t ret;
  sparse_matrix_t AM = fromEigenMatrix(A);
  sparse_matrix_t BM = fromEigenMatrix(B);

  matrix_descr descB;
  descB.type = SPARSE_MATRIX_TYPE_SYMMETRIC;
  descB.diag = SPARSE_DIAG_NON_UNIT;
  descB.mode = SPARSE_FILL_MODE_UPPER;

  sparse_matrix_t CM;
  ret = mkl_sparse_sypr(opt, AM, BM, descB, &CM, SPARSE_STAGE_FULL_MULT);
  if (ret != SPARSE_STATUS_SUCCESS) {
    std::cout << "Encounter sparse matrix multiplication Error: " << ret << std::endl;
    abort();
  }
  toEigenMatrix(CM, C, 1);

  mkl_sparse_destroy(CM);
  mkl_sparse_destroy(BM);
  mkl_sparse_destroy(AM);
#else
  if (transpose)
    C = A.transpose() * B * A;
  else
    C = A * B * A.transpose();
#endif
}

EIGEN_SUPPORT_INLINE void pgo::EigenSupport::aba(const SpMatD &A, const double *W, SpMatD &C, int transpose)
{
#if defined(MKL_INSPECTOR_EXECUTOR)
  sparse_operation_t opt = SPARSE_OPERATION_NON_TRANSPOSE;
  if (transpose)
    opt = SPARSE_OPERATION_TRANSPOSE;

  sparse_status_t ret;
  sparse_matrix_t AM = fromEigenMatrix(A);

  IDX diagSize = transpose ? A.rows() : A.cols();
  std::vector<int> rowStart(diagSize + 1);
  std::iota(rowStart.begin(), rowStart.end(), 0);

  sparse_matrix_t BM;
  ret = mkl_sparse_d_create_csr(&BM, SPARSE_INDEX_BASE_ZERO, (int)diagSize, (int)diagSize,
    rowStart.data(), rowStart.data() + 1,
    rowStart.data(), (double *)W);
  if (ret != SPARSE_STATUS_SUCCESS) {
    std::cout << "Encounter sparse matrix creation Error: " << ret << std::endl;
    abort();
  }

  matrix_descr descB;
  descB.type = SPARSE_MATRIX_TYPE_SYMMETRIC;
  descB.diag = SPARSE_DIAG_NON_UNIT;
  descB.mode = SPARSE_FILL_MODE_UPPER;

  sparse_matrix_t CM;
  ret = mkl_sparse_sypr(opt, AM, BM, descB, &CM, SPARSE_STAGE_FULL_MULT);
  if (ret != SPARSE_STATUS_SUCCESS) {
    std::cout << "Encounter sparse matrix multiplication Error: " << ret << std::endl;
    abort();
  }
  toEigenMatrix(CM, C, 1);

  mkl_sparse_destroy(CM);
  mkl_sparse_destroy(BM);
  mkl_sparse_destroy(AM);
#else
  abort();
#endif
}

EIGEN_SUPPORT_INLINE void pgo::EigenSupport::aTa(const SpMatD &A, SpMatD &C, int transpose)
{
#if defined(MKL_INSPECTOR_EXECUTOR)
  sparse_operation_t opt = SPARSE_OPERATION_NON_TRANSPOSE;
  if (transpose)
    opt = SPARSE_OPERATION_TRANSPOSE;

  sparse_matrix_t AM = fromEigenMatrix(A);
  sparse_matrix_t CM;
  sparse_status_t ret = mkl_sparse_syrk(opt, AM, &CM);

  if (ret != SPARSE_STATUS_SUCCESS) {
    std::cout << "Encounter sparse matrix multiplication Error: " << ret << std::endl;
    abort();
  }
  toEigenMatrix(CM, C, 1);

  mkl_sparse_destroy(CM);
  mkl_sparse_destroy(AM);
#else
  if (transpose) {
    C = A * A.transpose();
  }
  else {
    C = A.transpose() * A;
  }
#endif
}

EIGEN_SUPPORT_INLINE void pgo::EigenSupport::testSparseRoutine(int dim1, int dim2)
{
  std::cout << "========================================\n";
  std::cout << "Testing dim " << dim1 << 'x' << dim2 << "..." << std::endl;
  MXd A1 = MXd::Random(dim1, dim2);

  MXd B1 = MXd::Random(dim2, dim2);
  B1 = ((B1.transpose() + B1) * 0.5).eval();

  MXd B2 = MXd::Random(dim1, dim1);
  B2 = ((B2.transpose() + B2) * 0.5).eval();

  VXd v1 = VXd::Random(dim2);
  VXd v2 = VXd::Random(dim1);

  SpMatD A1sp = createSparseMatrix(A1);
  SpMatD B1sp = createSparseMatrix(B1);
  SpMatD B2sp = createSparseMatrix(B2);

  // std::cout << B2 << '\n' << std::endl;
  // std::cout << MXd(B2sp) << std::endl;

  std::cout << "Testing Asp x Bsp" << std::endl;
  SpMatD Csp;
  mm(A1sp, B1sp, Csp);
  std::cout << (MXd(Csp) - A1 * B1).norm() << std::endl;

  std::cout << "Testing Asp^T x Bsp" << std::endl;
  Csp.setZero();
  mm(A1sp, B2sp, Csp, 1);
  std::cout << (MXd(Csp) - A1.transpose() * B2).norm() << std::endl;

  std::cout << "Testing Asp x v" << std::endl;
  VXd b1(dim1);
  mv(A1sp, v1, b1);
  std::cout << (A1 * v1 - b1).norm() << std::endl;

  mv(A1sp, v1, b1, 1.0, 1.0);
  std::cout << (A1 * v1 * 2 - b1).norm() << std::endl;

  std::cout << "Testing Asp^T x v" << std::endl;
  VXd b2(dim2);
  mv(A1sp, v2, b2, 1);
  std::cout << (A1.transpose() * v2 - b2).norm() << std::endl;

  mv(A1sp, v2, b2, 1.0, 1.0, 1);
  std::cout << (A1.transpose() * v2 * 2 - b2).norm() << std::endl;

  std::cout << "Testing Asp x B" << std::endl;
  MXd C1(dim1, dim2);
  mm(A1sp, B1, C1);
  std::cout << (A1 * B1 - C1).norm() << std::endl;

  mm(A1sp, B1, C1, 1.0, 1.0);
  std::cout << (A1 * B1 * 2 - C1).norm() << std::endl;

  std::cout << "Testing Asp^T x B" << std::endl;
  MXd C2(dim2, dim1);
  mm(A1sp, B2, C2, 1);
  std::cout << (A1.transpose() * B2 - C2).norm() << std::endl;

  mm(A1sp, B2, C2, 1.0, 1.0, 1);
  std::cout << (A1.transpose() * B2 * 2 - C2).norm() << std::endl;

  std::cout << "Testing A^T A" << std::endl;
  Csp.setZero();
  aTa(A1sp, Csp, 1);
  std::cout << (MXd(Csp) - A1.transpose() * A1).norm() << std::endl;

  std::cout << "Testing A A^T" << std::endl;
  Csp.setZero();
  aTa(A1sp, Csp, 0);
  std::cout << (MXd(Csp) - A1 * A1.transpose()).norm() << std::endl;

  SpMatD Z2 = B2sp.triangularView<Eigen::Upper>();

  std::cout << "Testing A^T B A" << std::endl;
  Csp.setZero();
  aba(A1sp, B2sp, Csp, 1);
  std::cout << (MXd(Csp) - A1.transpose() * B2 * A1).norm() << std::endl;

  std::cout << "Using upper B only" << std::endl;
  Csp.setZero();
  aba(A1sp, Z2, Csp, 1);
  std::cout << (MXd(Csp) - A1.transpose() * B2 * A1).norm() << std::endl;

  std::cout << "Testing A B A^T" << std::endl;
  Csp.setZero();
  aba(A1sp, B1sp, Csp, 0);
  std::cout << (MXd(Csp) - A1 * B1 * A1.transpose()).norm() << std::endl;

  std::cout << "Testing A^T diag(v) A" << std::endl;
  Csp.setZero();
  aba(A1sp, v2.data(), Csp, 1);
  std::cout << (MXd(Csp) - A1.transpose() * v2.asDiagonal() * A1).norm() << std::endl;

  std::cout << "Testing A diag(v) A^T" << std::endl;
  Csp.setZero();
  aba(A1sp, v1.data(), Csp, 0);
  std::cout << (MXd(Csp) - A1 * v1.asDiagonal() * A1.transpose()).norm() << std::endl;

  std::cout << "Testing vT A v" << std::endl;

  MXd D = MXd::Random(dim1, dim1);
  SpMatD Dsp = createSparseMatrix(D);

  std::cout << vTMv(Dsp, v2, b1, 0) << "?=" << v2.dot(D * v2) << std::endl;
  std::cout << vTMv(Dsp, v2, b1, 1) << "?=" << v2.dot(D.transpose() * v2) << std::endl;
}

EIGEN_SUPPORT_INLINE void pgo::EigenSupport::transposeMapping(const SpMatD &A, const SpMatD &AT, SpMatI &mapping)
{
  std::vector<TripletI> entries;

  for (Eigen::Index outeri = 0; outeri < AT.outerSize(); outeri++) {
    for (SpMatD::InnerIterator it(AT, outeri); it; ++it) {
      int transposedRow = (int)it.col();
      int transposedCol = (int)it.row();

      const SpMatD::StorageIndex *rowStart = A.innerIndexPtr() + A.outerIndexPtr()[transposedRow];
      const SpMatD::StorageIndex *rowEnd = A.innerIndexPtr() + A.outerIndexPtr()[transposedRow + 1];

      auto itt = std::lower_bound(rowStart, rowEnd, transposedCol);

      if (itt != rowEnd && *itt == transposedCol) {
        Eigen::Index offset = itt - A.innerIndexPtr();
        entries.emplace_back((SpMatD::StorageIndex)it.row(), (SpMatD::StorageIndex)it.col(), offset);
      }
      else {
        throw std::domain_error("Different sparse matrix topology");
      }
    }
  }

  mapping.resize(AT.rows(), AT.cols());
  mapping.setFromTriplets(entries.begin(), entries.end());
}

EIGEN_SUPPORT_INLINE void pgo::EigenSupport::transposeTransfer(const SpMatD &A, SpMatD &AT, const SpMatI &mapping, int parallel)
{
  Eigen::Index nnz = AT.nonZeros();
  if (parallel) {
    tbb::parallel_for((Eigen::Index)0, nnz, [&](Eigen::Index i) {
      Eigen::Index newOffset = mapping.valuePtr()[i];
      AT.valuePtr()[i] = A.valuePtr()[newOffset]; },
      tbb::static_partitioner());
  }
  else {
    for (Eigen::Index i = 0; i < nnz; i++) {
      Eigen::Index newOffset = mapping.valuePtr()[i];
      AT.valuePtr()[i] = A.valuePtr()[newOffset];
    }
  }
}

EIGEN_SUPPORT_INLINE void pgo::EigenSupport::small2Big(const SpMatD &Asmall, const SpMatD &Abig, const std::vector<int> &dofs, SpMatI &mapping)
{
  // std::vector<TripletI> entries;
  tbb::concurrent_vector<TripletI> entries;
  entries.reserve(Asmall.nonZeros());

  // for (Eigen::Index outeri = 0; outeri < Asmall.outerSize(); outeri++) {
  tbb::parallel_for((IDX)0, Asmall.outerSize(), [&](IDX outeri) {
    for (SpMatD::InnerIterator it(Asmall, outeri); it; ++it) {
      Eigen::Index small_row = it.row();
      Eigen::Index small_col = it.col();

      Eigen::Index big_row = dofs[small_row];
      Eigen::Index big_col = dofs[small_col];

      const SpMatD::StorageIndex *rowStart = Abig.innerIndexPtr() + Abig.outerIndexPtr()[big_row];
      const SpMatD::StorageIndex *rowEnd = Abig.innerIndexPtr() + Abig.outerIndexPtr()[big_row + 1];

      auto itt = std::lower_bound(rowStart, rowEnd, big_col);

      if (itt != rowEnd && *itt == big_col) {
        Eigen::Index offset = itt - Abig.innerIndexPtr();
        entries.emplace_back((SpMatD::StorageIndex)small_row, (SpMatD::StorageIndex)small_col, offset);
      }
      else {
        throw std::domain_error("Different sparse matrix topology");
      }
    } });

  mapping.resize(Asmall.rows(), Asmall.cols());
  mapping.setFromTriplets(entries.begin(), entries.end());
}

EIGEN_SUPPORT_INLINE void pgo::EigenSupport::small2Big(const SpMatD &Asmall, const SpMatD &Abig, const std::vector<int> &dofsRow, const std::vector<int> &dofsCol, SpMatI &mapping)
{
  std::vector<TripletI> entries;
  entries.reserve(Asmall.nonZeros());

  for (Eigen::Index outeri = 0; outeri < Asmall.outerSize(); outeri++) {
    for (SpMatD::InnerIterator it(Asmall, outeri); it; ++it) {
      Eigen::Index small_row = it.row();
      Eigen::Index small_col = it.col();

      Eigen::Index big_row = dofsRow[small_row];
      Eigen::Index big_col = dofsCol[small_col];

      const SpMatD::StorageIndex *rowStart = Abig.innerIndexPtr() + Abig.outerIndexPtr()[big_row];
      const SpMatD::StorageIndex *rowEnd = Abig.innerIndexPtr() + Abig.outerIndexPtr()[big_row + 1];

      auto itt = std::lower_bound(rowStart, rowEnd, big_col);

      if (itt != rowEnd && *itt == big_col) {
        Eigen::Index offset = itt - Abig.innerIndexPtr();
        entries.emplace_back((SpMatD::StorageIndex)small_row, (SpMatD::StorageIndex)small_col, offset);
      }
      else {
        throw std::domain_error("Different sparse matrix topology");
      }
    }
  }

  mapping.resize(Asmall.rows(), Asmall.cols());
  mapping.setFromTriplets(entries.begin(), entries.end());
}

EIGEN_SUPPORT_INLINE void pgo::EigenSupport::small2Big(const SpMatD &Asmall, const SpMatD &Abig, int rowStart, int colStart, SpMatI &mapping)
{
  std::vector<TripletI> entries;
  entries.reserve(Asmall.nonZeros());

  for (Eigen::Index outeri = 0; outeri < Asmall.outerSize(); outeri++) {
    for (SpMatD::InnerIterator it(Asmall, outeri); it; ++it) {
      Eigen::Index row = it.row() + rowStart;
      Eigen::Index col = it.col() + colStart;
      const SpMatD::StorageIndex *rowStart = Abig.innerIndexPtr() + Abig.outerIndexPtr()[row];
      const SpMatD::StorageIndex *rowEnd = Abig.innerIndexPtr() + Abig.outerIndexPtr()[row + 1];

      auto itt = std::lower_bound(rowStart, rowEnd, col);

      if (itt != rowEnd && *itt == col) {
        Eigen::Index offset = itt - Abig.innerIndexPtr();
        entries.emplace_back((SpMatD::StorageIndex)it.row(), (SpMatD::StorageIndex)it.col(), offset);
      }
      else {
        throw std::domain_error("Different sparse matrix topology");
      }
    }
  }

  mapping.resize(Asmall.rows(), Asmall.cols());
  mapping.setFromTriplets(entries.begin(), entries.end());
}

EIGEN_SUPPORT_INLINE void pgo::EigenSupport::big2Small(const SpMatD &Abig, const SpMatD &Asmall, int rowStart, int colStart, SpMatI &mapping, int allowDifferentTopo)
{
  std::vector<TripletI> entries;
  entries.reserve(Abig.nonZeros());

  for (Eigen::Index outeri = 0; outeri < Abig.outerSize(); outeri++) {
    for (SpMatD::InnerIterator it(Abig, outeri); it; ++it) {
      Eigen::Index row = it.row() - rowStart;
      Eigen::Index col = it.col() - colStart;

      if (row >= 0 && col >= 0 &&
        row < Asmall.rows() && col < Asmall.cols()) {
        const SpMatD::StorageIndex *rowStart = Asmall.innerIndexPtr() + Asmall.outerIndexPtr()[row];
        const SpMatD::StorageIndex *rowEnd = Asmall.innerIndexPtr() + Asmall.outerIndexPtr()[row + 1];

        auto itt = std::lower_bound(rowStart, rowEnd, col);
        if (itt != rowEnd && *itt == col) {
          Eigen::Index offset = itt - Asmall.innerIndexPtr();
          entries.emplace_back((SpMatD::StorageIndex)it.row(), (SpMatD::StorageIndex)it.col(), offset);
        }
        else if (allowDifferentTopo) {
          entries.emplace_back((SpMatD::StorageIndex)it.row(), (SpMatD::StorageIndex)it.col(), -1);
        }
        else {
          throw std::domain_error("Different sparse matrix topology");
        }
      }
      else {
        entries.emplace_back((SpMatD::StorageIndex)it.row(), (SpMatD::StorageIndex)it.col(), -1);
      }
    }
  }

  mapping.resize(Abig.rows(), Abig.cols());
  mapping.setFromTriplets(entries.begin(), entries.end());
}

EIGEN_SUPPORT_INLINE void pgo::EigenSupport::big2Small(const SpMatD &Abig, const SpMatD &Asmall, const std::vector<int> &dofsRow, const std::vector<int> &dofsCol, SpMatI &mapping, int allowDifferentTopo)
{
  std::vector<TripletI> entries;
  entries.reserve(Abig.nonZeros());

  std::vector<int> rowMapping(Abig.rows()), colMapping(Abig.cols());
  for (IDX i = 0; i < Abig.rows(); i++) {
    auto iter = std::find(dofsRow.begin(), dofsRow.end(), (int)i);
    if (iter != dofsRow.end() && *iter == i)
      rowMapping[i] = (int)(iter - dofsRow.begin());
    else
      rowMapping[i] = -1;
  }

  for (IDX i = 0; i < Abig.cols(); i++) {
    auto iter = std::find(dofsCol.begin(), dofsCol.end(), (int)i);
    if (iter != dofsCol.end() && *iter == i)
      colMapping[i] = (int)(iter - dofsCol.begin());
    else
      colMapping[i] = -1;
  }

  for (Eigen::Index outeri = 0; outeri < Abig.outerSize(); outeri++) {
    for (SpMatD::InnerIterator it(Abig, outeri); it; ++it) {
      Eigen::Index rowBig = it.row();
      Eigen::Index colBig = it.col();

      Eigen::Index rowSmall = rowMapping[rowBig];
      Eigen::Index colSmall = colMapping[colBig];

      if (rowSmall >= 0 && colSmall >= 0 &&
        rowSmall < Asmall.rows() && colSmall < Asmall.cols()) {
        const SpMatD::StorageIndex *rowStart = Asmall.innerIndexPtr() + Asmall.outerIndexPtr()[rowSmall];
        const SpMatD::StorageIndex *rowEnd = Asmall.innerIndexPtr() + Asmall.outerIndexPtr()[rowSmall + 1];

        auto itt = std::lower_bound(rowStart, rowEnd, colSmall);
        if (itt != rowEnd && *itt == colSmall) {
          Eigen::Index offset = itt - Asmall.innerIndexPtr();
          entries.emplace_back((SpMatD::StorageIndex)rowBig, (SpMatD::StorageIndex)colBig, offset);
        }
        else if (allowDifferentTopo) {
          entries.emplace_back((SpMatD::StorageIndex)it.row(), (SpMatD::StorageIndex)it.col(), -1);
        }
        else {
          throw std::domain_error("Different sparse matrix topology");
        }
      }
      else {
        entries.emplace_back((SpMatD::StorageIndex)rowBig, (SpMatD::StorageIndex)colBig, -1);
      }
    }
  }

  mapping.resize(Abig.rows(), Abig.cols());
  mapping.setFromTriplets(entries.begin(), entries.end());
}

EIGEN_SUPPORT_INLINE void pgo::EigenSupport::removeRowsCols(const SpMatD &Abig, const SpMatD &Asmall, const std::vector<int> &removedDofs, SpMatI &mapping)
{
  std::vector<int> dofMappingsRow(Abig.rows());
  std::vector<int> dofMappingsCol(Abig.cols());
  for (Eigen::Index i = 0, j = 0; i < Abig.rows(); i++) {
    if (std::binary_search(removedDofs.begin(), removedDofs.end(), (int)i) == false) {
      dofMappingsRow[i] = (int)j++;
    }
    else {
      dofMappingsRow[i] = -1;
    }
  }

  for (Eigen::Index i = 0, j = 0; i < Abig.cols(); i++) {
    if (std::binary_search(removedDofs.begin(), removedDofs.end(), (int)i) == false) {
      dofMappingsCol[i] = (int)j++;
    }
    else {
      dofMappingsCol[i] = -1;
    }
  }

  std::vector<TripletI> entries;
  entries.reserve(Abig.nonZeros());

  for (Eigen::Index outeri = 0; outeri < Abig.outerSize(); outeri++) {
    for (SpMatD::InnerIterator it(Abig, outeri); it; ++it) {
      Eigen::Index row = dofMappingsRow[it.row()];
      Eigen::Index col = dofMappingsCol[it.col()];

      if (row >= 0 && col >= 0) {
        const SpMatD::StorageIndex *rowStart = Asmall.innerIndexPtr() + Asmall.outerIndexPtr()[row];
        const SpMatD::StorageIndex *rowEnd = Asmall.innerIndexPtr() + Asmall.outerIndexPtr()[row + 1];

        auto itt = std::lower_bound(rowStart, rowEnd, col);

        if (itt != rowEnd && *itt == col) {
          Eigen::Index offset = itt - Asmall.innerIndexPtr();
          entries.emplace_back((SpMatD::StorageIndex)it.row(),
            (SpMatD::StorageIndex)it.col(), offset);
        }
        else {
          throw std::domain_error("Different sparse matrix topology");
        }
      }
      else {
        entries.emplace_back((SpMatI::StorageIndex)it.row(),
          (SpMatI::StorageIndex)it.col(), -1);
      }
    }
  }

  mapping.resize(Abig.rows(), Abig.cols());
  mapping.setFromTriplets(entries.begin(), entries.end());
}

EIGEN_SUPPORT_INLINE void pgo::EigenSupport::removeRowsCols(const SpMatD &Abig, const std::vector<int> &removedDofs, SpMatD &mat)
{
  removeRowsCols(Abig, removedDofs, removedDofs, mat);
}

EIGEN_SUPPORT_INLINE void pgo::EigenSupport::removeRowsCols(const SpMatD &Abig, const std::vector<int> &removedRowDofs, const std::vector<int> &removedColDofs, SpMatD &mat)
{
  std::vector<int> dofMappingsRow(Abig.rows());
  std::vector<int> dofMappingsCol(Abig.cols());

  int newNumRows = 0;
  for (Eigen::Index i = 0; i < Abig.rows(); i++) {
    if (std::binary_search(removedRowDofs.begin(), removedRowDofs.end(), (int)i) == false) {
      dofMappingsRow[i] = (int)newNumRows++;
    }
    else {
      dofMappingsRow[i] = -1;
    }
  }

  int newNumCols = 0;
  for (Eigen::Index i = 0; i < Abig.cols(); i++) {
    if (std::binary_search(removedColDofs.begin(), removedColDofs.end(), (int)i) == false) {
      dofMappingsCol[i] = (int)newNumCols++;
    }
    else {
      dofMappingsCol[i] = -1;
    }
  }

  std::vector<TripletD> entries;
  entries.reserve(Abig.nonZeros());

  for (Eigen::Index outeri = 0; outeri < Abig.outerSize(); outeri++) {
    for (SpMatD::InnerIterator it(Abig, outeri); it; ++it) {
      Eigen::Index row = dofMappingsRow[it.row()];
      Eigen::Index col = dofMappingsCol[it.col()];

      if (row >= 0 && col >= 0) {
        entries.emplace_back((SpMatI::StorageIndex)row, (SpMatI::StorageIndex)col, it.value());
      }
    }
  }

  mat.resize(newNumRows, newNumCols);
  mat.setFromTriplets(entries.begin(), entries.end());
}

EIGEN_SUPPORT_INLINE void pgo::EigenSupport::selectRowsCols(const SpMatD &Abig, const std::vector<int> &selectedDofs, SpMatD &mat)
{
  selectRowsCols(Abig, selectedDofs, selectedDofs, mat);
}

EIGEN_SUPPORT_INLINE void pgo::EigenSupport::selectRowsCols(const SpMatD &Abig, const std::vector<int> &selectedRowDofs, const std::vector<int> &selectedColDofs, SpMatD &mat)
{
  std::vector<int> dofMappingsRow(Abig.rows());
  std::vector<int> dofMappingsCol(Abig.cols());

  int newNumRows = 0;
  for (Eigen::Index i = 0; i < Abig.rows(); i++) {
    if (std::binary_search(selectedRowDofs.begin(), selectedRowDofs.end(), (int)i)) {
      dofMappingsRow[i] = (int)newNumRows++;
    }
    else {
      dofMappingsRow[i] = -1;
    }
  }

  int newNumCols = 0;
  for (Eigen::Index i = 0; i < Abig.cols(); i++) {
    if (std::binary_search(selectedColDofs.begin(), selectedColDofs.end(), (int)i)) {
      dofMappingsCol[i] = (int)newNumCols++;
    }
    else {
      dofMappingsCol[i] = -1;
    }
  }

  std::vector<TripletD> entries;
  entries.reserve(Abig.nonZeros());

  for (Eigen::Index outeri = 0; outeri < Abig.outerSize(); outeri++) {
    for (SpMatD::InnerIterator it(Abig, outeri); it; ++it) {
      Eigen::Index row = dofMappingsRow[it.row()];
      Eigen::Index col = dofMappingsCol[it.col()];

      if (row >= 0 && col >= 0) {
        entries.emplace_back((SpMatD::StorageIndex)row, (SpMatD::StorageIndex)col, it.value());
      }
    }
  }

  mat.resize(newNumRows, newNumCols);
  mat.setFromTriplets(entries.begin(), entries.end());
}

EIGEN_SUPPORT_INLINE void pgo::EigenSupport::addSmallToBig(double coeff, const SpMatD &Asmall, SpMatD &Abig, double beta, const SpMatI &mapping, int parallel)
{
  Eigen::Index nnz = Asmall.nonZeros();
  if (parallel) {
    tbb::parallel_for((Eigen::Index)0, nnz, [&](Eigen::Index i) {
      Eigen::Index newOffset = mapping.valuePtr()[i];
      Abig.valuePtr()[newOffset] = Abig.valuePtr()[newOffset] * beta + coeff * Asmall.valuePtr()[i]; },
      tbb::static_partitioner());
  }
  else {
    for (Eigen::Index i = 0; i < nnz; i++) {
      Eigen::Index newOffset = mapping.valuePtr()[i];
      Abig.valuePtr()[newOffset] = Abig.valuePtr()[newOffset] * beta + coeff * Asmall.valuePtr()[i];
    }
  }
}

EIGEN_SUPPORT_INLINE void pgo::EigenSupport::addBigToSmall(double coeff, const SpMatD &Abig, SpMatD &Asmall, double beta, const SpMatI &mapping, int parallel)
{
  Eigen::Index nnz = Abig.nonZeros();
  if (parallel) {
    tbb::parallel_for((Eigen::Index)0, nnz, [&](Eigen::Index i) {
      Eigen::Index newOffset = mapping.valuePtr()[i];
      if (newOffset >= 0) {
        Asmall.valuePtr()[newOffset] = Asmall.valuePtr()[newOffset] * beta + coeff * Abig.valuePtr()[i];
      } },
      tbb::static_partitioner());
  }
  else {
    for (Eigen::Index i = 0; i < nnz; i++) {
      Eigen::Index newOffset = mapping.valuePtr()[i];
      if (newOffset >= 0) {
        Asmall.valuePtr()[newOffset] = Asmall.valuePtr()[newOffset] * beta + coeff * Abig.valuePtr()[i];
      }
    }
  }
}

EIGEN_SUPPORT_INLINE void pgo::EigenSupport::transferSmallToBig(const SpMatD &Asmall, SpMatD &Abig, const SpMatI &mapping, int parallel)
{
  Eigen::Index nnz = Asmall.nonZeros();
  if (parallel) {
    tbb::parallel_for((Eigen::Index)0, nnz, [&](Eigen::Index i) {
      Eigen::Index newOffset = mapping.valuePtr()[i];
      Abig.valuePtr()[newOffset] = Asmall.valuePtr()[i]; },
      tbb::static_partitioner());
  }
  else {
    for (Eigen::Index i = 0; i < nnz; i++) {
      Eigen::Index newOffset = mapping.valuePtr()[i];
      Abig.valuePtr()[newOffset] = Asmall.valuePtr()[i];
    }
  }
}

EIGEN_SUPPORT_INLINE void pgo::EigenSupport::transferSmallToBig(ConstRefVecXd Asmall, RefVecXd Abig, const std::vector<int> &mapping, int parallel)
{
  if (parallel) {
    tbb::parallel_for((IDX)0, Asmall.size(), [&](IDX i) {
      Eigen::Index newOffset = mapping[i];
      Abig[newOffset] = Asmall[i]; });
  }
  else {
    for (Eigen::Index i = 0; i < Asmall.size(); i++) {
      Eigen::Index newOffset = mapping[i];
      Abig[newOffset] = Asmall[i];
    }
  }
}

EIGEN_SUPPORT_INLINE void pgo::EigenSupport::transferBigToSmall(const SpMatD &Abig, SpMatD &Asmall, const SpMatI &mapping, int parallel)
{
  Eigen::Index nnz = Abig.nonZeros();
  if (parallel) {
    tbb::parallel_for((Eigen::Index)0, nnz, [&](Eigen::Index i) {
      Eigen::Index newOffset = mapping.valuePtr()[i];
      if (newOffset >= 0) {
        Asmall.valuePtr()[newOffset] = Abig.valuePtr()[i];
      } },
      tbb::static_partitioner());
  }
  else {
    for (Eigen::Index i = 0; i < nnz; i++) {
      Eigen::Index newOffset = mapping.valuePtr()[i];
      if (newOffset >= 0) {
        Asmall.valuePtr()[newOffset] = Abig.valuePtr()[i];
      }
    }
  }
}

EIGEN_SUPPORT_INLINE void pgo::EigenSupport::transferBigToSmall(ConstRefVecXd Abig, RefVecXd Asmall, const std::vector<int> &mapping, int parallel)
{
  if (parallel) {
    tbb::parallel_for((IDX)0, Abig.size(), [&](IDX i) {
      Eigen::Index newOffset = mapping[i];
      if (newOffset >= 0) {
        Asmall[newOffset] = Abig[i];
      } });
  }
  else {
    for (Eigen::Index i = 0; i < Abig.size(); i++) {
      Eigen::Index newOffset = mapping[i];
      if (newOffset >= 0) {
        Asmall[newOffset] = Abig[i];
      }
    }
  }
}

EIGEN_SUPPORT_INLINE void pgo::EigenSupport::removeRows(const SpMatD &Abig, const SpMatD &Asmall, const std::vector<int> &removedDofs, SpMatI &mapping)
{
  std::vector<int> dofMappingsRow(Abig.rows());
  for (Eigen::Index i = 0, j = 0; i < Abig.rows(); i++) {
    if (std::binary_search(removedDofs.begin(), removedDofs.end(), (int)i) == false) {
      dofMappingsRow[i] = (int)j++;
    }
    else {
      dofMappingsRow[i] = -1;
    }
  }

  std::vector<TripletI> entries;
  entries.reserve(Abig.nonZeros());

  for (Eigen::Index outeri = 0; outeri < Abig.outerSize(); outeri++) {
    for (SpMatD::InnerIterator it(Abig, outeri); it; ++it) {
      Eigen::Index row = dofMappingsRow[it.row()];
      Eigen::Index col = it.col();

      if (row >= 0) {
        const SpMatD::StorageIndex *rowStart = Asmall.innerIndexPtr() + Asmall.outerIndexPtr()[row];
        const SpMatD::StorageIndex *rowEnd = Asmall.innerIndexPtr() + Asmall.outerIndexPtr()[row + 1];

        auto itt = std::lower_bound(rowStart, rowEnd, col);

        if (itt != rowEnd && *itt == col) {
          Eigen::Index offset = itt - Asmall.innerIndexPtr();
          entries.emplace_back((SpMatD::StorageIndex)it.row(),
            (SpMatD::StorageIndex)it.col(), offset);
        }
        else {
          throw std::domain_error("Different sparse matrix topology");
        }
      }
      else {
        entries.emplace_back((SpMatI::StorageIndex)it.row(),
          (SpMatI::StorageIndex)it.col(), -1);
      }
    }
  }

  mapping.resize(Abig.rows(), Abig.cols());
  mapping.setFromTriplets(entries.begin(), entries.end());
}

EIGEN_SUPPORT_INLINE void pgo::EigenSupport::removeRows(int dofAll, const std::vector<int> &removedDofs, std::vector<int> &big2smallMapping, std::vector<int> &small2bigMapping)
{
  int inc = 0;
  small2bigMapping.resize(dofAll - removedDofs.size());
  big2smallMapping.resize(dofAll);

  for (int i = 0; i < dofAll; i++) {
    if (std::binary_search(removedDofs.begin(), removedDofs.end(), i)) {
      big2smallMapping[i] = -1;
    }
    else {
      big2smallMapping[i] = inc;
      small2bigMapping[inc] = i;
      inc++;
    }
  }
}

EIGEN_SUPPORT_INLINE std::ptrdiff_t pgo::EigenSupport::findEntryOffset(const SpMatD &A, int row, int col)
{
  const SpMatD::StorageIndex *rowBegin = A.innerIndexPtr() + A.outerIndexPtr()[row];
  const SpMatD::StorageIndex *rowEnd = A.innerIndexPtr() + A.outerIndexPtr()[row + 1];

  auto it = std::lower_bound(rowBegin, rowEnd, col);
  if (it == rowEnd || *it != col) {
    return -1;
    // std::cerr << "Cannot find the corresponding index!";
    // abort();
  }

  return it - A.innerIndexPtr();
}

EIGEN_SUPPORT_INLINE void pgo::EigenSupport::buildEntryMap(const SpMatD &A, EntryMap &entryMap)
{
  entryMap.reserve(A.nonZeros());
  for (Eigen::Index row = 0; row < A.outerSize(); row++) {
    int rowStart = A.outerIndexPtr()[row];
    int rowEnd = A.outerIndexPtr()[row + 1];

    for (int j = rowStart; j < rowEnd; j++) {
      int col = A.innerIndexPtr()[j];

      entryMap.emplace(std::make_pair((int)row, col), j);
    }
  }
}

// dense tilities
EIGEN_SUPPORT_INLINE pgo::EigenSupport::M3Xd pgo::EigenSupport::toM3Xd(const double *data, size_t size)
{
  return Eigen::Map<const M3Xd>(data, 3, size / 3);
}

EIGEN_SUPPORT_INLINE pgo::EigenSupport::MX3d pgo::EigenSupport::toMX3d(const double *data, size_t size)
{
  return Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>>(data, size / 3, 3);
}

EIGEN_SUPPORT_INLINE pgo::EigenSupport::VXd pgo::EigenSupport::toVXd(const double *data, size_t size)
{
  return Eigen::Map<const VXd>(data, size);
}

EIGEN_SUPPORT_INLINE void pgo::EigenSupport::rotationOnlySVD(const M3d &A, M3d *UOut, M3d *VOut, V3d *SOut)
{
  Eigen::JacobiSVD<M3d> svd(A, Eigen::ComputeFullU | Eigen::ComputeFullV);
  M3d U = svd.matrixU();
  M3d V = svd.matrixV();
  V3d S = svd.singularValues();

  if (U.determinant() < 0.0) {
    U.col(2) *= -1.0;
    S(2) *= -1.0;
  }
  if (V.determinant() < 0.0) {
    V.col(2) *= -1.0;
    S(2) *= -1.0;
  }

  if (UOut)
    *UOut = U;

  if (VOut)
    *VOut = V;

  if (SOut)
    *SOut = S;
}

EIGEN_SUPPORT_INLINE void pgo::EigenSupport::rigidFit(const M3Xd &p0, const M3Xd &p1, M3d &R, V3d &t, const double *weights)
{
  V3d p0Center(0, 0, 0), p1Center(0, 0, 0);

  if (p0.cols() != p1.cols()) {
    std::cerr << "# vertices mismatched." << std::endl;
    return;
  }

  int numPoints = static_cast<int>(p0.cols());
  double wAll = 0.0;
  for (int i = 0; i < numPoints; ++i) {
    double w = 1;
    if (weights)
      w = weights[i];

    p0Center += p0.col(i) * w;
    p1Center += p1.col(i) * w;
    wAll += w;
  }

  double invWAll = 1.0 / wAll;

  p0Center *= invWAll;
  p1Center *= invWAll;

  t.noalias() = p1Center - p0Center;

  // compute Apq
  // Aqq don't need to be evaluate, since we only care about the rotation
  M3d Apq = M3d::Zero();

  for (int i = 0; i < numPoints; i++) {
    double w = 1;
    if (weights)
      w = weights[i];

    V3d pi(p1.col(i)), qi(p0.col(i));

    pi -= p1Center;
    qi -= p0Center;

    M3d Ai;
    Ai(0, 0) = pi[0] * qi[0];
    Ai(1, 0) = pi[1] * qi[0];
    Ai(2, 0) = pi[2] * qi[0];

    Ai(0, 1) = pi[0] * qi[1];
    Ai(1, 1) = pi[1] * qi[1];
    Ai(2, 1) = pi[2] * qi[1];

    Ai(0, 2) = pi[0] * qi[2];
    Ai(1, 2) = pi[1] * qi[2];
    Ai(2, 2) = pi[2] * qi[2];

    Apq += Ai * w;
  }

  Apq *= invWAll;

  Eigen::JacobiSVD<M3d> svd(Apq, Eigen::ComputeFullU | Eigen::ComputeFullV);
  M3d U = svd.matrixU();
  M3d V = svd.matrixV();
  V3d S = svd.singularValues();

  R = U * V.transpose();

  if (R.determinant() < 0.0) {
    S(0) = 1.0;
    S(1) = 1.0;
    S(2) = -1.0;

    R = U * S.asDiagonal() * V.transpose();
  }

  t += p0Center - R * p0Center;
}

EIGEN_SUPPORT_INLINE int pgo::EigenSupport::saveSparseMatrix(const char *filename, const SpMatD &A)
{
  FILE *file = fopen(filename, "wb");
  if (!file) {
    std::cerr << "Cannot save to file " << filename << std::endl;
    return 1;
  }

  IDX nrow = A.rows(), ncol = A.cols();
  size_t ret = fwrite(&nrow, sizeof(IDX), 1, file);
  if (ret != 1)
    return 1;

  ret = fwrite(&ncol, sizeof(IDX), 1, file);
  if (ret != 1)
    return 1;

  IDX nnz = A.nonZeros();
  ret = fwrite(&nnz, sizeof(IDX), 1, file);
  if (ret != 1)
    return 1;

  for (IDX i = 0; i < A.outerSize(); i++) {
    for (SpMatD::InnerIterator it(A, i); it; ++it) {
      IDX r = it.row(), c = it.col();
      double v = it.value();

      ret = fwrite(&r, sizeof(IDX), 1, file);
      if (ret != 1)
        return 1;

      ret = fwrite(&c, sizeof(IDX), 1, file);
      if (ret != 1)
        return 1;

      ret = fwrite(&v, sizeof(double), 1, file);
      if (ret != 1)
        return 1;
    }
  }

  fclose(file);
  return 0;
}

EIGEN_SUPPORT_INLINE int pgo::EigenSupport::loadSparseMatrix(const char *filename, SpMatD &A)
{
  FILE *file = fopen(filename, "rb");
  if (!file) {
    std::cerr << "Cannot load from file " << filename << std::endl;
    return 1;
  }

  IDX nrow, ncol;
  size_t ret = fread(&nrow, sizeof(IDX), 1, file);
  if (ret != 1)
    return 1;

  ret = fread(&ncol, sizeof(IDX), 1, file);
  if (ret != 1)
    return 1;

  IDX nnz;
  ret = fread(&nnz, sizeof(IDX), 1, file);
  if (ret != 1)
    return 1;

  std::cout << "# rows:" << nrow << "; # cols:" << ncol << "; # nnz:" << nnz << std::endl;

  std::vector<TripletD> entries;
  entries.reserve(nnz);

  for (IDX i = 0; i < nnz; i++) {
    IDX r, c;
    double v;

    ret = fread(&r, sizeof(IDX), 1, file);
    if (ret != 1)
      return 1;

    ret = fread(&c, sizeof(IDX), 1, file);
    if (ret != 1)
      return 1;

    ret = fread(&v, sizeof(double), 1, file);
    if (ret != 1)
      return 1;

    entries.emplace_back((SpMatD::StorageIndex)r, (SpMatD::StorageIndex)c, v);
  }
  fclose(file);

  A.resize(nrow, ncol);
  A.setFromTriplets(entries.begin(), entries.end());

  return 0;
}

namespace pgo
{
namespace EigenSupport
{
struct SymbolicMmData
{
  std::vector<std::vector<std::pair<std::ptrdiff_t, std::ptrdiff_t>>> mulPairs;
  int transpose;
};
}  // namespace EigenSupport
}  // namespace pgo

EIGEN_SUPPORT_INLINE void pgo::EigenSupport::symbolicMm(const SpMatD &A, const SpMatD &B, SpMatD &C, SymbolicMmData **dat, int transpose)
{
#if defined(MKL_INSPECTOR_EXECUTOR)
  sparse_matrix_t AM, BM, CM;
  AM = fromEigenMatrix(A);
  BM = fromEigenMatrix(B);

  sparse_operation_t opt = SPARSE_OPERATION_NON_TRANSPOSE;
  if (transpose)
    opt = SPARSE_OPERATION_TRANSPOSE;

  matrix_descr desc;
  desc.diag = SPARSE_DIAG_NON_UNIT;
  desc.type = SPARSE_MATRIX_TYPE_GENERAL;

  sparse_status_t ret = mkl_sparse_sp2m(opt, desc, AM,
    SPARSE_OPERATION_NON_TRANSPOSE, desc, BM,
    SPARSE_STAGE_FULL_MULT, &CM);

  if (ret != SPARSE_STATUS_SUCCESS)
    abort();

  toEigenMatrix(CM, C, 0);
#else
  if (transpose) {
    C = A.transpose() * B;
  }
  else {
    C = A * B;
  }

#endif

  *dat = new SymbolicMmData;
  (*dat)->mulPairs.resize(C.nonZeros());

  SpMatD AT;
  if (transpose) {
    AT = A.transpose();
    AT.makeCompressed();
  }

  tbb::parallel_for(0, (int)C.nonZeros(), [&](int entryi) {
    auto iter = std::upper_bound(C.outerIndexPtr(), C.outerIndexPtr() + C.outerSize() + 1, entryi);
    if (*iter <= entryi || iter == C.outerIndexPtr())
      abort();

    int rowID = (int)(iter - C.outerIndexPtr() - 1);
    int colID = C.innerIndexPtr()[entryi];

    const SpMatD *AMat = &A;
    if (transpose) {
      AMat = &AT;
    }
    for (int j = AMat->outerIndexPtr()[rowID]; j < AMat->outerIndexPtr()[rowID + 1]; j++) {
      int rowA = rowID;
      int colA = AMat->innerIndexPtr()[j];

      int rowB = colA;
      int colB = colID;

      std::ptrdiff_t offsetA = j;
      if (transpose) {
        offsetA = findEntryOffset(A, colA, rowA);
        if (offsetA < 0)
          abort();
      }

      std::ptrdiff_t offsetB = findEntryOffset(B, rowB, colB);
      if (offsetB >= 0) {
        (*dat)->mulPairs[entryi].emplace_back(offsetA, offsetB);
      }
    } });

  (*dat)->transpose = transpose;
}

EIGEN_SUPPORT_INLINE void pgo::EigenSupport::destroySymbolicMmData(SymbolicMmData *dat)
{
  delete dat;
}

EIGEN_SUPPORT_INLINE void pgo::EigenSupport::mm(const SpMatD &A, const SpMatD &B, const SymbolicMmData *dat, SpMatD &C, int transpose)
{
  if (transpose != dat->transpose) {
    std::cerr << "Opt mismatched." << std::endl;
    return;
  }

  if (transpose) {
    if (A.cols() != C.rows() || B.cols() != C.cols()) {
      std::cerr << "Size mismatched." << std::endl;
      return;
    }
  }
  else {
    if (A.rows() != C.rows() || B.cols() != C.cols()) {
      std::cerr << "Size mismatched." << std::endl;
      return;
    }
  }

  tbb::parallel_for(0, (int)C.nonZeros(), [&](int entryi) {
    C.valuePtr()[entryi] = 0;

    for (size_t pri = 0; pri < dat->mulPairs[entryi].size(); pri++) {
      const auto &pr = dat->mulPairs[entryi][pri];

      C.valuePtr()[entryi] += A.valuePtr()[pr.first] * B.valuePtr()[pr.second];
    } });
}

#ifdef __GNUC__
#  pragma GCC diagnostic pop
#endif
