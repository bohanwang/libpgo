/*
author: Bohan Wang
copyright to USC, MIT
*/

#include "smoothRSEnergy.h"

#include "polarDecompositionDerivatives.h"
#include "pgoLogging.h"

#include <tbb/parallel_for.h>
#include <tbb/partitioner.h>
#include <tbb/enumerable_thread_specific.h>

#include <mkl.h>

#include <numeric>

using namespace pgo;
using namespace pgo::PredefinedPotentialEnergies;

namespace pgo::PredefinedPotentialEnergies
{
struct SmoothRSEnergyBufferBlock
{
  SmoothRSEnergyBufferBlock(int n3_, int nele_):
    n3(n3_), nele(nele_), u(n3), F(nele * 9), S(nele * 9), R(nele * 9), temp(nele * 9),
    LTLS(nele * 9), LTLR(nele * 9), dEdF(nele * 9)
  {
  }

  int n3, nele;
  ES::VXd u, F, S, R, temp, LTLS, LTLR, dEdF;
};

struct SparseCSR
{
  void init(const sparse_matrix_t &mat, int buildCopy)
  {
    sparse_status_t ret = mkl_sparse_d_export_csr(mat, &baseIndex, &numRows, &numCols, &rowStarts, &rowEnds, &colIndices, &values);
    PGO_ALOG(ret == SPARSE_STATUS_SUCCESS);

    if (buildCopy) {
      ret = mkl_sparse_d_create_csr(&m, baseIndex, numRows, numCols, rowStarts, rowEnds, colIndices, values);
      PGO_ALOG(ret == SPARSE_STATUS_SUCCESS);
    }
  }

  ~SparseCSR()
  {
    if (m) {
      mkl_sparse_destroy(m);
    }
  }

  sparse_index_base_t baseIndex;
  int numRows = 0, numCols = 0;
  int *rowStarts = nullptr, *rowEnds = nullptr, *colIndices = nullptr;
  double *values = nullptr;

  sparse_matrix_t m = nullptr;
};

struct SmoothRSEnergyBuffer
{
  SmoothRSEnergyBuffer(int n3, int nele):
    blockTLS(n3, nele) {}

  tbb::enumerable_thread_specific<SmoothRSEnergyBufferBlock> blockTLS;
  ES::SpMatD dSdFSp, dRdFSp, d2EdF2_2nd, dFG, LTLdFG, h1, h2, h3;

  sparse_matrix_t LTL9, G, dF;
  matrix_descr desc;
  // sparse_matrix_t dFG, LTLdFG, h1, h2, dSdF, dRdF;
  // SparseCSR dFGCopy, LTLdFGCopy;
};
}  // namespace pgo::PredefinedPotentialEnergies

template<int blockDim, typename MatrixType>
void SmoothRSEnergy::computeBlockSparsePattern(const ES::SpMatD &LTL, std::unordered_map<std::pair<int, int>, MatrixType, ES::IntPairHash, ES::IntPairEqual> &blocks) const
{
  blocks.reserve(LTL.outerSize());
  for (ES::IDX rowi = 0; rowi < LTL.outerSize(); rowi++) {
    for (ES::SpMatD::InnerIterator it(LTL, rowi); it; ++it) {
      int blockRowID = (int)it.row() / blockDim;
      int blockColID = (int)it.col() / blockDim;

      int localRowID = (int)it.row() % blockDim;
      int localColID = (int)it.col() % blockDim;

      auto iter = blocks.find(std::pair<int, int>(blockRowID, blockColID));
      if (iter != blocks.end()) {
        iter->second(localRowID, localColID) = it.value();
      }
      else {
        MatrixType mat;
        mat.setZero();
        mat(localRowID, localColID) = it.value();

        blocks.emplace(std::pair<int, int>(blockRowID, blockColID), mat);
      }
    }
  }
}

inline sparse_matrix_t fromEigenMatrix(const ES::SpMatD &A)
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

inline void toEigenMatrix(sparse_matrix_t CM, ES::SpMatD &C, int onlyUpper = 0)
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

  std::vector<ES::TripletD> entries;
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

SmoothRSEnergy::SmoothRSEnergy(const Mesh::TetMeshGeo &tetMesh_, const ES::SpMatD &G_, const ES::SpMatD &L, const double coeffs[2], int isPosition, int doff):
  tetMesh(tetMesh_), G(G_), coeffR(coeffs[0]), coeffS(coeffs[1]), inputIsPosition(isPosition), dofOffsets(doff)
{
  SPDLOG_LOGGER_INFO(Logging::lgr(), "Creating smooth RS energy...");

  ES::SpMatD L9;
  ES::expandN(L, L9, 9);
  ES::mm(L9, L9, LTL9, 1);

  restPositions.resize(tetMesh.numVertices() * 3);
  for (int i = 0; i < tetMesh.numVertices(); i++) {
    restPositions.segment<3>(i * 3) = tetMesh.pos(i);
  }

  evaluationBuf = std::make_shared<SmoothRSEnergyBuffer>(tetMesh.numVertices() * 3, tetMesh.numTets());

  std::vector<ES::TripletD> entries;
  for (int ei = 0; ei < tetMesh.numTets(); ei++) {
    for (int i = 0; i < 9; i++) {
      for (int j = 0; j < 9; j++) {
        int row = ei * 9 + i;
        int col = ei * 9 + j;

        entries.emplace_back(row, col, rand() / 1.0 / RAND_MAX);
      }
    }
  }
  dF.resize(LTL9.rows(), LTL9.cols());
  dF.setFromTriplets(entries.begin(), entries.end());

  dFblocks.resize(tetMesh.numTets());
  for (int ei = 0; ei < tetMesh.numTets(); ei++) {
    ES::M9i idx;
    for (int i = 0; i < 9; i++) {
      for (int j = 0; j < 9; j++) {
        int row = ei * 9 + i;
        int col = ei * 9 + j;

        idx(i, j) = (int)ES::findEntryOffset(dF, row, col);
        PGO_ALOG(idx(i, j) >= 0);
      }
    }

    dFblocks[ei] = idx;
  }

  SPDLOG_LOGGER_INFO(Logging::lgr(), "Computing evaluation buffer for smooth RS energy...");

  evaluationBuf->LTL9 = fromEigenMatrix(LTL9);
  evaluationBuf->G = fromEigenMatrix(G);
  evaluationBuf->dF = fromEigenMatrix(dF);

  evaluationBuf->desc.type = SPARSE_MATRIX_TYPE_GENERAL;
  evaluationBuf->desc.diag = SPARSE_DIAG_NON_UNIT;

  evaluationBuf->dSdFSp = dF;
  evaluationBuf->dRdFSp = dF;
  evaluationBuf->d2EdF2_2nd = dF;

  // dS/dF * G
  ES::symbolicMm(dF, G, evaluationBuf->dFG, &dFGMapping, 0);

  // LTL (dS/dF * G)
  ES::symbolicMm(LTL9, evaluationBuf->dFG, evaluationBuf->LTLdFG, &LTLdFGMapping, 0);

  // (dS/dF * G)^T * LTL (dS/dF * G)
  ES::symbolicMm(evaluationBuf->dFG, evaluationBuf->LTLdFG, evaluationBuf->h1, &h1Mapping, 1);

  // G^T d2F G
  ES::symbolicMm(G, evaluationBuf->dFG, evaluationBuf->h2, &h2Mapping, 1);

  hessTemplate = evaluationBuf->h1 + evaluationBuf->h2;

  evaluationBuf->h3 = evaluationBuf->h1;

  ES::small2Big(evaluationBuf->h1, hessTemplate, 0, 0, h1ToHess);
  ES::small2Big(evaluationBuf->h2, hessTemplate, 0, 0, h2ToHess);

  SPDLOG_LOGGER_INFO(Logging::lgr(), "Done.");

#if 0
  // symbolic multiplication
  matrix_descr spDesc;
  spDesc.type = SPARSE_MATRIX_TYPE_GENERAL;
  spDesc.diag = SPARSE_DIAG_NON_UNIT;

  mkl_sparse_copy(evaluationBuf->dF, spDesc, &evaluationBuf->dSdF);
  mkl_sparse_copy(evaluationBuf->dF, spDesc, &evaluationBuf->dRdF);

  // dS/dF * G
  sparse_status_t ret;
  ret = mkl_sparse_sp2m(SPARSE_OPERATION_NON_TRANSPOSE, spDesc, evaluationBuf->dF,
    SPARSE_OPERATION_NON_TRANSPOSE, spDesc, evaluationBuf->G,
    SPARSE_STAGE_FULL_MULT, &evaluationBuf->dFG);
  PGO_ALOG(ret == SPARSE_STATUS_SUCCESS);

  evaluationBuf->dFGCopy.init(evaluationBuf->dFG, 1);

  // LTL (dS/dF * G)
  ret = mkl_sparse_sp2m(SPARSE_OPERATION_NON_TRANSPOSE, spDesc, evaluationBuf->LTL9,
    SPARSE_OPERATION_NON_TRANSPOSE, spDesc, evaluationBuf->dFG,
    SPARSE_STAGE_FULL_MULT, &evaluationBuf->LTLdFG);
  PGO_ALOG(ret == SPARSE_STATUS_SUCCESS);

  evaluationBuf->LTLdFGCopy.init(evaluationBuf->LTLdFG, 1);

  // (dS/dF * G)^T * LTL (dS/dF * G)
  ret = mkl_sparse_sp2m(SPARSE_OPERATION_TRANSPOSE, spDesc, evaluationBuf->dFGCopy.m,
    SPARSE_OPERATION_NON_TRANSPOSE, spDesc, evaluationBuf->LTLdFGCopy.m,
    SPARSE_STAGE_FULL_MULT, &evaluationBuf->h1);
  PGO_ALOG(ret == SPARSE_STATUS_SUCCESS);

  ret = mkl_sparse_sp2m(SPARSE_OPERATION_TRANSPOSE, spDesc, evaluationBuf->G,
    SPARSE_OPERATION_NON_TRANSPOSE, spDesc, evaluationBuf->dFG,
    SPARSE_STAGE_FULL_MULT, &evaluationBuf->h2);
  PGO_ALOG(ret == SPARSE_STATUS_SUCCESS);

  sparse_index_base_t baseIndex;
  int numRows, numCols;
  int *rowStarts, *rowEnds, *colIndices;
  double *values;
  mkl_sparse_d_export_csr(evaluationBuf->h1, &baseIndex, &numRows, &numCols, &rowStarts, &rowEnds, &colIndices, &values);

  entries.clear();
  for (int row = 0; row < numRows; row++) {
    for (int j = rowStarts[row]; j < rowStarts[row + 1]; j++) {
      int col = colIndices[j];
      entries.emplace_back(row, col, 1.0);
    }
  }

  mkl_sparse_d_export_csr(evaluationBuf->h2, &baseIndex, &numRows, &numCols, &rowStarts, &rowEnds, &colIndices, &values);

  for (int row = 0; row < numRows; row++) {
    for (int j = rowStarts[row]; j < rowStarts[row + 1]; j++) {
      int col = colIndices[j];
      entries.emplace_back(row, col, 1.0);
    }
  }

  hess.resize(numRows, numCols);
  hess.setFromTriplets(entries.begin(), entries.end());
  ES::BuildEntryMap(hess, hessMap);

  evaluationBuf->desc = spDesc;

  evaluationBuf->dSdFSp = dF;
  evaluationBuf->dRdFSp = dF;
  evaluationBuf->d2EdF2_2nd = dF;
#endif

#if 0
  // dS/dF * G
  ES::SpMatD dFG;
  ES::mm(dF, G, dFG);

  // (dS/dF * G)^T * LTL (dS/dF * G)
  ES::SpMatD LTLdFG;
  ES::mm(LTL9, dFG, LTLdFG);
  ES::mm(dFG, LTLdFG, h1, 1);

  ES::SpMatD d2FG;
  ES::mm(dF, G, d2FG);
  ES::mm(G, d2FG, h2, 1);

  hess = h1 + h2;
#endif

  allDOFs.resize(tetMesh.numVertices() * 3);
  std::iota(allDOFs.begin(), allDOFs.end(), dofOffsets);
}

ES::M3d SmoothRSEnergy::getMat(const ES::VXd &matVec, int ei) const
{
  ES::M3d F;
  const double *FVecLocal = matVec.data() + ei * 9;

  F << FVecLocal[0], FVecLocal[3], FVecLocal[6],
    FVecLocal[1], FVecLocal[4], FVecLocal[7],
    FVecLocal[2], FVecLocal[5], FVecLocal[8];

  return F;
}

ES::V9d SmoothRSEnergy::packMat(const ES::M3d &m) const
{
  return ES::Mp<const ES::V9d>(m.data());
}

void SmoothRSEnergy::getu(const ES::ConstRefVecXd x, ES::VXd &u) const
{
  int n3 = tetMesh.numVertices() * 3;

  u = x.segment(dofOffsets, n3);
  if (inputIsPosition) {
    u -= restPositions;
  }
}

// 1/2 (L S(G u + I))^2 + 1/2 (L R(G u + I))^2
double SmoothRSEnergy::func(EigenSupport::ConstRefVecXd x) const
{
  int nele = tetMesh.numTets();
  int n3 = tetMesh.numVertices() * 3;

  auto &buf = evaluationBuf->blockTLS.local();

  ES::VXd &u = buf.u;
  getu(x, u);

  ES::VXd &F = buf.F, &S = buf.S, &R = buf.R;
  computeFSR(u, F, S, R);

  double energy = 0;
  ES::VXd &temp = buf.temp;

  sparse_operation_t opt = SPARSE_OPERATION_NON_TRANSPOSE;
  double e1, e2;
  int ret = mkl_sparse_d_dotmv(opt, 1.0, evaluationBuf->LTL9, evaluationBuf->desc, S.data(), 0.0, temp.data(), &e1);
  PGO_ALOG(ret == SPARSE_STATUS_SUCCESS);

  ret = mkl_sparse_d_dotmv(opt, 1.0, evaluationBuf->LTL9, evaluationBuf->desc, R.data(), 0.0, temp.data(), &e2);
  PGO_ALOG(ret == SPARSE_STATUS_SUCCESS);

  energy = e1 * coeffS + e2 * coeffR;
  return energy * 0.5;
}

// dE = d( 1/2 (L S(G u + I))^2 + 1/2 (L R(G u + I))^2 )
// = LTL S * dS/dF * G + LTL R * dR/dF * G
// = (LTL S)^T * dS/dF * G + (LTL R)^T * dR/dF * G
// transpose = GT (dS/dF)^T LTL S + GT (dR/dF)^T LTL R
//          = GT ((dS/dF)^T LTL S + (dR/dF)^T LTL R)

void SmoothRSEnergy::gradient(EigenSupport::ConstRefVecXd x, EigenSupport::RefVecXd grad) const
{
  // std::cout << "grad start." << std::endl;
  int nele = tetMesh.numTets();
  int n3 = tetMesh.numVertices() * 3;

  auto &buf = evaluationBuf->blockTLS.local();

  ES::VXd &u = buf.u;
  getu(x, u);

  ES::VXd &F = buf.F, &S = buf.S, &R = buf.R;
  computeFSR(u, F, S, R);

  // LTL S and LTL R
  ES::VXd &LTLS = buf.LTLS, &LTLR = buf.LTLR;

  sparse_operation_t opt = SPARSE_OPERATION_NON_TRANSPOSE;
  int ret = mkl_sparse_d_mv(opt, 1.0, evaluationBuf->LTL9, evaluationBuf->desc, S.data(), 0.0, LTLS.data());
  PGO_ALOG(ret == SPARSE_STATUS_SUCCESS);

  ret = mkl_sparse_d_mv(opt, 1.0, evaluationBuf->LTL9, evaluationBuf->desc, R.data(), 0.0, LTLR.data());
  PGO_ALOG(ret == SPARSE_STATUS_SUCCESS);

  ES::VXd &dEdF = buf.dEdF;

  // for (int ele = 0; ele < tetMesh.numTets(); ele++) {
  tbb::parallel_for(
    0, nele, [&](int ele) {
      ES::M3d Flocal = getMat(F, ele);
      ES::M3d Slocal = getMat(S, ele);
      ES::M3d Rlocal = getMat(R, ele);

      ES::M3d SInv = Slocal.fullPivHouseholderQr().inverse();

      ES::M3d dRdFi[9], dSdFi[9];
      NonlinearOptimization::PolarDecompositionDerivatives::computeFirstOrderDerivatives(Flocal, Rlocal, Slocal, SInv, dRdFi, dSdFi);

      ES::M9d dSdF, dRdF;
      for (int dof = 0; dof < 9; dof++) {
        dSdF.col(dof) = ES::Mp<const ES::V9d>(dSdFi[dof].data());
        dRdF.col(dof) = ES::Mp<const ES::V9d>(dRdFi[dof].data());
      }

      // d S6 / dF = d W S9 / dF = W * d S9 / dF
      // (d S6 / dF)^T = (d S9 / d F)^T * W^T

      ES::M9d dSdFT = dSdF.transpose();
      ES::M9d dRdFT = dRdF.transpose();

      // (ds/df)^T * LTL S + (dr/df)^T * LTL R
      dEdF.segment<9>(ele * 9) = dSdFT * LTLS.segment<9>(ele * 9) * coeffS + dRdFT * LTLR.segment<9>(ele * 9) * coeffR;
    },
    tbb::static_partitioner());

  opt = SPARSE_OPERATION_TRANSPOSE;
  ret = mkl_sparse_d_mv(opt, 1.0, evaluationBuf->G, evaluationBuf->desc, dEdF.data(), 0.0, grad.data());
  PGO_ALOG(ret == SPARSE_STATUS_SUCCESS);

  // std::cout << "grad done." << std::endl;
}

// grad^T = GT ((dS/dF)^T LTL S + (dR/dF)^T LTL R)
// ( S^T LTL dS/dF + R^T LTL dRdF ) * G
// first term: S
// h = (dS/dF)^T LTL dS/dF + S^T LTL d2S/dF2

void SmoothRSEnergy::hessian(EigenSupport::ConstRefVecXd x, EigenSupport::SpMatD &hess) const
{
  int nele = tetMesh.numTets();
  int n3 = tetMesh.numVertices() * 3;

  auto &buf = evaluationBuf->blockTLS.local();

  ES::VXd &u = buf.u;
  getu(x, u);

  ES::VXd &F = buf.F, &S = buf.S, &R = buf.R;
  computeFSR(u, F, S, R);

  // LTL S and LTL R
  ES::VXd &LTLS = buf.LTLS, &LTLR = buf.LTLR;
  sparse_operation_t opt = SPARSE_OPERATION_NON_TRANSPOSE;
  int ret = mkl_sparse_d_mv(opt, 1.0, evaluationBuf->LTL9, evaluationBuf->desc, S.data(), 0.0, LTLS.data());
  PGO_ALOG(ret == SPARSE_STATUS_SUCCESS);

  ret = mkl_sparse_d_mv(opt, 1.0, evaluationBuf->LTL9, evaluationBuf->desc, R.data(), 0.0, LTLR.data());
  PGO_ALOG(ret == SPARSE_STATUS_SUCCESS);

  ES::SpMatD &dSdFSp = evaluationBuf->dSdFSp;
  ES::SpMatD &dRdFSp = evaluationBuf->dRdFSp;
  ES::SpMatD &d2EdF2_2nd = evaluationBuf->d2EdF2_2nd;

  // compute d2EdF2
  // ds/df dr/df and second term
  // for (int ele = 0; ele < tetMesh.numTets(); ele++) {
  tbb::parallel_for(
    0, nele, [&](int ele) {
      ES::M3d Flocal = getMat(F, ele);
      ES::M3d Slocal = getMat(S, ele);
      ES::M3d Rlocal = getMat(R, ele);

      ES::M3d SInv = Slocal.fullPivHouseholderQr().inverse();

      ES::M3d dRdFi[9], dSdFi[9];
      ES::M3d dFdFi[9], dAdFi[9];
      ES::M3d dsqrtA_dA[9];
      ES::M9d MInv;
      NonlinearOptimization::PolarDecompositionDerivatives::computeFirstOrderDerivatives(Flocal, Rlocal, Slocal, SInv, dRdFi, dSdFi,
        dFdFi, dAdFi, dsqrtA_dA, &MInv);

      ES::M3d d2RdF2i[9][9];
      ES::M3d d2SdF2i[9][9];
      NonlinearOptimization::PolarDecompositionDerivatives::computeSecondOrderDerivatives(Flocal, Rlocal, Slocal, SInv, MInv, dRdFi, dSdFi, dFdFi, dAdFi, dsqrtA_dA,
        d2SdF2i, d2RdF2i);

      ES::M9d dSdF, dRdF;
      for (int dof = 0; dof < 9; dof++) {
        dSdF.col(dof) = ES::Mp<const ES::V9d>(dSdFi[dof].data());
        dRdF.col(dof) = ES::Mp<const ES::V9d>(dRdFi[dof].data());
      }

      for (int i = 0; i < 9; i++) {
        for (int j = 0; j < 9; j++) {
          int offset = dFblocks[ele](i, j);
          dSdFSp.valuePtr()[offset] = dSdF(i, j);
          dRdFSp.valuePtr()[offset] = dRdF(i, j);
        }
      }

      ES::M9d d2EdF2;
      for (int i = 0; i < 9; i++) {
        for (int j = 0; j < 9; j++) {
          ES::V9d dSdFidFj = ES::Mp<ES::V9d>(d2SdF2i[i][j].data());
          ES::V9d dRdFidFj = ES::Mp<ES::V9d>(d2RdF2i[i][j].data());

          d2EdF2(i, j) = dSdFidFj.dot(LTLS.segment<9>(ele * 9)) * coeffS + dRdFidFj.dot(LTLR.segment<9>(ele * 9)) * coeffR;
        }
      }

      for (int i = 0; i < 9; i++) {
        for (int j = 0; j < 9; j++) {
          int offset = dFblocks[ele](i, j);
          d2EdF2_2nd.valuePtr()[offset] = d2EdF2(i, j);
        }
      }
    },
    tbb::static_partitioner());

  memset(hess.valuePtr(), 0, sizeof(double) * hess.nonZeros());

  // dS/dF G
  ES::mm(dSdFSp, G, dFGMapping, evaluationBuf->dFG, 0);
  // LTL(dS/dF G)
  ES::mm(LTL9, evaluationBuf->dFG, LTLdFGMapping, evaluationBuf->LTLdFG, 0);
  // (dS/dF G)^T LTL (dS/dF G)
  ES::mm(evaluationBuf->dFG, evaluationBuf->LTLdFG, h1Mapping, evaluationBuf->h1, 1);

  // dR/dF G
  ES::mm(dRdFSp, G, dFGMapping, evaluationBuf->dFG, 0);
  // LTL(dR/dF G)
  ES::mm(LTL9, evaluationBuf->dFG, LTLdFGMapping, evaluationBuf->LTLdFG, 0);
  // (dR/dF G)^T LTL (dR/dF G)
  ES::mm(evaluationBuf->dFG, evaluationBuf->LTLdFG, h1Mapping, evaluationBuf->h3, 1);

  // (dS/dF G)^T LTL (dS/dF G) * coeff0 + (dR/dF G)^T LTL (dR/dF G) * coeff1
  cblas_daxpy(evaluationBuf->h1.nonZeros(), coeffR, evaluationBuf->h3.valuePtr(), coeffS, evaluationBuf->h1.valuePtr(), 1);

  // d2E/dF2 G
  ES::mm(d2EdF2_2nd, G, dFGMapping, evaluationBuf->dFG, 0);
  ES::mm(G, evaluationBuf->dFG, h2Mapping, evaluationBuf->h2, 1);

  ES::addSmallToBig(1.0, evaluationBuf->h1, hess, 0.0, h1ToHess, 1);
  ES::addSmallToBig(1.0, evaluationBuf->h2, hess, 1.0, h2ToHess, 1);

#if 0

  sparse_index_base_t baseIndex;
  int numRows, numCols;
  int *rowStarts, *rowEnds, *colIndices;
  double *values;

  sparse_index_base_t baseIndex_dF;
  int numRows_dF, numCols_dF;
  int *rowStarts_dF, *rowEnds_dF, *colIndices_dF;
  double *values_dF;
  ret = mkl_sparse_d_export_csr(evaluationBuf->dF, &baseIndex_dF, &numRows_dF, &numCols_dF, &rowStarts_dF, &rowEnds_dF, &colIndices_dF, &values_dF);
  PGO_ALOG(ret == SPARSE_STATUS_SUCCESS);

  //ES::SpMatD dSdFG, dRdFG;
  //ES::mm(dSdFSp, G, dSdFG);
  //ES::mm(dRdFSp, G, dRdFG);

  //// (dS/dF G)^T LTL (dS/dF G)
  //ES::SpMatD h1, h2;
  //ES::mm(LTL9, dSdFG, h1);
  //ES::mm(dSdFG, h1, h2, 1);

  //// (dR/dF G)^T LTL (dR/dF G)
  //ES::SpMatD h3, h4;
  //ES::mm(LTL9, dRdFG, h3);
  //ES::mm(dRdFG, h3, h4, 1);

  //h4.makeCompressed();
  //h2.makeCompressed();

  //ES::SpMatD h22 = h2;
  //h22.makeCompressed();

  //// (dS/dF G)^T LTL (dS/dF G) * coeff0 + (dR/dF G)^T LTL (dR/dF G) * coeff1
  //cblas_dscal(h2.nonZeros(), coeffS, h2.valuePtr(), 1);
  //cblas_daxpy(h2.nonZeros(), coeffR, h4.valuePtr(), 1, h2.valuePtr(), 1);

  // dS
  memcpy(values_dF, dSdFSp.valuePtr(), dSdFSp.nonZeros() * sizeof(double));

  // dS/dF * G
  ret = mkl_sparse_sp2m(SPARSE_OPERATION_NON_TRANSPOSE, evaluationBuf->desc, evaluationBuf->dF,
    SPARSE_OPERATION_NON_TRANSPOSE, evaluationBuf->desc, evaluationBuf->G,
    SPARSE_STAGE_FINALIZE_MULT, &evaluationBuf->dFG);
  PGO_ALOG(ret == SPARSE_STATUS_SUCCESS);

  // LTL (dS/dF * G)
  ret = mkl_sparse_sp2m(SPARSE_OPERATION_NON_TRANSPOSE, evaluationBuf->desc, evaluationBuf->LTL9,
    SPARSE_OPERATION_NON_TRANSPOSE, evaluationBuf->desc, evaluationBuf->dFG,
    SPARSE_STAGE_FINALIZE_MULT, &evaluationBuf->LTLdFG);
  PGO_ALOG(ret == SPARSE_STATUS_SUCCESS);

  // (dS/dF * G)^T * LTL (dS/dF * G)
  SparseCSR m1, m2, m3;
  m1.init(evaluationBuf->dFG, 1);
  m2.init(evaluationBuf->LTLdFG, 1);

  ret = mkl_sparse_sp2m(SPARSE_OPERATION_TRANSPOSE, evaluationBuf->desc, m1.m,
    SPARSE_OPERATION_NON_TRANSPOSE, evaluationBuf->desc, m2.m,
    SPARSE_STAGE_FULL_MULT, &m3.m);
  PGO_ALOG(ret == SPARSE_STATUS_SUCCESS);

  //ret = mkl_sparse_sp2m(SPARSE_OPERATION_TRANSPOSE, evaluationBuf->desc, evaluationBuf->dFG,
  //  SPARSE_OPERATION_NON_TRANSPOSE, evaluationBuf->desc, evaluationBuf->LTLdFG,
  //  SPARSE_STAGE_FINALIZE_MULT, &evaluationBuf->h1);
  //PGO_ALOG(ret == SPARSE_STATUS_SUCCESS);

  //ES::SpMatD hh2;
  //toEigenMatrix(m3.m, hh2);
  //LGI << (hh2 - h22).norm();
  //LGI << h22.norm();

  //ES::SpMatD dSdFG1;
  //toEigenMatrix(evaluationBuf->dFG, dSdFG1);
  //LGI << (dSdFG1 - dSdFG).norm();

  //ES::SpMatD dLTLdFG1;
  //toEigenMatrix(evaluationBuf->LTLdFG, dLTLdFG1);
  //LGI << (h1 - dLTLdFG1).norm();
  //LGI << h1.norm();

  // ES::SpMatD hh2, hh3;
  // toEigenMatrix(m3.m, hh2);
  // LGI << (h2 - hh2).norm();
  // ES::mm(dSdFG1, dLTLdFG1, hh3, 1);
  // LGI << (h2 - hh3).norm();
  // LGI << h2.norm();

  ret = mkl_sparse_d_export_csr(m3.m, &baseIndex, &numRows, &numCols, &rowStarts, &rowEnds, &colIndices, &values);
  PGO_ALOG(ret == SPARSE_STATUS_SUCCESS);

  for (int row = 0; row < numRows; row++) {
    for (int j = rowStarts[row]; j < rowStarts[row + 1]; j++) {
      int col = colIndices[j];

      auto iter = hessMap.find(std::make_pair(row, col));
      PGO_ALOG(iter != hessMap.end());
      std::ptrdiff_t offset = iter->second;
      hess.valuePtr()[offset] += values[j] * coeffS;
    }
  }

  // dR
  memcpy(values_dF, dRdFSp.valuePtr(), dRdFSp.nonZeros() * sizeof(double));

  // dR/dF * G
  ret = mkl_sparse_sp2m(SPARSE_OPERATION_NON_TRANSPOSE, evaluationBuf->desc, evaluationBuf->dF,
    SPARSE_OPERATION_NON_TRANSPOSE, evaluationBuf->desc, evaluationBuf->G,
    SPARSE_STAGE_FINALIZE_MULT, &evaluationBuf->dFG);
  PGO_ALOG(ret == SPARSE_STATUS_SUCCESS);

  // LTL (dR/dF * G)
  ret = mkl_sparse_sp2m(SPARSE_OPERATION_NON_TRANSPOSE, evaluationBuf->desc, evaluationBuf->LTL9,
    SPARSE_OPERATION_NON_TRANSPOSE, evaluationBuf->desc, evaluationBuf->dFG,
    SPARSE_STAGE_FINALIZE_MULT, &evaluationBuf->LTLdFG);
  PGO_ALOG(ret == SPARSE_STATUS_SUCCESS);

  mkl_sparse_destroy(m1.m);
  mkl_sparse_destroy(m2.m);
  mkl_sparse_destroy(m3.m);

  m1.init(evaluationBuf->dFG, 1);
  m2.init(evaluationBuf->LTLdFG, 1);

  // (dR/dF * G)^T * LTL (dR/dF * G)
  ret = mkl_sparse_sp2m(SPARSE_OPERATION_TRANSPOSE, evaluationBuf->desc, m1.m,
    SPARSE_OPERATION_NON_TRANSPOSE, evaluationBuf->desc, m2.m,
    SPARSE_STAGE_FULL_MULT, &m3.m);
  PGO_ALOG(ret == SPARSE_STATUS_SUCCESS);

  //ret = mkl_sparse_sp2m(SPARSE_OPERATION_TRANSPOSE, evaluationBuf->desc, evaluationBuf->dFG,
  //  SPARSE_OPERATION_NON_TRANSPOSE, evaluationBuf->desc, evaluationBuf->LTLdFG,
  //  SPARSE_STAGE_FINALIZE_MULT, &evaluationBuf->h1);
  //PGO_ALOG(ret == SPARSE_STATUS_SUCCESS);

  //ES::SpMatD hh4;
  //toEigenMatrix(m3.m, hh4);
  //LGI << (hh4 - h4).norm();
  //LGI << h4.norm();

  ret = mkl_sparse_d_export_csr(m3.m, &baseIndex, &numRows, &numCols, &rowStarts, &rowEnds, &colIndices, &values);
  PGO_ALOG(ret == SPARSE_STATUS_SUCCESS);

  for (int row = 0; row < numRows; row++) {
    for (int j = rowStarts[row]; j < rowStarts[row + 1]; j++) {
      int col = colIndices[j];

      auto iter = hessMap.find(std::make_pair(row, col));
      PGO_ALOG(iter != hessMap.end());
      std::ptrdiff_t offset = iter->second;
      hess.valuePtr()[offset] += values[j] * coeffR;
    }
  }

  //ES::SpMatD diff = h2 - hess;
  //double err = diff.norm();
  //LGI << err;

  // d2F
  memcpy(values_dF, d2EdF2_2nd.valuePtr(), d2EdF2_2nd.nonZeros() * sizeof(double));

  ret = mkl_sparse_sp2m(SPARSE_OPERATION_NON_TRANSPOSE, evaluationBuf->desc, evaluationBuf->dF,
    SPARSE_OPERATION_NON_TRANSPOSE, evaluationBuf->desc, evaluationBuf->G,
    SPARSE_STAGE_FINALIZE_MULT, &evaluationBuf->dFG);
  PGO_ALOG(ret == SPARSE_STATUS_SUCCESS);

  ret = mkl_sparse_sp2m(SPARSE_OPERATION_TRANSPOSE, evaluationBuf->desc, evaluationBuf->G,
    SPARSE_OPERATION_NON_TRANSPOSE, evaluationBuf->desc, evaluationBuf->dFG,
    SPARSE_STAGE_FINALIZE_MULT, &evaluationBuf->h2);
  PGO_ALOG(ret == SPARSE_STATUS_SUCCESS);

  ret = mkl_sparse_d_export_csr(evaluationBuf->h2, &baseIndex, &numRows, &numCols, &rowStarts, &rowEnds, &colIndices, &values);
  PGO_ALOG(ret == SPARSE_STATUS_SUCCESS);

  for (int row = 0; row < numRows; row++) {
    for (int j = rowStarts[row]; j < rowStarts[row + 1]; j++) {
      int col = colIndices[j];

      auto iter = hessMap.find(std::make_pair(row, col));
      PGO_ALOG(iter != hessMap.end());
      std::ptrdiff_t offset = iter->second;
      hess.valuePtr()[offset] += values[j] * coeffR;
    }
  }

  //ES::SpMatD h5, h6;
  //ES::mm(d2EdF2_2nd, G, h5);
  //ES::mm(G, h5, h6, 1);

  //ES::SpMatD hess1 = h2 + h6;
  //hess1.makeCompressed();

  //LGI << (hess - hess1).norm();
  //LGI << hess1.norm();
#endif

#if 0
  ES::SpMatD dSdFG, dRdFG;
  ES::mm(dSdFSp, G, dSdFG);
  ES::mm(dRdFSp, G, dRdFG);

  // (dS/dF G)^T LTL (dS/dF G)
  ES::SpMatD h1, h2;
  ES::mm(LTL9, dSdFG, h1);
  ES::mm(dSdFG, h1, h2, 1);

  // (dR/dF G)^T LTL (dR/dF G)
  ES::SpMatD h3, h4;
  ES::mm(LTL9, dRdFG, h3);
  ES::mm(dRdFG, h3, h4, 1);

  h4.makeCompressed();
  h2.makeCompressed();

  // (dS/dF G)^T LTL (dS/dF G) * coeff0 + (dR/dF G)^T LTL (dR/dF G) * coeff1
  cblas_dscal(h2.nonZeros(), coeffS, h2.valuePtr(), 1);
  cblas_daxpy(h2.nonZeros(), coeffR, h4.valuePtr(), 1, h2.valuePtr(), 1);

  ES::SpMatD h5, h6;
  ES::mm(d2EdF2_2nd, G, h5);
  ES::mm(G, h5, h6, 1);

  ES::SpMatD hess1 = h2 + h6;
  hess1.makeCompressed();

  // sanity check
  for (int outeri = 0; outeri < hess1.outerSize(); outeri++) {
    PGO_ALOG(hess1.outerIndexPtr()[outeri] == hess.outerIndexPtr()[outeri]);
  }

  PGO_ALOG(hess.nonZeros() == hess1.nonZeros());
  for (int ni = 0; ni < hess1.nonZeros(); ni++) {
    PGO_ALOG(hess1.innerIndexPtr()[ni] == hess.innerIndexPtr()[ni]);
  }
  memcpy(hess.valuePtr(), hess1.valuePtr(), sizeof(double) * hess1.nonZeros());
#endif
}

void SmoothRSEnergy::computeFSR(const ES::VXd &u, ES::VXd &F, ES::VXd &S, ES::VXd &R) const
{
  int nele = tetMesh.numTets();
  ES::mv(G, u, F);

  // for (int ele = 0; ele < tetMesh.numTets(); ele++) {
  tbb::parallel_for(
    0, nele, [&](int ele) {
      Eigen::Map<ES::M3d>(F.data() + ele * 9) += ES::M3d::Identity();

      ES::M3d Flocal = getMat(F, ele);
      ES::M3d Slocal, Rlocal;

      Slocal = NonlinearOptimization::PolarDecompositionDerivatives::computeS(Flocal);
      Rlocal = NonlinearOptimization::PolarDecompositionDerivatives::computeR(Flocal, Slocal);

      S.segment<9>(ele * 9) = packMat(Slocal);
      R.segment<9>(ele * 9) = packMat(Rlocal);
    },
    tbb::static_partitioner());
}

void SmoothRSEnergy::convertRowMajorFG(const ES::SpMatD &GRowMajor, ES::SpMatD &GColMajor)
{
  int nele = (int)GRowMajor.rows() / 9;
  int indexMapping[9] = { 0, 3, 6, 1, 4, 7, 2, 5, 8 };

  std::vector<ES::TripletD> entries;
  for (int ei = 0; ei < nele; ei++) {
    for (int rowi = 0; rowi < 9; rowi++) {
      int srcRowID = ei * 9 + rowi;
      int tgtRowID = ei * 9 + indexMapping[rowi];

      for (ES::SpMatD::InnerIterator it(GRowMajor, srcRowID); it; ++it) {
        entries.emplace_back(tgtRowID, (int)it.col(), it.value());
      }
    }
  }

  GColMajor.resize(GRowMajor.rows(), GRowMajor.cols());
  GColMajor.setFromTriplets(entries.begin(), entries.end());
}