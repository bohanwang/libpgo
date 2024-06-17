/*
author: Bohan Wang
copyright to USC
*/

#include "svdDerivatives.h"

using namespace pgo;

namespace pgo::NonlinearOptimization::SVDDerivatives
{
template<typename RealScalar>
void makeJacobi(const RealScalar &A11, const RealScalar &A12, const RealScalar &A22, RealScalar &c, RealScalar &s);
template<typename MatrixType, typename RealScalar, typename Index>
void real2x2JacobiSVD(const MatrixType &matrix, Index p, Index q,
  Eigen::JacobiRotation<RealScalar> &j_left, Eigen::JacobiRotation<RealScalar> &j_right);
template<typename DerivedA, typename DerivedS, typename DerivedU, typename DerivedV>
void unorderedSquareMatrixSVD(const DerivedA &A, DerivedS &S, DerivedU *U = nullptr, DerivedV *V = nullptr);
template<typename DerivedA, typename DerivedS, typename DerivedU, typename DerivedV>
void unorderedSquareMatrixSVDDerivatices(const DerivedA &F, const DerivedS &S, const DerivedU &U, const DerivedV &V,
  const DerivedA &dF, DerivedS &dS, DerivedU *dU = nullptr, DerivedV *dV = nullptr,
  const DerivedA *d2F = nullptr, DerivedS *d2S = nullptr, DerivedU *d2U = nullptr, DerivedV *d2V = nullptr);
}  // namespace pgo::NonlinearOptimization::SVDDerivatives

template<typename RealScalar>
void NonlinearOptimization::SVDDerivatives::makeJacobi(const RealScalar &A11, const RealScalar &A12, const RealScalar &A22, RealScalar &c, RealScalar &s)
{
  if (std::abs(A12) < std::numeric_limits<RealScalar>::min()) {
    c = RealScalar(1);
    s = RealScalar(0);
  }
  else {
    RealScalar rho = (A11 - A22) / A12;
    RealScalar delta = std::sqrt(rho * rho + RealScalar(4));
    RealScalar t;
    if (rho > RealScalar(0)) {
      t = (rho - delta) * RealScalar(0.5);
    }
    else {
      t = (rho + delta) * RealScalar(0.5);
    }

    c = 1 / sqrt(1 + t * t);
    s = c * t;
  }
}

template<typename MatrixType, typename RealScalar, typename Index>
void NonlinearOptimization::SVDDerivatives::real2x2JacobiSVD(const MatrixType &matrix, Index p, Index q,
  Eigen::JacobiRotation<RealScalar> &j_left, Eigen::JacobiRotation<RealScalar> &j_right)
{
  Eigen::Matrix<RealScalar, 2, 2> m;
  m << RealScalar(matrix.coeff(p, p)), RealScalar(matrix.coeff(p, q)),
    RealScalar(matrix.coeff(q, p)), RealScalar(matrix.coeff(q, q));

  Eigen::JacobiRotation<RealScalar> rot1;
  RealScalar t = m.coeff(0, 0) + m.coeff(1, 1);
  RealScalar d = m.coeff(1, 0) - m.coeff(0, 1);

  RealScalar z = std::sqrt(t * t + d * d);
  rot1.s() = d / z;
  rot1.c() = t / z;

#if 1
  if (std::abs(d) < (std::numeric_limits<RealScalar>::min)()) {
    rot1.s() = RealScalar(0);
    rot1.c() = RealScalar(1);
  }
  else {
    // If d!=0, then t/d cannot overflow because the magnitude of the
    // entries forming d are not too small compared to the ones forming t.
    RealScalar u = t / d;
    RealScalar tmp = std::sqrt(RealScalar(1) + u * u);
    rot1.s() = RealScalar(1) / tmp;
    rot1.c() = u / tmp;
  }
#endif

  m.applyOnTheLeft(0, 1, rot1);
  // makeJacobi(m(0, 0), m(0, 1), m(1, 1), j_right->c(), j_right->s());
  j_right.makeJacobi(m, 0, 1);
  j_left = rot1 * j_right.transpose();
}

template<typename DerivedA, typename DerivedS, typename DerivedU, typename DerivedV>
void NonlinearOptimization::SVDDerivatives::unorderedSquareMatrixSVD(const DerivedA &A, DerivedS &S, DerivedU *U, DerivedV *V)
{
  static_assert(DerivedA::ColsAtCompileTime == DerivedA::RowsAtCompileTime &&
    DerivedA::RowsAtCompileTime > 0 && DerivedA::ColsAtCompileTime > 0);

  // It is copied from Eigen
  using RealScalar = typename DerivedA::Scalar;

  // currently we stop when we reach precision 2*epsilon as the last bit of precision can require an unreasonable number of iterations,
  // only worsening the precision of U and V as we accumulate more rotations
  const RealScalar precision = RealScalar(2.0) * Eigen::NumTraits<RealScalar>::epsilon();

  // limit for denormal numbers to be considered zero in order to avoid infinite loops (see bug 286)
  constexpr RealScalar considerAsZero = std::numeric_limits<RealScalar>::min();

  constexpr ES::IDX dim = DerivedA::RowsAtCompileTime;

  if (U && (U->rows() == 0 || U->cols() == 0)) {
    U->resize(dim, dim);
  }

  if (V && (V->rows() == 0 || V->cols() == 0)) {
    V->resize(dim, dim);
  }

  if (S.size() == 0)
    S.resize(dim);

  // Scaling factor to reduce over/under-flows
  RealScalar scale = A.cwiseAbs().maxCoeff();
  if (scale == RealScalar(0))
    scale = RealScalar(1);

  /*** step 1. The R-SVD step: we use a QR decomposition to reduce to the case of a square matrix */
  // we omit this step as we assume the matrix to be a square matrix
  Eigen::Matrix<typename DerivedA::Scalar, DerivedA::ColsAtCompileTime, DerivedA::ColsAtCompileTime, DerivedA::Options> workMatrix;

  workMatrix.noalias() = A / scale;

  if (U)
    U->setIdentity();

  if (V)
    V->setIdentity();

  /*** step 2. The main Jacobi SVD iteration. ***/
  RealScalar maxDiagEntry = workMatrix.cwiseAbs().diagonal().maxCoeff();

  bool finished = false;
  while (!finished) {
    finished = true;

    // do a sweep: for all index pairs (p,q), perform SVD of the corresponding 2x2 sub-matrix
    for (ES::IDX p = 1; p < dim; ++p) {
      for (ES::IDX q = 0; q < p; ++q) {
        // if this 2x2 sub-matrix is not diagonal already...
        // notice that this comparison will evaluate to false if any NaN is involved, ensuring that NaN's don't
        // keep us iterating forever. Similarly, small denormal numbers are considered zero.
        RealScalar threshold = std::max(considerAsZero, precision * maxDiagEntry);

        if (std::abs(workMatrix.coeff(p, q)) > threshold || std::abs(workMatrix.coeff(q, p)) > threshold) {
          finished = false;

          // perform SVD decomposition of 2x2 sub-matrix corresponding to indices p,q to make it diagonal
          // the complex to real operation returns true if the updated 2x2 block is not already diagonal
          Eigen::JacobiRotation<RealScalar> j_left, j_right;
          real2x2JacobiSVD(workMatrix, p, q, j_left, j_right);

          // accumulate resulting Jacobi rotations
          workMatrix.applyOnTheLeft(p, q, j_left);
          if (U)
            U->applyOnTheRight(p, q, j_left.transpose());

          workMatrix.applyOnTheRight(p, q, j_right);
          if (V)
            V->applyOnTheRight(p, q, j_right);

          // keep track of the largest diagonal coefficient
          maxDiagEntry = std::max(maxDiagEntry, std::max(std::abs(workMatrix.coeff(p, p)), abs(workMatrix.coeff(q, q))));
        }
      }
    }
  }

  /*** step 3. The work matrix is now diagonal ***/
  for (ES::IDX i = 0; i < dim; ++i) {
    RealScalar a = (RealScalar)workMatrix.coeff(i, i);
    S.coeffRef(i) = a;
  }

  S *= scale;
}

template<typename DerivedA, typename DerivedS, typename DerivedU, typename DerivedV>
void NonlinearOptimization::SVDDerivatives::unorderedSquareMatrixSVDDerivatices(const DerivedA &, const DerivedS &S, const DerivedU &U, const DerivedV &V,
  const DerivedA &dF, DerivedS &dS, DerivedU *dU, DerivedV *dV,
  const DerivedA *d2F, DerivedS *d2S, DerivedU *d2U, DerivedV *d2V)
{
  static_assert(DerivedA::ColsAtCompileTime == DerivedA::RowsAtCompileTime &&
    DerivedA::RowsAtCompileTime > 0 && DerivedA::ColsAtCompileTime > 0);

  static_assert(DerivedS::RowsAtCompileTime > 0 && DerivedS::ColsAtCompileTime == 1);

  using RealScalar = typename DerivedA::Scalar;
  constexpr ES::IDX dim = DerivedA::RowsAtCompileTime;

  // F = U S VT
  // dF = dU S VT + U dS VT + U S dVT
  // UT dF V = UTdU S + dS + S (VTdV)T

  // both are skew symmetric <-- UT U = I
  // dY = U^T dU
  // dZ = V^T dV

  // skew symmetric mulitplying with diagonal matrix == > matrix with zero diag
  // so dS is the diag of dA
  if (dS.size() == 0)
    dS.resize(dim);

  Eigen::Matrix<RealScalar, dim, dim> dA = U.transpose() * dF * V;
  for (ES::IDX i = 0; i < dim; i++)
    dS(i) = dA(i, i);

  DerivedS S2 = S.cwiseProduct(S);
  Eigen::Matrix<RealScalar, dim, dim> dUlocal, dVlocal;

  if (dU || (d2F && (d2S || d2U || d2V))) {
    // UT dF V = UTdU S + dS + S (VTdV)T
    // K = Ibar x dA = dY S - S dZ
    // KT = Ibar x dAT = -S dY + dZ S

    // K S = dY S^2 - S dZ S      ----- 1
    // S KT = -S^2 dY + S dZ S    ----- 2
    // 1 + 2 =>
    // K S + S KT = dY S^2 - S^2 dY
    Eigen::Matrix<RealScalar, dim, dim> dY;
    Eigen::Matrix<RealScalar, dim, dim> W = dA * S.asDiagonal() + S.asDiagonal() * dA.transpose();

    dY.setZero();
    for (ES::IDX row = 0; row < dim; row++) {
      for (ES::IDX col = 0; col < dim; col++) {
        if (row == col)
          continue;

        double coeff = S2(col) - S2(row);
        if (std::abs(coeff) < 1e-10)
          dY(row, col) = 0;
        else {
          double rhs = W(row, col);
          dY(row, col) = rhs / coeff;
        }
      }
    }

    dUlocal = U * dY;
  }

  if (dV || (d2F && (d2S || d2U || d2V))) {
    // UT dF V = UTdU S + dS + S (VTdV)T
    // K = Ibar x dA = dY S - S dZ
    // KT = Ibar x dAT = -S dY + dZ S

    // S K  = S dY S - S^2 dZ     ----- 1
    // KT S = -S dY S + dZ S^2    ----- 2
    // 1 + 2 =>
    // S K + KT S= dZ S^2 - S^2 dZ
    Eigen::Matrix<RealScalar, dim, dim> dZ;
    Eigen::Matrix<RealScalar, dim, dim> W = S.asDiagonal() * dA + dA.transpose() * S.asDiagonal();

    dZ.setZero();
    for (ES::IDX row = 0; row < dim; row++) {
      for (ES::IDX col = 0; col < dim; col++) {
        if (row == col)
          continue;

        double coeff = S2(col) - S2(row);        
        if (std::abs(coeff) < 1e-10)
          dZ(row, col) = 0;
        else {
          double rhs = W(row, col);
          dZ(row, col) = rhs / coeff;
        }
      }
    }

    dVlocal = V * dZ;
  }

  if (dU)
    *dU = dUlocal;

  if (dV)
    *dV = dVlocal;

  Eigen::Matrix<RealScalar, dim, 1> d2Slocal;
  if (d2F && (d2S || d2U || d2V)) {
    // UT dF V = UTdU S + dS + S (VTdV)T
    // I x (UT dF V) = dS
    // I x (dUT dF V + UT d2F V + UT dF dV) = d2S
    Eigen::Matrix<RealScalar, dim, dim> dUTdFV = dUlocal.transpose() * dF * V;
    Eigen::Matrix<RealScalar, dim, dim> UTd2FV = U.transpose() * *d2F * V;
    Eigen::Matrix<RealScalar, dim, dim> UTdFdV = U.transpose() * dF * dVlocal;

    for (ES::IDX i = 0; i < dim; i++)
      d2Slocal(i) = dUTdFV(i, i) + UTd2FV(i, i) + UTdFdV(i, i);
  }

  if (d2S)
    *d2S = d2Slocal;

  if (d2F && d2U) {
    throw std::runtime_error("Not supported yet");
  }

  if (d2F && d2V) {
    throw std::runtime_error("Not supported yet");
  }
}

void NonlinearOptimization::SVDDerivatives::unorderedSquareMatrixSVD3(const ES::M3d &F, ES::V3d &S, ES::M3d *U, ES::M3d *V)
{
  unorderedSquareMatrixSVD(F, S, U, V);
}

void NonlinearOptimization::SVDDerivatives::unorderedSquareMatrixSVD3Derivatices(const ES::M3d &F, const ES::V3d &S, const ES::M3d &U, const ES::M3d &V,
  const ES::M3d &dF, ES::V3d &dS, ES::M3d *dU, ES::M3d *dV,
  const ES::M3d *d2F, ES::V3d *d2S, ES::M3d *d2U, ES::M3d *d2V)
{
  unorderedSquareMatrixSVDDerivatices(F, S, U, V, dF, dS, dU, dV, d2F, d2S, d2U, d2V);
}

void NonlinearOptimization::SVDDerivatives::unorderedSquareMatrixSVD2(const ES::M2d &F, ES::V2d &S, ES::M2d *U, ES::M2d *V)
{
  unorderedSquareMatrixSVD(F, S, U, V);
}

void NonlinearOptimization::SVDDerivatives::unorderedSquareMatrixSVD2Derivatices(const ES::M2d &F, const ES::V2d &S, const ES::M2d &U, const ES::M2d &V,
  const ES::M2d &dF, ES::V2d &dS, ES::M2d *dU, ES::M2d *dV,
  const ES::M2d *d2F, ES::V2d *d2S, ES::M2d *d2U, ES::M2d *d2V)
{
  unorderedSquareMatrixSVDDerivatices(F, S, U, V, dF, dS, dU, dV, d2F, d2S, d2U, d2V);
}

#if 0

template<typename MatType1, typename MatType2, typename MatType3, typename MatType4>
void firstOrder(const MatType1 &dA, const MatType2 &U, const MatType3 &S, const MatType4 &V, MatType3 &dS, MatType2 *dU, MatType4 *dV)
{
  typedef MatType1::Scalar RealScalar;
  typedef MatType1::Index IDX;
  typedef Eigen::Matrix<RealScalar, MatType1::ColsAtCompileTime, MatType1::ColsAtCompileTime, MatType1::Options> SmallSquareMat;

  RealScalar eps = RealScalar(1e-5);
  // both are skew symmetric
  // dY = U^T dU
  // dZ = V^T dV

  // dP = U^T dA V = dY S + dS + S dZ^T \in R^{k x k}
  SmallSquareMat dP;
  dP = U.transpose() * dA * V;

  // dS = I * dP, because dY S and S dZ^T has zero diagnoal
  // * is element-wise product
  IDX minDOF = std::min(dA.rows(), dA.cols());

  if (dS.size() < minDOF)
    dS.resize(minDOF);

  for (IDX i = 0; i < minDOF; i++)
    dS[i] = dP(i, i);

  // dY S + S dZ^T = Ibar * dP

  // dY S - S dZ = Ibar * dP = W (1)
  // (dY S)_ij - (S dZ)_ij = W_ij
  // dY_ij * S_j - S_i * dZ_ij = W_ij

  // transpose above equation
  // -S dY + dZ S = Ibar * dP^T (2)
  // -S_i dY_ij + dZ_ij * S_j = W^T_ij

  // i.e. we solve
  // dY_ij * S_j - S_i * dZ_ij = W_ij
  // -S_i dY_ij + dZ_ij * S_j = W^T_ij

  // for any i,j where i < j.

  if (dU || dV) {
    SmallSquareMat dY(dP.rows(), dP.cols());
    SmallSquareMat dZ(dP.rows(), dP.cols());

    dY.setZero(dP.rows(), dP.cols());
    dZ.setZero(dP.rows(), dP.cols());
    // std::vector<int> equalSingularValues(dP.rows() * dP.rows(), 0);

    for (ES::IDX i = 0; i < dP.rows(); i++) {
      for (ES::IDX j = i + 1; j < dP.rows(); j++) {
        Eigen::Matrix<RealScalar, 2, 2> sys;
        Eigen::Matrix<RealScalar, 2, 1> rhs, x;

        if (fabs(S[j] - S[i]) < eps) {
          // equalSingularValues[i * dP.rows() + j] = 1;

          sys << S[j] + eps, -S[i],
            -S[i], S[j] + eps;
        }
        else {
          sys << S[j], -S[i], -S[i], S[j];
        }

        rhs[0] = dP(i, j);
        rhs[1] = dP(j, i);

        x = sys.fullPivHouseholderQr().solve(rhs);

        dY(i, j) = x[0];
        dY(j, i) = -x[0];

        dZ(i, j) = x[1];
        dZ(j, i) = -x[1];
      }
    }
#  if 0
    bool hasEqualSingularValues = false;
    for (ES::IDX i = 0; i < dP.rows(); i++) {
      for (ES::IDX j = i + 1; j < dP.rows(); j++) {
        if (equalSingularValues[i * dP.rows() + j]) {
          hasEqualSingularValues = true;
          break;
        }
      }
      if (hasEqualSingularValues)
        break;
    }

    if (hasEqualSingularValues) {
      ES::MXd sys(dP.rows() * dP.rows() * 2, dP.rows() * dP.rows() * 2);
      ES::VXd rhs(dP.rows() * dP.rows() * 2);
      ES::VXd x(dP.rows() * dP.rows() * 2);
      std::vector<int> equalSingularValuesID(dP.rows(), 0);

      for (ES::IDX i = 0; i < dP.rows(); i++) {
        // check whether there are any similar singular values
        int count = 0;
        for (ES::IDX j = i + 1; j < dP.rows(); j++) {
          if (equalSingularValues[i * dP.rows() + j]) {
            equalSingularValuesID[count + 1] = (int)j;
            count++;
          }
        }

        if (count) {
          equalSingularValuesID[0] = (int)i;
          int totalValues = count + 1;

          for (int c = 0; c < totalValues * (totalValues - 1); c++) {
            for (int r = 0; r < totalValues * (totalValues - 1); r++) {
              sys(r, c) = 0;
            }
          }

          count = 0;
          for (int i1 = 0; i1 < totalValues; i1++) {
            for (int i2 = i1 + 1; i2 < totalValues; i2++) {
              int s1 = equalSingularValuesID[i1];
              int s2 = equalSingularValuesID[i2];

              sys(count * 2, count * 2) = S[s2];
              sys(count * 2, count * 2 + 1) = -S[s1];

              sys(count * 2 + 1, count * 2) = -S[s1];
              sys(count * 2 + 1, count * 2 + 1) = S[s2];

              rhs[count * 2] = Ibar_dP(s1, s2);
              rhs[count * 2 + 1] = Ibar_dP(s2, s1);

              count++;
            }
          }

          // ATA x - AT b = 0
          ES::MXd ATA = sys.block(0, 0, count * 2, count * 2).transpose() * sys.block(0, 0, count * 2, count * 2);
          ES::VXd ATb = sys.block(0, 0, count * 2, count * 2).transpose() * rhs.head(count * 2);
          x.head(count * 2) = ATA.fullPivLu().solve(ATb);

          ES::VXd z1 = sys.block(0, 0, count * 2, count * 2) * x.head(count * 2);

          count = 0;
          for (int i1 = 0; i1 < totalValues; i1++) {
            for (int i2 = i1 + 1; i2 < totalValues; i2++) {
              int s1 = equalSingularValuesID[i1];
              int s2 = equalSingularValuesID[i2];

              dY(s1, s2) = x[count * 2];
              dY(s2, s1) = -x[count * 2];

              dZ(s1, s2) = x[count * 2 + 1];
              dZ(s2, s1) = -x[count * 2 + 1];

              count++;

              equalSingularValues[s1 * dP.rows() + s2] = 0;
            }
          }
        }
      }
    }
#  endif

    if (dU) {
      *dU = U * dY;
    }

    if (dV) {
      *dV = V * dZ;
    }
  }
}


template<typename MatType1, typename MatType2, typename MatType3, typename MatType4>
void secondOrder(const MatType1 &dA, const MatType1 &d2A,
  const MatType2 &U, const MatType3 &S, const MatType4 &V,
  const MatType2 &dU, const MatType4 &dV, MatType3 &d2S)
{
  typedef MatType1::Scalar RealScalar;
  typedef MatType1::Index IDX;
  typedef Eigen::Matrix<RealScalar, MatType1::ColsAtCompileTime, MatType1::ColsAtCompileTime, MatType1::Options> SmallSquareMat;

  IDX minDOF = std::min(dA.rows(), dA.cols());

  if (d2S.size() < minDOF) {
    d2S.resize(minDOF);
  }

  // both algorithm work fine
  // the first one requires more computation
  // the second one does not
#  if 0
  // dU^T U + U^T dU = 0
  // d2U^T U + dU^T dU + dU^T dU + U^T d2U = 0
  // d2U^T U + U^T d2U = - 2 dU^T dU
  // (U^T d2U)^T + U^T d2U = -2 dU^T dU
  SmallSquareMat dUTdU = dU.transpose() * dU;
  SmallSquareMat dVTdV = dV.transpose() * dV;

  ES::VXd diagUTd2U = dUTdU.diagonal() * -1.0;
  ES::VXd diagVTd2V = dVTdV.diagonal() * -1.0;

  // dA = dU S VT + U dS VT + U S dVT
  // d2A = d2U S VT + dU dS VT + dU S dVT
  //     + dU dS VT + U d2S VT + U dS dVT
  //     + dU S dVT + U dS dVT + U S d2VT
  // UT d2A V = (UT d2U) S + d2S + S (d2VT V)
  //          + 2 UT dU dS + 2 dS dVT V + 2 UT dU S dVT V
  // I * (U^T d2A V) = I * ( (UT d2U) S + d2S + S (d2VT V) + 2 UT dU S dVT V ) + I * (  2 UT dU dS + 2 dS dVT V )
  // I * (U^T d2A V) = I * ( (UT d2U) S + d2S + S (d2VT V) + 2 UT dU S dVT V ) + 0
  // d2S = I * (U^T d2A V - 2 UT dU S dVT V - (UT d2U) S - S (d2VT V))
  SmallSquareMat dP = U.transpose() * d2A * V;
  dP -= U.transpose() * dU * S.asDiagonal() * dV.transpose() * V * 2;

  for (IDX i = 0; i < minDOF; i++)
    d2S[i] = dP(i, i) - diagUTd2U[i] * S[i] - diagVTd2V[i] * S[i];
#  else
  // * is element-wise product
  // I *(UT dA V) = dS;
  // d2S = I * (dU^T dA V + UT d2A V + UT dA dV)
  SmallSquareMat dP1 = dU.transpose() * dA * V + U.transpose() * d2A * V + U.transpose() * dA * dV;
  for (IDX i = 0; i < minDOF; i++)
    d2S[i] = dP1(i, i);
#  endif
}

void firstOrder1()
{
  ES::M3d A = ES::V3d(1, 2, 3).asDiagonal();
  ES::M3d dA = ES::M3d::Zero();
  dA(0, 1) = 1.0;

  ES::M3d d2A = ES::M3d::Zero();

  ES::M3d U, V, dU, dV;
  ES::V3d S, dS, d2S;
  unorderedSVD(A, S, &U, &V);

  firstOrder(dA, U, S, V, dS, &dU, &dV);
  secondOrder(dA, d2A, U, S, V, dU, dV, d2S);

  for (ES::IDX i = 0; i < S.size(); i++) {
    ES::IDX m = A.rows(), n = A.cols();
    ES::IDX p = m + n;

    ES::M3d B11 = U.transpose() * dA * V;
    Eigen::SelfAdjointEigenSolver<ES::M3d> eig;
    eig.compute((B11 + B11.transpose()) * 0.5, Eigen::ComputeEigenvectors);

    ES::V3d eigval = eig.eigenvalues();
    ES::M3d eigvec = eig.eigenvectors();

    double b11 = U.col(i).dot(dA * V.col(i));
    /*
    ES::MXd Q1 = ES::MXd::Zero(p * 3, p * 3);

    Q1.block(0, 0, m, m) = ES::MXd::Identity(m, m) * S[i];
    Q1.block(m, m, n, n) = ES::MXd::Identity(n, n) * S[i];

    Q1.block(p, p, m, m) = ES::MXd::Identity(m, m) * S[i];
    Q1.block(p + m, p + m, n, n) = ES::MXd::Identity(n, n) * -S[i];

    Q1.block(p * 2, p * 2, m, m) = ES::MXd::Identity(m, m) * -S[i];
    Q1.block(p * 2 + m, p * 2 + m, n, n) = ES::MXd::Identity(n, n) * -S[i];

    Q1.block(0, m, m, n) = A;
    Q1.block(m, 0, n, m) = A.transpose();

    Q1.block(p, p + m, m, n) = A;
    Q1.block(p + m, p, n, m) = A.transpose();

    Q1.block(p * 2, p * 2 + m, m, n) = A;
    Q1.block(p * 2 + m, p * 2, n, m) = A.transpose();

    Q1.block(p, m, m, n) = dA;
    Q1.block(p + m, 0, n, m) = dA.transpose();

    Q1.block(p * 2, p + m, m, n) = dA;
    Q1.block(2 * p + m, p, n, m) = dA.transpose();

    Q1.block(p * 2, m, m, n) = ES::MXd::Zero(m, n);
    Q1.block(p * 2 + m, 0, n, m) = ES::MXd::Zero(n, m);

    Q1.block(p * 2, p, m, m) = ES::MXd::Identity(m, m) * dS[i];
    Q1.block(p * 2 + m, p + m, n, n) = ES::MXd::Identity(n, n) * -dS[i];

    ES::MXd Q2 = ES::MXd::Zero(p * 3, p * 3);
    Q2.block(p, 0, m, m) = ES::MXd::Identity(m, m);
    Q2.block(p + m, m, n, n) = ES::MXd::Identity(n, n);

    Q2.block(p * 2, 0, m, m) = ES::MXd::Identity(m, m);
    Q2.block(p * 2 + m, m, n, n) = ES::MXd::Identity(n, n);

    std::cout << "Q1:\n"
              << Q1 << "\n\n";
    std::cout << "Q2:\n"
              << Q2 << "\n\n";
   */
    ES::MXd Q1 = ES::MXd::Zero(p * 2, p * 2);
    Q1.block(0, 0, m, m) = ES::MXd::Identity(m, m) * -S[i];
    Q1.block(m, m, n, n) = ES::MXd::Identity(n, n) * -S[i];
    Q1.block(p, p, m, m) = ES::MXd::Identity(m, m) * -S[i];
    Q1.block(p + m, p + m, n, n) = ES::MXd::Identity(n, n) * -S[i];

    Q1.block(0, m, m, n) = A;
    Q1.block(m, 0, n, m) = A.transpose();

    Q1.block(p, p + m, m, n) = A;
    Q1.block(p + m, p, n, m) = A.transpose();

    Q1.block(p + m, 0, n, m) = dA.transpose();
    Q1.block(p, m, m, n) = dA;

    ES::MXd Q2 = ES::MXd::Zero(p * 2, p * 2);
    Q2.block(p, 0, m, m) = ES::MXd::Identity(m, m);
    Q2.block(p + m, m, n, n) = ES::MXd::Identity(n, n);

    std::cout << "Q1:\n"
              << Q1 << "\n\n";
    std::cout << "Q2:\n"
              << Q2 << "\n\n";

    ES::VXd x(p * 2);
    x.segment(0, m) = U.col(i);
    x.segment(m, n) = V.col(i);
    x.segment(p, m) = dU.col(i);
    x.segment(p + m, n) = dV.col(i);

    x.normalize();

    std::cout << std::setprecision(15) << (Q1 * x - Q2 * x * dS[i]) << std::endl;
    std::cout << std::setprecision(15) << "x:" << x.transpose() << std::endl;
    std::cout << "det Q1=" << Q1.determinant() << std::endl;

    Eigen::JacobiSVD<ES::MXd> svd;
    svd.compute(Q1);
    std::cout << "Q1 S: " << svd.singularValues().transpose() << std::endl;

    ES::MXd Q11 = Q1;
    ES::MXd Q21 = Q2;
    int ilo, ihi;
    ES::VXd lscale(Q11.rows()), rscale(Q21.rows());

    ES::VXd a1(Q11.rows()), a2(Q11.rows()), b1(Q11.rows());
    double abnrm, bbnrm;
    ES::VXd rconde(Q11.rows()), rcondv(Q11.rows());

    ES::MXd retv = Q11;
    int rr = LAPACKE_dggevx(LAPACK_COL_MAJOR, 'B', 'N', 'N', 'E', Q11.rows(), Q11.data(), Q11.rows(), Q21.data(), Q21.rows(), a1.data(), a2.data(), b1.data(), retv.data(), retv.rows(), retv.data(), retv.rows(),
      &ilo, &ihi, lscale.data(), rscale.data(), &abnrm, &bbnrm, rconde.data(), rcondv.data());

    std::cout << rconde.transpose() << std::endl;

    int ret = LAPACKE_dggbal(LAPACK_COL_MAJOR, 'B', Q11.rows(), Q11.data(), Q11.rows(), Q21.data(), Q21.rows(), &ilo, &ihi, lscale.data(), rscale.data());

    std::cout << "Q1:\n"
              << Q11 << "\n\n";
    std::cout << "Q2:\n"
              << Q21 << "\n\n";

    ES::VXd tau(Q2.rows());
    ret = LAPACKE_dgeqrf(LAPACK_COL_MAJOR, Q21.rows(), Q21.cols(), Q21.data(), Q21.rows(), tau.data());

    LAPACKE_dormqr(LAPACK_COL_MAJOR, 'L', 'T', Q11.rows(), Q11.rows(), Q11.cols(), Q21.data(), Q21.rows(), tau.data(), Q11.data(), Q11.rows());

    std::cout << "Q1:\n"
              << Q11 << "\n\n";
    std::cout << "Q2:\n"
              << Q21 << "\n\n";

    ES::MXd O = ES::MXd::Identity(Q11.rows(), Q11.cols());
    ES::MXd Z = ES::MXd::Identity(Q11.rows(), Q11.cols());
    ret = LAPACKE_dgghrd(LAPACK_COL_MAJOR, 'I', 'I', Q11.rows(), 1, Q11.rows(), Q11.data(), Q11.rows(), Q21.data(), Q21.rows(), O.data(), O.rows(), Z.data(), Z.rows());

    std::cout << "Q1:\n"
              << Q11 << "\n\n";
    std::cout << "Q2:\n"
              << Q21 << "\n\n";

    ES::VXd v1(Q11.rows()), v2(Q11.rows()), v3(Q11.rows());
    ret = LAPACKE_dhgeqz(LAPACK_COL_MAJOR, 'E', 'N', 'N', O.rows(), 1, O.rows(), Q11.data(), Q11.rows(), Q21.data(), Q21.rows(), v1.data(), v2.data(), v3.data(), O.data(), O.size(), Z.data(), Z.size());

    Eigen::GeneralizedEigenSolver<ES::MXd> eig1;
    eig1.compute(Q1, Q2, true);

    Eigen::MatrixXcd eigvec1 = eig1.eigenvectors();
    Eigen::VectorXcd eigval1 = eig1.alphas();
    Eigen::VectorXcd eigval2 = eig1.betas();

    std::cout << "eval:\n"
              << eigval1.transpose() << "\n"
              << eigval2.transpose() << "\n";

    std::cout << "evec:\n"
              << eigvec1 << "\n\n";
  }
}
#endif