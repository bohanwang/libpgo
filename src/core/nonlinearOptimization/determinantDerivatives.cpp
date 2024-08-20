/*
author: Bohan Wang
copyright to USC
*/

#include "determinantDerivatives.h"

#include "EigenSupport.h"

namespace pgo
{
namespace ES = EigenSupport;

namespace NonlinearOptimization
{
namespace Determinant
{
template<int power>
inline double mypow(double val)
{
  double ret = 1;
  for (int i = 0; i < power; i++)
    ret *= val;

  return ret;
}

template<>
inline double mypow<0>(double)
{
  return 1;
}

template<>
inline double mypow<1>(double val)
{
  return val;
}

template<>
inline double mypow<2>(double val)
{
  return val * val;
}

template<>
inline double mypow<3>(double val)
{
  return val * val * val;
}

template<>
inline double mypow<4>(double val)
{
  return val * val * val * val;
}

#define MY_POWER(val, power) (mypow<power>(val))

namespace Dim3
{
inline const double &elt(const double A[9], int row, int col)
{
  return A[col * 3 + row];
}

inline const double &elt(const double A[81], int row0, int col0, int row1, int col1)
{
  int row = col0 * 3 + row0;
  int col = col1 * 3 + row1;

  return A[col * 9 + row];
}

inline double &elt(double A[9], int row, int col)
{
  return A[col * 3 + row];
}

inline double &elt(double A[81], int row0, int col0, int row1, int col1)
{
  int row = col0 * 3 + row0;
  int col = col1 * 3 + row1;

  return A[col * 9 + row];
}

using Map3 = Eigen::Map<ES::M3d>;
using MapC3 = Eigen::Map<const ES::M3d>;

const double dFdSData[9 * 6] = {
  1, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 1, 0, 1, 0, 0, 0, 0, 0,
  0, 0, 1, 0, 0, 0, 1, 0, 0,
  0, 0, 0, 0, 1, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 1, 0, 1, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 1
};

typedef Eigen::Matrix<double, 9, 6> M9x6d;
typedef Eigen::Matrix<double, 6, 9> M6x9d;
typedef Eigen::Matrix<double, 6, 6> M6d;
typedef Eigen::Matrix<double, 6, 1> V6d;

typedef Eigen::Map<M9x6d> Map9x6d;
typedef Eigen::Map<const M9x6d> MapC9x6d;

M9x6d dFdS = MapC9x6d(dFdSData);
M6x9d dSdF = dFdS.transpose();
}  // namespace Dim3
}  // namespace Determinant

double Determinant::Dim3::det(const double M[9])
{
  return (MapC3(M)).determinant();
}

void Determinant::Dim3::ddetA_dA(const double A[9], double ddetA_dAOut[9])
{
  elt(ddetA_dAOut, 0, 0) = -(elt(A, 1, 2) * elt(A, 2, 1)) + elt(A, 1, 1) * elt(A, 2, 2);
  elt(ddetA_dAOut, 0, 1) = elt(A, 1, 2) * elt(A, 2, 0) - elt(A, 1, 0) * elt(A, 2, 2);
  elt(ddetA_dAOut, 0, 2) = -(elt(A, 1, 1) * elt(A, 2, 0)) + elt(A, 1, 0) * elt(A, 2, 1);

  elt(ddetA_dAOut, 1, 0) = elt(A, 0, 2) * elt(A, 2, 1) - elt(A, 0, 1) * elt(A, 2, 2);
  elt(ddetA_dAOut, 1, 1) = -(elt(A, 0, 2) * elt(A, 2, 0)) + elt(A, 0, 0) * elt(A, 2, 2);
  elt(ddetA_dAOut, 1, 2) = elt(A, 0, 1) * elt(A, 2, 0) - elt(A, 0, 0) * elt(A, 2, 1);

  elt(ddetA_dAOut, 2, 0) = -(elt(A, 0, 2) * elt(A, 1, 1)) + elt(A, 0, 1) * elt(A, 1, 2);
  elt(ddetA_dAOut, 2, 1) = elt(A, 0, 2) * elt(A, 1, 0) - elt(A, 0, 0) * elt(A, 1, 2);
  elt(ddetA_dAOut, 2, 2) = -(elt(A, 0, 1) * elt(A, 1, 0)) + elt(A, 0, 0) * elt(A, 1, 1);
}

void Determinant::Dim3::d2detA_dA2(const double A[9], double ret[81])
{
  elt(ret, 0, 0, 0, 0) = 0;
  elt(ret, 0, 0, 0, 1) = 0;
  elt(ret, 0, 0, 0, 2) = 0;

  elt(ret, 0, 0, 1, 0) = 0;
  elt(ret, 0, 0, 1, 1) = elt(A, 2, 2);
  elt(ret, 0, 0, 1, 2) = -elt(A, 2, 1);

  elt(ret, 0, 0, 2, 0) = 0;
  elt(ret, 0, 0, 2, 1) = -elt(A, 1, 2);
  elt(ret, 0, 0, 2, 2) = elt(A, 1, 1);

  elt(ret, 0, 1, 0, 0) = 0;
  elt(ret, 0, 1, 0, 1) = 0;
  elt(ret, 0, 1, 0, 2) = 0;

  elt(ret, 0, 1, 1, 0) = -elt(A, 2, 2);
  elt(ret, 0, 1, 1, 1) = 0;
  elt(ret, 0, 1, 1, 2) = elt(A, 2, 0);

  elt(ret, 0, 1, 2, 0) = elt(A, 1, 2);
  elt(ret, 0, 1, 2, 1) = 0;
  elt(ret, 0, 1, 2, 2) = -elt(A, 1, 0);

  elt(ret, 0, 2, 0, 0) = 0;
  elt(ret, 0, 2, 0, 1) = 0;
  elt(ret, 0, 2, 0, 2) = 0;

  elt(ret, 0, 2, 1, 0) = elt(A, 2, 1);
  elt(ret, 0, 2, 1, 1) = -elt(A, 2, 0);
  elt(ret, 0, 2, 1, 2) = 0;

  elt(ret, 0, 2, 2, 0) = -elt(A, 1, 1);
  elt(ret, 0, 2, 2, 1) = elt(A, 1, 0);
  elt(ret, 0, 2, 2, 2) = 0;

  elt(ret, 1, 0, 0, 0) = 0;
  elt(ret, 1, 0, 0, 1) = -elt(A, 2, 2);
  elt(ret, 1, 0, 0, 2) = elt(A, 2, 1);

  elt(ret, 1, 0, 1, 0) = 0;
  elt(ret, 1, 0, 1, 1) = 0;
  elt(ret, 1, 0, 1, 2) = 0;

  elt(ret, 1, 0, 2, 0) = 0;
  elt(ret, 1, 0, 2, 1) = elt(A, 0, 2);
  elt(ret, 1, 0, 2, 2) = -elt(A, 0, 1);

  elt(ret, 1, 1, 0, 0) = elt(A, 2, 2);
  elt(ret, 1, 1, 0, 1) = 0;
  elt(ret, 1, 1, 0, 2) = -elt(A, 2, 0);

  elt(ret, 1, 1, 1, 0) = 0;
  elt(ret, 1, 1, 1, 1) = 0;
  elt(ret, 1, 1, 1, 2) = 0;

  elt(ret, 1, 1, 2, 0) = -elt(A, 0, 2);
  elt(ret, 1, 1, 2, 1) = 0;
  elt(ret, 1, 1, 2, 2) = elt(A, 0, 0);

  elt(ret, 1, 2, 0, 0) = -elt(A, 2, 1);
  elt(ret, 1, 2, 0, 1) = elt(A, 2, 0);
  elt(ret, 1, 2, 0, 2) = 0;

  elt(ret, 1, 2, 1, 0) = 0;
  elt(ret, 1, 2, 1, 1) = 0;
  elt(ret, 1, 2, 1, 2) = 0;

  elt(ret, 1, 2, 2, 0) = elt(A, 0, 1);
  elt(ret, 1, 2, 2, 1) = -elt(A, 0, 0);
  elt(ret, 1, 2, 2, 2) = 0;

  elt(ret, 2, 0, 0, 0) = 0;
  elt(ret, 2, 0, 0, 1) = elt(A, 1, 2);
  elt(ret, 2, 0, 0, 2) = -elt(A, 1, 1);

  elt(ret, 2, 0, 1, 0) = 0;
  elt(ret, 2, 0, 1, 1) = -elt(A, 0, 2);
  elt(ret, 2, 0, 1, 2) = elt(A, 0, 1);

  elt(ret, 2, 0, 2, 0) = 0;
  elt(ret, 2, 0, 2, 1) = 0;
  elt(ret, 2, 0, 2, 2) = 0;

  elt(ret, 2, 1, 0, 0) = -elt(A, 1, 2);
  elt(ret, 2, 1, 0, 1) = 0;
  elt(ret, 2, 1, 0, 2) = elt(A, 1, 0);

  elt(ret, 2, 1, 1, 0) = elt(A, 0, 2);
  elt(ret, 2, 1, 1, 1) = 0;
  elt(ret, 2, 1, 1, 2) = -elt(A, 0, 0);

  elt(ret, 2, 1, 2, 0) = 0;
  elt(ret, 2, 1, 2, 1) = 0;
  elt(ret, 2, 1, 2, 2) = 0;

  elt(ret, 2, 2, 0, 0) = elt(A, 1, 1);
  elt(ret, 2, 2, 0, 1) = -elt(A, 1, 0);
  elt(ret, 2, 2, 0, 2) = 0;

  elt(ret, 2, 2, 1, 0) = -elt(A, 0, 1);
  elt(ret, 2, 2, 1, 1) = elt(A, 0, 0);
  elt(ret, 2, 2, 1, 2) = 0;

  elt(ret, 2, 2, 2, 0) = 0;
  elt(ret, 2, 2, 2, 1) = 0;
  elt(ret, 2, 2, 2, 2) = 0;
}

void Determinant::Dim3::ddetA_dA_sym(const double param[6], double ret[6])
{
  //ddetA_dA[0] = -MY_POWER(param[4], 2) + param[3] * param[5];
  //ddetA_dA[1] = 2 * param[2] * param[4] - 2 * param[1] * param[5];
  //ddetA_dA[2] = -2 * param[2] * param[3] + 2 * param[1] * param[4];
  //ddetA_dA[3] = -MY_POWER(param[2], 2) + param[0] * param[5];
  //ddetA_dA[4] = 2 * param[1] * param[2] - 2 * param[0] * param[4];
  //ddetA_dA[5] = -MY_POWER(param[1], 2) + param[0] * param[3];

  double F[9] = {
    param[0], param[1], param[2],
    param[1], param[3], param[4],
    param[2], param[4], param[5]
  };

  ES::V9d ddetF_dF;
  ddetA_dA(F, ddetF_dF.data());

  (Eigen::Map<V6d>(ret)) = dSdF * ddetF_dF;
}

void Determinant::Dim3::ddetA_dA_diag(const double param[3], double ret[3])
{
  ret[0] = param[1] * param[2];
  ret[1] = param[0] * param[2];
  ret[2] = param[0] * param[1];
}

void Determinant::Dim3::d2detA_dA2_sym(const double param[6], double ret[36])
{
  double F[9] = {
    param[0], param[1], param[2],
    param[1], param[3], param[4],
    param[2], param[4], param[5]
  };

  ES::M9d d2detF_dF2;
  d2detA_dA2(F, d2detF_dF2.data());

  (Eigen::Map<M6d>(ret)) = dSdF * d2detF_dF2 * dFdS;
}

void Determinant::Dim3::d2detA_dA2_diag(const double param[3], double ret[9])
{
  elt(ret, 0, 0) = 0;
  elt(ret, 0, 1) = param[2];
  elt(ret, 0, 2) = param[1];

  elt(ret, 1, 0) = param[2];
  elt(ret, 1, 1) = 0;
  elt(ret, 1, 2) = param[0];

  elt(ret, 2, 0) = param[1];
  elt(ret, 2, 1) = param[0];
  elt(ret, 2, 2) = 0;
}

namespace Determinant
{
namespace Dim2
{
using M2d = Eigen::Matrix2d;
using Map2 = Eigen::Map<M2d>;
using MapC2 = Eigen::Map<const M2d>;

inline const double &elt(const double A[4], int row, int col)
{
  return A[col * 2 + row];
}

inline const double &elt(const double A[16], int row0, int col0, int row1, int col1)
{
  int row = col0 * 2 + row0;
  int col = col1 * 2 + row1;

  return A[col * 4 + row];
}

inline double &elt(double A[4], int row, int col)
{
  return A[col * 2 + row];
}

inline double &elt(double A[16], int row0, int col0, int row1, int col1)
{
  int row = col0 * 2 + row0;
  int col = col1 * 2 + row1;

  return A[col * 4 + row];
}

const double dFdSData[4 * 3] = {
  1, 0, 0, 0,
  0, 1, 1, 0,
  0, 0, 0, 1
};

typedef Eigen::Matrix<double, 4, 3> M4x3d;
typedef Eigen::Matrix<double, 3, 4> M3x4d;

typedef Eigen::Map<M4x3d> Map4x3d;
typedef Eigen::Map<const M4x3d> MapC4x3d;

M4x3d dFdS = MapC4x3d(dFdSData);
M3x4d dSdF = dFdS.transpose();
}  // namespace Dim2
}  // namespace Determinant

double Determinant::Dim2::det(const double M[4])
{
  return (MapC2(M)).determinant();
}

void Determinant::Dim2::ddetA_dA(const double A[4], double ret[4])
{
  elt(ret, 0, 0) = elt(A, 1, 1);
  elt(ret, 0, 1) = -elt(A, 1, 0);
  elt(ret, 1, 0) = -elt(A, 0, 1);
  elt(ret, 1, 1) = elt(A, 0, 0);
}

void Determinant::Dim2::d2detA_dA2(const double[4], double ret[16])
{
  elt(ret, 0, 0, 0, 0) = 0;
  elt(ret, 0, 0, 0, 1) = 0;
  elt(ret, 0, 0, 1, 0) = 0;
  elt(ret, 0, 0, 1, 1) = 1;

  elt(ret, 0, 1, 0, 0) = 0;
  elt(ret, 0, 1, 0, 1) = 0;
  elt(ret, 0, 1, 1, 0) = -1;
  elt(ret, 0, 1, 1, 1) = 0;

  elt(ret, 1, 0, 0, 0) = 0;
  elt(ret, 1, 0, 0, 1) = -1;
  elt(ret, 1, 0, 1, 0) = 0;
  elt(ret, 1, 0, 1, 1) = 0;

  elt(ret, 1, 1, 0, 0) = 1;
  elt(ret, 1, 1, 0, 1) = 0;
  elt(ret, 1, 1, 1, 0) = 0;
  elt(ret, 1, 1, 1, 1) = 0;
}

void Determinant::Dim2::ddetA_dA_sym(const double param[3], double ret[3])
{
  double F[4] = {
    param[0], param[1],
    param[1], param[2]
  };

  ES::V4d ddetF_dF;
  ddetA_dA(F, ddetF_dF.data());

  (Eigen::Map<ES::V3d>(ret)) = dSdF * ddetF_dF;
}

void Determinant::Dim2::ddetA_dA_diag(const double param[2], double ret[2])
{
  ret[0] = param[1];
  ret[1] = param[0];
}

void Determinant::Dim2::d2detA_dA2_sym(const double param[3], double ret[9])
{
  double F[4] = {
    param[0], param[1],
    param[1], param[2]
  };

  ES::M4d d2detF_dF2;
  d2detA_dA2(F, d2detF_dF2.data());

  (Eigen::Map<ES::M3d>(ret)) = dSdF * d2detF_dF2 * dFdS;
}

void Determinant::Dim2::d2detA_dA2_diag(const double[2], double ret[4])
{
  elt(ret, 0, 0) = 0;
  elt(ret, 0, 1) = 1;

  elt(ret, 1, 0) = 1;
  elt(ret, 1, 1) = 0;
}

}  // namespace NonlinearOptimization
}  // namespace pgo
