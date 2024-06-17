/*
author: Bohan Wang
copyright to USC
*/

#pragma once

namespace pgo
{
namespace NonlinearOptimization
{
namespace Determinant
{
namespace Dim3
{
// here if it is a matrix, then it will be column major
double det(const double M[9]);
// A = { A00, A10, A20, A01, A11, A21, A02, A12, A22 }
void ddetA_dA(const double A[9], double ddetA_dAOut[9]);
// A = { A00, A01, A02, A11, A12, A22 }
void ddetA_dA_sym(const double A[6], double ddetA_dAOut[6]);
// A = { A00, A11, A22 }
void ddetA_dA_diag(const double A[3], double ddetA_dAOut[3]);

// follow the same def above
void d2detA_dA2(const double A[9], double d2detA_dA2Out[81]);
void d2detA_dA2_sym(const double A[6], double d2detA_dA2Out[36]);
void d2detA_dA2_diag(const double A[3], double d2detA_dA2Out[9]);
}  // namespace Dim3

namespace Dim2
{
// here if it is a matrix, then it will be column major
double det(const double M[4]);
// A = { A00, A10, A01, A11 }
void ddetA_dA(const double A[4], double ddetA_dAOut[4]);
// A = { A00, A01, A11 }
void ddetA_dA_sym(const double A[3], double ddetA_dAOut[3]);
// A = { A00, A11 }
void ddetA_dA_diag(const double A[2], double ddetA_dAOut[2]);

// follow the same def above
void d2detA_dA2(const double A[4], double d2detA_dA2Out[16]);
void d2detA_dA2_sym(const double A[3], double d2detA_dA2Out[9]);
void d2detA_dA2_diag(const double A[2], double d2detA_dA2Out[4]);
}  // namespace Dim2

}  // namespace Determinant
}  // namespace NonlinearOptimization
}  // namespace pgo