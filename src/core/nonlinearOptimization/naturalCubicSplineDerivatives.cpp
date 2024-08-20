/*
author: Bohan Wang
copyright to USC
*/

#include "naturalCubicSplineDerivatives.h"

#include "pgoLogging.h"

using namespace pgo::NonlinearOptimization;

NaturalCubicSpline2DWithParameterDerivatives::NaturalCubicSpline2DWithParameterDerivatives(int numPoints, const double *xValues, const double *yValues)
{
  n = numPoints;
  xNodeValue = ES::Mp<const ES::VXd>(xValues, n);

  if (yValues)
    yNodeValue = ES::Mp<const ES::VXd>(yValues, n);
  else
    yNodeValue.setZero(n);

  A.setZero(4 * n - 4, 4 * n - 4);
  rhs.setZero(4 * n - 4);

  int inc = 0;
  // f_i(t) = a_i t^3 + b_i t^2 + c_i t + d_i

  // f_i(0) = y_i, i = 0...n-2
  // d_i = y_i
  for (int i = 0; i < n - 1; i++) {
    A(inc, 4 * i + 3) = 1.0;
    rhs[inc] = yNodeValue[inc];
    inc++;
  }

  // f_{n-2}(1) = y_{n-1}
  // a_n-2 + b_n-2 + c_n-2 + d_n-2 = y_{n-1}
  A(inc, 4 * (n - 2)) = 1.0;
  A(inc, 4 * (n - 2) + 1) = 1.0;
  A(inc, 4 * (n - 2) + 2) = 1.0;
  A(inc, 4 * (n - 2) + 3) = 1.0;
  rhs[inc] = yNodeValue[inc];
  inc++;

  // f_i(1)=f_{i+1}(0), i = 0...n-3
  // a_i + b_i + c_i + d_i = d_i+1
  for (int i = 0; i < n - 2; i++) {
    A(inc, 4 * i) = 1.0;
    A(inc, 4 * i + 1) = 1.0;
    A(inc, 4 * i + 2) = 1.0;
    A(inc, 4 * i + 3) = 1.0;
    A(inc, 4 * (i + 1) + 3) = -1.0;
    rhs[inc] = 0.0;
    inc++;
  }

  // f_i'(t) = 3a_i t^2 + 2 b_i t + c_i

  // f_i'(1)=f_{i+1}'(0), i = 0...n-3
  // 3 a_i + 2 b_i + c_i = c_i+1
  for (int i = 0; i < n - 2; i++) {
    A(inc, 4 * i) = 3.0;
    A(inc, 4 * i + 1) = 2.0;
    A(inc, 4 * i + 2) = 1.0;
    A(inc, 4 * (i + 1) + 2) = -1.0;
    rhs[inc] = 0.0;
    inc++;
  }

  // f_i''(t) = 6a_i t + 2 b_i

  // f_i''(1)=f_{i+1}''(0), i = 0...n-3
  // 6 a_i + 2 b_i  = 2 b_i+1
  for (int i = 0; i < n - 2; i++) {
    A(inc, 4 * i) = 6.0;
    A(inc, 4 * i + 1) = 2.0;
    A(inc, 4 * (i + 1) + 1) = -2.0;
    rhs[inc] = 0.0;
    inc++;
  }

  // f_0''(0) = 0
  // 2 b_0 = 0
  A(inc, 1) = 2.0;
  rhs[inc] = 0.0;
  inc++;

  // f_n-2''(1) = 0
  // 6 a_{n-2} + 2 b_{n-2} = 0
  A(inc, 4 * (n - 2)) = 6.0;
  A(inc, 4 * (n - 2) + 1) = 2.0;

  rhs[inc] = 0.0;
  inc++;

  PGO_ALOG(inc == 4 * n - 4);

  AInv = A.fullPivHouseholderQr().inverse();
}

double NaturalCubicSpline2DWithParameterDerivatives::y(double x, const double *yValues) const
{
  const double *yValuePtr = yNodeValue.data();
  if (yValues) {
    yValuePtr = yValues;
  }

  // if x is left outside of the spline
  if (x < xNodeValue[0]) {
    ES::V4d coeff = AInv.topRows<4>().leftCols(n) * ES::Mp<const ES::VXd>(yValuePtr, n);

    // (x - x0) / (x1 - x0) * (t1 - t0) = t
    // dt/dx = (t1 - t0) / (x1 - x0)

    // f_i'(t) = 3a_i t^2 + 2 b_i t + c_i
    double k = coeff[2] / (xNodeValue[1] - xNodeValue[0]);

    // k x0 + b = y0
    double b = yValuePtr[0] - k * xNodeValue[0];

    return k * x + b;
  }
  else if (x > xNodeValue[n - 1]) {
    ES::V4d coeff = AInv.bottomRows<4>().leftCols(n) * ES::Mp<const ES::VXd>(yValuePtr, n);

    double k = (3 * coeff[0] + 2 * coeff[1] + coeff[2]) / (xNodeValue[n - 1] - xNodeValue[n - 2]);

    // k x0 + b = y0
    double b = yValuePtr[n - 1] - k * xNodeValue[n - 1];

    return k * x + b;
  }
  else {
    auto iter = std::upper_bound(xNodeValue.data(), xNodeValue.data() + xNodeValue.size(), x);
    std::ptrdiff_t i = iter - xNodeValue.data();

    if (i == xNodeValue.size())
      i--;

    ES::V4d coeff = AInv.middleRows<4>((i - 1) * 4).leftCols(n) * ES::Mp<const ES::VXd>(yValuePtr, n);

    double x0 = xNodeValue[i - 1];
    double x1 = xNodeValue[i];
    double t = (x - x0) / (x1 - x0);

    return coeff[0] * t * t * t + coeff[1] * t * t + coeff[2] * t + coeff[3];
  }
}

double NaturalCubicSpline2DWithParameterDerivatives::dy_dx(double x, const double *yValues) const
{
  const double *yValuePtr = yNodeValue.data();
  if (yValues) {
    yValuePtr = yValues;
  }

  // if x is left outside of the spline
  if (x < xNodeValue[0]) {
    ES::V4d coeff = AInv.topRows<4>().leftCols(n) * ES::Mp<const ES::VXd>(yValuePtr, n);

    // (x - x0) / (x1 - x0) * (t1 - t0) = t
    // dt/dx = (t1 - t0) / (x1 - x0)

    // f_i'(t) = 3a_i t^2 + 2 b_i t + c_i
    double k = coeff[2] / (xNodeValue[1] - xNodeValue[0]);
    return k;
  }
  else if (x > xNodeValue[n - 1]) {
    ES::V4d coeff = AInv.bottomRows<4>().leftCols(n) * ES::Mp<const ES::VXd>(yValuePtr, n);

    double k = (3 * coeff[0] + 2 * coeff[1] + coeff[2]) / (xNodeValue[n - 1] - xNodeValue[n - 2]);
    return k;
  }
  else {
    auto iter = std::upper_bound(xNodeValue.data(), xNodeValue.data() + xNodeValue.size(), x);
    std::ptrdiff_t i = iter - xNodeValue.data();

    if (i == xNodeValue.size())
      i--;

    ES::V4d coeff = AInv.middleRows<4>((i - 1) * 4).leftCols(n) * ES::Mp<const ES::VXd>(yValuePtr, n);

    double x0 = xNodeValue[i - 1];
    double x1 = xNodeValue[i];
    double t = (x - x0) / (x1 - x0);

    return (3.0 * coeff[0] * t * t + 2.0 * coeff[1] * t + coeff[2]) / (x1 - x0);
  }
}

double NaturalCubicSpline2DWithParameterDerivatives::d2y_dx2(double x, const double *yValues) const
{
  const double *yValuePtr = yNodeValue.data();
  if (yValues) {
    yValuePtr = yValues;
  }

  // if x is left outside of the spline
  if (x < xNodeValue[0]) {
    return 0;
  }
  else if (x > xNodeValue[n - 1]) {
    return 0;
  }
  else {
    auto iter = std::upper_bound(xNodeValue.data(), xNodeValue.data() + xNodeValue.size(), x);
    std::ptrdiff_t i = iter - xNodeValue.data();

    if (i == xNodeValue.size())
      i--;

    ES::V4d coeff = AInv.middleRows<4>((i - 1) * 4).leftCols(n) * ES::Mp<const ES::VXd>(yValuePtr, n);

    double x0 = xNodeValue[i - 1];
    double x1 = xNodeValue[i];
    double t = (x - x0) / (x1 - x0);

    return (6.0 * coeff[0] * t + 2.0 * coeff[1]) / ((x1 - x0) * (x1 - x0));
  }
}

void NaturalCubicSpline2DWithParameterDerivatives::dy_dparam(double x, const double *yValues, double *dparam) const
{
  // if x is left outside of the spline
  if (x < xNodeValue[0]) {
    ES::V4d dk_dabcd(0.0, 0.0, 1.0, 0.0);
    ES::V4d coeff = dk_dabcd * (x - xNodeValue[0]) / (xNodeValue[1] - xNodeValue[0]);

    (ES::Mp<ES::VXd>(dparam, n)) = AInv.block(0, 0, 4, n).transpose() * coeff;
    dparam[0] += 1;
  }
  else if (x > xNodeValue[n - 1]) {
    double x0 = xNodeValue[n - 2];
    double x1 = xNodeValue[n - 1];

    ES::V4d dk_dabcd(3.0, 2.0, 1.0, 0.0);
    ES::V4d coeff = dk_dabcd * (x - x1) / (x1 - x0);

    (ES::Mp<ES::VXd>(dparam, n)) = AInv.block((n - 2) * 4, 0, 4, n).transpose() * coeff;
    dparam[n - 1] += 1;

    // k x0 + b = y0
    // double b = yValuePtr[n - 1] - k * xNodeValue[n - 1];
  }
  else {
    auto iter = std::upper_bound(xNodeValue.data(), xNodeValue.data() + xNodeValue.size(), x);
    std::ptrdiff_t i = iter - xNodeValue.data();

    if (i == xNodeValue.size())
      i--;

    double x0 = xNodeValue[i - 1];
    double x1 = xNodeValue[i];
    double t = (x - x0) / (x1 - x0);

    ES::V4d coeff(t * t * t, t * t, t, 1);

    (ES::Mp<ES::VXd>(dparam, n)) = AInv.block((i - 1) * 4, 0, 4, n).transpose() * coeff;
  }
}

void NaturalCubicSpline2DWithParameterDerivatives::d2y_dparam2(double x, const double *yValues, double *dparam) const
{
  (ES::Mp<ES::MXd>(dparam, n, n)).setZero();
}

void NaturalCubicSpline2DWithParameterDerivatives::d2y_dparam_dx(double x, const double *yValues, double *dparam) const
{
  // if x is left outside of the spline
  if (x < xNodeValue[0]) {
    ES::V4d dk_dabcd(0.0, 0.0, 1.0, 0.0);
    ES::V4d coeff = dk_dabcd / (xNodeValue[1] - xNodeValue[0]);

    (ES::Mp<ES::VXd>(dparam, n)) = AInv.block(0, 0, 4, n).transpose() * coeff;
  }
  else if (x > xNodeValue[n - 1]) {
    double x0 = xNodeValue[n - 2];
    double x1 = xNodeValue[n - 1];

    ES::V4d dk_dabcd(3.0, 2.0, 1.0, 0.0);
    ES::V4d coeff = dk_dabcd / (x1 - x0);

    (ES::Mp<ES::VXd>(dparam, n)) = AInv.block((n - 2) * 4, 0, 4, n).transpose() * coeff;
  }
  else {
    auto iter = std::upper_bound(xNodeValue.data(), xNodeValue.data() + xNodeValue.size(), x);
    std::ptrdiff_t i = iter - xNodeValue.data();

    if (i == xNodeValue.size())
      i--;

    double x0 = xNodeValue[i - 1];
    double x1 = xNodeValue[i];
    double t = (x - x0) / (x1 - x0);

    ES::V4d coeff(3 * t * t, 2 * t, 1, 0);
    coeff /= (x1 - x0);

    (ES::Mp<ES::VXd>(dparam, n)) = AInv.block((i - 1) * 4, 0, 4, n).transpose() * coeff;
  }
}

double NaturalCubicSpline2DWithParameterDerivatives::dy_dparam(double x, const double *, int parami) const
{
  // if x is left outside of the spline
  if (x < xNodeValue[0]) {
    ES::V4d dk_dabcd(0.0, 0.0, 1.0, 0.0);
    ES::V4d coeff = dk_dabcd * (x - xNodeValue[0]) / (xNodeValue[1] - xNodeValue[0]);

    double val = AInv.block<4, 1>(0, parami).dot(coeff);

    if (parami == 0)
      val += 1;

    return val;
  }
  else if (x > xNodeValue[n - 1]) {
    double x0 = xNodeValue[n - 2];
    double x1 = xNodeValue[n - 1];

    ES::V4d dk_dabcd(3.0, 2.0, 1.0, 0.0);
    ES::V4d coeff = dk_dabcd * (x - x1) / (x1 - x0);

    double val = AInv.block<4, 1>((n - 2) * 4, parami).dot(coeff);

    if (parami == n - 1)
      val += 1;

    return val;
  }
  else {
    auto iter = std::upper_bound(xNodeValue.data(), xNodeValue.data() + xNodeValue.size(), x);
    std::ptrdiff_t i = iter - xNodeValue.data();

    if (i == xNodeValue.size())
      i--;

    double x0 = xNodeValue[i - 1];
    double x1 = xNodeValue[i];
    double t = (x - x0) / (x1 - x0);

    ES::V4d coeff(t * t * t, t * t, t, 1);

    return AInv.block<4, 1>((i - 1) * 4, parami).dot(coeff);
  }
}

double NaturalCubicSpline2DWithParameterDerivatives::d2y_dparam2(double, const double *, int, int) const
{
  return 0;
}

double NaturalCubicSpline2DWithParameterDerivatives::d2y_dparam_dx(double x, const double *, int parami) const
{
  // if x is left outside of the spline
  if (x < xNodeValue[0]) {
    ES::V4d dk_dabcd(0.0, 0.0, 1.0, 0.0);
    ES::V4d coeff = dk_dabcd / (xNodeValue[1] - xNodeValue[0]);
    double val = AInv.block<4, 1>(0, parami).dot(coeff);
    return val;
  }
  else if (x > xNodeValue[n - 1]) {
    double x0 = xNodeValue[n - 2];
    double x1 = xNodeValue[n - 1];

    ES::V4d dk_dabcd(3.0, 2.0, 1.0, 0.0);
    ES::V4d coeff = dk_dabcd / (x1 - x0);

    double val = AInv.block<4, 1>((n - 2) * 4, parami).dot(coeff);
    return val;
  }
  else {
    auto iter = std::upper_bound(xNodeValue.data(), xNodeValue.data() + xNodeValue.size(), x);
    std::ptrdiff_t offsetIdx = iter - xNodeValue.data();

    if (offsetIdx == xNodeValue.size())
      offsetIdx--;

    double x0 = xNodeValue[offsetIdx - 1];
    double x1 = xNodeValue[offsetIdx];
    double t = (x - x0) / (x1 - x0);

    ES::V4d coeff(3 * t * t, 2 * t, 1, 0);
    coeff /= (x1 - x0);

    return AInv.block<4, 1>((offsetIdx - 1) * 4, parami).dot(coeff);
  }
}

NaturalCubicSpline2DAsNonlinearConstraints::NaturalCubicSpline2DAsNonlinearConstraints(int np, int xoff, int yoff, int poff, int na, double y0Prime, double ynPrime):
  ConstraintFunctions(na), numPoints(np), xDOFOffset(xoff), yDOFOffset(yoff), paramDOFOffset(poff), nAll(na)
{
  endDeriv[0] = y0Prime;
  endDeriv[1] = ynPrime;

  // create jacobian
  int si = 0;
  int n = numPoints - 1;
  std::vector<ES::TripletD> entries;

  // ===========================
  // f0 (x0) - y0 = 0
  // a0 x0^3 + b0 * x0^2 + c0 * x0 + d0

  // dx0 = a0 * 3 x0^2 + b0 * 2 x0 + c0
  entries.emplace_back(si, xDOFOffset, 1.0);

  // da0 = x0^3, db0 = x0^2, dc0 = x0, dd0 = 1
  for (int j = 0; j < 4; j++)
    entries.emplace_back(si, paramDOFOffset + j, 1.0);

  // dy0 = -1
  entries.emplace_back(si, yDOFOffset, 1.0);

  si++;

  // ===========================
  // fn-1 (xn) - yn = 0
  // dxn
  entries.emplace_back(si, xDOFOffset + n, 1.0);

  // dparam_n-1
  for (int j = 0; j < 4; j++)
    entries.emplace_back(si, paramDOFOffset + (n - 1) * 4 + j, 1.0);

  // dyn
  entries.emplace_back(si, yDOFOffset + n, 1.0);

  si++;

  // ===========================
  // fi (x_i+1) - f_{i+1}(x_i+1) = 0
  for (int i = 0; i < numPoints - 2; i++) {
    // dx_{i+1}
    entries.emplace_back(si, xDOFOffset + i + 1, 1.0);

    // dparam_i
    for (int j = 0; j < 4; j++)
      entries.emplace_back(si, paramDOFOffset + i * 4 + j, 1.0);

    for (int j = 0; j < 4; j++)
      entries.emplace_back(si, paramDOFOffset + i * 4 + 4 + j, 1.0);

    si++;
  }

  // ==============================
  // fi (x_i+1)  = y_i+1
  for (int i = 0; i < n - 1; i++) {
    // dx_{i+1}
    entries.emplace_back(si, xDOFOffset + i + 1, 1.0);

    // dparam_i
    for (int j = 0; j < 4; j++)
      entries.emplace_back(si, paramDOFOffset + i * 4 + j, 1.0);

    // dy_{i+1}
    entries.emplace_back(si, yDOFOffset + i + 1, 1.0);

    si++;
  }

  // ==============================
  // f0'(x0) = y0'
  // 3 a0 x0^2 + 2 b0 * x0 + c0 - y0' = 0

  // dx0
  entries.emplace_back(si, xDOFOffset, 1.0);

  // dparam0
  for (int j = 0; j < 4; j++)
    entries.emplace_back(si, paramDOFOffset + j, 1.0);

  si++;

  // ==============================
  // fn-1'(xn) = yn'

  // dxn
  entries.emplace_back(si, xDOFOffset + n, 1.0);

  // dparamn
  for (int j = 0; j < 4; j++)
    entries.emplace_back(si, paramDOFOffset + (n - 1) * 4 + j, 1.0);

  si++;

  // ==============================
  // fi' (x_i+1) - f'_{i+1}(x_i+1) = 0
  for (int i = 0; i < numPoints - 2; i++) {
    // dx_{i+1}
    entries.emplace_back(si, xDOFOffset + i + 1, 1.0);

    // dparam_i
    for (int j = 0; j < 4; j++)
      entries.emplace_back(si, paramDOFOffset + i * 4 + j, 1.0);

    for (int j = 0; j < 4; j++)
      entries.emplace_back(si, paramDOFOffset + i * 4 + 4 + j, 1.0);

    si++;
  }

  // fi'' (x_i+1) - f''_{i+1}(x_i+1) = 0
  // fi'' = 6 a_i x + 2 b_i
  for (int i = 0; i < numPoints - 2; i++) {
    // dx_{i+1}
    entries.emplace_back(si, xDOFOffset + i + 1, 1.0);

    // dparam_i
    for (int j = 0; j < 4; j++)
      entries.emplace_back(si, paramDOFOffset + i * 4 + j, 1.0);

    for (int j = 0; j < 4; j++)
      entries.emplace_back(si, paramDOFOffset + i * 4 + 4 + j, 1.0);

    si++;
  }

  PGO_ALOG(si == 4 * numPoints - 4);

  jacobianTemplate.resize(4 * numPoints - 4, nAll);
  jacobianTemplate.setFromTriplets(entries.begin(), entries.end());

  // ==================================================
  // hessian
  si = 0;
  entries.clear();

  // ===========================
  // f0 (x0) - y0 = 0
  // dx0 dx0
  entries.emplace_back(xDOFOffset, xDOFOffset, 1.0);
  // dx0 da0 = 3 x0^2, dx0 db0 = 2 x0, dx0 dc0 = 1, dx0 dd0 = 0
  for (int j = 0; j < 4; j++) {
    entries.emplace_back(xDOFOffset, paramDOFOffset + j, 1.0);
    entries.emplace_back(paramDOFOffset + j, xDOFOffset, 1.0);
  }

  // ===========================
  // fn-1 (xn) - yn = 0
  // an-1 xn^3 + bn-1 * xn^2 + cn-1 * xn + dn-1 - yn = 0
  // dxn dxn
  entries.emplace_back(xDOFOffset + n, xDOFOffset + n, 1.0);
  // dxn dparam_n-1
  for (int j = 0; j < 4; j++) {
    entries.emplace_back(xDOFOffset + n, paramDOFOffset + (n - 1) * 4 + j, 1.0);
    entries.emplace_back(paramDOFOffset + (n - 1) * 4 + j, xDOFOffset + n, 1.0);
  }

  // fi (x_i+1) - f_{i+1}(x_i+1) = 0
  for (int i = 0; i < numPoints - 2; i++) {
    // dx_{i+1} dx_{i+1}
    entries.emplace_back(xDOFOffset + i + 1, xDOFOffset + i + 1, 1.0);

    // dx_{i+1} dparam_i
    for (int j = 0; j < 4; j++) {
      entries.emplace_back(xDOFOffset + i + 1, paramDOFOffset + i * 4 + j, 1.0);
      entries.emplace_back(paramDOFOffset + i * 4 + j, xDOFOffset + i + 1, 1.0);

      entries.emplace_back(xDOFOffset + i + 1, paramDOFOffset + 4 + i * 4 + j, 1.0);
      entries.emplace_back(paramDOFOffset + 4 + i * 4 + j, xDOFOffset + i + 1, 1.0);
    }
  }

  lambdahTemplate.resize(nAll, nAll);
  lambdahTemplate.setFromTriplets(entries.begin(), entries.end());
}

void NaturalCubicSpline2DAsNonlinearConstraints::solveParams(ES::RefVecXd x) const
{
  int si = 0;
  int n = numPoints - 1;

  std::vector<ES::TripletD> entries;
  ES::VXd rhs(4 * n);

  // f0 (x0) = y0
  // x0^3 a0 + x0^2 b0 + x0 c0 + d0 = y0
  double x0 = getxi(x, 0);
  entries.emplace_back(si, 0, x0 * x0 * x0);
  entries.emplace_back(si, 1, x0 * x0);
  entries.emplace_back(si, 2, x0);
  entries.emplace_back(si, 3, 1.0);

  rhs[si] = getyi(x, 0);

  si++;

  // fn-1 (xn) = yn
  double xn = getxi(x, n);
  entries.emplace_back(si, (n - 1) * 4, xn * xn * xn);
  entries.emplace_back(si, (n - 1) * 4 + 1, xn * xn);
  entries.emplace_back(si, (n - 1) * 4 + 2, xn);
  entries.emplace_back(si, (n - 1) * 4 + 3, 1.0);

  rhs[si] = getyi(x, n);

  si++;

  // fi (x_i+1) - f_{i+1}(x_i+1) = 0
  for (int i = 0; i < n - 1; i++) {
    double xi1 = getxi(x, i + 1);

    entries.emplace_back(si, i * 4, xi1 * xi1 * xi1);
    entries.emplace_back(si, i * 4 + 1, xi1 * xi1);
    entries.emplace_back(si, i * 4 + 2, xi1);
    entries.emplace_back(si, i * 4 + 3, 1.0);

    entries.emplace_back(si, i * 4 + 4, -xi1 * xi1 * xi1);
    entries.emplace_back(si, i * 4 + 4 + 1, -xi1 * xi1);
    entries.emplace_back(si, i * 4 + 4 + 2, -xi1);
    entries.emplace_back(si, i * 4 + 4 + 3, -1.0);

    rhs[si] = 0;

    si++;
  }

  // fi (x_i+1)  = y_i+1
  for (int i = 0; i < numPoints - 2; i++) {
    double xi1 = getxi(x, i + 1);

    entries.emplace_back(si, i * 4, xi1 * xi1 * xi1);
    entries.emplace_back(si, i * 4 + 1, xi1 * xi1);
    entries.emplace_back(si, i * 4 + 2, xi1);
    entries.emplace_back(si, i * 4 + 3, 1.0);

    rhs[si] = getyi(x, i + 1);

    si++;
  }

  // f0'(x0) = y0'
  entries.emplace_back(si, 0, 3 * x0 * x0);
  entries.emplace_back(si, 1, 2 * x0);
  entries.emplace_back(si, 2, 1.0);

  rhs[si] = endDeriv[0];

  si++;

  // fn-1'(xn) = yn'
  entries.emplace_back(si, (n - 1) * 4, 3 * xn * xn);
  entries.emplace_back(si, (n - 1) * 4 + 1, 2 * xn);
  entries.emplace_back(si, (n - 1) * 4 + 2, 1.0);

  rhs[si] = endDeriv[1];

  si++;

  // fi' (x_i+1) - f'_{i+1}(x_i+1) = 0
  for (int i = 0; i < numPoints - 2; i++) {
    double xi1 = getxi(x, i + 1);

    entries.emplace_back(si, i * 4, 3.0 * xi1 * xi1);
    entries.emplace_back(si, i * 4 + 1, 2.0 * xi1);
    entries.emplace_back(si, i * 4 + 2, 1.0);

    entries.emplace_back(si, i * 4 + 4, -3.0 * xi1 * xi1);
    entries.emplace_back(si, i * 4 + 4 + 1, -2.0 * xi1);
    entries.emplace_back(si, i * 4 + 4 + 2, -1.0);

    rhs[si] = 0;

    si++;
  }

  // fi'' (x_i+1) - f''_{i+1}(x_i+1) = 0
  for (int i = 0; i < numPoints - 2; i++) {
    double xi1 = getxi(x, i + 1);

    entries.emplace_back(si, i * 4, 6.0 * xi1);
    entries.emplace_back(si, i * 4 + 1, 2.0);

    entries.emplace_back(si, i * 4 + 4, -6.0 * xi1);
    entries.emplace_back(si, i * 4 + 4 + 1, -2.0);

    rhs[si] = 0;

    si++;
  }

  PGO_ALOG(si == 4 * numPoints - 4);

  ES::SpMatD sys(4 * n, 4 * n);
  sys.setFromTriplets(entries.begin(), entries.end());

  ES::LUSolver solver;
  solver.analyzePattern(sys);
  solver.factorize(sys);

  x.segment(paramDOFOffset, 4 * n) = solver.solve(rhs);
}

void NaturalCubicSpline2DAsNonlinearConstraints::func(ES::ConstRefVecXd x, ES::RefVecXd g) const
{
  int si = 0;

  // f0 (x0) = y0
  g[si++] = compute_fi(getParami(x, 0), getxi(x, 0)) - getyi(x, 0);

  // fn-1 (xn) = yn
  g[si++] = compute_fi(getParami(x, numPoints - 2), getxi(x, numPoints - 1)) - getyi(x, numPoints - 1);

  // fi (x_i+1) - f_{i+1}(x_i+1) = 0
  for (int i = 0; i < numPoints - 2; i++) {
    g[si++] = compute_fi(getParami(x, i), getxi(x, i + 1)) - compute_fi(getParami(x, i + 1), getxi(x, i + 1));
  }

  // fi (x_i+1)  = y_i+1
  for (int i = 0; i < numPoints - 2; i++) {
    g[si++] = compute_fi(getParami(x, i), getxi(x, i + 1)) - getyi(x, i + 1);
  }

  // f0'(x0) = y0'
  g[si++] = compute_fi_1stDeriv(getParami(x, 0), getxi(x, 0)) - endDeriv[0];

  // fn-1'(xn) = yn'
  g[si++] = compute_fi_1stDeriv(getParami(x, numPoints - 2), getxi(x, numPoints - 1)) - endDeriv[1];

  // fi' (x_i+1) - f'_{i+1}(x_i+1) = 0
  for (int i = 0; i < numPoints - 2; i++) {
    g[si++] = compute_fi_1stDeriv(getParami(x, i), getxi(x, i + 1)) - compute_fi_1stDeriv(getParami(x, i + 1), getxi(x, i + 1));
  }

  // fi'' (x_i+1) - f''_{i+1}(x_i+1) = 0
  for (int i = 0; i < numPoints - 2; i++) {
    g[si++] = compute_fi_2ndDeriv(getParami(x, i), getxi(x, i + 1)) - compute_fi_2ndDeriv(getParami(x, i + 1), getxi(x, i + 1));
  }

  PGO_ALOG(si == 4 * numPoints - 4);
}

void NaturalCubicSpline2DAsNonlinearConstraints::jacobian(ES::ConstRefVecXd x, ES::SpMatD &jac) const
{
  int si = 0;
  int n = numPoints - 1;
  double *jacVals = jac.valuePtr();
  std::ptrdiff_t offset = -1;

  // ===========================
  // f0 (x0) - y0 = 0
  // a0 x0^3 + b0 * x0^2 + c0 * x0 + d0

  // dx0 = a0 * 3 x0^2 + b0 * 2 x0 + c0
  double x0 = getxi(x, 0);
  offset = ES::findEntryOffset(jac, si, xDOFOffset);
  PGO_ALOG(offset >= 0);

  jacVals[offset] = compute_fi_1stDeriv(getParami(x, 0), x0);

  // da0 = x0^3, db0 = x0^2, dc0 = x0, dd0 = 1
  offset = ES::findEntryOffset(jac, si, paramDOFOffset);
  PGO_ALOG(offset >= 0);

  jacVals[offset] = x0 * x0 * x0;
  jacVals[offset + 1] = x0 * x0;
  jacVals[offset + 2] = x0;
  jacVals[offset + 3] = 1.0;

  // dy0 = -1
  offset = ES::findEntryOffset(jac, si, yDOFOffset);
  PGO_ALOG(offset >= 0);

  jacVals[offset] = -1.0;

  si++;

  // ===========================
  // fn-1 (xn) - yn = 0
  // an-1 xn^3 + bn-1 * xn^2 + cn-1 * xn + dn-1 - yn = 0
  double xn = getxi(x, n);
  // dxn
  offset = ES::findEntryOffset(jac, si, xDOFOffset + n);
  PGO_ALOG(offset >= 0);

  jacVals[offset] = compute_fi_1stDeriv(getParami(x, n - 1), xn);

  // dparam_n-1
  offset = ES::findEntryOffset(jac, si, paramDOFOffset + (n - 1) * 4);
  PGO_ALOG(offset >= 0);

  jacVals[offset] = xn * xn * xn;
  jacVals[offset + 1] = xn * xn;
  jacVals[offset + 2] = xn;
  jacVals[offset + 3] = 1.0;

  // dyn
  offset = ES::findEntryOffset(jac, si, yDOFOffset + n);
  PGO_ALOG(offset >= 0);

  jacVals[offset] = -1.0;

  si++;

  // fi (x_i+1) - f_{i+1}(x_i+1) = 0
  for (int i = 0; i < numPoints - 2; i++) {
    double xi1 = getxi(x, i + 1);

    // dx_{i+1}
    offset = ES::findEntryOffset(jac, si, xDOFOffset + i + 1);
    PGO_ALOG(offset >= 0);

    jacVals[offset] = compute_fi_1stDeriv(getParami(x, i), xi1) - compute_fi_1stDeriv(getParami(x, i + 1), xi1);

    // dparam_i
    offset = ES::findEntryOffset(jac, si, paramDOFOffset + i * 4);
    PGO_ALOG(offset >= 0);

    jacVals[offset] = xi1 * xi1 * xi1;
    jacVals[offset + 1] = xi1 * xi1;
    jacVals[offset + 2] = xi1;
    jacVals[offset + 3] = 1.0;

    // dparam_{i + 1}
    offset = ES::findEntryOffset(jac, si, paramDOFOffset + i * 4 + 4);
    PGO_ALOG(offset >= 0);

    jacVals[offset] = -xi1 * xi1 * xi1;
    jacVals[offset + 1] = -xi1 * xi1;
    jacVals[offset + 2] = -xi1;
    jacVals[offset + 3] = -1.0;

    si++;
  }

  // ==============================
  // fi (x_i+1)  = y_i+1
  for (int i = 0; i < n - 1; i++) {
    double xi1 = getxi(x, i + 1);
    // dx_{i+1}
    offset = ES::findEntryOffset(jac, si, xDOFOffset + i + 1);
    PGO_ALOG(offset >= 0);

    jacVals[offset] = compute_fi_1stDeriv(getParami(x, i), xi1);

    // dparam_i
    offset = ES::findEntryOffset(jac, si, paramDOFOffset + i * 4);
    PGO_ALOG(offset >= 0);

    jacVals[offset] = xi1 * xi1 * xi1;
    jacVals[offset + 1] = xi1 * xi1;
    jacVals[offset + 2] = xi1;
    jacVals[offset + 3] = 1.0;

    // dy_{i+1}
    offset = ES::findEntryOffset(jac, si, yDOFOffset + i + 1);
    PGO_ALOG(offset >= 0);

    jacVals[offset] = -1.0;

    si++;
  }

  // ==============================
  // f0'(x0) = y0'
  // 3 a0 x0^2 + 2 b0 * x0 + c0 - y0' = 0

  // dx0
  offset = ES::findEntryOffset(jac, si, xDOFOffset);
  PGO_ALOG(offset >= 0);

  jacVals[offset] = compute_fi_2ndDeriv(getParami(x, 0), x0);

  // dparam0
  offset = ES::findEntryOffset(jac, si, paramDOFOffset);
  PGO_ALOG(offset >= 0);

  jacVals[offset] = 3 * x0 * x0;
  jacVals[offset + 1] = 2 * x0;
  jacVals[offset + 2] = 1.0;
  jacVals[offset + 3] = 0;

  si++;

  // ==============================
  // fn-1'(xn) = yn'

  // dxn
  offset = ES::findEntryOffset(jac, si, xDOFOffset + n);
  PGO_ALOG(offset >= 0);

  jacVals[offset] = compute_fi_2ndDeriv(getParami(x, n - 1), xn);

  // dparamn
  offset = ES::findEntryOffset(jac, si, paramDOFOffset + (n - 1) * 4);
  PGO_ALOG(offset >= 0);

  jacVals[offset] = 3 * xn * xn;
  jacVals[offset + 1] = 2 * xn;
  jacVals[offset + 2] = 1.0;
  jacVals[offset + 3] = 0;

  si++;

  // fi' (x_i+1) - f'_{i+1}(x_i+1) = 0
  for (int i = 0; i < numPoints - 2; i++) {
    double xi1 = getxi(x, i + 1);

    // dx_{i+1}
    offset = ES::findEntryOffset(jac, si, xDOFOffset + i + 1);
    PGO_ALOG(offset >= 0);

    jacVals[offset] = compute_fi_2ndDeriv(getParami(x, i), xi1) - compute_fi_2ndDeriv(getParami(x, i + 1), xi1);

    // dparam_i
    offset = ES::findEntryOffset(jac, si, paramDOFOffset + i * 4);
    PGO_ALOG(offset >= 0);

    jacVals[offset] = 3 * xi1 * xi1;
    jacVals[offset + 1] = 2 * xi1;
    jacVals[offset + 2] = 1;
    jacVals[offset + 3] = 0;

    // dparam_{i + 1}
    offset = ES::findEntryOffset(jac, si, paramDOFOffset + i * 4 + 4);
    PGO_ALOG(offset >= 0);

    jacVals[offset] = -3 * xi1 * xi1;
    jacVals[offset + 1] = -2 * xi1;
    jacVals[offset + 2] = -1.0;
    jacVals[offset + 3] = 0;

    si++;
  }

  // fi'' (x_i+1) - f''_{i+1}(x_i+1) = 0
  // fi'' = 6 a_i x + 2 b_i
  for (int i = 0; i < numPoints - 2; i++) {
    double xi1 = getxi(x, i + 1);
    ES::V4d parami = getParami(x, i);
    ES::V4d parami1 = getParami(x, i + 1);

    // dx_{i+1}
    offset = ES::findEntryOffset(jac, si, xDOFOffset + i + 1);
    PGO_ALOG(offset >= 0);

    jacVals[offset] = 6 * parami[0] - 6 * parami1[0];

    // dparam_i
    offset = ES::findEntryOffset(jac, si, paramDOFOffset + i * 4);
    PGO_ALOG(offset >= 0);

    jacVals[offset] = 6.0 * xi1;
    jacVals[offset + 1] = 2;
    jacVals[offset + 2] = 0;
    jacVals[offset + 3] = 0;

    // dparam_{i + 1}
    offset = ES::findEntryOffset(jac, si, paramDOFOffset + i * 4 + 4);
    PGO_ALOG(offset >= 0);

    jacVals[offset] = -6 * xi1;
    jacVals[offset + 1] = -2;
    jacVals[offset + 2] = 0;
    jacVals[offset + 3] = 0;

    si++;
  }

  PGO_ALOG(si == 4 * numPoints - 4);
}

void NaturalCubicSpline2DAsNonlinearConstraints::hessian(ES::ConstRefVecXd x, ES::ConstRefVecXd lambda, ES::SpMatD &hess) const
{
  int si = 0;
  int n = numPoints - 1;
  std::ptrdiff_t offset = -1;

  memset(hess.valuePtr(), 0, sizeof(double) * hess.nonZeros());

  // ===========================
  // f0 (x0) - y0 = 0
  // a0 x0^3 + b0 * x0^2 + c0 * x0 + d0

  // dx0 = a0 * 3 * x0^2 + b0 * 2 * x0 + c0

  // dx0 dx0 = a0 * 6 x0 + b0 * 2
  double x0 = getxi(x, 0);
  hess.coeffRef(xDOFOffset, xDOFOffset) += compute_fi_2ndDeriv(getParami(x, 0), x0) * lambda[si];

  // dx0 da0 = 3 x0^2, dx0 db0 = 2 x0, dx0 dc0 = 1, dx0 dd0 = 0
  hess.coeffRef(xDOFOffset, paramDOFOffset) += 3 * x0 * x0 * lambda[si];
  hess.coeffRef(xDOFOffset, paramDOFOffset + 1) += 2 * x0 * lambda[si];
  hess.coeffRef(xDOFOffset, paramDOFOffset + 2) += 1 * lambda[si];
  hess.coeffRef(xDOFOffset, paramDOFOffset + 3) += 0;

  hess.coeffRef(paramDOFOffset, xDOFOffset) += 3 * x0 * x0 * lambda[si];
  hess.coeffRef(paramDOFOffset + 1, xDOFOffset) += 2 * x0 * lambda[si];
  hess.coeffRef(paramDOFOffset + 2, xDOFOffset) += 1 * lambda[si];
  hess.coeffRef(paramDOFOffset + 3, xDOFOffset) += 0;

  PGO_ALOG(hess.isCompressed() == true);

  // dx0 dy0 = 0
  // dparam0 dy0 = 0
  // dparam0 dparam0 = 0

  si++;

  // ===========================
  // fn-1 (xn) - yn = 0
  // an-1 xn^3 + bn-1 * xn^2 + cn-1 * xn + dn-1 - yn = 0
  double xn = getxi(x, n);

  // dxn dxn
  hess.coeffRef(xDOFOffset + n, xDOFOffset + n) += compute_fi_2ndDeriv(getParami(x, n - 1), xn) * lambda[si];

  // dxn dparam_n-1
  hess.coeffRef(xDOFOffset + n, paramDOFOffset + (n - 1) * 4) += 3 * xn * xn * lambda[si];
  hess.coeffRef(xDOFOffset + n, paramDOFOffset + (n - 1) * 4 + 1) += 2 * xn * lambda[si];
  hess.coeffRef(xDOFOffset + n, paramDOFOffset + (n - 1) * 4 + 2) += 1;
  hess.coeffRef(xDOFOffset + n, paramDOFOffset + (n - 1) * 4 + 3) += 0;

  hess.coeffRef(paramDOFOffset + (n - 1) * 4, xDOFOffset + n) += 3 * xn * xn * lambda[si];
  hess.coeffRef(paramDOFOffset + (n - 1) * 4 + 1, xDOFOffset + n) += 2 * xn * lambda[si];
  hess.coeffRef(paramDOFOffset + (n - 1) * 4 + 2, xDOFOffset + n) += 1 * lambda[si];
  hess.coeffRef(paramDOFOffset + (n - 1) * 4 + 3, xDOFOffset + n) += 0;

  PGO_ALOG(hess.isCompressed() == true);

  si++;

  // fi (x_i+1) - f_{i+1}(x_i+1) = 0
  for (int i = 0; i < numPoints - 2; i++) {
    double xi1 = getxi(x, i + 1);

    // dx_{i+1} dx_{i+1}
    hess.coeffRef(xDOFOffset + i + 1, xDOFOffset + i + 1) += (compute_fi_2ndDeriv(getParami(x, i), xi1) - compute_fi_2ndDeriv(getParami(x, i + 1), xi1)) * lambda[si];

    // dx_{i+1} dparam_i
    hess.coeffRef(xDOFOffset + i + 1, paramDOFOffset + i * 4) += 3 * xi1 * xi1 * lambda[si];
    hess.coeffRef(xDOFOffset + i + 1, paramDOFOffset + i * 4 + 1) += 2 * xi1 * lambda[si];
    hess.coeffRef(xDOFOffset + i + 1, paramDOFOffset + i * 4 + 2) += 1 * lambda[si];
    hess.coeffRef(xDOFOffset + i + 1, paramDOFOffset + i * 4 + 3) += 0;

    hess.coeffRef(paramDOFOffset + i * 4, xDOFOffset + i + 1) += 3 * xi1 * xi1 * lambda[si];
    hess.coeffRef(paramDOFOffset + i * 4 + 1, xDOFOffset + i + 1) += 2 * xi1 * lambda[si];
    hess.coeffRef(paramDOFOffset + i * 4 + 2, xDOFOffset + i + 1) += 1 * lambda[si];
    hess.coeffRef(paramDOFOffset + i * 4 + 3, xDOFOffset + i + 1) += 0;

    // dx_{i+1} dparam_{i + 1}
    hess.coeffRef(xDOFOffset + i + 1, paramDOFOffset + i * 4 + 4) += -3 * xi1 * xi1 * lambda[si];
    hess.coeffRef(xDOFOffset + i + 1, paramDOFOffset + i * 4 + 4 + 1) += -2 * xi1 * lambda[si];
    hess.coeffRef(xDOFOffset + i + 1, paramDOFOffset + i * 4 + 4 + 2) += -1 * lambda[si];
    hess.coeffRef(xDOFOffset + i + 1, paramDOFOffset + i * 4 + 4 + 3) += 0;

    hess.coeffRef(paramDOFOffset + i * 4 + 4, xDOFOffset + i + 1) += -3 * xi1 * xi1 * lambda[si];
    hess.coeffRef(paramDOFOffset + i * 4 + 4 + 1, xDOFOffset + i + 1) += -2 * xi1 * lambda[si];
    hess.coeffRef(paramDOFOffset + i * 4 + 4 + 2, xDOFOffset + i + 1) += -1 * lambda[si];
    hess.coeffRef(paramDOFOffset + i * 4 + 4 + 3, xDOFOffset + i + 1) += 0;

    si++;
  }

  PGO_ALOG(hess.isCompressed() == true);

  // ==============================
  // fi (x_i+1)  = y_i+1
  for (int i = 0; i < n - 1; i++) {
    double xi1 = getxi(x, i + 1);

    // dx_{i+1} dx_{i+1}
    hess.coeffRef(xDOFOffset + i + 1, xDOFOffset + i + 1) += compute_fi_2ndDeriv(getParami(x, i), xi1) * lambda[si];

    // dx_{i+1} dparam_i
    hess.coeffRef(xDOFOffset + i + 1, paramDOFOffset + i * 4) += 3 * xi1 * xi1 * lambda[si];
    hess.coeffRef(xDOFOffset + i + 1, paramDOFOffset + i * 4 + 1) += 2 * xi1 * lambda[si];
    hess.coeffRef(xDOFOffset + i + 1, paramDOFOffset + i * 4 + 2) += 1 * lambda[si];
    hess.coeffRef(xDOFOffset + i + 1, paramDOFOffset + i * 4 + 3) += 0;

    hess.coeffRef(paramDOFOffset + i * 4, xDOFOffset + i + 1) += 3 * xi1 * xi1 * lambda[si];
    hess.coeffRef(paramDOFOffset + i * 4 + 1, xDOFOffset + i + 1) += 2 * xi1 * lambda[si];
    hess.coeffRef(paramDOFOffset + i * 4 + 2, xDOFOffset + i + 1) += 1 * lambda[si];
    hess.coeffRef(paramDOFOffset + i * 4 + 3, xDOFOffset + i + 1) += 0;

    si++;
  }

  PGO_ALOG(hess.isCompressed() == true);

  // ==============================
  // f0'(x0) = y0'
  // 3 a0 x0^2 + 2 b0 * x0 + c0 - y0' = 0

  // dx0 = 6 a0 x0 + 2 b0
  // dx0 dx0 = 6 a0
  ES::V4d param0 = getParami(x, 0);
  hess.coeffRef(xDOFOffset, xDOFOffset) += 6 * param0[0] * lambda[si];

  // dx0 dparam0
  hess.coeffRef(xDOFOffset, paramDOFOffset) += 6 * x0 * lambda[si];
  hess.coeffRef(xDOFOffset, paramDOFOffset + 1) += 2 * lambda[si];
  hess.coeffRef(xDOFOffset, paramDOFOffset + 2) += 0;
  hess.coeffRef(xDOFOffset, paramDOFOffset + 3) += 0;

  hess.coeffRef(paramDOFOffset, xDOFOffset) += 6 * x0 * lambda[si];
  hess.coeffRef(paramDOFOffset + 1, xDOFOffset) += 2 * lambda[si];
  hess.coeffRef(paramDOFOffset + 2, xDOFOffset) += 0;
  hess.coeffRef(paramDOFOffset + 3, xDOFOffset) += 0;

  si++;

  PGO_ALOG(hess.isCompressed() == true);

  // ==============================
  // fn-1'(xn) = yn'

  // dxn dxn
  ES::V4d paramn = getParami(x, n - 1);
  hess.coeffRef(xDOFOffset + n, xDOFOffset + n) += 6 * paramn[0] * lambda[si];

  // dx0 dparam0
  hess.coeffRef(xDOFOffset + n, paramDOFOffset + (n - 1) * 4) += 6 * xn * lambda[si];
  hess.coeffRef(xDOFOffset + n, paramDOFOffset + (n - 1) * 4 + 1) += 2 * lambda[si];
  hess.coeffRef(xDOFOffset + n, paramDOFOffset + (n - 1) * 4 + 2) += 0;
  hess.coeffRef(xDOFOffset + n, paramDOFOffset + (n - 1) * 4 + 3) += 0;

  hess.coeffRef(paramDOFOffset + (n - 1) * 4, xDOFOffset + n) += 6 * xn * lambda[si];
  hess.coeffRef(paramDOFOffset + (n - 1) * 4 + 1, xDOFOffset + n) += 2 * lambda[si];
  hess.coeffRef(paramDOFOffset + (n - 1) * 4 + 2, xDOFOffset + n) += 0;
  hess.coeffRef(paramDOFOffset + (n - 1) * 4 + 3, xDOFOffset + n) += 0;

  si++;

  PGO_ALOG(hess.isCompressed() == true);

  // fi' (x_i+1) - f'_{i+1}(x_i+1) = 0
  for (int i = 0; i < numPoints - 2; i++) {
    double xi1 = getxi(x, i + 1);
    ES::V4d parami = getParami(x, i);
    ES::V4d parami1 = getParami(x, i + 1);

    // dx_{i+1} dx_{i+1}
    hess.coeffRef(xDOFOffset + i + 1, xDOFOffset + i + 1) += (6 * parami[0] - 6 * parami1[0]) * lambda[si];

    // dx_{i+1} dparam_i
    hess.coeffRef(xDOFOffset + i + 1, paramDOFOffset + i * 4) += 6 * xi1 * lambda[si];
    hess.coeffRef(xDOFOffset + i + 1, paramDOFOffset + i * 4 + 1) += 2 * lambda[si];
    hess.coeffRef(xDOFOffset + i + 1, paramDOFOffset + i * 4 + 2) += 0;
    hess.coeffRef(xDOFOffset + i + 1, paramDOFOffset + i * 4 + 3) += 0;

    hess.coeffRef(paramDOFOffset + i * 4, xDOFOffset + i + 1) += 6 * xi1 * lambda[si];
    hess.coeffRef(paramDOFOffset + i * 4 + 1, xDOFOffset + i + 1) += 2 * lambda[si];
    hess.coeffRef(paramDOFOffset + i * 4 + 2, xDOFOffset + i + 1) += 0;
    hess.coeffRef(paramDOFOffset + i * 4 + 3, xDOFOffset + i + 1) += 0;

    // dx_{i+1} dparam_{i + 1}
    hess.coeffRef(xDOFOffset + i + 1, paramDOFOffset + i * 4 + 4) += -6 * xi1 * lambda[si];
    hess.coeffRef(xDOFOffset + i + 1, paramDOFOffset + i * 4 + 4 + 1) += -2 * lambda[si];
    hess.coeffRef(xDOFOffset + i + 1, paramDOFOffset + i * 4 + 4 + 2) += 0;
    hess.coeffRef(xDOFOffset + i + 1, paramDOFOffset + i * 4 + 4 + 3) += 0;

    hess.coeffRef(paramDOFOffset + i * 4 + 4, xDOFOffset + i + 1) += -6 * xi1 * lambda[si];
    hess.coeffRef(paramDOFOffset + i * 4 + 4 + 1, xDOFOffset + i + 1) += -2 * lambda[si];
    hess.coeffRef(paramDOFOffset + i * 4 + 4 + 2, xDOFOffset + i + 1) += 0;
    hess.coeffRef(paramDOFOffset + i * 4 + 4 + 3, xDOFOffset + i + 1) += 0;

    si++;
  }

  PGO_ALOG(hess.isCompressed() == true);

  // fi'' (x_i+1) - f''_{i+1}(x_i+1) = 0
  // fi'' = 6 a_i x + 2 b_i

  // dx = 6 a_i
  for (int i = 0; i < numPoints - 2; i++) {
    // dx_{i+1} dx_{i+1} = 0

    // dparam_i
    // dx_{i+1} dparam_i
    hess.coeffRef(xDOFOffset + i + 1, paramDOFOffset + i * 4) += 6 * lambda[si];

    hess.coeffRef(paramDOFOffset + i * 4, xDOFOffset + i + 1) += 6 * lambda[si];

    // dx_{i+1} dparam_{i + 1}
    hess.coeffRef(xDOFOffset + i + 1, paramDOFOffset + i * 4 + 4) += -6 * lambda[si];

    hess.coeffRef(paramDOFOffset + i * 4 + 4, xDOFOffset + i + 1) += -6 * lambda[si];

    si++;
  }

  PGO_ALOG(si == 4 * numPoints - 4);

  PGO_ALOG(hess.isCompressed() == true);
}