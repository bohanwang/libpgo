/*
author: Bohan Wang
copyright to USC
*/

#include "finiteDifference.h"

#include "potentialEnergy.h"
#include "constraintFunctions.h"
#include "EigenSupport.h"

#include <random>
#include <iostream>
#include <cmath>
#include <numeric>
#include <fstream>

using namespace pgo;
using namespace pgo::NonlinearOptimization;
namespace ES = pgo::EigenSupport;

const double coeffs5[5] = {
  1.0 / 12.0,
  -2.0 / 3.0,
  0,
  2.0 / 3.0,
  -1.0 / 12.0
};

const double coeffs3[3] = {
  -1.0 / 2.0,
  0,
  1.0 / 2.0
};

FiniteDifference::FiniteDifference(METHOD m, double eps_):
  method(m), eps(eps_)
{
  if (method == M_FIVE_POINT) {
    coeffs = coeffs5;
    left = -2;
    right = 2;
  }
  else if (method == M_THREE_POINT) {
    coeffs = coeffs3;
    left = -1;
    right = 1;
  }
}

void FiniteDifference::testEnergy(std::shared_ptr<const PotentialEnergy> energy, bool testGradient, bool testHessian, double range, const double *x, int dofRange,
  double *gradRelError, double *hessRelError) const
{
  int nAll = energy->getNumDOFs();
  ES::VXd xcur(nAll);

  std::vector<int> dofs;
  if (dofRange >= 0) {
    std::default_random_engine eng;
    // std::random_device rd;
    // std::mt19937 eng(rd());
    std::uniform_int_distribution<int> distrib(0, nAll - 1);

    for (int i = 0; i < dofRange; i++) {
      dofs.push_back(distrib(eng));
      // dofs.push_back(i);
    }
  }
  else {
    dofs.resize(nAll);
    std::iota(dofs.begin(), dofs.end(), 0);
  }

  if (x == nullptr) {
    std::default_random_engine eng;
    std::uniform_real_distribution<double> distrib(-range, range);
    for (int i = 0; i < nAll; i++) {
      xcur[i] = distrib(eng);
    }
  }
  else {
    xcur = Eigen::Map<const ES::VXd>(x, nAll);
  }

  if (testGradient) {
    ES::VXd xtemp = xcur, grad = ES::VXd::Zero(nAll), grad1 = ES::VXd::Zero(nAll);
    int inc = 0;
    for (int dofi : dofs) {
      xtemp = xcur;

      double sum = 0;
      for (int i = left; i <= right; i++) {
        double c = coeffs[i - left];
        xtemp[dofi] = xcur[dofi] + eps * i;

        double eng = energy->func(xtemp);
        sum += eng * c;
      }

      grad[dofi] = sum / eps;

      if (inc % 100 == 0 || inc == (int)dofs.size() - 1)
        std::cout << inc << ' ' << std::flush;
      inc++;
    }
    std::cout << std::endl;

    energy->gradient(xcur, grad1);

    double absSquaredNorm = 0, squaredNormFD = 0, squaredNormExact = 0;
    for (int dofi : dofs) {
      absSquaredNorm += (grad1[dofi] - grad[dofi]) * (grad1[dofi] - grad[dofi]);
      squaredNormExact += grad1[dofi] * grad1[dofi];
      squaredNormFD += grad[dofi] * grad[dofi];

      if (std::abs(grad[dofi] - grad1[dofi]) > 1e-6)
        std::cout << dofi << ',' << grad[dofi] << "------" << grad1[dofi] << "=====" << (grad[dofi] - grad1[dofi]) << ',' << absSquaredNorm << '\n';
    }

    double relError = sqrt(absSquaredNorm / squaredNormFD);
    std::cout << std::endl;
    std::cout << "Gradient:\n"
              << "  ||g_fd||=" << sqrt(squaredNormFD) << "\n"
              << "  ||g||=" << sqrt(squaredNormExact) << "\n"
              << "  rel err: " << relError << "\n"
              << "  abs error: " << sqrt(absSquaredNorm) << std::endl;

    if (gradRelError) {
      *gradRelError = relError;
    }
  }

  if (testHessian) {
    ES::VXd xtemp = xcur, grad(nAll), gradSum(nAll);
    ES::SpMatD hess;
    energy->createHessian(hess);

    memset(hess.valuePtr(), 0, sizeof(double) * hess.nonZeros());
    energy->hessian(xcur, hess);

    double exactNorm = 0;
    for (int dofi : dofs) {
      exactNorm += hess.col(dofi).squaredNorm();
    }
    exactNorm = sqrt(exactNorm);

    double norm = 0;
    double absErr = 0;
    int inc = 0;
    for (int dofi : dofs) {
      xtemp.noalias() = xcur;

      memset(gradSum.data(), 0, sizeof(double) * gradSum.size());
      for (int i = left; i <= right; i++) {
        double c = coeffs[i - left];
        xtemp[dofi] = xcur[dofi] + eps * i;

        memset(grad.data(), 0, sizeof(double) * grad.size());
        energy->gradient(xtemp, grad);
        gradSum += grad * c;
      }
      gradSum /= eps;

      norm += gradSum.squaredNorm();
      absErr += (hess.col(dofi) - gradSum).squaredNorm();

      if (inc % 100 == 0 || inc == (int)dofs.size() - 1)
        std::cout << inc << ' ' << std::flush;
      inc++;
    }

    std::cout << std::endl;
    std::cout << "Hessian:\n"
              << "  ||h_fd||=" << sqrt(norm) << "\n"
              << "  ||h||=" << exactNorm << "\n"
              << "  rel err: " << sqrt(absErr / norm) << "\n"
              << "  abs error: " << sqrt(absErr) << std::endl;

    if (hessRelError) {
      *hessRelError = sqrt(absErr / norm);
    }
  }
}

void FiniteDifference::testEnergy(EvalFunc evalFunc, int n, bool testGradient, bool testHessian, double range, const double *x,
  double *gradRelError, double *hessRelError) const
{
  int nAll = n;
  ES::VXd xcur(nAll);

  if (x == nullptr) {
    std::default_random_engine eng;
    std::uniform_real_distribution<double> distrib(-range, range);
    // std::cout << "x: ";
    for (int i = 0; i < nAll; i++) {
      xcur[i] = distrib(eng);
      // std::cout << xcur[i] << ' ';
    }
    // std::cout << std::endl;
  }
  else {
    xcur = Eigen::Map<const ES::VXd>(x, nAll);
  }

  if (testGradient) {
    ES::VXd xtemp(nAll), grad(nAll), grad1(nAll);
    for (int dofi = 0; dofi < nAll; dofi++) {
      xtemp = xcur;

      double sum = 0;
      for (int i = left; i <= right; i++) {
        double c = coeffs[i - left];
        xtemp[dofi] = xcur[dofi] + eps * i;

        double eng;
        evalFunc(xtemp.data(), &eng, nullptr, nullptr);
        sum += eng * c;
      }

      grad[dofi] = sum / eps;
    }

    evalFunc(xcur.data(), nullptr, grad1.data(), nullptr);

    /*std::cout << grad << "\n\n\n"
              << grad1 << std::endl;*/

    double nrm = grad.norm();
    double dif = (grad - grad1).norm();
    double err = nrm < 1e-30 ? dif : dif / nrm;

    // for (int i = 0; i < grad.size(); i++) {
    //   std::cout << grad[i] << ',' << grad1[i] << ',' << fabs(grad1[i] - grad[i]) << std::endl;
    // }

    if (err > 1e-6 && nrm > 1e-9) {
      std::cout
        << "Gradient:\n"
        << "  ||g_fd||=" << grad.norm() << "\n"
        << "  ||g||=" << grad1.norm() << "\n"
        << "  rel err: " << err << "\n"
        << "  abs error: " << dif << std::endl;
    }
    else {
      // std::cout << "Grad passed. Error: " << err << std::endl;
    }

    if (gradRelError)
      *gradRelError = err;
  }

  if (testHessian) {
    ES::VXd xtemp(nAll), grad(nAll);
    ES::MXd hess(nAll, nAll), hess1(nAll, nAll);

    for (int dofi = 0; dofi < nAll; dofi++) {
      xtemp = xcur;

      hess.col(dofi) = ES::VXd::Zero(nAll);
      for (int i = left; i <= right; i++) {
        double c = coeffs[i - left];
        xtemp[dofi] = xcur[dofi] + eps * i;

        evalFunc(xtemp.data(), nullptr, grad.data(), nullptr);
        hess.col(dofi) += grad * c;
      }
      hess.col(dofi) /= eps;
    }

    evalFunc(xcur.data(), nullptr, nullptr, hess1.data());

    double nrm = hess.norm();
    double dif = (hess - hess1).norm();
    double err = nrm < 1e-30 ? dif : dif / nrm;

    if (err > 1e-6 && nrm > 1e-9) {
      std::cout
        << "Hessian:\n"
        << "  ||h_fd||=" << hess.norm() << "\n"
        << "  ||h||=" << hess1.norm() << "\n"
        << "  rel err: " << err << "\n"
        << "  abs error: " << dif << std::endl;

      std::ofstream("temp.txt") << hess << "\n\n\n"
                                << hess1 << "\n\n\n"
                                << (hess - hess1) << std::endl;
    }
    else {
      // std::cout << "Hess passed. Error: " << err << std::endl;
    }

    if (hessRelError)
      *hessRelError = err;
  }
}

void FiniteDifference::testConstraints(std::shared_ptr<const ConstraintFunctions> cfunc, double range, const double *x, const double *lambda) const
{
  std::default_random_engine eng;

  ES::SpMatD jac;
  cfunc->createJacobian(jac);

  int nAll = (int)jac.cols();
  int cAll = (int)jac.rows();
  ES::VXd xcur(nAll);

  if (x == nullptr) {
    std::uniform_real_distribution<double> distrib(-range, range);

    for (int i = 0; i < nAll; i++) {
      xcur[i] = distrib(eng);
    }
  }
  else {
    xcur = Eigen::Map<const ES::VXd>(x, nAll);
  }

  ES::VXd lambdaCur(cAll);
  if (lambda == nullptr) {
    std::uniform_real_distribution<double> distrib(-range, range);

    for (int i = 0; i < cAll; i++) {
      lambdaCur[i] = distrib(eng);
    }
  }
  else {
    lambdaCur = Eigen::Map<const ES::VXd>(lambda, cAll);
  }

  ES::VXd xtemp(nAll), grad(cAll), gradSum(cAll);
  cfunc->jacobian(xcur, jac);

  double norm = 0;
  double absErr = 0;

  // ES::MXd jacDense0 = jac;
  // ES::MXd jacDense1 = jacDense0;

  for (int dofi = 0; dofi < nAll; dofi++) {
    xtemp = xcur;

    gradSum = ES::VXd::Zero(cAll);
    for (int i = left; i <= right; i++) {
      double c = coeffs[i - left];
      xtemp[dofi] = xcur[dofi] + eps * i;

      cfunc->func(xtemp, grad);
      gradSum += grad * c;
    }
    gradSum /= eps;

    // jacDense1.col(dofi) = gradSum;

    norm += gradSum.squaredNorm();
    absErr += (jac.col(dofi) - gradSum).squaredNorm();
  }

  std::cout << "Jacobian:\n"
            << "  ||J_fd||=" << sqrt(norm) << "\n"
            << "  ||J||=" << jac.norm() << "\n";

  if (norm < 1e-30)
    norm = 1.0;

  std::cout
    << "  rel err: " << sqrt(absErr / norm) << "\n"
    << "  abs error: " << sqrt(absErr) << std::endl;

  // std::ofstream("F:/fd.txt") << jacDense0 << "\n\n" << jacDense1 << "\n\n" << (jacDense0 - jacDense1) << std::endl;

  // hessian
  grad.resize(nAll);
  gradSum.resize(nAll);

  ES::SpMatD hess;
  cfunc->createHessian(hess);
  cfunc->hessian(xcur, lambdaCur, hess);

  norm = 0;
  absErr = 0;

  // here we assume our function is lambda^T J
  // then the derivative with respect to x is lambda^T H

  for (int dofi = 0; dofi < nAll; dofi++) {
    xtemp = xcur;

    gradSum.setZero(nAll);

    for (int i = left; i <= right; i++) {
      double c = coeffs[i - left];
      xtemp[dofi] = xcur[dofi] + eps * i;

      cfunc->jacobian(xtemp, jac);
      ES::mv(jac, lambdaCur, grad, 1);

      gradSum += grad * c;
    }
    gradSum /= eps;

    norm += gradSum.squaredNorm();
    absErr += (hess.col(dofi) - gradSum).squaredNorm();
  }

  std::cout << "Hessian:\n"
            << "  ||h_fd||=" << sqrt(norm) << "\n"
            << "  ||h||=" << hess.norm() << "\n";

  if (norm < 1e-30)
    norm = 1.0;

  std::cout
    << "  rel err: " << sqrt(absErr / norm) << "\n"
    << "  abs error: " << sqrt(absErr) << std::endl;
}

void FiniteDifference::testVecFunc(EvalVecFunc evalFunc, int m, int n, double range, const double *x, double *err) const
{
  std::default_random_engine eng;

  ES::VXd xcur(n);
  if (x == nullptr) {
    std::uniform_real_distribution<double> distrib(-range, range);

    for (int i = 0; i < n; i++) {
      xcur[i] = distrib(eng);
    }
  }
  else {
    xcur = Eigen::Map<const ES::VXd>(x, n);
  }

  ES::VXd xtemp(n), g(m);
  ES::MXd jacExact(m, n), jacFD(m, n);
  evalFunc(xcur.data(), nullptr, jacExact.data());

  for (int dofi = 0; dofi < n; dofi++) {
    xtemp.noalias() = xcur;

    jacFD.col(dofi).setZero();
    for (int i = left; i <= right; i++) {
      double c = coeffs[i - left];
      xtemp[dofi] = xcur[dofi] + eps * i;

      g.setZero();
      evalFunc(xtemp.data(), g.data(), nullptr);
      jacFD.col(dofi) += g * c;
    }
    jacFD.col(dofi) /= eps;
  }

  // std::cout << "Jacobian:\n"
  //           << "  ||J_fd||=" << jacFD.norm() << "\n"
  //           << "  ||J||=" << jacExact.norm() << std::endl;

  // std::cout << "Difference in last column:\n";
  // std::cout << jacFD.rightCols(1) - jacExact.rightCols(1) << std::endl;

  double norm = jacFD.norm();
  double absErr = (jacExact - jacFD).norm();
  if (norm < 1e-30)
    norm = 1.0;

  // std::cout
  //   << "  rel err: " << absErr / norm << "\n"
  //   << "  abs error: " << absErr << std::endl;

  if (err)
    *err = absErr / norm;
}

void FiniteDifference::randomSeq(double range, double *x, int n)
{
  std::default_random_engine eng;
  std::uniform_real_distribution<double> distrib(-range, range);

  for (int i = 0; i < n; i++) {
    x[i] = distrib(eng);
  }
}

void FiniteDifference::gradient(const double *x, double *grad, int n, std::vector<int> *dofs, std::function<double(const double *)> energyFunc)
{
  ES::VXd xtemp(n), xcur = Eigen::Map<const ES::VXd>(x, n);
  ES::Mp<ES::VXd> gradMap(grad, n);
  gradMap.setZero();

  std::vector<int> usedDOFs;
  if (dofs) {
    usedDOFs = *dofs;
    std::sort(usedDOFs.begin(), usedDOFs.end());
  }
  else {
    usedDOFs.resize(n);
    std::iota(usedDOFs.begin(), usedDOFs.end(), 0);
  }

  for (int dofi : usedDOFs) {
    xtemp.noalias() = xcur;

    double sum = 0;
    for (int i = left; i <= right; i++) {
      double c = coeffs[i - left];
      xtemp[dofi] = xcur[dofi] + eps * i;

      double eng = energyFunc(xtemp.data());
      sum += eng * c;
    }

    gradMap(dofi) = sum / eps;

    if (dofi % 100 == 0)
      std::cout << dofi << ' ' << std::flush;
  }
  std::cout << std::endl;
}

void FiniteDifference::hessian(const double *x, double *hess, int n, std::vector<int> *dofs, std::function<void(const double *, double *)> gradFunc)
{
  ES::VXd xtemp(n), xcur = Eigen::Map<const ES::VXd>(x, n);
  ES::VXd grad = ES::VXd::Zero(n), gradSum(n);
  ES::Mp<ES::MXd> hessMap(hess, n, n);
  hessMap.setZero();

  std::vector<int> usedDOFs;
  if (dofs) {
    usedDOFs = *dofs;
    std::sort(usedDOFs.begin(), usedDOFs.end());
  }
  else {
    usedDOFs.resize(n);
    std::iota(usedDOFs.begin(), usedDOFs.end(), 0);
  }

  for (int dofi : usedDOFs) {
    xtemp.noalias() = xcur;

    memset(gradSum.data(), 0, sizeof(double) * gradSum.size());
    for (int i = left; i <= right; i++) {
      double c = coeffs[i - left];
      xtemp[dofi] = xcur[dofi] + eps * i;

      memset(grad.data(), 0, sizeof(double) * grad.size());
      gradFunc(xtemp.data(), grad.data());

      gradSum += grad * c;
    }
    hessMap.col(dofi) = gradSum / eps;

    if (dofi % 100 == 0)
      std::cout << dofi << ' ' << std::flush;
  }

  for (int i = 0; i < n; i++) {
    if (std::binary_search(usedDOFs.begin(), usedDOFs.end(), i))
      continue;

    hessMap.row(i).setZero();
  }
}
