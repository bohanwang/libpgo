/*
author: Bohan Wang
copyright to USC
*/

#include "naturalCubicSplineFitting.h"
#include "naturalCubicSplineDerivatives.h"
#include "constraintFunctionsAssember.h"
#include "linearConstraintFunctions.h"
#include "potentialEnergy.h"
#include "knitroOptimizer.h"
#include "knitroProblem.h"

#include "pgoLogging.h"

#include <numeric>

using namespace pgo;
using namespace pgo::NonlinearOptimization;

namespace pgo::NonlinearOptimization
{
class FittingEnergy : public PotentialEnergy
{
public:
  FittingEnergy(const ES::VXd &x, const ES::VXd &y, int n):
    xVals(x), yVals(y), numControlPoints(n)
  {
    allDOFs.resize(n + n + (n - 1) * 4);
    std::iota(allDOFs.begin(), allDOFs.end(), 0);

    std::vector<ES::TripletD> entries;
    for (int i = 0; i < numControlPoints - 1; i++) {
      for (int r = 0; r < 4; r++) {
        for (int c = 0; c < 4; c++) {
          entries.emplace_back(numControlPoints * 2 + i * 4 + r, numControlPoints * 2 + i * 4 + c, 1.0);
        }
      }
    }

    hessTemplate.resize(allDOFs.size(), allDOFs.size());
    hessTemplate.setFromTriplets(entries.begin(), entries.end());
  }

  virtual double func(ES::ConstRefVecXd x) const override
  {
    double energy = 0;

    for (ES::IDX i = 0; i < xVals.size(); i++) {
      int j = 0;
      for (; j < numControlPoints - 1; j++) {
        if (x[j] <= xVals[i] && xVals[i] <= x[j + 1])
          break;
      }
      PGO_ALOG(j >= 0 && j < numControlPoints - 1);

      double a = x[numControlPoints * 2 + j * 4];
      double b = x[numControlPoints * 2 + j * 4 + 1];
      double c = x[numControlPoints * 2 + j * 4 + 2];
      double d = x[numControlPoints * 2 + j * 4 + 3];

      double y = a * xVals[i] * xVals[i] * xVals[i] + b * xVals[i] * xVals[i] + c * xVals[i] + d;
      energy += (y - yVals[i]) * (y - yVals[i]);
    }

    return energy * 0.5;
  }

  virtual void gradient(ES::ConstRefVecXd x, ES::RefVecXd grad) const override
  {
    grad.setZero();

    for (ES::IDX i = 0; i < xVals.size(); i++) {
      int j = 0;
      for (; j < numControlPoints - 1; j++) {
        if (x[j] <= xVals[i] && xVals[i] <= x[j + 1])
          break;
      }
      PGO_ALOG(j >= 0 && j < numControlPoints - 1);

      double a = x[numControlPoints * 2 + j * 4];
      double b = x[numControlPoints * 2 + j * 4 + 1];
      double c = x[numControlPoints * 2 + j * 4 + 2];
      double d = x[numControlPoints * 2 + j * 4 + 3];

      double y = a * xVals[i] * xVals[i] * xVals[i] + b * xVals[i] * xVals[i] + c * xVals[i] + d;
      double diff = y - yVals[i];

      grad[numControlPoints * 2 + j * 4] += diff * xVals[i] * xVals[i] * xVals[i];
      grad[numControlPoints * 2 + j * 4 + 1] += diff * xVals[i] * xVals[i];
      grad[numControlPoints * 2 + j * 4 + 2] += diff * xVals[i];
      grad[numControlPoints * 2 + j * 4 + 3] += diff;
    }
  }

  virtual void hessian(ES::ConstRefVecXd x, ES::SpMatD &hess) const override
  {
    memset(hess.valuePtr(), 0, sizeof(double) * hess.nonZeros());

    for (ES::IDX i = 0; i < xVals.size(); i++) {
      int j = 0;
      for (; j < numControlPoints - 1; j++) {
        if (x[j] <= xVals[i] && xVals[i] <= x[j + 1])
          break;
      }
      PGO_ALOG(j >= 0 && j < numControlPoints - 1);

      ES::V4d coeff(xVals[i] * xVals[i] * xVals[i], xVals[i] * xVals[i], xVals[i], 1.0);
      ES::M4d ccT;
      ES::tensorProduct(ccT, coeff, coeff);

      for (int row = 0; row < 4; row++) {
        for (int col = 0; col < 4; col++) {
          int gRow = numControlPoints * 2 + j * 4 + row;
          int gCol = numControlPoints * 2 + j * 4 + col;
          hess.coeffRef(gRow, gCol) += ccT(row, col);
        }
      }
    }

    PGO_ALOG(hess.isCompressed() == true);
  }

  virtual void createHessian(EigenSupport::SpMatD &hess) const override { hess = hessTemplate; }

  virtual void getDOFs(std::vector<int> &dofs) const override { dofs = allDOFs; }
  virtual int getNumDOFs() const override { return (int)allDOFs.size(); }

  virtual int isQuadratic() const override { return 0; }
  virtual int hasHessianVector() const override { return 0; }

protected:
  const ES::VXd &xVals;
  const ES::VXd &yVals;

  int numControlPoints;

  ES::SpMatD hessTemplate;
  std::vector<int> allDOFs;
};
}  // namespace pgo::NonlinearOptimization

NaturalCubicSplineFitting::NaturalCubicSplineFitting(const ES::VXd &x, const ES::VXd &y, int n, double ld, double rd, double minxd):
  xVals(x), yVals(y), numControlPoints(n), leftDerivative(ld), rightDerivative(rd), minXDist(minxd)
{
  nAll = numControlPoints + numControlPoints + (numControlPoints - 1) * 4;
  splineC = std::make_shared<NaturalCubicSpline2DAsNonlinearConstraints>(numControlPoints,
    0, numControlPoints, numControlPoints * 2, nAll, leftDerivative, rightDerivative);
}

void NaturalCubicSplineFitting::buildSpline(const double *cpX, const double *cpY)
{
  controlPointsX = ES::Mp<const ES::VXd>(cpX, numControlPoints);
  controlPointsY = ES::Mp<const ES::VXd>(cpY, numControlPoints);

  buildSpline(cpX, cpY, parameters);
}

void NaturalCubicSplineFitting::buildSpline(const double *cpX, const double *cpY, ES::VXd &param) const
{
  ES::VXd xparam(nAll);
  for (int i = 0; i < numControlPoints; i++) {
    xparam[i] = cpX[i];
    xparam[i + numControlPoints] = cpY[i];
  }

  splineC->solveParams(xparam);
  param = xparam.segment(numControlPoints * 2, (numControlPoints - 1) * 4);
}

int NaturalCubicSplineFitting::fit(const char *solverConfigFilename)
{
  std::shared_ptr<FittingEnergy> energy = std::make_shared<FittingEnergy>(xVals, yVals, numControlPoints);

  ES::SpMatD C(numControlPoints - 1, nAll);
  ES::VXd d(numControlPoints - 1);
  d.setZero();

  std::vector<ES::TripletD> entries;
  for (int i = 0; i < numControlPoints - 1; i++) {
    entries.emplace_back(i, i, -1.0);
    entries.emplace_back(i, i + 1, 1.0);
  }
  C.setFromTriplets(entries.begin(), entries.end());
  std::shared_ptr<LinearConstraintFunctions> sortedXC = std::make_shared<LinearConstraintFunctions>(C, d);

  std::shared_ptr<ConstraintFunctionsAssembler> constraints = std::make_shared<ConstraintFunctionsAssembler>(nAll);
  constraints->addConstraint(sortedXC);
  constraints->addConstraint(splineC);
  constraints->init();

  double xLeft = xVals[0];
  double xRight = xVals[xVals.size() - 1];
  double gap = (xRight - xLeft) / (numControlPoints - 1);

  ES::VXd xparam(nAll);
  for (int i = 0; i < numControlPoints; i++) {
    double xcur = xLeft + i * gap;

    int closestX = -1;
    double dist = 1e100;
    for (int j = 0; j < (int)xVals.size(); j++) {
      if (std::abs(xVals[j] - xcur) < dist) {
        dist = std::abs(xVals[j] - xcur);
        closestX = j;
      }
    }

    xparam[i] = xVals[closestX];
    xparam[i + numControlPoints] = yVals[closestX];
  }

  // xparam.segment(numControlPoints * 2, (numControlPoints - 1) * 4).setConstant(1.0);
  splineC->solveParams(xparam);

  ES::VXd xlow(nAll);
  ES::VXd xhi(nAll);

  // for x
  xlow.segment(0, numControlPoints).setConstant(xLeft);
  xhi.segment(0, numControlPoints).setConstant(xRight);

  xlow[0] = xLeft;
  xhi[0] = xLeft;

  xlow[numControlPoints - 1] = xRight;
  xhi[numControlPoints - 1] = xRight;

  // for y
  xlow.segment(numControlPoints, numControlPoints).setConstant(-1e20);
  xhi.segment(numControlPoints, numControlPoints).setConstant(1e20);

  xlow[numControlPoints] = yVals[0];
  xhi[numControlPoints] = yVals[0];

  xlow[numControlPoints + numControlPoints - 1] = yVals[yVals.size() - 1];
  xhi[numControlPoints + numControlPoints - 1] = yVals[yVals.size() - 1];

  // for param
  xlow.segment(numControlPoints * 2, (numControlPoints - 1) * 4).setConstant(-1e20);
  xhi.segment(numControlPoints * 2, (numControlPoints - 1) * 4).setConstant(1e20);

  ES::VXd clow(constraints->getNumConstraints());
  ES::VXd chi(constraints->getNumConstraints());

  // x sorted
  clow.head(C.rows()).setConstant(minXDist);
  chi.head(C.rows()).setConstant(1e20);

  // spline
  clow.tail(splineC->getNumConstraints()).setConstant(0);
  chi.tail(splineC->getNumConstraints()).setConstant(0);

  std::unique_ptr<KnitroProblem> problem = std::make_unique<KnitroProblem>(energy, constraints);
  problem->setInit(xparam);
  problem->setRange(xlow, xhi);
  problem->setConstraintsRange(clow, chi);

  KnitroOptimizer solver(problem.get());

  if (solverConfigFilename) {
    solver.setConfigFile(solverConfigFilename);
  }

  solver.setMaxIter(1000);

  solver.setFeasTol(1e-8);
  solver.setOptTol(1e-6);
  solver.setVerbose(3);
  solver.enableMultiEvaluation(0);
  solver.emphasisFeasibility(1);

  solver.init();

  int solverRet = solver.solve();

  xparam = Eigen::Map<const ES::VXd>(solver.getx(), problem->getn());

  controlPointsX = xparam.segment(0, numControlPoints);
  controlPointsY = xparam.segment(numControlPoints, numControlPoints);

  parameters = xparam.segment(numControlPoints * 2, (numControlPoints - 1) * 4);

  return solverRet;
}

double NaturalCubicSplineFitting::func(double x) const
{
  PGO_ALOG(controlPointsX.size() && controlPointsY.size() && parameters.size());

  int j = 0;
  for (; j < numControlPoints - 1; j++) {
    if (controlPointsX[j] <= x && x <= controlPointsX[j + 1])
      break;
  }
  PGO_ALOG(j >= 0 && j < numControlPoints - 1);

  double a = parameters[j * 4];
  double b = parameters[j * 4 + 1];
  double c = parameters[j * 4 + 2];
  double d = parameters[j * 4 + 3];

  return a * x * x * x + b * x * x + c * x + d;
}
