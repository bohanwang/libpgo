/*
author: Bohan Wang
copyright to USC
*/

#include "lagrangian.h"
#include "EigenSupport.h"

#include <numeric>

using namespace pgo;
using namespace pgo::NonlinearOptimization;

namespace ES = pgo::EigenSupport;

Lagrangian::Lagrangian(PotentialEnergy_const_p eng, ConstraintFunctions_const_p cnt):
  energy(eng), constraints(cnt)
{
  grad.resize(eng->getNumDOFs());
  eng->createHessian(energyHessian);

  if (constraints) {
    std::vector<ES::TripletD> entries;

    // put energy hessian in
    for (Eigen::Index i = 0; i < energyHessian.outerSize(); i++) {
      for (ES::SpMatD::InnerIterator it(energyHessian, i); it; ++it) {
        entries.emplace_back((int)it.row(), (int)it.col(), 1.0);
      }
    }

    if (constraints->isLinear() == false) {
      // put constraints hessian in
      constraints->createHessian(constraintHessian);
      for (Eigen::Index i = 0; i < constraintHessian.outerSize(); i++) {
        for (ES::SpMatD::InnerIterator it(constraintHessian, i); it; ++it) {
          entries.emplace_back((int)it.row(), (int)it.col(), 1.0);
        }
      }
    }

    // put constraints jacobian in
    constraints->createJacobian(constraintJac);
    // [ H J' ]
    // [ J 0  ]

    for (Eigen::Index i = 0; i < constraintJac.outerSize(); i++) {
      for (ES::SpMatD::InnerIterator it(constraintJac, i); it; ++it) {
        entries.emplace_back((int)it.row() + energy->getNumDOFs(), (int)it.col(), 1.0);
        entries.emplace_back((int)it.col(), (int)it.row() + energy->getNumDOFs(), 1.0);
      }
    }

    hess.resize(energyHessian.rows() + constraintJac.rows(), energyHessian.rows() + constraintJac.rows());
    hess.setFromTriplets(entries.begin(), entries.end());

    // compute mapping matrix
    ES::small2Big(energyHessian, hess, 0, 0, eMapping);

    if (constraints->isLinear() == false) {
      ES::small2Big(constraintHessian, hess, 0, 0, cHessMapping);
    }

    constraintJacT = constraintJac.transpose();
    ES::transposeMapping(constraintJac, constraintJacT, cJacTMapping);

    ES::small2Big(constraintJacT, hess, 0, energy->getNumDOFs(), cJacMappingTR);
    ES::small2Big(constraintJac, hess, energy->getNumDOFs(), 0, cJacMappingBL);

    allDOFs.resize(energyHessian.rows() + constraintJac.rows());
    std::iota(allDOFs.begin(), allDOFs.end(), 0);

    g.resize(constraints->getNumConstraints());
  }
  else {
    hess = energyHessian;
    ES::small2Big(energyHessian, hess, 0, 0, eMapping);

    allDOFs.resize(energyHessian.rows());
    std::iota(allDOFs.begin(), allDOFs.end(), 0);
  }
}

double Lagrangian::func(ES::ConstRefVecXd x) const
{
  double energyAll = energy->func(x.head(energy->getNumDOFs()));

  if (constraints) {
    // energy
    constraints->func(x.head(energy->getNumDOFs()), g);
    // lambda^T C
    energyAll += x.tail(constraints->getNumConstraints()).dot(g);
  }

  return energyAll;
}

void Lagrangian::gradient(ES::ConstRefVecXd x, ES::RefVecXd gradOut) const
{
  memset(gradOut.data(), 0, sizeof(double) * gradOut.size());

  if (constraints) {
    // J^T lambda
    constraints->jacobian(x.head(energy->getNumDOFs()), constraintJac);
    ES::mv(constraintJac, x.tail(constraints->getNumConstraints()), gradOut.head(energy->getNumDOFs()), 1);

    constraints->func(x.head(energy->getNumDOFs()), gradOut.tail(constraints->getNumConstraints()));
  }

  // dE / dx
  energy->gradient(x.head(energy->getNumDOFs()), grad);
  gradOut.head(energy->getNumDOFs()) += grad;
}

void Lagrangian::hessian(ES::ConstRefVecXd x, ES::SpMatD &hessOut) const
{
  memset(hessOut.valuePtr(), 0, sizeof(double) * hessOut.nonZeros());

  energy->hessian(x.head(energy->getNumDOFs()), energyHessian);
  ES::addSmallToBig(1.0, energyHessian, hessOut, 0.0, eMapping, 1);

  if (constraints) {
    if (constraints->isLinear() == false) {
      constraints->hessian(x.head(energy->getNumDOFs()), x.tail(constraints->getNumConstraints()), constraintHessian);
      ES::addSmallToBig(1.0, constraintHessian, hessOut, 1.0, cHessMapping, 1);
    }

    constraints->jacobian(x.head(energy->getNumDOFs()), constraintJac);
    ES::addSmallToBig(1.0, constraintJac, hessOut, 1.0, cJacMappingBL, 1);

    ES::transposeTransfer(constraintJac, constraintJacT, cJacTMapping);
    ES::addSmallToBig(1.0, constraintJacT, hessOut, 1.0, cJacMappingTR, 1);
  }
}
