/*
author: Bohan Wang
copyright to USC, MIT
*/

#include "constraintFunctionsAssember.h"
#include "pgoLogging.h"

#include <numeric>

using namespace pgo;
using namespace pgo::NonlinearOptimization;
namespace ES = pgo::EigenSupport;

namespace pgo::NonlinearOptimization
{
struct ConstraintFunctionsAssemblerBuffer
{
  std::vector<ES::SpMatD> constraints1D_H, constraints1D_J;
  std::vector<ES::SpMatD> constraintsND_H, constraintsND_J;

  ES::SpMatD jacBuf0, jacBuf1, hess;
  ES::VXd xleft, xright;
};

}  // namespace pgo::NonlinearOptimization

ConstraintFunctionsAssembler::ConstraintFunctionsAssembler(int na):
  ConstraintFunctions(na)
{
  buf = std::make_shared<ConstraintFunctionsAssemblerBuffer>();
}

ConstraintFunctionsAssembler::~ConstraintFunctionsAssembler()
{
}

void ConstraintFunctionsAssembler::addConstraint(ConstraintFunction_p constt)
{
  constraints1D.push_back(constt);
}

void ConstraintFunctionsAssembler::addConstraint(std::shared_ptr<ConstraintFunctions> constts)
{
  constraintsND.push_back(constts);
}

void ConstraintFunctionsAssembler::init()
{
  std::vector<ES::TripletD> entriesj, entriesh;

  ES::IDX rowOffset = 0;
  for (size_t i = 0; i < constraints1D.size(); i++) {
    ES::SpMatD jac, hess;
    const auto &dofs = constraints1D[i]->getDOFs();

    constraints1D[i]->createHessian(hess);
    constraints1D[i]->createJacobian(jac);

    PGO_ALOG(jac.rows() == 1);
    for (Eigen::Index outeri = 0; outeri < jac.outerSize(); outeri++) {
      for (ES::SpMatD::InnerIterator it(jac, outeri); it; ++it) {
        entriesj.emplace_back((int)i, dofs[it.col()], 1.0);
      }
    }

    PGO_ALOG(hess.rows() == hess.cols());
    for (Eigen::Index outeri = 0; outeri < hess.outerSize(); outeri++) {
      for (ES::SpMatD::InnerIterator it(hess, outeri); it; ++it) {
        entriesh.emplace_back(dofs[it.row()], dofs[it.col()], 1.0);
      }
    }

    buf->constraints1D_H.push_back(hess);
    buf->constraints1D_J.push_back(jac);
  }

  rowOffset = (ES::IDX)constraints1D.size();

  for (size_t i = 0; i < constraintsND.size(); i++) {
    ES::SpMatD jac, hess;
    constraintsND[i]->createHessian(hess);
    constraintsND[i]->createJacobian(jac);

    for (Eigen::Index outeri = 0; outeri < jac.outerSize(); outeri++) {
      for (ES::SpMatD::InnerIterator it(jac, outeri); it; ++it) {
        entriesj.emplace_back(rowOffset + it.row(), it.col(), 1.0);
      }
    }

    for (Eigen::Index outeri = 0; outeri < hess.outerSize(); outeri++) {
      for (ES::SpMatD::InnerIterator it(hess, outeri); it; ++it) {
        entriesh.emplace_back(it.row(), it.col(), 1.0);
      }
    }

    buf->constraintsND_H.push_back(hess);
    buf->constraintsND_J.push_back(jac);

    constraintNDOffsets.push_back(rowOffset);
    rowOffset += jac.rows();
  }
  constraintNDOffsets.push_back(rowOffset);

  jacobianTemplate.resize(rowOffset, nAll);
  jacobianTemplate.setFromTriplets(entriesj.begin(), entriesj.end());

  lambdahTemplate.resize(nAll, nAll);
  lambdahTemplate.setFromTriplets(entriesh.begin(), entriesh.end());

  for (size_t i = 0; i < constraints1D.size(); i++) {
    ES::SpMatI mapping;
    const auto &dofs = constraints1D[i]->getDOFs();

    std::vector<int> dofsRow;
    dofsRow.push_back((int)i);

    ES::small2Big(buf->constraints1D_J[i], jacobianTemplate, dofsRow, dofs, mapping);
    jacobian1DMappings.push_back(mapping);

    ES::SpMatI mappingh;
    ES::small2Big(buf->constraints1D_H[i], lambdahTemplate, dofs, mappingh);
    hessian1DMappings.push_back(mappingh);
  }

  for (size_t i = 0; i < constraintsND.size(); i++) {
    ES::SpMatI mapping;

    std::vector<int> dofsRow;
    for (ES::IDX row = constraintNDOffsets[i]; row < constraintNDOffsets[i + 1]; row++)
      dofsRow.push_back((int)row);

    std::vector<int> dofsCol(nAll);
    std::iota(dofsCol.begin(), dofsCol.end(), 0);

    ES::small2Big(buf->constraintsND_J[i], jacobianTemplate, dofsRow, dofsCol, mapping);
    jacobianNDMappings.push_back(mapping);

    ES::SpMatI mappingh;
    ES::small2Big(buf->constraintsND_H[i], lambdahTemplate, dofsCol, mappingh);
    hessianNDMappings.push_back(mappingh);
  }

  buf->jacBuf0 = jacobianTemplate;
  buf->jacBuf1 = jacobianTemplate;

  buf->xleft.resize(nAll);
  buf->xright.resize(nAll);

  buf->hess = lambdahTemplate;

  hasHessianVectorRoutine = true;
  for (const auto &c : constraintsND) {
    if (c->hasHessianVector() == false) {
      hasHessianVectorRoutine = false;
    }
  }
}

void ConstraintFunctionsAssembler::func(ES::ConstRefVecXd x, ES::RefVecXd g) const
{
  for (size_t i = 0; i < constraints1D.size(); i++) {
    g[i] = constraints1D[i]->func(x);
  }

  for (size_t i = 0; i < constraintsND.size(); i++) {
    ES::IDX dim = constraintNDOffsets[i + 1] - constraintNDOffsets[i];
    ES::IDX start = constraintNDOffsets[i];
    constraintsND[i]->func(x, g.segment(start, dim));
  }
}

void ConstraintFunctionsAssembler::jacobian(ES::ConstRefVecXd x, ES::SpMatD &jac) const
{
  for (size_t i = 0; i < constraints1D.size(); i++) {
    constraints1D[i]->jacobian(x, buf->constraints1D_J[i]);
    ES::addSmallToBig(1.0, buf->constraints1D_J[i], jac, 0.0, jacobian1DMappings[i]);
  }

  for (size_t i = 0; i < constraintsND.size(); i++) {
    constraintsND[i]->jacobian(x, buf->constraintsND_J[i]);
    ES::addSmallToBig(1.0, buf->constraintsND_J[i], jac, 0.0, jacobianNDMappings[i]);
  }
}

void ConstraintFunctionsAssembler::hessian(ES::ConstRefVecXd x, ES::ConstRefVecXd lambda, ES::SpMatD &hess) const
{
  for (size_t i = 0; i < constraints1D.size(); i++) {
    if (buf->constraints1D_H[i].nonZeros() == 0)
      continue;

    constraints1D[i]->hessian(x, buf->constraints1D_H[i]);
    ES::addSmallToBig(lambda[i], buf->constraints1D_H[i], hess, 1.0, hessian1DMappings[i]);
  }

  for (size_t i = 0; i < constraintsND.size(); i++) {
    if (buf->constraintsND_H[i].nonZeros() == 0)
      continue;

    ES::IDX dim = constraintNDOffsets[i + 1] - constraintNDOffsets[i];
    ES::IDX start = constraintNDOffsets[i];

    constraintsND[i]->hessian(x, lambda.segment(start, dim), buf->constraintsND_H[i]);
    ES::addSmallToBig(1.0, buf->constraintsND_H[i], hess, 1.0, hessianNDMappings[i]);
  }
}

void ConstraintFunctionsAssembler::hessianVector(ES::ConstRefVecXd x, ES::ConstRefVecXd lambda, ES::ConstRefVecXd vec, ES::RefVecXd hessVec) const
{
  if (hasHessianVectorRoutine) {
    memset(buf->hess.valuePtr(), 0, sizeof(double) * buf->hess.nonZeros());
    for (size_t i = 0; i < constraints1D.size(); i++) {
      constraints1D[i]->hessian(x, buf->constraints1D_H[i]);
      ES::addSmallToBig(lambda[i], buf->constraints1D_H[i], buf->hess, 1.0, hessian1DMappings[i]);
    }

    ES::mv(buf->hess, vec, hessVec);

    ES::VXd temp(hessVec.size());
    for (size_t i = 0; i < constraints1D.size(); i++) {
      constraintsND[i]->hessianVector(x, lambda, vec, temp);
      hessVec += temp;
    }
  }
}