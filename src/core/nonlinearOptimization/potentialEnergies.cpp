/*
author: Bohan Wang
copyright to USC
*/

#include "potentialEnergies.h"
#include "pgoLogging.h"
#include "EigenSupport.h"

#include <tbb/concurrent_vector.h>
#include <tbb/parallel_for.h>

#include <numeric>
#include <iostream>

using namespace pgo::NonlinearOptimization;
namespace ES = pgo::EigenSupport;

namespace pgo::NonlinearOptimization
{
class PotentialEnergiesBuffer
{
public:
  std::vector<double> energyBuffers;

  std::vector<ES::VXd> xlocals;
  std::vector<ES::VXd> vecs;

  std::vector<ES::VXd> gradients;
  std::vector<ES::VXd> hessianVectors;
  std::vector<ES::SpMatD> hessianMatrices;
};

PotentialEnergies::PotentialEnergies(int n):
  nAll(n)
{
  buffer = std::make_shared<PotentialEnergiesBuffer>();
}

PotentialEnergies::PotentialEnergies(int n, std::shared_ptr<PotentialEnergiesBuffer> buf):
  nAll(n), buffer(buf)
{
}

PotentialEnergies::~PotentialEnergies()
{
}
}  // namespace pgo::NonlinearOptimization

void PotentialEnergies::init()
{
  // std::vector<ES::TripletD> entries;
  tbb::concurrent_vector<ES::TripletD> entries;
  for (auto energy : potentialEnergies) {
    ES::SpMatD h;
    energy->createHessian(h);

    std::vector<int> dofs;
    energy->getDOFs(dofs);

    /*
    for (Eigen::Index outeri = 0; outeri < h.outerSize(); outeri++) {
      for (ES::SpMatD::InnerIterator it(h, outeri); it; ++it) {
        entries.emplace_back(
          (ES::SpMatD::StorageIndex)dofs[it.row()],
          (ES::SpMatD::StorageIndex)dofs[it.col()],
          1.0);
      }
    }
    */

    tbb::parallel_for((ES::IDX)0, h.outerSize(), [&](ES::IDX outeri) {
      for (ES::SpMatD::InnerIterator it(h, outeri); it; ++it) {
        entries.emplace_back(
          (ES::SpMatD::StorageIndex)dofs[it.row()],
          (ES::SpMatD::StorageIndex)dofs[it.col()],
          1.0);
      }
    });

    buffer->hessianMatrices.push_back(h);
  }

  hessianAll.resize(nAll, nAll);
  hessianAll.setFromTriplets(entries.begin(), entries.end());

  for (size_t i = 0; i < potentialEnergies.size(); i++) {
    std::vector<int> dofs;
    potentialEnergies[i]->getDOFs(dofs);

    ES::SpMatI mapping;
    if (buffer->hessianMatrices[i].nonZeros()) {
      ES::small2Big(buffer->hessianMatrices[i], hessianAll, dofs, mapping);
    }

    hessianMatrixMappings.push_back(mapping);

    buffer->xlocals.push_back(ES::VXd::Zero(dofs.size()));
    buffer->vecs.push_back(ES::VXd::Zero(dofs.size()));
    buffer->gradients.push_back(ES::VXd::Zero(dofs.size()));
    buffer->hessianVectors.push_back(ES::VXd::Zero(dofs.size()));
    energyDOFs.emplace_back(std::move(dofs));
  }

  allDOFs.resize(nAll);
  std::iota(allDOFs.begin(), allDOFs.end(), 0);

  isQuadraticEnergy = 0;
  for (size_t i = 0; i < potentialEnergies.size(); i++) {
    if (potentialEnergies[i]->isQuadratic() == 0) {
      isQuadraticEnergy = 0;
      break;
    }
  }

  hasHessianVectorProduct = 1;
  for (size_t i = 0; i < potentialEnergies.size(); i++) {
    if (potentialEnergies[i]->hasHessianVector() == 0) {
      hasHessianVectorProduct = 0;
      break;
    }
  }
}

void PotentialEnergies::setDOFs(const std::vector<int> &dofs)
{
  PGO_ALOG(dofs.size() == allDOFs.size());
  allDOFs = dofs;
}

void PotentialEnergies::mapx(ES::ConstRefVecXd x, const std::vector<int> &dofs, ES::RefVecXd xlocal) const
{
  for (int i = 0; i < (int)dofs.size(); i++) {
    xlocal(i) = x(dofs[i]);
  }
}

double PotentialEnergies::func(ES::ConstRefVecXd x) const
{
  double energyAll = 0;
  for (size_t i = 0; i < potentialEnergies.size(); i++) {
    if (energyCoeffs[i] == 0) {
      continue;
    }

    mapx(x, energyDOFs[i], buffer->xlocals[i]);
    double eng = potentialEnergies[i]->func(buffer->xlocals[i]) * energyCoeffs[i];
    energyAll += eng;
    // std::cout << eng << ',';
  }

  // std::cout << std::endl;
  return energyAll;
}

void PotentialEnergies::gradient(ES::ConstRefVecXd x, ES::RefVecXd grad) const
{
  grad.setZero();

  for (size_t i = 0; i < potentialEnergies.size(); i++) {
    if (energyCoeffs[i] == 0) {
      continue;
    }

    mapx(x, energyDOFs[i], buffer->xlocals[i]);
    potentialEnergies[i]->gradient(buffer->xlocals[i], buffer->gradients[i]);

    for (Eigen::Index j = 0; j < buffer->gradients[i].size(); j++)
      grad[energyDOFs[i][j]] += buffer->gradients[i][j] * energyCoeffs[i];
  }
}

void PotentialEnergies::hessian(ES::ConstRefVecXd x, ES::SpMatD &hess) const
{
  std::memset(hess.valuePtr(), 0, sizeof(double) * hess.nonZeros());
  for (size_t i = 0; i < potentialEnergies.size(); i++) {
    if (energyCoeffs[i] == 0) {
      continue;
    }

    if (buffer->hessianMatrices[i].nonZeros() == 0)
      continue;

    mapx(x, energyDOFs[i], buffer->xlocals[i]);
    potentialEnergies[i]->hessian(buffer->xlocals[i], buffer->hessianMatrices[i]);

    ES::addSmallToBig(energyCoeffs[i], buffer->hessianMatrices[i], hess, 1.0, hessianMatrixMappings[i]);
  }
}

void PotentialEnergies::printEnergy(EigenSupport::ConstRefVecXd x) const
{
  for (size_t i = 0; i < potentialEnergies.size(); i++) {
    mapx(x, energyDOFs[i], buffer->xlocals[i]);
    double energy = potentialEnergies[i]->func(buffer->xlocals[i]);
    std::cout << "Energy " << i << ": " << energy << ',' << energy * energyCoeffs[i] << std::endl;
  }
}

void PotentialEnergies::hessianVector(EigenSupport::ConstRefVecXd x, EigenSupport::ConstRefVecXd vec, EigenSupport::RefVecXd hessVec) const
{
  hessVec.setZero();

  if (hasHessianVectorProduct) {
    for (size_t i = 0; i < potentialEnergies.size(); i++) {
      if (energyCoeffs[i] == 0) {
        continue;
      }

      mapx(x, energyDOFs[i], buffer->xlocals[i]);
      mapx(vec, energyDOFs[i], buffer->vecs[i]);
      potentialEnergies[i]->hessianVector(buffer->xlocals[i], buffer->vecs[i], buffer->hessianVectors[i]);

      for (Eigen::Index j = 0; j < buffer->hessianVectors[i].size(); j++)
        hessVec[energyDOFs[i][j]] += buffer->hessianVectors[i][j] * energyCoeffs[i];
    }
  }
}