/*
copyright to Bohan Wang
*/

#pragma once

#include "EigenDef.h"
#include "ipc/core/surfaceIPCPairs.h"

#include <array>
#include <vector>

namespace pgo::Contact::CIPC
{

// Lagged, regularized Coulomb friction (IPC technical supplement, sections 8-9).
// update() freezes closest points, tangent spaces, and barrier normal forces.
// The reference positions remain x_n throughout all lagging iterations of a step.
class SurfaceIPCFriction
{
public:
  void clear();
  void update(EigenSupport::ConstRefVecXd referencePositions,
    EigenSupport::ConstRefVecXd laggedPositions,
    const std::vector<PTPair> &ptPairs, const std::vector<EEPair> &eePairs,
    double dhat, double kappa, double epsEE, double frictionCoeff,
    double frictionEpsV, double timestep);

  double computeEnergy(EigenSupport::ConstRefVecXd positions) const;
  void addGradient(EigenSupport::ConstRefVecXd positions, EigenSupport::RefVecXd gradient) const;
  void computeHessian(EigenSupport::ConstRefVecXd positions, EigenSupport::SpMatD &hessian) const;
  std::size_t numPairs() const { return contacts_.size(); }

private:
  struct Contact
  {
    std::array<int, 4> vertices;
    std::array<double, 4> coefficients;
    EigenSupport::M3d tangentProjector;
    double frictionNormalForce;
  };

  void validatePositions(EigenSupport::ConstRefVecXd positions) const;
  EigenSupport::V3d slip(const Contact &contact, EigenSupport::ConstRefVecXd positions) const;

  EigenSupport::VXd referencePositions_;
  std::vector<Contact> contacts_;
  double smoothingLength_ = 1.0;
};

}  // namespace pgo::Contact::CIPC
