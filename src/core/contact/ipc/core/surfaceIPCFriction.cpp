/*
copyright to Bohan Wang
*/

#include "surfaceIPCFriction.h"
#include "ipc/geometry/ipcBarrier.h"
#include "ipc/geometry/ipcDistancePrimitives.h"

#include <cmath>
#include <stdexcept>

namespace pgo::Contact::CIPC
{
using namespace EigenSupport;

void SurfaceIPCFriction::clear()
{
  contacts_.clear();
  referencePositions_.resize(0);
}

void SurfaceIPCFriction::update(ConstRefVecXd reference, ConstRefVecXd lagged,
  const std::vector<PTPair> &ptPairs, const std::vector<EEPair> &eePairs,
  double dhat, double kappa, double epsEE, double mu, double epsV, double timestep)
{
  if (reference.size() != lagged.size() || reference.size() % 3 != 0 ||
    !reference.allFinite() || !lagged.allFinite())
    throw std::invalid_argument("IPC friction requires matching finite surface position vectors.");
  if (!std::isfinite(dhat) || dhat <= 0 || !std::isfinite(kappa) || kappa <= 0 ||
    !std::isfinite(epsEE) || epsEE < 0 || !std::isfinite(mu) || mu < 0 ||
    !std::isfinite(epsV) || epsV <= 0 || !std::isfinite(timestep) || timestep <= 0 ||
    !std::isfinite(epsV * timestep) || epsV * timestep <= 0)
    throw std::invalid_argument("Invalid IPC friction parameters (mu >= 0, epsv and timestep > 0 required).");

  clear();
  referencePositions_ = reference;
  smoothingLength_ = epsV * timestep;
  if (mu == 0)
    return;

  auto vertex = [&](int i) -> V3d {
    if (i < 0 || i >= lagged.size() / 3)
      throw std::invalid_argument("IPC friction pair has an out-of-range vertex index.");
    return lagged.segment<3>(3 * i);
  };
  auto addContact = [&](const std::array<int, 4> &vertices, double d2,
                      const V12d &distanceGradient, const V3d &separation, double weight) {
    // The envelope theorem gives grad_i(d^2) = 2 * c_i * separation,
    // including the point/edge boundary cases selected by the distance kernel.
    // This recovers the signed closest-point barycentric coefficients without
    // a second (potentially inconsistent) closest-feature classification.
    const double separationSquared = separation.squaredNorm();
    if (!distanceGradient.allFinite() || !std::isfinite(separationSquared) || separationSquared <= 0)
      throw std::runtime_error("IPC friction encountered an invalid closest-point separation.");
    Contact contact;
    contact.vertices = vertices;
    for (int i = 0; i < 4; ++i)
      contact.coefficients[i] = distanceGradient.segment<3>(3 * i).dot(separation) / (2 * separationSquared);
    const V3d normal = separation / std::sqrt(separationSquared);
    contact.tangentProjector = M3d::Identity() - normal * normal.transpose();
    // Differentiate this repository's normalized barrier with respect to
    // unsigned distance, including its quadrature weight and EE mollifier.
    contact.frictionNormalForce = mu * (-2 * std::sqrt(d2) * weight * kappa * barrier::dbds(d2, dhat * dhat));
    if (!std::isfinite(contact.frictionNormalForce) || contact.frictionNormalForce < 0)
      throw std::runtime_error("IPC friction encountered an invalid barrier normal force.");
    if (contact.frictionNormalForce > 0)
      contacts_.push_back(contact);
  };

  for (const PTPair &pair : ptPairs) {
    const V3d p = vertex(pair.p), t0 = vertex(pair.t0), t1 = vertex(pair.t1), t2 = vertex(pair.t2);
    const double d2 = distance::computePTSqDist(p, t0, t1, t2);
    if (!std::isfinite(d2) || d2 <= 0)
      throw std::runtime_error("IPC friction requires strictly separated lagged contacts.");
    if (d2 >= dhat * dhat)
      continue;
    const V12d g = distance::computePTSqDistGrad(p, t0, t1, t2);
    addContact({ pair.p, pair.t0, pair.t1, pair.t2 }, d2, g, g.head<3>() / 2, pair.weight);
  }
  for (const EEPair &pair : eePairs) {
    const V3d a0 = vertex(pair.ea0), a1 = vertex(pair.ea1), b0 = vertex(pair.eb0), b1 = vertex(pair.eb1);
    const double d2 = distance::computeEESqDist(a0, a1, b0, b1);
    if (!std::isfinite(d2) || d2 <= 0)
      throw std::runtime_error("IPC friction requires strictly separated lagged contacts.");
    if (d2 >= dhat * dhat)
      continue;
    const double mollifier = epsEE > 0 ? distance::eeMollifier(a0, a1, b0, b1, epsEE) : 1.0;
    const V12d g = distance::computeEESqDistGrad(a0, a1, b0, b1);
    addContact({ pair.ea0, pair.ea1, pair.eb0, pair.eb1 }, d2, g,
      (g.segment<3>(0) + g.segment<3>(3)) / 2, pair.weight * mollifier);
  }
}

void SurfaceIPCFriction::validatePositions(ConstRefVecXd positions) const
{
  if (positions.size() != referencePositions_.size() || !positions.allFinite())
    throw std::invalid_argument("IPC friction positions must match the finite reference state.");
}

V3d SurfaceIPCFriction::slip(const Contact &contact, ConstRefVecXd positions) const
{
  V3d relative = V3d::Zero();
  for (int i = 0; i < 4; ++i) {
    const int offset = 3 * contact.vertices[i];
    relative += contact.coefficients[i] * (positions.segment<3>(offset) - referencePositions_.segment<3>(offset));
  }
  return contact.tangentProjector * relative;
}

double SurfaceIPCFriction::computeEnergy(ConstRefVecXd positions) const
{
  validatePositions(positions);
  double energy = 0;
  for (const Contact &contact : contacts_) {
    const double length = slip(contact, positions).norm();
    const double t = length / smoothingLength_;
    // C2 potential with a C1 force, identical to the IPC f0 regularization.
    const double f0 = length < smoothingLength_ ?
      smoothingLength_ * (t * t * (1 - t / 3) + 1.0 / 3) :
      length;
    energy += contact.frictionNormalForce * f0;
  }
  return energy;
}

void SurfaceIPCFriction::addGradient(ConstRefVecXd positions, RefVecXd gradient) const
{
  validatePositions(positions);
  if (gradient.size() != positions.size())
    throw std::invalid_argument("IPC friction gradient size must match the surface state.");
  for (const Contact &contact : contacts_) {
    const V3d r = slip(contact, positions);
    const double length = r.norm();
    const double f1OverLength = length < smoothingLength_ ?
      (2 - length / smoothingLength_) / smoothingLength_ :
      1 / length;
    const V3d force = contact.frictionNormalForce * f1OverLength * (contact.tangentProjector * r);
    for (int i = 0; i < 4; ++i)
      gradient.segment<3>(3 * contact.vertices[i]) += contact.coefficients[i] * force;
  }
}

void SurfaceIPCFriction::computeHessian(ConstRefVecXd positions, SpMatD &hessian) const
{
  validatePositions(positions);
  std::vector<TripletD> triplets;
  triplets.reserve(contacts_.size() * 144);
  for (const Contact &contact : contacts_) {
    const V3d r = slip(contact, positions);
    const double length = r.norm();
    M3d H;
    if (length == 0) {
      // Analytic zero-slip limit; never divide by zero or discard its stiffness.
      H = (2 / smoothingLength_) * M3d::Identity();
    }
    else {
      const V3d direction = r / length;
      const M3d outer = direction * direction.transpose();
      if (length < smoothingLength_) {
        const double t = length / smoothingLength_;
        H = ((2 - t) * M3d::Identity() - t * outer) / smoothingLength_;
      }
      else {
        H = (M3d::Identity() - outer) / length;
      }
    }
    const M3d block = contact.frictionNormalForce * contact.tangentProjector * H * contact.tangentProjector;
    for (int i = 0; i < 4; ++i) {
      for (int j = 0; j < 4; ++j) {
        const double weight = contact.coefficients[i] * contact.coefficients[j];
        if (weight == 0)
          continue;
        for (int r = 0; r < 3; ++r)
          for (int c = 0; c < 3; ++c)
            triplets.emplace_back(3 * contact.vertices[i] + r, 3 * contact.vertices[j] + c, weight * block(r, c));
      }
    }
  }
  hessian.resize(positions.size(), positions.size());
  hessian.setFromTriplets(triplets.begin(), triplets.end());
  hessian.makeCompressed();
}

}  // namespace pgo::Contact::CIPC
