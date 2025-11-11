/*
author: Bohan Wang
copyright to USC,MIT,NUS
*/

#pragma once
namespace pgo
{
namespace SolidDeformationModel
{
class InvariantBasedMaterial
{
public:
  InvariantBasedMaterial() {}
  virtual ~InvariantBasedMaterial() {}

  virtual double compute_psi(const double invariants[3]) const = 0;
  virtual void compute_dpsi_dI(const double invariants[3], double gradient[3]) const = 0;
  // invariants is a 3-vector, hessian is a 3x3 symmetric matrix, unrolled into a 6-vector, 
  // in the following order: (11, 12, 13, 22, 23, 33).
  virtual void compute_d2psi_dI2(const double invariants[3], double hessian[6]) const = 0;
};

}  // namespace SolidDeformationModel
}  // namespace pgo