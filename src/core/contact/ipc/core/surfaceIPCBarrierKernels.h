#pragma once

#include "EigenDef.h"

namespace pgo
{
namespace Contact
{
namespace CIPC
{
namespace barrier_kernels
{

struct LocalContribution
{
  double energy;
  EigenSupport::V12d gradient;
  EigenSupport::M12d hessian;
  bool active;
};

struct DynamicPointContribution
{
  double energy = 0.0;
  EigenSupport::V3d gradient = EigenSupport::V3d::Zero();
  EigenSupport::M3d hessian = EigenSupport::M3d::Zero();
  bool active = false;
};

struct DynamicEdgeContribution
{
  double energy = 0.0;
  EigenSupport::V6d gradient = EigenSupport::V6d::Zero();
  EigenSupport::M6d hessian = EigenSupport::M6d::Zero();
  bool active = false;
};

struct DynamicTriangleContribution
{
  double energy = 0.0;
  EigenSupport::V9d gradient = EigenSupport::V9d::Zero();
  EigenSupport::M9d hessian = EigenSupport::M9d::Zero();
  bool active = false;
};

LocalContribution pointTriangle(
  const EigenSupport::V3d &p,
  const EigenSupport::V3d &t0,
  const EigenSupport::V3d &t1,
  const EigenSupport::V3d &t2,
  double weight,
  double dhat2,
  double kappa,
  bool needGradient,
  bool needHessian);

LocalContribution edgeEdge(
  const EigenSupport::V3d &ea0,
  const EigenSupport::V3d &ea1,
  const EigenSupport::V3d &eb0,
  const EigenSupport::V3d &eb1,
  double weight,
  double dhat2,
  double kappa,
  double epsEe,
  bool needGradient,
  bool needHessian);

DynamicPointContribution pointStaticTriangle(
  const EigenSupport::V3d &p,
  const EigenSupport::V3d &t0,
  const EigenSupport::V3d &t1,
  const EigenSupport::V3d &t2,
  double weight,
  double dhat2,
  double kappa,
  bool needGradient,
  bool needHessian);

DynamicTriangleContribution staticPointTriangle(
  const EigenSupport::V3d &p,
  const EigenSupport::V3d &t0,
  const EigenSupport::V3d &t1,
  const EigenSupport::V3d &t2,
  double weight,
  double dhat2,
  double kappa,
  bool needGradient,
  bool needHessian);

DynamicEdgeContribution edgeStaticEdge(
  const EigenSupport::V3d &ea0,
  const EigenSupport::V3d &ea1,
  const EigenSupport::V3d &eb0,
  const EigenSupport::V3d &eb1,
  double weight,
  double dhat2,
  double kappa,
  double epsEe,
  bool needGradient,
  bool needHessian);

}  // namespace barrier_kernels
}  // namespace CIPC
}  // namespace Contact
}  // namespace pgo
