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
  double energy = 0.0;
  EigenSupport::V12d gradient = EigenSupport::V12d::Zero();
  EigenSupport::M12d hessian = EigenSupport::M12d::Zero();
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

}  // namespace barrier_kernels
}  // namespace CIPC
}  // namespace Contact
}  // namespace pgo
