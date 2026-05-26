#include "ipc/core/surfaceIPCBarrierKernels.h"

#include "ipc/geometry/ipcBarrier.h"
#include "ipc/geometry/ipcDistancePrimitives.h"
#include "ipc/geometry/ipcHessianProjection.h"

namespace pgo
{
namespace Contact
{
namespace IPC
{
namespace barrier_kernels
{

namespace
{
LocalContribution inactiveLocalContribution()
{
  LocalContribution c;
  c.energy = 0.0;
  c.gradient.setZero();
  c.hessian.setZero();
  c.active = false;
  return c;
}
}  // namespace

LocalContribution pointTriangle(
  const EigenSupport::V3d &p,
  const EigenSupport::V3d &t0,
  const EigenSupport::V3d &t1,
  const EigenSupport::V3d &t2,
  double weight,
  double dhat2,
  double kappa,
  bool needGradient,
  bool needHessian)
{
  // --- fetch distance / gradient / hessian ---
  double d2;
  EigenSupport::V12d gd2;
  EigenSupport::M12d Hd2;

  if (needGradient && needHessian) {
    auto all = distance::computePTSqDistAll(p, t0, t1, t2);
    d2 = all.d2;
    if (d2 >= dhat2 || d2 <= 0.0)
      return inactiveLocalContribution();
    gd2 = all.grad;
    Hd2 = all.hess;
  }
  else {
    d2 = distance::computePTSqDist(p, t0, t1, t2);
    if (d2 >= dhat2 || d2 <= 0.0)
      return inactiveLocalContribution();
    if (needGradient || needHessian)
      gd2 = distance::computePTSqDistGrad(p, t0, t1, t2);
    if (needHessian)
      Hd2 = distance::computePTSqDistHess(p, t0, t1, t2);
  }

  // --- assemble ---
  LocalContribution c;
  c.active = true;
  double wk = weight * kappa;
  c.energy = wk * barrier::b(d2, dhat2);

  if (needGradient)
    c.gradient = wk * barrier::dbds(d2, dhat2) * gd2;
  else
    c.gradient.setZero();

  if (needHessian) {
    double db = barrier::dbds(d2, dhat2);
    double d2b = barrier::d2bds2(d2, dhat2);
    c.hessian = wk * (d2b * gd2 * gd2.transpose() + db * Hd2);
    c.hessian = projectToPSD(c.hessian);
  }
  else {
    c.hessian.setZero();
  }

  return c;
}

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
  bool needHessian)
{
  // --- fetch distance / gradient / hessian ---
  double d2;
  EigenSupport::V12d gd2;
  EigenSupport::M12d Hd2;

  if (needGradient && needHessian) {
    auto all = distance::computeEESqDistAll(ea0, ea1, eb0, eb1);
    d2 = all.d2;
    if (d2 >= dhat2 || d2 <= 0.0)
      return inactiveLocalContribution();
    gd2 = all.grad;
    Hd2 = all.hess;
  }
  else {
    d2 = distance::computeEESqDist(ea0, ea1, eb0, eb1);
    if (d2 >= dhat2 || d2 <= 0.0)
      return inactiveLocalContribution();
    if (needGradient || needHessian)
      gd2 = distance::computeEESqDistGrad(ea0, ea1, eb0, eb1);
    if (needHessian)
      Hd2 = distance::computeEESqDistHess(ea0, ea1, eb0, eb1);
  }

  // --- mollifier (only when epsEe > 0) ---
  double m = 1.0;
  EigenSupport::V12d gm;
  EigenSupport::M12d Hm;
  if (epsEe > 0.0) {
    m = distance::eeMollifier(ea0, ea1, eb0, eb1, epsEe);
    if (needGradient || needHessian)
      gm = distance::eeMollifierGrad(ea0, ea1, eb0, eb1, epsEe);
    if (needHessian)
      Hm = distance::eeMollifierHess(ea0, ea1, eb0, eb1, epsEe);
  }

  // --- assemble ---
  LocalContribution c;
  c.active = true;
  double wk = weight * kappa;
  double bv = barrier::b(d2, dhat2);
  c.energy = wk * m * bv;

  if (needGradient) {
    double dbv = barrier::dbds(d2, dhat2);
    if (epsEe > 0.0)
      c.gradient = wk * (gm * bv + m * dbv * gd2);
    else
      c.gradient = wk * dbv * gd2;
  }
  else {
    c.gradient.setZero();
  }

  if (needHessian) {
    double gp = barrier::dbds(d2, dhat2);
    double gpp = barrier::d2bds2(d2, dhat2);
    EigenSupport::V12d gb = gp * gd2;
    if (epsEe > 0.0)
      c.hessian = wk * (bv * Hm + gm * gb.transpose() + gb * gm.transpose()
        + m * (gpp * gd2 * gd2.transpose() + gp * Hd2));
    else
      c.hessian = wk * (gpp * gd2 * gd2.transpose() + gp * Hd2);
    c.hessian = projectToPSD(c.hessian);
  }
  else {
    c.hessian.setZero();
  }

  return c;
}

DynamicPointContribution pointStaticTriangle(
  const EigenSupport::V3d &p,
  const EigenSupport::V3d &t0,
  const EigenSupport::V3d &t1,
  const EigenSupport::V3d &t2,
  double weight,
  double dhat2,
  double kappa,
  bool needGradient,
  bool needHessian)
{
  DynamicPointContribution c;
  const LocalContribution local =
    pointTriangle(p, t0, t1, t2, weight, dhat2, kappa, needGradient, needHessian);
  c.active = local.active;
  if (!c.active)
    return c;

  c.energy = local.energy;
  if (needGradient)
    c.gradient = local.gradient.head<3>();
  if (needHessian)
    c.hessian = local.hessian.block<3, 3>(0, 0);
  return c;
}

DynamicTriangleContribution staticPointTriangle(
  const EigenSupport::V3d &p,
  const EigenSupport::V3d &t0,
  const EigenSupport::V3d &t1,
  const EigenSupport::V3d &t2,
  double weight,
  double dhat2,
  double kappa,
  bool needGradient,
  bool needHessian)
{
  DynamicTriangleContribution c;
  const LocalContribution local =
    pointTriangle(p, t0, t1, t2, weight, dhat2, kappa, needGradient, needHessian);
  c.active = local.active;
  if (!c.active)
    return c;

  c.energy = local.energy;
  if (needGradient)
    c.gradient = local.gradient.segment<9>(3);
  if (needHessian)
    c.hessian = local.hessian.block<9, 9>(3, 3);
  return c;
}

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
  bool needHessian)
{
  DynamicEdgeContribution c;
  const LocalContribution local =
    edgeEdge(ea0, ea1, eb0, eb1, weight, dhat2, kappa, epsEe, needGradient, needHessian);
  c.active = local.active;
  if (!c.active)
    return c;

  c.energy = local.energy;
  if (needGradient)
    c.gradient = local.gradient.head<6>();
  if (needHessian)
    c.hessian = local.hessian.block<6, 6>(0, 0);
  return c;
}

}  // namespace barrier_kernels
}  // namespace IPC
}  // namespace Contact
}  // namespace pgo
