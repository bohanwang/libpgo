#include "ipc/core/surfaceIPCBarrierKernels.h"

#include "ipc/geometry/ipcBarrier.h"
#include "ipc/geometry/ipcDistancePrimitives.h"
#include "ipc/geometry/ipcHessianProjection.h"

namespace pgo
{
namespace Contact
{
namespace CIPC
{
namespace barrier_kernels
{

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
  LocalContribution c;

  double d2 = distance::computePTSqDist(p, t0, t1, t2);
  if (d2 >= dhat2 || d2 <= 0.0)
    return c;

  c.active = true;
  double wk = weight * kappa;
  c.energy = wk * barrier::b(d2, dhat2);

  if (needGradient || needHessian) {
    EigenSupport::V12d gd2 = distance::computePTSqDistGrad(p, t0, t1, t2);
    double db = barrier::dbds(d2, dhat2);

    if (needGradient) {
      c.gradient = wk * db * gd2;
    }

    if (needHessian) {
      EigenSupport::M12d Hd2 = distance::computePTSqDistHess(p, t0, t1, t2);
      double d2b = barrier::d2bds2(d2, dhat2);
      c.hessian = wk * (d2b * gd2 * gd2.transpose() + db * Hd2);
      c.hessian = projectToPSD(c.hessian);
    }
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
  LocalContribution c;

  double d2 = distance::computeEESqDist(ea0, ea1, eb0, eb1);
  if (d2 >= dhat2 || d2 <= 0.0)
    return c;

  c.active = true;
  double wk = weight * kappa;

  double m = 1.0;
  EigenSupport::V12d gm = EigenSupport::V12d::Zero();
  EigenSupport::M12d Hm = EigenSupport::M12d::Zero();
  if (epsEe > 0.0) {
    m = distance::eeMollifier(ea0, ea1, eb0, eb1, epsEe);
    if (needGradient || needHessian)
      gm = distance::eeMollifierGrad(ea0, ea1, eb0, eb1, epsEe);
    if (needHessian)
      Hm = distance::eeMollifierHess(ea0, ea1, eb0, eb1, epsEe);
  }

  double bv = barrier::b(d2, dhat2);
  c.energy = wk * m * bv;

  if (needGradient || needHessian) {
    EigenSupport::V12d gd2 = distance::computeEESqDistGrad(ea0, ea1, eb0, eb1);

    if (needGradient) {
      if (epsEe > 0.0) {
        double dbv = barrier::dbds(d2, dhat2);
        c.gradient = wk * (gm * bv + m * dbv * gd2);
      } else {
        c.gradient = wk * barrier::dbds(d2, dhat2) * gd2;
      }
    }

    if (needHessian) {
      EigenSupport::M12d Hd2 = distance::computeEESqDistHess(ea0, ea1, eb0, eb1);
      double gp = barrier::dbds(d2, dhat2);
      double gpp = barrier::d2bds2(d2, dhat2);
      EigenSupport::V12d gb = barrier::dbds(d2, dhat2) * gd2;

      if (epsEe > 0.0) {
        c.hessian = wk * (bv * Hm + gm * gb.transpose() + gb * gm.transpose()
          + m * (gpp * gd2 * gd2.transpose() + gp * Hd2));
      } else {
        c.hessian = wk * (gpp * gd2 * gd2.transpose() + gp * Hd2);
      }
      c.hessian = projectToPSD(c.hessian);
    }
  }

  return c;
}

}  // namespace barrier_kernels
}  // namespace CIPC
}  // namespace Contact
}  // namespace pgo
