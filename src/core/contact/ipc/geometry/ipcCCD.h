/*
copyright to Bohan Wang
*/

#pragma once

#include "EigenDef.h"

namespace pgo
{
namespace Contact
{
namespace IPC
{

namespace ccd
{

double pointTriangleCCD(const EigenSupport::V3d &p, const EigenSupport::V3d &t0,
  const EigenSupport::V3d &t1, const EigenSupport::V3d &t2,
  const EigenSupport::V3d &dp, const EigenSupport::V3d &dt0,
  const EigenSupport::V3d &dt1, const EigenSupport::V3d &dt2,
  double thickness = 0.0,
  double tMax = 1.0);

double edgeEdgeCCD(const EigenSupport::V3d &ea0, const EigenSupport::V3d &ea1,
  const EigenSupport::V3d &eb0, const EigenSupport::V3d &eb1,
  const EigenSupport::V3d &dea0, const EigenSupport::V3d &dea1,
  const EigenSupport::V3d &deb0, const EigenSupport::V3d &deb1,
  double thickness = 0.0,
  double tMax = 1.0);

}  // namespace ccd

}  // namespace IPC
}  // namespace Contact
}  // namespace pgo
