#pragma once

#include "EigenSupport.h"

#include <cstdint>
#include <vector>

namespace pgo::Contact::CIPC
{

struct IPCExternalSurface
{
  EigenSupport::MXd vertices;
  EigenSupport::MXi triangles;
};

struct IPCCollisionSurfaceData
{
  EigenSupport::MXd vertices;
  EigenSupport::MXi triangles;
  EigenSupport::SpMatD displacementMap;
  // This is an ownership/role mask, not a property inferred from W. A
  // deformable vertex may have locally zero mapping rows; an external vertex
  // must have zero mapping rows.
  std::vector<uint8_t> vertexIsDeformable;
};

IPCCollisionSurfaceData assembleIPCCollisionSurface(
  const EigenSupport::MXd &deformableVertices,
  const EigenSupport::MXi &deformableTriangles,
  const EigenSupport::SpMatD &deformableDisplacementMap,
  const std::vector<IPCExternalSurface> &externalSurfaces);

}  // namespace pgo::Contact::CIPC
