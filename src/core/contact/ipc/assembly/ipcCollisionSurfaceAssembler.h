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

// Where one external surface landed inside the combined collision surface.
struct IPCExternalSurfaceRange
{
  int firstVertex = 0;
  int numVertices = 0;
  int firstTriangle = 0;
  int numTriangles = 0;
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
  // One entry per external surface, in input order; deformable vertices and
  // triangles come first, so every range starts after them.
  std::vector<IPCExternalSurfaceRange> externalRanges;
};

IPCCollisionSurfaceData assembleIPCCollisionSurface(
  const EigenSupport::MXd &deformableVertices,
  const EigenSupport::MXi &deformableTriangles,
  const EigenSupport::SpMatD &deformableDisplacementMap,
  const std::vector<IPCExternalSurface> &externalSurfaces);

}  // namespace pgo::Contact::CIPC
