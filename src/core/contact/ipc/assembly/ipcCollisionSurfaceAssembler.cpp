#include "ipcCollisionSurfaceAssembler.h"

#include <algorithm>
#include <limits>
#include <stdexcept>

namespace pgo::Contact::CIPC
{

IPCCollisionSurfaceData assembleIPCCollisionSurface(
  const EigenSupport::MXd &deformableVertices,
  const EigenSupport::MXi &deformableTriangles,
  const EigenSupport::SpMatD &deformableDisplacementMap,
  const std::vector<IPCExternalSurface> &externalSurfaces)
{
  namespace ES = pgo::EigenSupport;

  if (deformableVertices.cols() != 3 || deformableTriangles.cols() != 3)
    throw std::invalid_argument("Collision-surface assembler requires triangle mesh matrices with three columns.");
  if (deformableDisplacementMap.rows() != deformableVertices.rows() * 3)
    throw std::invalid_argument("Deformable collision displacement-map row count is invalid.");

  std::int64_t totalVertices = deformableVertices.rows();
  std::int64_t totalTriangles = deformableTriangles.rows();
  for (const IPCExternalSurface &surface : externalSurfaces) {
    if (surface.vertices.cols() != 3 || surface.triangles.cols() != 3)
      throw std::invalid_argument("External collision surfaces must use V/F matrices with three columns.");
    if (!surface.vertices.allFinite())
      throw std::invalid_argument("External collision surface contains non-finite vertices.");
    if (surface.triangles.size() > 0 &&
      (surface.triangles.minCoeff() < 0 || surface.triangles.maxCoeff() >= surface.vertices.rows()))
      throw std::invalid_argument("External collision surface contains an out-of-range face index.");
    totalVertices += surface.vertices.rows();
    totalTriangles += surface.triangles.rows();
  }
  if (totalVertices > std::numeric_limits<int>::max() || totalTriangles > std::numeric_limits<int>::max())
    throw std::invalid_argument("Combined collision surface exceeds supported 32-bit index counts.");

  IPCCollisionSurfaceData data;
  data.vertices.resize(static_cast<Eigen::Index>(totalVertices), 3);
  data.triangles.resize(static_cast<Eigen::Index>(totalTriangles), 3);
  data.vertices.topRows(deformableVertices.rows()) = deformableVertices;
  data.triangles.topRows(deformableTriangles.rows()) = deformableTriangles;
  data.vertexIsDeformable.assign(static_cast<std::size_t>(totalVertices), uint8_t{ 0 });
  std::fill_n(
    data.vertexIsDeformable.begin(), static_cast<std::size_t>(deformableVertices.rows()), uint8_t{ 1 });

  Eigen::Index vertexOffset = deformableVertices.rows();
  Eigen::Index triangleOffset = deformableTriangles.rows();
  for (const IPCExternalSurface &surface : externalSurfaces) {
    data.vertices.middleRows(vertexOffset, surface.vertices.rows()) = surface.vertices;
    data.triangles.middleRows(triangleOffset, surface.triangles.rows()) =
      surface.triangles.array() + static_cast<int>(vertexOffset);
    vertexOffset += surface.vertices.rows();
    triangleOffset += surface.triangles.rows();
  }

  std::vector<ES::TripletD> triplets;
  triplets.reserve(static_cast<std::size_t>(deformableDisplacementMap.nonZeros()));
  for (Eigen::Index outer = 0; outer < deformableDisplacementMap.outerSize(); ++outer) {
    for (ES::SpMatD::InnerIterator it(deformableDisplacementMap, outer); it; ++it)
      triplets.emplace_back(it.row(), it.col(), it.value());
  }
  data.displacementMap.resize(
    static_cast<Eigen::Index>(totalVertices) * 3, deformableDisplacementMap.cols());
  data.displacementMap.setFromTriplets(triplets.begin(), triplets.end());

  if (!data.vertices.allFinite())
    throw std::invalid_argument("Combined collision surface contains non-finite vertices.");
  if (data.triangles.size() > 0 &&
    (data.triangles.minCoeff() < 0 || data.triangles.maxCoeff() >= data.vertices.rows()))
    throw std::invalid_argument("Combined collision surface contains an out-of-range face index.");
  if (data.displacementMap.rows() != data.vertices.rows() * 3 ||
    data.displacementMap.cols() != deformableDisplacementMap.cols())
    throw std::invalid_argument("Combined collision displacement map has invalid dimensions.");
  return data;
}

}  // namespace pgo::Contact::CIPC
