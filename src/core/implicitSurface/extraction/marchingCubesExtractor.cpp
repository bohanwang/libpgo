#include "extraction/marchingCubesExtractor.h"

#include "EigenDef.h"
#include "libiglInterface.h"

#include <stdexcept>

namespace pgo::ImplicitSurface {

void extractMarchingCubes(const DenseGrid &field, const MarchingCubesOptions &options,
  Mesh::TriMeshGeo &outMesh)
{
  const EigenSupport::V3d &bmin = field.gridSpec().bmin;
  const EigenSupport::V3d &bmax = field.gridSpec().bmax;
  const int resolution = field.gridSpec().resolution;

  // Map the DenseGrid data as an Eigen vector for the libigl interface.
  Eigen::Map<const EigenSupport::VXd> fieldMap(field.data(), field.size());

  libiglInterface::computeMarchingCubes(bmin, bmax, resolution, fieldMap, outMesh, options.isoOffset);
}

}  // namespace pgo::ImplicitSurface
