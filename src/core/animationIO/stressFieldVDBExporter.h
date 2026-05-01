#pragma once

#include "EigenDef.h"
#include "tetMesh.h"

#include <memory>
#include <string>
#include <vector>

namespace pgo
{
namespace AnimationIO
{
// Builds a per-frame OpenVDB scalar field from a tet-mesh simulation.
//
// Inputs (all written by tools/runSim/runIPCSim):
//   - .veg volumetric mesh (rest configuration)
//   - per-frame deformation files: deform{frame:04d}.u
//       each .u stores a (3*numVertices) x 3 matrix; columns are [u, uvel, uacc]
//   - per-frame von Mises stress JSON: von_mises{frame:04d}.json
//       fields: "frame", "values" (numElements scalars, one per tet)
//
// Output: a sequence of .vdb files, one per frame, each containing a single
// FloatGrid named "von_mises". Per-element stress values are splatted into the
// 3D field at the deformed configuration.
class StressFieldVDBExporter
{
public:
  StressFieldVDBExporter() = default;
  ~StressFieldVDBExporter() = default;

  int loadTetMesh(const char *vegFilename);

  // Loads the deformation sequence written by runIPCSim.
  // `folderPath` is the "states" subfolder; `pattern` is a fmt-style pattern
  // for the per-frame file name (e.g. "deform{:04d}.u"). Frames in
  // [frameStart, frameEnd) are loaded.
  int loadDeformationSequence(const char *folderPath, const char *pattern,
    int frameStart, int frameEnd);

  // Loads the von Mises stress sequence written by runIPCSim.
  // `folderPath` is the "stress" subfolder; `pattern` is a fmt-style pattern
  // for the per-frame file name (e.g. "von_mises{:04d}.json").
  int loadVonMisesSequence(const char *folderPath, const char *pattern,
    int frameStart, int frameEnd);

  // Writes one .vdb file per frame to `outputDirectory`, using
  // "{prefix}{frame:04d}.vdb" as the file name. `voxelSize` is in world units;
  // when `voxelSize <= 0` it is auto-derived from the rest mesh's average tet
  // edge length / 2. Each tet element splats its stress into voxels whose
  // centers lie inside the deformed tet.
  int exportAnimationVDB(const char *outputDirectory, const char *prefix,
    double voxelSize) const;

  int numFrames() const { return static_cast<int>(displacements.size()); }

private:
  std::shared_ptr<VolumetricMeshes::TetMesh> tetMesh;
  std::vector<EigenSupport::VXd> displacements;
  std::vector<EigenSupport::VXd> stresses;
  int frameStart = 0;
};
}  // namespace AnimationIO
}  // namespace pgo
