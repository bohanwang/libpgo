#pragma once

#include "field/denseGrid.h"

namespace pgo::ImplicitSurface {

// Offset a mesh unsigned-distance field by -0.5*thickness to form a shell:
//   out[i] = meshDistance[i] - 0.5 * thickness
// Negative values are inside the shell, positive outside.
void thickenMeshShell(const DenseGrid &meshDistance, double thickness, DenseGrid &out);

}  // namespace pgo::ImplicitSurface
