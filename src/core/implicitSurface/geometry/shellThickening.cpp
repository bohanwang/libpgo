#include "geometry/shellThickening.h"

#include <stdexcept>

namespace pgo::ImplicitSurface {

void thickenMeshShell(const DenseGrid &meshDistance, double thickness, DenseGrid &out)
{
  if (&out == &meshDistance)
    throw std::runtime_error("thickenMeshShell: out must not alias meshDistance");
  if (out.gridSpec() != meshDistance.gridSpec())
    throw std::runtime_error("thickenMeshShell: grid spec mismatch");

  const double halfThickness = 0.5 * thickness;
  const int total = meshDistance.size();
  for (int i = 0; i < total; ++i)
    out[i] = meshDistance[i] - halfThickness;
}

}  // namespace pgo::ImplicitSurface
