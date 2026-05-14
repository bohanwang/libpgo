#include "field/denseGrid.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace pgo::ImplicitSurface {

void validateGridSpec(const GridSpec &spec)
{
  if (spec.resolution < 2)
    throw std::runtime_error("GridSpec resolution must be at least 2");

  for (int axis = 0; axis < 3; ++axis) {
    if (!std::isfinite(spec.bmin[axis]) || !std::isfinite(spec.bmax[axis]))
      throw std::runtime_error("GridSpec bounds must be finite");
    if (spec.bmax[axis] <= spec.bmin[axis])
      throw std::runtime_error("GridSpec: bmax[" + std::to_string(axis) + "] must be greater than bmin[" + std::to_string(axis) + "]");
  }
}

bool GridSpec::operator==(const GridSpec &other) const
{
  return resolution == other.resolution &&
         bmin == other.bmin &&
         bmax == other.bmax;
}

int linearIndex(int x, int y, int z, int resolution)
{
  return z * resolution * resolution + y * resolution + x;
}

DenseGrid::DenseGrid(const GridSpec &spec)
  : spec_(spec)
{
  validateGridSpec(spec_);

  const std::size_t total = static_cast<std::size_t>(spec_.resolution) *
                            static_cast<std::size_t>(spec_.resolution) *
                            static_cast<std::size_t>(spec_.resolution);
  values_.resize(total, 0.0);
}

double &DenseGrid::at(int x, int y, int z)
{
  return values_[linearIndex(x, y, z, spec_.resolution)];
}

const double &DenseGrid::at(int x, int y, int z) const
{
  return values_[linearIndex(x, y, z, spec_.resolution)];
}

void DenseGrid::fill(double value)
{
  std::fill(values_.begin(), values_.end(), value);
}

void DenseGrid::setZero()
{
  fill(0.0);
}

}  // namespace pgo::ImplicitSurface
