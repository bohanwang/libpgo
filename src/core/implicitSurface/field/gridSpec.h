#pragma once

#include "EigenSupport.h"

namespace pgo::ImplicitSurface {

struct GridSpec {
  EigenSupport::V3d bmin;
  EigenSupport::V3d bmax;
  int resolution = 0;

  bool operator==(const GridSpec &other) const;
};

void validateGridSpec(const GridSpec &spec);
int linearIndex(int x, int y, int z, int resolution);

}  // namespace pgo::ImplicitSurface
