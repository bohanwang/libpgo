#pragma once

#include "EigenSupport.h"
#include "triMeshGeo.h"

namespace pgo::Tools::FBMSUnionSurface
{

struct ThicknessSearchResult
{
  double thickness = 0.0;
  double volume = 0.0;
  Mesh::TriMeshGeo mesh;
  bool budgetExceededAtMin = false;
  bool budgetUnreachedAtMax = false;
  int iterations = 0;
};

}  // namespace pgo::Tools::FBMSUnionSurface
