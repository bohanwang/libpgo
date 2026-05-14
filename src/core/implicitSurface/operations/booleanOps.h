#pragma once

#include "field/denseGrid.h"

namespace pgo::ImplicitSurface {

enum class BooleanOp { Union, Intersection, Difference };

// Apply boolean op: Union=min(a,b), Intersection=max(a,b), Difference=max(a,-b).
// All grids must have matching GridSpec. out must not alias a or b.
void applyBoolean(const DenseGrid &a, const DenseGrid &b, BooleanOp op, DenseGrid &out);

}  // namespace pgo::ImplicitSurface
