/*
author: Bohan Wang
copyright to USC
*/

#pragma once

#include "EigenDef.h"

namespace pgo
{
namespace NonlinearOptimization
{
// Optional capability for energies whose line-search evaluations need internal
// state frozen over the step segment [x, x+dx] (e.g. IPC's swept active-set
// superset). The solver detects this via dynamic_cast and only invokes the
// lifecycle when an energy advertises it; energies that do not need it (the vast
// majority) simply do not inherit this interface.
class LineSearchAwareEnergy
{
public:
  virtual ~LineSearchAwareEnergy() = default;

  // Freeze internal state valid for all evaluations along the segment [x, x+dx].
  virtual void beginLineSearch(EigenSupport::ConstRefVecXd x, EigenSupport::ConstRefVecXd dx) const = 0;
  virtual void endLineSearch() const = 0;

  // Largest line-search alpha for which the frozen state stays valid. IPC builds
  // its swept superset for alpha in [0, 1], so the default is 1.0. Line-search
  // methods that may probe beyond this must not use the frozen state.
  virtual double maxValidLineSearchAlpha() const { return 1.0; }
};
}  // namespace NonlinearOptimization
}  // namespace pgo
