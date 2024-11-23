/*
author: Bohan Wang
copyright to USC
*/

#pragma once

#include "EigenSupport.h"

#include <iostream>
#include <iomanip>

namespace pgo
{
namespace ES = EigenSupport;

namespace NonlinearOptimization
{
namespace SVDDerivatives
{
void unorderedSquareMatrixSVD3(const EigenSupport::M3d &F, EigenSupport::V3d &S, EigenSupport::M3d *U = nullptr, EigenSupport::M3d *V = nullptr);
void unorderedSquareMatrixSVD3Derivatices(const EigenSupport::M3d &F, const EigenSupport::V3d &S, const EigenSupport::M3d &U, const EigenSupport::M3d &V,
  const EigenSupport::M3d &dF, EigenSupport::V3d &dS, EigenSupport::M3d *dU = nullptr, EigenSupport::M3d *dV = nullptr,
  const EigenSupport::M3d *d2F = nullptr, EigenSupport::V3d *d2S = nullptr, EigenSupport::M3d *d2U = nullptr, EigenSupport::M3d *d2V = nullptr);

void unorderedSquareMatrixSVD2(const EigenSupport::M2d &F, EigenSupport::V2d &S, EigenSupport::M2d *U = nullptr, EigenSupport::M2d *V = nullptr);
void unorderedSquareMatrixSVD2Derivatices(const EigenSupport::M2d &F, const EigenSupport::V2d &S, const EigenSupport::M2d &U, const EigenSupport::M2d &V,
  const EigenSupport::M2d &dF, EigenSupport::V2d &dS, EigenSupport::M2d *dU = nullptr, EigenSupport::M2d *dV = nullptr,
  const EigenSupport::M2d *d2F = nullptr, EigenSupport::V2d *d2S = nullptr, EigenSupport::M2d *d2U = nullptr, EigenSupport::M2d *d2V = nullptr);

}  // namespace SVDDerivatives
}  // namespace NonlinearOptimization
}  // namespace pgo