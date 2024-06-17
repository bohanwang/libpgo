/*
author: Bohan Wang
copyright to USC
*/

#pragma once

#include "EigenSupport.h"

namespace pgo
{
namespace ES = EigenSupport;
namespace NonlinearOptimization
{
namespace PolarDecompositionDerivatives
{
// function
EigenSupport::M3d computeS(const EigenSupport::M3d &F);
EigenSupport::M3d computeR(const EigenSupport::M3d &F, const EigenSupport::M3d &S, EigenSupport::M3d *SInvOut = nullptr);

// 1st order derivatives
EigenSupport::M3d compute_dS_dF(const EigenSupport::M3d &F, int ei, EigenSupport::M3d *ROut = nullptr, EigenSupport::M3d *SInvOut = nullptr, EigenSupport::M3d *dFdFOut = nullptr);
EigenSupport::M3d compute_dR_dF(const EigenSupport::M3d &F, int ei);
void computeFirstOrderDerivatives(const EigenSupport::M3d &F, const EigenSupport::M3d &R, const EigenSupport::M3d &S, const EigenSupport::M3d &SInv,
  EigenSupport::M3d dRdF[9], EigenSupport::M3d dSdF[9],
  EigenSupport::M3d dFdFOut[9] = nullptr, EigenSupport::M3d dAdFOut[9] = nullptr, EigenSupport::M3d dsqrtA_dAOut[9] = nullptr, EigenSupport::M9d *MInvOut = nullptr);

// 2nd order derivatives
void computeSecondOrderDerivatives(const EigenSupport::M3d &F, const EigenSupport::M3d &R, const EigenSupport::M3d &S, const EigenSupport::M3d &SInv, const EigenSupport::M9d &MInv,
  const EigenSupport::M3d dRdF[9], const EigenSupport::M3d dSdF[9], const EigenSupport::M3d dFdF[9], const EigenSupport::M3d dAdF[9], const EigenSupport::M3d dsqrtA_dA[9],
  EigenSupport::M3d d2S_dF2[9][9], EigenSupport::M3d d2R_dF2[9][9]);

// compute everything
void compute(const EigenSupport::M3d &F, EigenSupport::M3d &R, EigenSupport::M3d &S, EigenSupport::M3d dSdF[9], EigenSupport::M3d dRdF[9],
  EigenSupport::M3d d2SdF2[9][9], EigenSupport::M3d d2RdF2[9][9]);
void fdTest(const EigenSupport::M3d &F);
}  // namespace PolarDecompositionDerivatives
}  // namespace NonlinearOptimization
}  // namespace pgo