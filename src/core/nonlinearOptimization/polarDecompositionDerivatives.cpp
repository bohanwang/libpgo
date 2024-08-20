/*
author: Bohan Wang
copyright to USC
*/

#include "polarDecompositionDerivatives.h"

#include <iostream>

namespace pgo
{
namespace NonlinearOptimization
{
namespace PolarDecompositionDerivatives
{
ES::M9d KroneckerProduct3(const ES::M3d &A, const ES::M3d &B)
{
  ES::M9d ret;

  for (int col = 0; col < 3; col++) {
    for (int row = 0; row < 3; row++) {
      ret.block<3, 3>(row * 3, col * 3) = A(row, col) * B;
    }
  }

  return ret;
}

ES::M9d KroneckerSum3(const ES::M3d &A, const ES::M3d &B)
{
  ES::M3d I3 = ES::M3d::Identity();
  ES::M9d AI = KroneckerProduct3(A, I3);
  ES::M9d IB = KroneckerProduct3(I3, B);
  return AI + IB;
}

ES::M3d compute_dsqrtA_dAi(const ES::M9d &MInv, int idx)
{
  ES::V9d dAi = ES::V9d::Zero();
  dAi[idx] = 1.0;

  ES::V9d dSqrtA = MInv * dAi;
  return Eigen::Map<const ES::M3d>(dSqrtA.data());
}

ES::M3d computeS(const ES::M3d &F)
{
  ES::M3d A = F.transpose() * F;
  Eigen::SelfAdjointEigenSolver<ES::M3d> eig(A);
  ES::M3d V = eig.eigenvectors();
  ES::V3d D = eig.eigenvalues();
  ES::V3d DSqrt = D.array().sqrt();

  return V * DSqrt.asDiagonal() * V.transpose();
}

ES::M3d computeR(const ES::M3d &F, const ES::M3d &S, ES::M3d *SInvOut)
{
  ES::M3d SInv = S.fullPivLu().inverse();
  ES::M3d R = F * SInv;

  if (SInvOut)
    *SInvOut = SInv;

  return R;
}

inline int toIndex(int row, int col)
{
  return col * 3 + row;
};

ES::M3d compute_dF_dF(int ele)
{
  ES::M3d dFdF = ES::M3d::Zero();
  dFdF.data()[ele] = 1.0;
  return dFdF;
}

ES::M3d compute_dA_dF(const ES::M3d &F, const ES::M3d &dFdF)
{
  return dFdF.transpose() * F + F.transpose() * dFdF;
}

ES::M3d compute_dS_dF(const ES::M3d &dAdF, const ES::M3d dsqrtA_dA[9])
{
  ES::M3d dSdF = ES::M3d::Zero();
  for (int j = 0; j < 9; j++) {
    dSdF += dsqrtA_dA[j] * dAdF.data()[j];
  }

  return dSdF;
}

ES::M3d compute_dR_dF(const ES::M3d &R, const ES::M3d &SInv, const ES::M3d &dFdF, const ES::M3d &dSdF)
{
  ES::M3d dRdFi_S = dFdF - R * dSdF;
  return dRdFi_S * SInv;
}

ES::M3d compute_dS_dF(const ES::M3d &F, int ele, ES::M3d *ROut, ES::M3d *SInvOut, ES::M3d *dFdFOut)
{
  ES::M3d S = computeS(F);
  ES::M9d M = KroneckerSum3(S, S);
  ES::M9d MInv = M.fullPivLu().inverse();
  ES::M3d SInv, R;
  R = computeR(F, S, &SInv);

  ES::M3d dFdF = compute_dF_dF(ele);
  ES::M3d dAdF = compute_dA_dF(F, dFdF);

  ES::M3d dsqrtA_dA[9];
  for (int i = 0; i < 9; i++)
    dsqrtA_dA[i] = compute_dsqrtA_dAi(MInv, i);

  ES::M3d dSdF = compute_dS_dF(dAdF, dsqrtA_dA);

  if (ROut)
    *ROut = R;

  if (SInvOut)
    *SInvOut = SInv;

  if (dFdFOut)
    *dFdFOut = dFdF;

  return dSdF;
}

ES::M3d compute_dR_dF(const ES::M3d &F, int ele)
{
  ES::M3d SInv, R, dFdF;
  ES::M3d dSdF = compute_dS_dF(F, ele, &R, &SInv, &dFdF);
  ES::M3d dRdF = compute_dR_dF(R, SInv, dFdF, dSdF);

  return dRdF;
}

void computeFirstOrderDerivatives(const ES::M3d &F, const ES::M3d &R, const ES::M3d &S, const ES::M3d &SInv,
  ES::M3d dRdF[9], ES::M3d dSdF[9], ES::M3d dFdFOut[9], ES::M3d dAdFOut[9], ES::M3d dsqrtA_dAOut[9], ES::M9d *MInvOut)
{
  ES::M9d M = KroneckerSum3(S, S);
  ES::M9d MInv = M.fullPivLu().inverse();

  ES::M3d dFdF[9], dAdF[9];
  ES::M3d dsqrtA_dA[9];
  for (int ei = 0; ei < 9; ei++) {
    dFdF[ei] = compute_dF_dF(ei);
    dAdF[ei] = compute_dA_dF(F, dFdF[ei]);
    dsqrtA_dA[ei] = compute_dsqrtA_dAi(MInv, ei);
  }

  for (int ei = 0; ei < 9; ei++) {
    dSdF[ei] = compute_dS_dF(dAdF[ei], dsqrtA_dA);
    dRdF[ei] = compute_dR_dF(R, SInv, dFdF[ei], dSdF[ei]);
  }

  if (MInvOut)
    *MInvOut = MInv;

  if (dFdFOut) {
    for (int ei = 0; ei < 9; ei++) {
      dFdFOut[ei] = dFdF[ei];
    }
  }

  if (dAdFOut) {
    for (int ei = 0; ei < 9; ei++) {
      dAdFOut[ei] = dAdF[ei];
    }
  }

  if (dsqrtA_dAOut) {
    for (int ei = 0; ei < 9; ei++) {
      dsqrtA_dAOut[ei] = dsqrtA_dA[ei];
    }
  }
}

ES::M3d compute_d2sqrtA_dAidAj(const ES::M3d dsqrtA_dA[9], const ES::M9d &MInv, int ei, int ej)
{
  ES::M3d Zij = dsqrtA_dA[ei] * dsqrtA_dA[ej];
  ES::M3d Zji = dsqrtA_dA[ej] * dsqrtA_dA[ei];
  ES::M3d rhs = -Zij - Zji;

  ES::V9d rhs9 = Eigen::Map<const ES::V9d>(rhs.data());
  ES::V9d x9 = MInv * rhs9;

  return Eigen::Map<const ES::M3d>(x9.data());
}

ES::M3d compute_d2S_dF2(const ES::M3d &F, int ei, int ej, const ES::M3d dsqrtA_dA[9], const ES::M9d &MInv, const ES::M3d dFdF[9])
{
  ES::M3d d2sqrtAdA2[9][9];
  for (int i = 0; i < 9; i++) {
    for (int j = 0; j < 9; j++) {
      d2sqrtAdA2[i][j] = compute_d2sqrtA_dAidAj(dsqrtA_dA, MInv, i, j);
    }
  }

  // dA / dFi
  ES::M3d dAdFi = dFdF[ei].transpose() * F + F.transpose() * dFdF[ei];
  ES::M3d dAdFj = dFdF[ej].transpose() * F + F.transpose() * dFdF[ej];

  ES::M3d ret = ES::M3d::Zero();
  for (int k = 0; k < 9; k++) {
    for (int l = 0; l < 9; l++) {
      ret += d2sqrtAdA2[k][l] * dAdFi.data()[k] * dAdFj.data()[l];
    }
  }

  ES::M3d d2AdFidFj = dFdF[ei].transpose() * dFdF[ej] + dFdF[ej].transpose() * dFdF[ei];
  for (int k = 0; k < 9; k++) {
    ret += dsqrtA_dA[k] * d2AdFidFj.data()[k];
  }

  return ret;
}

ES::M3d compute_d2R_dF2([[maybe_unused]] const ES::M3d &F, int ei, int ej, const ES::M3d &R, const ES::M3d &SInv, const ES::M3d dRdF[9], const ES::M3d dSdF[9], const ES::M3d &d2SdF2)
{
  ES::M3d ret = ES::M3d::Zero() - R * d2SdF2 - dRdF[ei] * dSdF[ej] - dRdF[ej] * dSdF[ei];
  return ret * SInv;
}

void computeSecondOrderDerivatives([[maybe_unused]] const ES::M3d &F, const ES::M3d &R, [[maybe_unused]] const ES::M3d &S, const ES::M3d &SInv, const ES::M9d &MInv,
  const ES::M3d dRdF[9], const ES::M3d dSdF[9], const ES::M3d dFdF[9], const ES::M3d dAdF[9], const ES::M3d dsqrtA_dA[9],
  ES::M3d d2S_dF2[9][9], ES::M3d d2R_dF2[9][9])
{
  ES::M3d d2sqrtAdA2[9][9];
  for (int i = 0; i < 9; i++) {
    for (int j = 0; j < 9; j++) {
      d2sqrtAdA2[i][j] = compute_d2sqrtA_dAidAj(dsqrtA_dA, MInv, i, j);
    }
  }

  for (int ei = 0; ei < 9; ei++) {
    for (int ej = 0; ej < 9; ej++) {
      // dA / dFi
      ES::M3d dAdFi = dAdF[ei];
      ES::M3d dAdFj = dAdF[ej];

      ES::M3d ret = ES::M3d::Zero();
      for (int k = 0; k < 9; k++) {
        for (int l = 0; l < 9; l++) {
          ret += d2sqrtAdA2[k][l] * dAdFi.data()[k] * dAdFj.data()[l];
        }
      }

      ES::M3d d2AdFidFj = dFdF[ei].transpose() * dFdF[ej] + dFdF[ej].transpose() * dFdF[ei];
      for (int k = 0; k < 9; k++) {
        ret += dsqrtA_dA[k] * d2AdFidFj.data()[k];
      }

      d2S_dF2[ei][ej] = ret;

      ES::M3d temp = ES::M3d::Zero() - R * d2S_dF2[ei][ej] - dRdF[ei] * dSdF[ej] - dRdF[ej] * dSdF[ei];
      d2R_dF2[ei][ej] = temp * SInv;
    }
  }
}

void compute(const ES::M3d &F, ES::M3d &R, ES::M3d &S, ES::M3d dSdF[9], ES::M3d dRdF[9],
  ES::M3d d2SdF2[9][9], ES::M3d d2RdF2[9][9])
{
  S = computeS(F);
  ES::M9d M = KroneckerSum3(S, S);
  ES::M9d MInv = M.fullPivLu().inverse();
  ES::M3d SInv;
  R = computeR(F, S, &SInv);

  if (dSdF == nullptr && dRdF == nullptr &&
    d2RdF2 == nullptr && d2SdF2 == nullptr)
    return;

  ES::M3d dFdF[9];
  for (int col = 0; col < 3; col++) {
    for (int row = 0; row < 3; row++) {
      int idx = toIndex(row, col);
      dFdF[idx] = compute_dF_dF(idx);
    }
  }

  ES::M3d dsqrtA_dA[9];
  for (int i = 0; i < 9; i++)
    dsqrtA_dA[i] = compute_dsqrtA_dAi(MInv, i);

  for (int col = 0; col < 3; col++) {
    for (int row = 0; row < 3; row++) {
      int idx = toIndex(row, col);
      dSdF[idx] = ES::M3d::Zero();

      // dA / dFi
      ES::M3d dAdFi = dFdF[idx].transpose() * F + F.transpose() * dFdF[idx];
      for (int j = 0; j < 9; j++) {
        dSdF[idx] += dsqrtA_dA[j] * dAdFi.data()[j];
      }
    }
  }

  for (int col = 0; col < 3; col++) {
    for (int row = 0; row < 3; row++) {
      int idx = toIndex(row, col);
      ES::M3d dRdFi_S = dFdF[idx] - R * dSdF[idx];
      dRdF[idx] = dRdFi_S * SInv;
    }
  }

  if (d2SdF2 == nullptr && d2RdF2 == nullptr)
    return;

  ES::M3d d2sqrtAdA2[9][9];
  for (int i = 0; i < 9; i++) {
    for (int j = 0; j < 9; j++) {
      ES::M3d Zij = dsqrtA_dA[i] * dsqrtA_dA[j];
      ES::M3d Zji = dsqrtA_dA[j] * dsqrtA_dA[i];
      ES::M3d rhs = -Zij - Zji;

      ES::V9d rhs9 = Eigen::Map<const ES::V9d>(rhs.data());
      ES::V9d x9 = MInv * rhs9;

      d2sqrtAdA2[i][j] = Eigen::Map<const ES::M3d>(x9.data());
    }
  }

  for (int i = 0; i < 9; i++) {
    for (int j = 0; j < 9; j++) {
      // dA / dFi
      ES::M3d dAdFi = dFdF[i].transpose() * F + F.transpose() * dFdF[i];
      ES::M3d dAdFj = dFdF[j].transpose() * F + F.transpose() * dFdF[j];

      ES::M3d ret = ES::M3d::Zero();
      for (int k = 0; k < 9; k++) {
        for (int l = 0; l < 9; l++) {
          ret += d2sqrtAdA2[k][l] * dAdFi.data()[k] * dAdFj.data()[l];
        }
      }

      ES::M3d d2AdFidFj = dFdF[i].transpose() * dFdF[j] + dFdF[j].transpose() * dFdF[i];
      for (int k = 0; k < 9; k++) {
        ret += dsqrtA_dA[k] * d2AdFidFj.data()[k];
      }

      d2SdF2[i][j] = ret;
    }
  }

  for (int i = 0; i < 9; i++) {
    for (int j = 0; j < 9; j++) {
      ES::M3d ret = ES::M3d::Zero() - R * d2SdF2[i][j] - dRdF[i] * dSdF[j] - dRdF[j] * dSdF[i];
      d2RdF2[i][j] = ret * SInv;
    }
  }
}

void fdTest(const ES::M3d &F)
{
  std::cout << F << std::endl;

  ES::M3d Q, S, SInv;
  S = computeS(F);
  Q = computeR(F, S, &SInv);

  ES::M3d dFdF[9], dAdF[9], dsqrtAdA[9];
  ES::M3d dRdF[9], dSdF[9];
  ES::M9d MInv;
  computeFirstOrderDerivatives(F, Q, S, SInv, dRdF, dSdF, dFdF, dAdF, dsqrtAdA, &MInv);

  ES::M3d d2SdF2[9][9], d2RdF2[9][9];
  computeSecondOrderDerivatives(F, Q, S, SInv, MInv, dRdF, dSdF, dFdF, dAdF, dsqrtAdA, d2SdF2, d2RdF2);

  const double coeffs5[5] = {
    1.0 / 12.0,
    -2.0 / 3.0,
    0,
    2.0 / 3.0,
    -1.0 / 12.0
  };

  std::cout << "==================================================\n";
  std::cout << "dSdF\n";

  for (int ele = 0; ele < 9; ele++) {
    double eps = 1e-7;
    ES::M3d ret = ES::M3d::Zero();
    for (int i = -2; i <= 2; i++) {
      ES::M3d Fcur = F;
      Fcur.data()[ele] = F.data()[ele] + eps * i;

      ES::M3d Scur = computeS(Fcur);
      ret += Scur * coeffs5[i + 2];
    }
    ret /= eps;

    ES::M3d diff = ret - dSdF[ele];
    std::cout << ele << ": rel=" << diff.norm() / ret.norm() << "; abs=" << diff.norm() << std::endl;
  }

  std::cout << "==================================================\n";
  std::cout << "dRdF\n";

  for (int ele = 0; ele < 9; ele++) {
    double eps = 1e-7;
    ES::M3d ret = ES::M3d::Zero();
    for (int i = -2; i <= 2; i++) {
      ES::M3d Fcur = F;
      Fcur.data()[ele] = F.data()[ele] + eps * i;

      ES::M3d SCur = computeS(Fcur);
      ES::M3d Rcur = computeR(Fcur, SCur);
      ret += Rcur * coeffs5[i + 2];
    }
    ret /= eps;

    ES::M3d diff = ret - dRdF[ele];
    std::cout << ele << ": rel=" << diff.norm() / ret.norm() << "; abs=" << diff.norm() << std::endl;
  }

  std::cout << "==================================================\n";
  std::cout << "d2SdF2\n";

  for (int ei = 0; ei < 9; ei++) {
    for (int ej = 0; ej < 9; ej++) {
      double eps = 1e-7;
      ES::M3d ret = ES::M3d::Zero();
      for (int i = -2; i <= 2; i++) {
        ES::M3d Fcur = F;
        Fcur.data()[ej] = F.data()[ej] + eps * i;

        ES::M3d Rcur = compute_dS_dF(Fcur, ei);
        ret += Rcur * coeffs5[i + 2];
      }
      ret /= eps;

      ES::M3d diff = ret - d2SdF2[ei][ej];
      std::cout << ei << ',' << ej << ": rel=" << diff.norm() / ret.norm() << "; abs=" << diff.norm() << std::endl;
    }
  }

  std::cout << "==================================================\n";
  std::cout << "d2RdF2\n";

  for (int ei = 0; ei < 9; ei++) {
    for (int ej = 0; ej < 9; ej++) {
      double eps = 1e-7;
      ES::M3d ret = ES::M3d::Zero();
      for (int i = -2; i <= 2; i++) {
        ES::M3d Fcur = F;
        Fcur.data()[ej] = F.data()[ej] + eps * i;

        ES::M3d Scur = compute_dR_dF(Fcur, ei);
        ret += Scur * coeffs5[i + 2];
      }
      ret /= eps;

      ES::M3d diff = ret - d2RdF2[ei][ej];
      // std::cout << ret << '\n' << diff << '\n' << dSdF << std::endl;
      std::cout << ei << ',' << ej << ": rel=" << diff.norm() / ret.norm() << "; abs=" << diff.norm() << std::endl;
    }
  }
}
}  // namespace PolarDecompositionDerivatives
}  // namespace NonlinearOptimization
}  // namespace pgo
