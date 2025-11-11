#include "elasticModel3DMooneyRivlin.h"

using namespace pgo;
using namespace pgo::SolidDeformationModel;
namespace ES = pgo::EigenSupport;

namespace pgo::SolidDeformationModel
{
// Flatten a 3x3 into a 9x1 (column-major vec)
inline ES::V9d vec(const ES::M3d &A)
{
  ES::V9d v;
  int k = 0;
  for (int j = 0; j < 3; ++j)
    for (int i = 0; i < 3; ++i)
      v(k++) = A(i, j);
  return v;
}

inline double inner(const ES::M3d &A, const ES::M3d &B)
{
  return (A.cwiseProduct(B)).sum();
}

// Action of the Hessian of I2 on dF (linear operator form)
// I2 = 0.5 * (I1^2 - tr(C^2)), with C = F^T F, I1 = tr(C)
// ∇_F I2 = 2 (I1 F - F C)
inline ES::M3d H_I2_apply(const ES::M3d &F, const ES::M3d &C, double I1, const ES::M3d &dF)
{
  const double s = 2.0 * inner(F, dF);  // δI1 = 2 F : dF
  ES::M3d HI2 = ES::M3d::Zero();

  // δ{ ∇I2 } = 2[ (δI1)F + I1 δF - δF C - F(δF^T F + F^T δF) ]
  HI2 = 2.0 * (s * F + I1 * dF - dF * C - F * (dF.transpose() * F + F.transpose() * dF));
  return HI2;
}

// Action of the Hessian of J (i.e., δ(J F^{-T})) on dF
inline ES::M3d H_J_apply(const ES::M3d &Finv, const ES::M3d &FinvT, double J, const ES::M3d &dF)
{
  const double trTerm = inner(FinvT, dF);                  // tr(F^{-1} δF) = F^{-T} : δF
  const double dJ = J * trTerm;                            // δJ
  const ES::M3d dFinvT = -FinvT * dF.transpose() * FinvT;  // δ(F^{-T}) = -F^{-T} δF^T F^{-T}
  return dJ * FinvT + J * dFinvT;                          // δ(JF^{-T})
}

// Action of the Hessian of I1bar and I2bar on dF
inline ES::M3d H_I1bar_apply(const ES::M3d &F, const ES::M3d &Finv, const ES::M3d &FinvT,
  double J, double Jm23, double I1, const ES::M3d &dF)
{
  // Precomputations
  const ES::M3d gradI1 = 2.0 * F;
  const double s = inner(FinvT, dF);             // tr(F^{-1} δF)
  const double dJm23 = (-2.0 / 3.0) * Jm23 * s;  // δ( J^{-2/3} )
  const double dI1 = inner(gradI1, dF);          // δI1 = 2 F : δF
  const ES::M3d dFinvT = -FinvT * dF.transpose() * FinvT;

  // ∇I1bar = J^{-2/3} ∇I1 + I1 * (-2/3) J^{-2/3} F^{-T}
  // δ(∇I1bar) = δ(J^{-2/3}) ∇I1 + J^{-2/3} δ(∇I1)
  //           + δI1 * (-2/3) J^{-2/3} F^{-T}
  //           + I1 * (-2/3) δ(J^{-2/3}) F^{-T}
  //           + I1 * (-2/3) J^{-2/3} δ(F^{-T})
  return dJm23 * gradI1 + Jm23 * (2.0 * dF) + dI1 * (-2.0 / 3.0) * Jm23 * FinvT + I1 * (-2.0 / 3.0) * dJm23 * FinvT + I1 * (-2.0 / 3.0) * Jm23 * dFinvT;
}

inline ES::M3d H_I2bar_apply(const ES::M3d &F, const ES::M3d &C,
  const ES::M3d &Finv, const ES::M3d &FinvT,
  double J, double Jm43, double I1, double I2, const ES::M3d &dF)
{
  const ES::M3d gradI2 = 2.0 * (I1 * F - F * C);
  const double s = inner(FinvT, dF);             // tr(F^{-1} δF)
  const double dJm43 = (-4.0 / 3.0) * Jm43 * s;  // δ(J^{-4/3})
  const double dI2 = inner(gradI2, dF);
  const ES::M3d HI2 = H_I2_apply(F, C, I1, dF);
  const ES::M3d dFinvT = -FinvT * dF.transpose() * FinvT;

  // ∇I2bar = J^{-4/3} ∇I2 + I2 * (-4/3) J^{-4/3} F^{-T}
  return dJm43 * gradI2 + Jm43 * HI2 + dI2 * (-4.0 / 3.0) * Jm43 * FinvT + I2 * (-4.0 / 3.0) * dJm43 * FinvT + I2 * (-4.0 / 3.0) * Jm43 * dFinvT;
}
}  // namespace pgo::SolidDeformationModel

ElasticModel3DMooneyRivlin::ElasticModel3DMooneyRivlin(int N_, const double *Cpq_, int M_, const double *D_):
  N(N_), M(M_), Cpq(N + 1, N + 1), D(M)
{
  for (int p = 0; p <= N; ++p) {
    for (int q = 0; q <= N; ++q) {
      Cpq(p, q) = Cpq_[q * (N + 1) + p];
    }
  }

  for (int m = 0; m < M; ++m) {
    D(m) = D_[m];
  }
}

void ElasticModel3DMooneyRivlin::updateParameters(int N_, const double *Cpq_, int M_, const double *D_)
{
  N = N_;
  M = M_;
  Cpq.resize(N + 1, N + 1);
  D.resize(M);

  for (int p = 0; p <= N; ++p) {
    for (int q = 0; q <= N; ++q) {
      Cpq(p, q) = Cpq_[q * (N + 1) + p];
    }
  }

  for (int m = 0; m < M; ++m) {
    D(m) = D_[m];
  }
}

double ElasticModel3DMooneyRivlin::compute_psi(const double *param, const double _F[9],
  const double U[9], const double V[9], const double S[3]) const
{
  const ES::M3d F = ES::Mp<const ES::M3d>(_F);

  // Kinematics and invariants
  const ES::M3d C = F.transpose() * F;
  const double I1 = C.trace();
  const double I2 = 0.5 * (I1 * I1 - (C * C).trace());
  const double J = F.determinant();
  const ES::M3d Finv = F.fullPivHouseholderQr().inverse();

  const double Jm23 = std::pow(J * J, -1.0 / 3.0);
  const double Jm43 = std::pow(J * J * J * J, -1.0 / 3.0);
  const double I1bar = Jm23 * I1;
  const double I2bar = Jm43 * I2;

  double W = 0.0;
  const double f1 = I1bar - 3.0;
  const double f2 = I2bar - 3.0;

  for (int p = 0; p <= N; ++p) {
    for (int q = 0; q <= N; ++q) {
      if (p == 0 && q == 0)
        continue;

      const double c = Cpq(p, q);

      if (c == 0.0)
        continue;

      W += c * std::pow(f1, p) * std::pow(f2, q);
    }
  }

  for (int m = 1; m <= M; ++m) {
    const double invDm = 1.0 / D[m - 1];
    W += invDm * std::pow(J - 1.0, 2 * m);
  }

  return W;
}

void ElasticModel3DMooneyRivlin::compute_P(const double *param, const double _F[9],
  const double U[9], const double V[9], const double S[3], double P[9]) const
{
  const ES::M3d F = ES::Mp<const ES::M3d>(_F);

  // Kinematics and invariants
  const ES::M3d C = F.transpose() * F;
  const double I1 = C.trace();
  const double I2 = 0.5 * (I1 * I1 - (C * C).trace());
  const double J = F.determinant();
  const ES::M3d Finv = F.fullPivHouseholderQr().inverse();
  const ES::M3d FinvT = Finv.transpose();

  const double Jm23 = std::pow(J * J, -1.0 / 3.0);
  const double Jm43 = std::pow(J * J * J * J, -1.0 / 3.0);
  const double I1bar = Jm23 * I1;
  const double I2bar = Jm43 * I2;

  // Gradient building blocks
  const ES::M3d gradI1 = 2.0 * F;
  const ES::M3d gradI2 = 2.0 * (I1 * F - F * C);
  const ES::M3d gradJ = J * FinvT;  // ∇_F J
  const ES::M3d g1bar = Jm23 * gradI1 + I1 * (-2.0 / 3.0) * Jm23 * FinvT;
  const ES::M3d g2bar = Jm43 * gradI2 + I2 * (-4.0 / 3.0) * Jm43 * FinvT;

  const double f1 = I1bar - 3.0;
  const double f2 = I2bar - 3.0;

  // (b) First derivative dW/dF ---------------------------------------------
  ES::M3d dW_dF;
  dW_dF.setZero();

  for (int p = 0; p <= N; ++p) {
    for (int q = 0; q <= N; ++q) {
      if (p == 0 && q == 0)
        continue;

      const double c = Cpq(p, q);
      if (c == 0.0)
        continue;

      const double pow1 = (p > 0) ? std::pow(f1, p - 1) : 0.0;
      const double pow2 = (q > 0) ? std::pow(f2, q - 1) : 0.0;

      const double alpha = (p > 0) ? c * p * pow1 * std::pow(f2, q) : 0.0;  // d/dF contribution via I1bar
      const double beta = (q > 0) ? c * q * std::pow(f1, p) * pow2 : 0.0;   // via I2bar

      dW_dF += alpha * g1bar + beta * g2bar;
    }
  }

  // volumetric gradient: φ'(J) * ∇J
  double phiPrime = 0.0;
  for (int m = 1; m <= M; ++m) {
    const double invDm = 1.0 / D[m - 1];
    phiPrime += invDm * (2.0 * m) * std::pow(J - 1.0, 2 * m - 1);
  }

  dW_dF += phiPrime * gradJ;

  (ES::Mp<ES::M3d>(P)) = dW_dF;
}

void ElasticModel3DMooneyRivlin::compute_dPdF(const double *param, const double _F[9],
  const double U[9], const double V[9], const double S[3], double dPdF[81]) const
{
  const ES::M3d F = ES::Mp<const ES::M3d>(_F);

  // Kinematics and invariants
  const ES::M3d C = F.transpose() * F;
  const double I1 = C.trace();
  const double I2 = 0.5 * (I1 * I1 - (C * C).trace());
  const double J = F.determinant();
  const ES::M3d Finv = F.fullPivHouseholderQr().inverse();
  const ES::M3d FinvT = Finv.transpose();

  const double Jm23 = std::pow(J * J, -1.0 / 3.0);
  const double Jm43 = std::pow(J * J * J * J, -1.0 / 3.0);
  const double I1bar = Jm23 * I1;
  const double I2bar = Jm43 * I2;

  // Gradient building blocks
  const ES::M3d gradI1 = 2.0 * F;
  const ES::M3d gradI2 = 2.0 * (I1 * F - F * C);
  const ES::M3d gradJ = J * FinvT;  // ∇_F J
  const ES::M3d g1bar = Jm23 * gradI1 + I1 * (-2.0 / 3.0) * Jm23 * FinvT;
  const ES::M3d g2bar = Jm43 * gradI2 + I2 * (-4.0 / 3.0) * Jm43 * FinvT;

  const double f1 = I1bar - 3.0;
  const double f2 = I2bar - 3.0;

  // volumetric
  double phiPrime = 0.0, phiDoublePrime = 0.0;
  for (int m = 1; m <= M; ++m) {
    const double invDm = 1.0 / D[m - 1];
    phiPrime += invDm * (2.0 * m) * std::pow(J - 1.0, 2 * m - 1);
    phiDoublePrime += invDm * (2.0 * m) * (2.0 * m - 1.0) * std::pow(J - 1.0, 2 * m - 2);
  }

  auto H_apply = [&](const ES::M3d &dF) -> ES::M3d {
    ES::M3d out = ES::M3d::Zero();

    // Distortional part
    const double dI1bar = inner(g1bar, dF);
    const double dI2bar = inner(g2bar, dF);

    for (int p = 0; p <= N; ++p) {
      for (int q = 0; q <= N; ++q) {
        if (p == 0 && q == 0)
          continue;

        const double c = Cpq(p, q);
        if (c == 0.0)
          continue;

        const double pow_f1_p = std::pow(f1, p);
        const double pow_f2_q = std::pow(f2, q);
        const double pow_f1_p1 = (p > 0) ? std::pow(f1, p - 1) : 0.0;
        const double pow_f2_q1 = (q > 0) ? std::pow(f2, q - 1) : 0.0;
        const double pow_f1_p2 = (p > 1) ? std::pow(f1, p - 2) : 0.0;
        const double pow_f2_q2 = (q > 1) ? std::pow(f2, q - 2) : 0.0;

        const double alpha = (p > 0) ? c * p * pow_f1_p1 * pow_f2_q : 0.0;
        const double beta = (q > 0) ? c * q * pow_f1_p * pow_f2_q1 : 0.0;

        // Derivatives of alpha and beta w.r.t. F through I1bar, I2bar
        const double a11 = (p > 1) ? c * p * (p - 1) * pow_f1_p2 * pow_f2_q : 0.0;      // ∂α/∂I1bar
        const double a12 = (p > 0 && q > 0) ? c * p * q * pow_f1_p1 * pow_f2_q1 : 0.0;  // ∂α/∂I2bar
        const double b11 = a12;                                                         // ∂β/∂I1bar
        const double b12 = (q > 1) ? c * q * (q - 1) * pow_f1_p * pow_f2_q2 : 0.0;      // ∂β/∂I2bar

        // α * H(I1bar) + β * H(I2bar)
        out += alpha * H_I1bar_apply(F, Finv, FinvT, J, Jm23, I1, dF);
        out += beta * H_I2bar_apply(F, C, Finv, FinvT, J, Jm43, I1, I2, dF);

        // (∂α/∂F) g1bar + (∂β/∂F) g2bar
        out += (a11 * dI1bar + a12 * dI2bar) * g1bar;
        out += (b11 * dI1bar + b12 * dI2bar) * g2bar;
      }
    }

    // Volumetric part: φ''(J) (∇J:δF) ∇J + φ'(J) * H_J_apply[δF]
    out += phiDoublePrime * inner(gradJ, dF) * gradJ + phiPrime * H_J_apply(Finv, FinvT, J, dF);

    return out;
  };

  // Assemble 9x9 Hessian by applying H to each basis perturbation δF = E_ij
  ES::M9d d2W_dF2;

  d2W_dF2.setZero();
  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 3; ++i) {
      ES::M3d dF = ES::M3d::Zero();
      dF(i, j) = 1.0;

      const ES::V9d col = vec(H_apply(dF));

      const int k = i + 3 * j;  // column-major index of (i,j)
      d2W_dF2.col(k) = col;
    }
  }

  (ES::Mp<ES::M9d>(dPdF)) = d2W_dF2;
}