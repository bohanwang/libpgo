#include "elasticModel2DFundamentalFormsFabric.h"

namespace pgo
{
namespace ES = pgo::EigenSupport;

namespace SolidDeformationModel
{

inline ES::V2d normalizeWithMetric(const ES::V2d &v, const ES::M2d &A0)
{
  const double n = std::sqrt(v.dot(A0 * v));

  if (!(n > 0.0))
    throw std::runtime_error("Direction has nonpositive metric length.");

  return v / n;
}

inline ES::V3d symToMandel(const ES::M2d &S)
{
  ES::V3d v;

  v << S(0, 0), S(1, 1), std::sqrt(2.0) * S(0, 1);

  return v;
}

inline ES::M2d mandelToSym(const ES::V3d &v)
{
  ES::M2d S;

  S.setZero();
  S(0, 0) = v(0);
  S(1, 1) = v(1);
  S(0, 1) = S(1, 0) = v(2) / std::sqrt(2.0);

  return S;
}

inline ES::M4d pack_Mandel3To4(const ES::M3d &M3)
{
  // int index[4] = { 0, 2, 2, 1 };
  // ES::M4d M4;
  // M4.setZero();
  // for (int i = 0; i < 4; i++) {
  //   for (int j = 0; j < 4; j++) {
  //     M4(i, j) = M3(index[i], index[j]);
  //   }
  // }

  ES::M4d M4;
  Eigen::Matrix<double, 3, 4> S;
  S << 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0,
    0.0, 1 / std::sqrt(2.0), 1 / std::sqrt(2.0), 0.0;

  M4 = S.transpose() * M3 * S;

  return M4;
}

inline ES::M2d symmetrize(const ES::M2d &M)
{
  return 0.5 * (M + M.transpose());
}

// ===================== Smooth ramp (C1) =====================

inline double smoothRamp(double x, double eps)
{
  const double e = std::max(eps, 1e-16);  // guard
  const double srt = std::sqrt(x * x + e * e);
  return 0.5 * (x + srt);
}

inline double smoothRampGrad(double x, double eps)
{
  const double e = std::max(eps, 1e-16);  // guard
  const double srt = std::sqrt(x * x + e * e);
  return 0.5 * (1.0 + x / srt);
}

inline double smoothRampHess(double x, double eps)
{
  const double e = std::max(eps, 1e-16);  // guard
  const double srt = std::sqrt(x * x + e * e);
  return 0.5 * (e * e) / (srt * srt * srt);
}

double ElasticModel2DFundamentalFormsFabric::compute_psi_a(const double *param, const double a_[4], const double abar_[4]) const
{
  // Unpack parameters
  // Tiny isotropic matrix term (optional)
  double mu0 = param[0];

  // Warp (I4) exponential fiber law
  double k1_4 = param[1];
  double k2_4 = param[2];

  // Weft (I6) exponential fiber law
  double k1_6 = param[3];
  double k2_6 = param[4];

  // Shear trellising with saturation: (ks/2) * d8^2 / (1 + alpha d8^2)
  double ks = param[5];
  double alpha = param[6];

  // Bending (orthotropic)
  double kappa11 = param[7];
  double kappa22 = param[8];
  double kappa12 = param[9];

  double I8_0 = param[10];

  double h = param[11];  // shell thickness

  ES::M2d a = ES::Mp<const ES::M2d>(a_);
  ES::M2d abar = ES::Mp<const ES::M2d>(abar_);
  ES::M2d abarInv = abar.inverse();
  ES::M2d C = abarInv * a;

  ES::V2d e1, e2;
  e1 = normalizeWithMetric(warpDir, abar);
  e2 = normalizeWithMetric(weftDir, abar);

  double I1 = C.trace();
  double I4 = e1.dot(C * e1);
  double I6 = e2.dot(C * e2);
  double I8 = e1.dot(C * e2);

  // --- Membrane energy pieces (C1 fibers + shear saturation) ---
  // Isotropic matrix (tiny)
  const double psi_iso = 0.5 * mu0 * (I1 - 2.0);

  // Warp fiber (I4)
  const double r4 = smoothRamp(I4 - 1.0, eps4);  // s4, ds4, d2s4
  const double E4 = std::exp(k2_4 * r4 * r4);
  const double psi4 = (k1_4 / (2.0 * k2_4)) * (E4 - 1.0);

  // Weft fiber (I6)
  const double r6 = smoothRamp(I6 - 1.0, eps6);
  const double E6 = std::exp(k2_6 * r6 * r6);
  const double psi6 = (k1_6 / (2.0 * k2_6)) * (E6 - 1.0);

  // Shear trellising (I8) with saturation
  const double d8 = I8 - I8_0;
  const double denom = 1.0 + alpha * d8 * d8;
  const double psi8 = 0.5 * ks * (d8 * d8) / denom;

  // Total membrane energy
  const double psi_mem = psi_iso + psi4 + psi6 + psi8;

  return h * psi_mem;
}

double ElasticModel2DFundamentalFormsFabric::compute_psi_b(const double *param, const double b_[4], const double abar_[4], const double bbar_[4]) const
{
  // Unpack parameters
  // Tiny isotropic matrix term (optional)
  double mu0 = param[0];

  // Warp (I4) exponential fiber law
  double k1_4 = param[1];
  double k2_4 = param[2];

  // Weft (I6) exponential fiber law
  double k1_6 = param[3];
  double k2_6 = param[4];

  // Shear trellising with saturation: (ks/2) * d8^2 / (1 + alpha d8^2)
  double ks = param[5];
  double alpha = param[6];

  // Bending (orthotropic)
  double kappa11 = param[7];
  double kappa22 = param[8];
  double kappa12 = param[9];

  double I8_0 = param[10];

  double h = param[11];  // shell thickness

  ES::M2d b = ES::Mp<const ES::M2d>(b_);
  ES::M2d abar = ES::Mp<const ES::M2d>(abar_);
  ES::M2d bbar = ES::Mp<const ES::M2d>(bbar_);

  ES::M2d abarInv = abar.inverse();
  ES::M2d b_eff = abarInv * (b - bbar);

  // -------- Bending (correct Mandel scaling) --------
  // ψ_b = 0.5*(k11 b11^2 + k22 b22^2 + 2 k12 b12^2)
  // δψ_b = k11 b11 δb11 + k22 b22 δb22 + 2 k12 b12 δb12

  const double psi_bend =
    0.5 * (kappa11 * b_eff(0, 0) * b_eff(0, 0) + kappa22 * b_eff(1, 1) * b_eff(1, 1) + 2.0 * kappa12 * b_eff(0, 1) * b_eff(0, 1));

  return h * h * h / 12 * psi_bend;
}

void ElasticModel2DFundamentalFormsFabric::compute_dpsi_da(const double *param, const double a_[4], const double abar_[4], double da_[4]) const
{
  // Unpack parameters
  // Tiny isotropic matrix term (optional)
  double mu0 = param[0];

  // Warp (I4) exponential fiber law
  double k1_4 = param[1];
  double k2_4 = param[2];

  // Weft (I6) exponential fiber law
  double k1_6 = param[3];
  double k2_6 = param[4];

  // Shear trellising with saturation: (ks/2) * d8^2 / (1 + alpha d8^2)
  double ks = param[5];
  double alpha = param[6];

  // Bending (orthotropic)
  double kappa11 = param[7];
  double kappa22 = param[8];
  double kappa12 = param[9];

  double I8_0 = param[10];

  double h = param[11];  // shell thickness

  ES::M2d a = ES::Mp<const ES::M2d>(a_);
  ES::M2d abar = ES::Mp<const ES::M2d>(abar_);
  ES::M2d abarInv = abar.inverse();
  ES::M2d C = abarInv * a;

  ES::V2d e1, e2;
  e1 = normalizeWithMetric(warpDir, abar);
  e2 = normalizeWithMetric(weftDir, abar);

  double I4 = e1.dot(C * e1);
  double I6 = e2.dot(C * e2);
  double I8 = e1.dot(C * e2);

  ES::M2d dI1_da = symmetrize(abarInv);
  ES::M2d dI4_da = symmetrize(abarInv * (e1 * e1.transpose()));
  ES::M2d dI6_da = symmetrize(abarInv * (e2 * e2.transpose()));
  ES::M2d dI8_da = symmetrize(abarInv * (e2 * e1.transpose()));

  // --- Membrane energy pieces (C1 fibers + shear saturation) ---
  // Isotropic matrix (tiny)
  const double dpsi_dI1 = 0.5 * mu0;

  // Warp fiber (I4)
  const double r4 = smoothRamp(I4 - 1.0, eps4);  // s4, ds4, d2s4
  const double dr4 = smoothRampGrad(I4 - 1.0, eps4);
  const double E4 = std::exp(k2_4 * r4 * r4);
  const double dpsi_dI4 = k1_4 * r4 * E4 * dr4;

  // Weft fiber (I6)
  const double r6 = smoothRamp(I6 - 1.0, eps6);
  const double dr6 = smoothRampGrad(I6 - 1.0, eps6);
  const double E6 = std::exp(k2_6 * r6 * r6);
  const double dpsi_dI6 = k1_6 * r6 * E6 * dr6;

  // Shear trellising (I8) with saturation
  const double d8 = I8 - I8_0;
  const double denom = 1.0 + alpha * d8 * d8;
  const double dpsi_dI8 = ks * d8 / (denom * denom);

  // -------- Gradient wrt 'a' (2x2), Hessian Haa (3x3 Mandel) --------
  ES::M2d grad_a_mat = dpsi_dI1 * dI1_da + dpsi_dI4 * dI4_da + dpsi_dI6 * dI6_da + dpsi_dI8 * dI8_da;
  (ES::Mp<ES::M2d>(da_)) = grad_a_mat * h;
}

void ElasticModel2DFundamentalFormsFabric::compute_dpsi_db(const double *param, const double b_[4], const double abar_[4], const double bbar_[4], double db_[4]) const
{
  // Unpack parameters
  // Tiny isotropic matrix term (optional)
  double mu0 = param[0];

  // Warp (I4) exponential fiber law
  double k1_4 = param[1];
  double k2_4 = param[2];

  // Weft (I6) exponential fiber law
  double k1_6 = param[3];
  double k2_6 = param[4];

  // Shear trellising with saturation: (ks/2) * d8^2 / (1 + alpha d8^2)
  double ks = param[5];
  double alpha = param[6];

  // Bending (orthotropic)
  double kappa11 = param[7];
  double kappa22 = param[8];
  double kappa12 = param[9];

  double I8_0 = param[10];

  double h = param[11];  // shell thickness

  ES::M2d abar = ES::Mp<const ES::M2d>(abar_);
  ES::M2d bbar = ES::Mp<const ES::M2d>(bbar_);
  ES::M2d b = ES::Mp<const ES::M2d>(b_);

  ES::M2d abarInv = abar.inverse();
  ES::M2d b_eff = abarInv * (b - bbar);

  // -------- Bending (correct Mandel scaling) --------
  // ψ_b = 0.5*(k11 b11^2 + k22 b22^2 + 2 k12 b12^2)
  // δψ_b = k11 b11 δb11 + k22 b22 δb22 + 2 k12 b12 δb12
  ES::M2d grad_b_mat = ES::M2d::Zero();
  grad_b_mat(0, 0) = kappa11 * b_eff(0, 0);
  grad_b_mat(1, 1) = kappa22 * b_eff(1, 1);
  grad_b_mat(0, 1) = grad_b_mat(1, 0) = kappa12 * b_eff(0, 1);  // NOTE: not 2*kappa12

  grad_b_mat *= h * h * h / 12;
  (ES::Mp<ES::M2d>(db_)) = grad_b_mat;
}

void ElasticModel2DFundamentalFormsFabric::compute_d2psi_da2(const double *param, const double a_[4], const double abar_[4], double da2_[16]) const
{
  // Unpack parameters
  // Tiny isotropic matrix term (optional)
  // double mu0 = param[0];

  // Warp (I4) exponential fiber law
  double k1_4 = param[1];
  double k2_4 = param[2];

  // Weft (I6) exponential fiber law
  double k1_6 = param[3];
  double k2_6 = param[4];

  // Shear trellising with saturation: (ks/2) * d8^2 / (1 + alpha d8^2)
  double ks = param[5];
  double alpha = param[6];

  // Bending (orthotropic)
  double kappa11 = param[7];
  double kappa22 = param[8];
  double kappa12 = param[9];

  double I8_0 = param[10];

  double h = param[11];  // shell thickness

  ES::M2d a = ES::Mp<const ES::M2d>(a_);
  ES::M2d abar = ES::Mp<const ES::M2d>(abar_);
  ES::M2d abarInv = abar.inverse();
  ES::M2d C = abarInv * a;

  ES::V2d e1, e2;
  e1 = normalizeWithMetric(warpDir, abar);
  e2 = normalizeWithMetric(weftDir, abar);

  // double I1 = C.trace();
  double I4 = e1.dot(C * e1);
  double I6 = e2.dot(C * e2);
  double I8 = e1.dot(C * e2);

  ES::M2d dI1_da = symmetrize(abarInv);
  ES::M2d dI4_da = symmetrize(abarInv * (e1 * e1.transpose()));
  ES::M2d dI6_da = symmetrize(abarInv * (e2 * e2.transpose()));
  ES::M2d dI8_da = symmetrize(abarInv * (e2 * e1.transpose()));

  ES::V3d dI1_m = symToMandel(dI1_da);
  ES::V3d dI4_m = symToMandel(dI4_da);
  ES::V3d dI6_m = symToMandel(dI6_da);
  ES::V3d dI8_m = symToMandel(dI8_da);

  // --- Membrane energy pieces (C1 fibers + shear saturation) ---
  // Isotropic matrix (tiny)
  const double d2psi_dI1 = 0.0;

  // Warp fiber (I4)
  const double r4 = smoothRamp(I4 - 1.0, eps4);  // s4, ds4, d2s4
  const double dr4 = smoothRampGrad(I4 - 1.0, eps4);
  const double d2r4 = smoothRampHess(I4 - 1.0, eps4);

  const double E4 = std::exp(k2_4 * r4 * r4);
  const double d2psi_dI4 = k1_4 * ((dr4 * dr4) * E4 * (1.0 + 2.0 * k2_4 * r4 * r4) + r4 * E4 * d2r4);

  // Weft fiber (I6)
  const double r6 = smoothRamp(I6 - 1.0, eps6);
  const double dr6 = smoothRampGrad(I6 - 1.0, eps6);
  const double d2r6 = smoothRampHess(I6 - 1.0, eps6);

  const double E6 = std::exp(k2_6 * r6 * r6);
  const double d2psi_dI6 =
    k1_6 * ((dr6 * dr6) * E6 * (1.0 + 2.0 * k2_6 * r6 * r6) + r6 * E6 * d2r6);

  // Shear trellising (I8) with saturation
  const double d8 = I8 - I8_0;
  const double denom = 1.0 + alpha * d8 * d8;
  const double d2psi_dI8 = ks * (1.0 - 3.0 * alpha * d8 * d8) / std::pow(denom, 3);

  ES::M3d Haa = ES::M3d::Zero();
  Haa.noalias() += d2psi_dI1 * (dI1_m * dI1_m.transpose());  // zero anyway
  Haa.noalias() += d2psi_dI4 * (dI4_m * dI4_m.transpose());
  Haa.noalias() += d2psi_dI6 * (dI6_m * dI6_m.transpose());
  Haa.noalias() += d2psi_dI8 * (dI8_m * dI8_m.transpose());

  ES::M4d Haa4;
  Haa4 = pack_Mandel3To4(Haa) * h;

  (ES::Mp<ES::M4d>(da2_)) = Haa4;
}

void ElasticModel2DFundamentalFormsFabric::compute_d2psi_db2(const double *param, const double b_[4], const double abar_[4], const double bbar_[4], double db2_[16]) const
{
  // Unpack parameters
  // Tiny isotropic matrix term (optional)
  // double mu0 = param[0];

  // Warp (I4) exponential fiber law
  double k1_4 = param[1];
  double k2_4 = param[2];

  // Weft (I6) exponential fiber law
  double k1_6 = param[3];
  double k2_6 = param[4];

  // Shear trellising with saturation: (ks/2) * d8^2 / (1 + alpha d8^2)
  double ks = param[5];
  double alpha = param[6];

  // Bending (orthotropic)
  double kappa11 = param[7];
  double kappa22 = param[8];
  double kappa12 = param[9];

  double I8_0 = param[10];

  double h = param[11];  // shell thickness

  ES::M2d abar = ES::Mp<const ES::M2d>(abar_);
  ES::M2d bbar = ES::Mp<const ES::M2d>(bbar_);
  ES::M2d b = ES::Mp<const ES::M2d>(b_);

  ES::M2d abarInv = abar.inverse();
  ES::M2d b_eff = abarInv * (b - bbar);  ///???

  // In Mandel, the Hbb diag must be [k11, k22, k12] to match energy correctly.
  ES::M3d Hbb = ES::M3d::Zero();
  Hbb(0, 0) = kappa11;
  Hbb(1, 1) = kappa22;
  Hbb(2, 2) = kappa12;  // NOTE: not 2*kappa12

  ES::M4d Hbb4;
  Hbb4 = pack_Mandel3To4(Hbb);

  Hbb4 *= h * h * h / 12;

  (ES::Mp<ES::M4d>(db2_)) = Hbb4;
}

#if 0

void ElasticModelThinShellFabric::compute_dP_dparam(const double *param, int i, const double F[9],
  const double U[9], const double V[9], const double S[3], double *ret) const
{
  // Unpack parameters
  // Tiny isotropic matrix term (optional)
  double mu0 = param[0];

  // Warp (I4) exponential fiber law
  double k1_4 = param[1];
  double k2_4 = param[2];

  // Weft (I6) exponential fiber law
  double k1_6 = param[3];
  double k2_6 = param[4];

  // Shear trellising with saturation: (ks/2) * d8^2 / (1 + alpha d8^2)
  double ks = param[5];
  double alpha = param[6];

  // Bending (orthotropic)
  double kappa11 = param[7];
  double kappa22 = param[8];
  double kappa12 = param[9];

  double I8_0 = param[10];

  double h = param[11];  // shell thickness

  ES::M2d a = unpack_a(F);
  ES::M2d b = unpack_b(F);
  ES::M2d abar = unpack_abar(F);
  ES::M2d bbar = unpack_bbar(F);

  ES::M2d abarInv = abar.inverse();
  ES::M2d C = abarInv * a;
  ES::M2d b_eff = abarInv * (b - bbar);

  ES::V2d e1, e2;
  e1 = normalizeWithMetric(warpDir, abar);
  e2 = normalizeWithMetric(weftDir, abar);

  // double I1 = C.trace();
  double I4 = e1.dot(C * e1);
  double I6 = e2.dot(C * e2);
  double I8 = e1.dot(C * e2);

  ES::M2d dI1_da = symmetrize(abarInv);
  ES::M2d dI4_da = symmetrize(abarInv * (e1 * e1.transpose()));
  ES::M2d dI6_da = symmetrize(abarInv * (e2 * e2.transpose()));
  ES::M2d dI8_da = symmetrize(abarInv * (e2 * e1.transpose()));

  ES::V3d dI1_m = symToMandel(dI1_da);
  ES::V3d dI4_m = symToMandel(dI4_da);
  ES::V3d dI6_m = symToMandel(dI6_da);
  ES::V3d dI8_m = symToMandel(dI8_da);

  // --- Membrane energy pieces (C1 fibers + shear saturation) ---
  // Isotropic matrix (tiny)
  const double dpsi_dI1 = 0.5 * mu0;

  // Warp fiber (I4)
  const double r4 = smoothRamp(I4 - 1.0, eps4);  // s4, ds4, d2s4
  const double dr4 = smoothRampGrad(I4 - 1.0, eps4);
  const double E4 = std::exp(k2_4 * r4 * r4);
  const double dpsi_dI4 = k1_4 * r4 * E4 * dr4;

  // Weft fiber (I6)
  const auto r6 = smoothRamp(I6 - 1.0, eps6);
  const double dr6 = smoothRampGrad(I6 - 1.0, eps6);
  const double E6 = std::exp(k2_6 * r6 * r6);
  const double dpsi_dI6 = k1_6 * r6 * E6 * dr6;

  // Shear trellising (I8) with saturation
  const double d8 = I8 - I8_0;
  const double denom = 1.0 + alpha * d8 * d8;
  const double dpsi_dI8 = ks * d8 / (denom * denom);
  const double d2psi_dI8 = ks * (1.0 - 3.0 * alpha * d8 * d8) / std::pow(denom, 3);

  // -------- Gradient wrt 'a' (2x2), Hessian Haa (3x3 Mandel) --------
  const ES::M2d grad_a_mat =
    dpsi_dI1 * dI1_da + dpsi_dI4 * dI4_da + dpsi_dI6 * dI6_da + dpsi_dI8 * dI8_da;

  // -------- Bending (correct Mandel scaling) --------
  // ψ_b = 0.5*(k11 b11^2 + k22 b22^2 + 2 k12 b12^2)
  // δψ_b = k11 b11 δb11 + k22 b22 δb22 + 2 k12 b12 δb12
  ES::M2d grad_b_mat = ES::M2d::Zero();
  grad_b_mat(0, 0) = kappa11 * b_eff(0, 0);
  grad_b_mat(1, 1) = kappa22 * b_eff(1, 1);
  grad_b_mat(0, 1) = grad_b_mat(1, 0) = kappa12 * b_eff(0, 1);  // NOTE: not 2*kappa12

  // Helper to write a-block columns compactly: J_a_col = coeff * dIi_m
  // auto add_Ja = [&](ParamIndex idx, double coeff, const Vec3 &dIi_m) {
  //   const int j = static_cast<int>(idx);
  //   out.J_grad_params.block<3, 1>(0, j).noalias() += coeff * dIi_m;
  // };

  ES::V8d outGrad;
  outGrad.setZero();

  if (i == 0) {
    outGrad.head<4>() = ES::Mp<ES::V4d>(dI1_da.data()) * 0.5 * h;
  }
  else if (i == 1) {
    // k1_4: ∂(∂ψ/∂I4)/∂k1_4 = s4*exp(k2_4 s4^2) * ds4/dI4
    outGrad.head<4>() = ES::Mp<ES::V4d>(dI4_da.data()) * (r4 * E4 * dr4) * h;
  }
  else if (i == 2) {
    // k2_4: ∂(∂ψ/∂I4)/∂k2_4 = k1_4 * s4^3 * exp(k2_4 s4^2) * ds4/dI4
    outGrad.head<4>() = ES::Mp<ES::V4d>(dI4_da.data()) * (k1_4 * r4 * r4 * r4 * E4 * dr4) * h;
  }
  else if (i == 3) {
    // k1_6: ∂(∂ψ/∂I6)/∂k1_6 = s6*exp(k2_6 s6^2) * ds6/dI6
    outGrad.head<4>() = ES::Mp<ES::V4d>(dI6_da.data()) * (r6 * E6 * dr6) * h;
  }
  else if (i == 4) {
    // k2_6: ∂(∂ψ/∂I6)/∂k2_6 = k1_6 * s6^3 * exp(k2_6 s6^2) * ds6/dI6
    outGrad.head<4>() = ES::Mp<ES::V4d>(dI6_da.data()) * (k1_6 * r6 * r6 * r6 * E6 * dr6) * h;
  }
  else if (i == 5) {
    // ks:   ∂(∂ψ/∂I8)/∂ks = d8 / (1 + α d8^2)^2
    outGrad.head<4>() = ES::Mp<ES::V4d>(dI8_da.data()) * (d8 / (denom * denom)) * h;
  }
  else if (i == 6) {
    // alpha: ∂(∂ψ/∂I8)/∂alpha = -2 ks d8^3 / (1 + α d8^2)^3
    outGrad.head<4>() = ES::Mp<ES::V4d>(dI8_da.data()) * (-2.0 * ks * d8 * d8 * d8 / std::pow(denom, 3)) * h;
  }
  else if (i == 7) {
    // In Mandel: grad_b_M(0) = k11*b11; grad_b_M(1) = k22*b22; grad_b_M(2) = sqrt(2)*k12*b12
    outGrad.tail<4>()(0) = b_eff(0, 0) * h * h * h / 12;
  }
  else if (i == 8) {
    // Bending params affect only b-block
    outGrad.tail<4>()(3) = b_eff(1, 1) * h * h * h / 12;
  }
  else if (i == 9) {
    outGrad.tail<4>()(2) = b_eff(0, 1) * h * h * h / 12 + b_eff(1, 0) * h * h * h / 12;
  }
  else if (i == 10) {
    outGrad.head<4>() = ES::Mp<ES::V4d>(dI8_da.data()) * (-d2psi_dI8) * h;
  }
  else if (i == 11) {
    outGrad.head<4>() = ES::Mp<const ES::V4d>(grad_a_mat.data());
    outGrad.tail<4>() = ES::Mp<const ES::V4d>(grad_b_mat.data()) * 3 * h * h / 12;
  }
  else {
    throw std::runtime_error("Invalid material parameter index in ElasticModelThinShellFabric::compute_dP_dparam");
  }

  (ES::Mp<ES::V8d>(ret)) = outGrad;
}
#endif

}  // namespace SolidDeformationModel
}  // namespace pgo