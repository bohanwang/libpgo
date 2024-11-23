/*
author: Bohan Wang
copyright to USC,MIT
*/

#include "elasticModelStableNeoHookeanMaterial.h"

#include "EigenDef.h"

using namespace pgo::SolidDeformationModel;

namespace ES = pgo::EigenSupport;

ElasticModelStableNeoHookeanMaterial::ElasticModelStableNeoHookeanMaterial(double mu_, double lambda_):
  _mu(mu_), _lambda(lambda_ + mu_)
{
  _ratio = _mu / _lambda;
}

ElasticModelStableNeoHookeanMaterial::~ElasticModelStableNeoHookeanMaterial()
{
}

double ElasticModelStableNeoHookeanMaterial::compute_psi(const double * /*param*/, const double FIn[9], const double UIn[9], const double VIn[9], const double SIn[3]) const
{
  ES::M3d F = Eigen::Map<const ES::M3d>(FIn);
  ES::M3d U = Eigen::Map<const ES::M3d>(UIn);
  ES::M3d V = Eigen::Map<const ES::M3d>(VIn);

  const double Ic = F.squaredNorm();
  const double Jminus1 = F.determinant() - 1.0 - _ratio;
  return 0.5 * (_mu * (Ic - 3.0) + _lambda * Jminus1 * Jminus1) - _lambda * _ratio * _ratio * 0.5;
}

static ES::M3d PartialJpartialF(const ES::M3d &F)
{
  ES::M3d pJpF;

  pJpF.col(0) = F.col(1).cross(F.col(2));
  pJpF.col(1) = F.col(2).cross(F.col(0));
  pJpF.col(2) = F.col(0).cross(F.col(1));

  return pJpF;
}

void ElasticModelStableNeoHookeanMaterial::compute_P(const double * /*param*/, const double FIn[9], const double UIn[9], const double VIn[9], const double SIn[3], double Pout[9]) const
{
  ES::M3d F = Eigen::Map<const ES::M3d>(FIn);
  ES::M3d U = Eigen::Map<const ES::M3d>(UIn);
  ES::M3d V = Eigen::Map<const ES::M3d>(VIn);
  ES::V3d S(SIn[0], SIn[1], SIn[2]);

  const ES::M3d pJpF = PartialJpartialF(F);
  const double Jminus1 = F.determinant() - 1.0 - _ratio;

  (ES::Mp<ES::M3d>(Pout)) = _mu * F + _lambda * Jminus1 * pJpF;
}

static void BuildTwistAndFlipEigenvectors(const ES::M3d &U, const ES::M3d &V, ES::M9d &Q)
{
  static const double scale = 1.0 / std::sqrt(2.0);
  const ES::M3d sV = scale * V;

  using M3 = Eigen::Matrix<double, 3, 3, Eigen::ColMajor>;

  M3 A;
  A << sV(0, 2) * U(0, 1), sV(1, 2) * U(0, 1), sV(2, 2) * U(0, 1),
    sV(0, 2) * U(1, 1), sV(1, 2) * U(1, 1), sV(2, 2) * U(1, 1),
    sV(0, 2) * U(2, 1), sV(1, 2) * U(2, 1), sV(2, 2) * U(2, 1);

  M3 B;
  B << sV(0, 1) * U(0, 2), sV(1, 1) * U(0, 2), sV(2, 1) * U(0, 2),
    sV(0, 1) * U(1, 2), sV(1, 1) * U(1, 2), sV(2, 1) * U(1, 2),
    sV(0, 1) * U(2, 2), sV(1, 1) * U(2, 2), sV(2, 1) * U(2, 2);

  M3 C;
  C << sV(0, 2) * U(0, 0), sV(1, 2) * U(0, 0), sV(2, 2) * U(0, 0),
    sV(0, 2) * U(1, 0), sV(1, 2) * U(1, 0), sV(2, 2) * U(1, 0),
    sV(0, 2) * U(2, 0), sV(1, 2) * U(2, 0), sV(2, 2) * U(2, 0);

  M3 D;
  D << sV(0, 0) * U(0, 2), sV(1, 0) * U(0, 2), sV(2, 0) * U(0, 2),
    sV(0, 0) * U(1, 2), sV(1, 0) * U(1, 2), sV(2, 0) * U(1, 2),
    sV(0, 0) * U(2, 2), sV(1, 0) * U(2, 2), sV(2, 0) * U(2, 2);

  M3 E;
  E << sV(0, 1) * U(0, 0), sV(1, 1) * U(0, 0), sV(2, 1) * U(0, 0),
    sV(0, 1) * U(1, 0), sV(1, 1) * U(1, 0), sV(2, 1) * U(1, 0),
    sV(0, 1) * U(2, 0), sV(1, 1) * U(2, 0), sV(2, 1) * U(2, 0);

  M3 F;
  F << sV(0, 0) * U(0, 1), sV(1, 0) * U(0, 1), sV(2, 0) * U(0, 1),
    sV(0, 0) * U(1, 1), sV(1, 0) * U(1, 1), sV(2, 0) * U(1, 1),
    sV(0, 0) * U(2, 1), sV(1, 0) * U(2, 1), sV(2, 0) * U(2, 1);

  // Twist eigenvectors
  Eigen::Map<M3>(Q.data()) = B - A;
  Eigen::Map<M3>(Q.data() + 9) = D - C;
  Eigen::Map<M3>(Q.data() + 18) = F - E;

  // Flip eigenvectors
  Eigen::Map<M3>(Q.data() + 27) = A + B;
  Eigen::Map<M3>(Q.data() + 36) = C + D;
  Eigen::Map<M3>(Q.data() + 45) = E + F;
}

static ES::M9d ProjectHessianWithAnalyticalFormulasNew(const double &mu, const double &lambda, const ES::M3d &F, const ES::M3d &U, const ES::M3d &V, const ES::V3d &S)
{
  ES::V9d eigenvalues;
  ES::M9d eigenvectors;

  const double J = F.determinant();

  // Compute the twist and flip eigenvalues
  {
    // Twist eigenvalues
    eigenvalues.segment<3>(0) = S;
    // Flip eigenvalues
    eigenvalues.segment<3>(3) = -S;
    const double evScale = lambda * (J - 1.0) - mu;
    eigenvalues.segment<6>(0) *= evScale;
    eigenvalues.segment<6>(0).array() += mu;
  }

  // Compute the twist and flip eigenvectors
  BuildTwistAndFlipEigenvectors(U, V, eigenvectors);

  // Compute the remaining three eigenvalues and eigenvectors
  {
    ES::M3d A;
    const double s0s0 = S(0) * S(0);
    const double s1s1 = S(1) * S(1);
    const double s2s2 = S(2) * S(2);
    A(0, 0) = mu + lambda * s1s1 * s2s2;
    A(1, 1) = mu + lambda * s0s0 * s2s2;
    A(2, 2) = mu + lambda * s0s0 * s1s1;
    const double evScale = lambda * (2.0 * J - 1.0) - mu;
    A(0, 1) = evScale * S(2);
    A(1, 0) = A(0, 1);
    A(0, 2) = evScale * S(1);
    A(2, 0) = A(0, 2);
    A(1, 2) = evScale * S(0);
    A(2, 1) = A(1, 2);

    const Eigen::SelfAdjointEigenSolver<ES::M3d> Aeigs(A);
    eigenvalues.segment<3>(6) = Aeigs.eigenvalues();

    Eigen::Map<ES::M3d>(eigenvectors.data() + 54) = U * Aeigs.eigenvectors().col(0).asDiagonal() * V.transpose();
    Eigen::Map<ES::M3d>(eigenvectors.data() + 63) = U * Aeigs.eigenvectors().col(1).asDiagonal() * V.transpose();
    Eigen::Map<ES::M3d>(eigenvectors.data() + 72) = U * Aeigs.eigenvectors().col(2).asDiagonal() * V.transpose();
  }

  // Clamp the eigenvalues
  for (int i = 0; i < 9; i++) {
    if (eigenvalues(i) < 0.0) {
      eigenvalues(i) = 0.0;
    }
  }

  return eigenvectors * eigenvalues.asDiagonal() * eigenvectors.transpose();
}

static ES::V9d PartialJpartialFVec(const ES::M3d &F)
{
  ES::V9d pJpF;
  pJpF.segment<3>(0) = F.col(1).cross(F.col(2));
  pJpF.segment<3>(3) = F.col(2).cross(F.col(0));
  pJpF.segment<3>(6) = F.col(0).cross(F.col(1));
  return pJpF;
}

static ES::M3d CrossProductMatrix(const ES::M3d &F, const int col, const double &scale)
{
  ES::M3d cpm;
  cpm << 0, -scale * F(2, col), scale * F(1, col),
    scale * F(2, col), 0, -scale * F(0, col),
    -scale * F(1, col), scale * F(0, col), 0;
  return cpm;
}

static ES::M9d ComputeFJSecondDerivContribs(const double &lambda, const double &ratio, const ES::M3d &F)
{
  const double scale = lambda * (F.determinant() - 1.0 - ratio);

  const ES::M3d ahat = CrossProductMatrix(F, 0, scale);
  const ES::M3d bhat = CrossProductMatrix(F, 1, scale);
  const ES::M3d chat = CrossProductMatrix(F, 2, scale);

  ES::M9d FJ;
  FJ.block<3, 3>(0, 0).setZero();
  FJ.block<3, 3>(0, 3) = -chat;
  FJ.block<3, 3>(0, 6) = bhat;

  FJ.block<3, 3>(3, 0) = chat;
  FJ.block<3, 3>(3, 3).setZero();
  FJ.block<3, 3>(3, 6) = -ahat;

  FJ.block<3, 3>(6, 0) = -bhat;
  FJ.block<3, 3>(6, 3) = ahat;
  FJ.block<3, 3>(6, 6).setZero();

  return FJ;
}

void ElasticModelStableNeoHookeanMaterial::compute_dPdF(const double * /*param*/, const double FIn[9], const double UIn[9], const double VIn[9], const double SIn[3], double dPdFOut[81]) const
{
  ES::M3d F = Eigen::Map<const ES::M3d>(FIn);
  ES::M3d U = Eigen::Map<const ES::M3d>(UIn);
  ES::M3d V = Eigen::Map<const ES::M3d>(VIn);
  ES::V3d S(SIn[0], SIn[1], SIn[2]);

  ES::M9d dPdF;
  if (enforceSPD_) {
    const ES::V9d pjpf = PartialJpartialFVec(F);
    dPdF = _mu * ES::M9d::Identity() + _lambda * pjpf * pjpf.transpose() + ComputeFJSecondDerivContribs(_lambda, _ratio, F);
  }
  else {
    dPdF = ProjectHessianWithAnalyticalFormulasNew(_mu, _lambda, F, U, V, S);
  }

  (Eigen::Map<Eigen::Matrix<double, 9, 9>>(dPdFOut)) = dPdF;
}

void ElasticModelStableNeoHookeanMaterial::setMaterial(double mu_, double lambda_)
{
  _mu = mu_;
  _lambda = (lambda_ + mu_);
  _ratio = _mu / _lambda;
}
