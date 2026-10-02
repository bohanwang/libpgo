#include "cudaDeformationModelEnergy.h"
#include "cudaHexJacobian.h"

#include "cudaDeviceBuffer.h"
#include "simulationMesh.h"
#include "deformationModelManager.h"
#include <libuipc/math/make_spd.h>

#include <Eigen/Core>
#include <Eigen/LU>
#include <cuda_runtime.h>

#include <cmath>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace pgo::SolidDeformationModel
{
namespace
{
using CudaUtilities::DeviceBuffer;
using CudaUtilities::checkCuda;
using Matrix9 = Eigen::Matrix<double, 9, 9>;
constexpr int kBlockSize = 128;
constexpr double kDetEpsilon = 1e-8;
constexpr double kBoundaryTolerance = 1e-12;

template<int NumVertices>
struct QuadratureData
{
  int vertices[NumVertices];
  // Shape gradients, indexed by local vertex then material coordinate.
  double shapeGradient[3 * NumVertices];
  double volume;
  double mu;
  double lambda;
  double compression;
  double volumePenalty;
  bool stvkVolume;
};

template<int NumVertices>
__device__ void deformationGradient(const QuadratureData<NumVertices> &tet, const double *x, bool displacement, double F[9])
{
  double local[3 * NumVertices];
  for (int a = 0; a < NumVertices; ++a)
    for (int i = 0; i < 3; ++i)
      local[3 * a + i] = x[3 * tet.vertices[a] + i];
  for (int c = 0; c < 3; ++c) {
    for (int i = 0; i < 3; ++i) {
      double value = 0.0;
      for (int a = 1; a < NumVertices; ++a)
        value += (local[3 * a + i] - local[i]) * tet.shapeGradient[3 * a + c];
      F[3 * c + i] = value + (displacement && c == i ? 1.0 : 0.0);
    }
  }
}

__device__ double determinant(const double F[9])
{
  return F[0] * (F[4] * F[8] - F[5] * F[7]) - F[3] * (F[1] * F[8] - F[2] * F[7]) + F[6] * (F[1] * F[5] - F[2] * F[4]);
}

__device__ void cofactor(const double F[9], double C[9])
{
  for (int c = 0; c < 3; ++c) {
    const int a = (c + 1) % 3, b = (c + 2) % 3;
    for (int r = 0; r < 3; ++r) {
      const int i = (r + 1) % 3, j = (r + 2) % 3;
      C[3 * c + r] = F[3 * a + i] * F[3 * b + j] - F[3 * a + j] * F[3 * b + i];
    }
  }
}

__device__ double secondInvariant(const double G[9])
{
  return G[0] * G[4] + G[0] * G[8] + G[4] * G[8] - G[3] * G[1] - G[6] * G[2] - G[7] * G[5];
}

__device__ double determinantMagnitude(const double F[9])
{
  return fabs(F[0] * F[4] * F[8]) + fabs(F[0] * F[5] * F[7]) +
    fabs(F[3] * F[1] * F[8]) + fabs(F[3] * F[2] * F[7]) +
    fabs(F[6] * F[1] * F[5]) + fabs(F[6] * F[2] * F[4]);
}

__device__ double stableNeoEnergy(const double G[9], double mu, double lambda)
{
  double squaredNorm = 0.0;
  for (int i = 0; i < 9; ++i)
    squaredNorm += G[i] * G[i];
  // This expression retains quadratic energies near F=I. Near other rotations,
  // use C-I invariants instead: their quadratic terms avoid cancellation between
  // the large rotation-dependent terms of the displacement-gradient expression.
  double Q[9], normQ = 0.0;
  for (int c = 0; c < 3; ++c) {
    for (int r = 0; r < 3; ++r) {
      double value = G[3 * c + r] + G[3 * r + c];
      for (int k = 0; k < 3; ++k)
        value += G[3 * r + k] * G[3 * c + k];
      Q[3 * c + r] = value;
      normQ += value * value;
    }
  }
  const double nonlinearDet = secondInvariant(G) + determinant(G);
  const double deltaJ = G[0] + G[4] + G[8] + nonlinearDet;
  // ||C-I||_F < 1/2 keeps every eigenvalue of C above 1/2, so evaluating its
  // determinant from these invariants remains well-conditioned. Include
  // moderate contact strain: falling back too early reintroduces cancellation
  // at a few quadrature points even when the total body strain is small.
  if (normQ < 0.25 && deltaJ > -1.0) {
    const double higherInvariants = secondInvariant(Q) + determinant(Q);
    const double deltaDetC = Q[0] + Q[4] + Q[8] + higherInvariants;
    const double stableDeltaJ = deltaDetC / (sqrt(1.0 + deltaDetC) + 1.0);
    return 0.5 * (lambda + mu) * stableDeltaJ * stableDeltaJ - 0.5 * mu * higherInvariants;
  }
  return 0.5 * mu * squaredNorm - mu * nonlinearDet + 0.5 * lambda * deltaJ * deltaJ;
}

__device__ double stableNeoStress(const double G[9], double mu, double lambda,
  double P[9], double cofactorF[9], double *stressMagnitude = nullptr)
{
  double cofactorG[9];
  cofactor(G, cofactorG);
  const double trace = G[0] + G[4] + G[8];
  const double invariant2 = secondInvariant(G), invariant3 = determinant(G);
  const double deltaJ = trace + invariant2 + invariant3;
  const double traceMagnitude = fabs(G[0]) + fabs(G[4]) + fabs(G[8]);
  const double deltaMagnitude = traceMagnitude + determinantMagnitude(G) +
    fabs(G[0] * G[4]) + fabs(G[0] * G[8]) + fabs(G[4] * G[8]) +
    fabs(G[3] * G[1]) + fabs(G[6] * G[2]) + fabs(G[7] * G[5]);
  for (int c = 0; c < 3; ++c) {
    for (int r = 0; r < 3; ++r) {
      const int i = 3 * c + r;
      const double linearCofactor = (r == c ? trace : 0.0) - G[3 * r + c];
      cofactorF[i] = (r == c ? 1.0 : 0.0) + linearCofactor + cofactorG[i];
      P[i] = mu * (G[i] - linearCofactor - cofactorG[i]) + lambda * deltaJ * cofactorF[i];
      if (stressMagnitude) {
        const int a = (c + 1) % 3, b = (c + 2) % 3;
        const int j = (r + 1) % 3, k = (r + 2) % 3;
        const double cofactorMagnitude = fabs(G[3 * a + j] * G[3 * b + k]) + fabs(G[3 * a + k] * G[3 * b + j]);
        const double linearMagnitude = (r == c ? traceMagnitude : 0.0) + fabs(G[3 * r + c]);
        const double magnitudeF = (r == c ? 1.0 : 0.0) + linearMagnitude + cofactorMagnitude;
        stressMagnitude[i] = mu * (fabs(G[i]) + linearMagnitude + cofactorMagnitude) +
          lambda * deltaMagnitude * magnitudeF;
      }
    }
  }
  return deltaJ;
}

// Q = F^T F-I evaluated from G=F-I retains small strains near the rest pose.
__device__ void strainTensor(const double G[9], double Q[9])
{
  for (int c = 0; c < 3; ++c)
    for (int r = 0; r < 3; ++r) {
      double value = G[3 * c + r] + G[3 * r + c];
      for (int k = 0; k < 3; ++k)
        value += G[3 * r + k] * G[3 * c + k];
      Q[3 * c + r] = value;
    }
}

__device__ double determinantOffset(const double G[9])
{
  return G[0] + G[4] + G[8] + secondInvariant(G) + determinant(G);
}

__device__ double stvkVolumeEnergy(const double G[9], double mu, double lambda,
  double compression, double volume)
{
  double Q[9], norm = 0;
  strainTensor(G, Q);
  for (int i = 0; i < 9; ++i)
    norm += Q[i] * Q[i];
  const double trace = Q[0] + Q[4] + Q[8];
  const double deltaJ = determinantOffset(G);
  const double compressionStrain = fmin(deltaJ >= -1 ? deltaJ : -2 - deltaJ, 0.0);
  return lambda * trace * trace / 8 + mu * norm / 4 -
    compression * compressionStrain * compressionStrain * compressionStrain / 2592 + volume * deltaJ * deltaJ / 2;
}

__device__ void stvkVolumeStress(const double G[9], double mu, double lambda,
  double compression, double volume, double P[9], double C[9])
{
  double Q[9], F[9];
  strainTensor(G, Q);
  for (int i = 0; i < 9; ++i)
    F[i] = G[i] + (i % 4 == 0 ? 1.0 : 0.0);
  cofactor(F, C);
  const double trace = Q[0] + Q[4] + Q[8], deltaJ = determinantOffset(G), J = 1 + deltaJ;
  const double compressionStrain = fmin(J >= 0 ? deltaJ : -2 - deltaJ, 0.0);
  const double compressionDerivative = compression > 0 && J == 0 ? nan("") :
                                                                   -compression * compressionStrain * compressionStrain * copysign(1.0, J) / 864;
  for (int c = 0; c < 3; ++c)
    for (int r = 0; r < 3; ++r) {
      double value = 0.5 * lambda * trace * F[3 * c + r];
      for (int k = 0; k < 3; ++k)
        value += mu * F[3 * k + r] * Q[3 * c + k];
      P[3 * c + r] = value + (compressionDerivative + volume * deltaJ) * C[3 * c + r];
    }
}

__device__ double stvkVolumeEnergyChange(const double G[9], const double F[9], const double D[9],
  double mu, double lambda, double compression, double volume, double &magnitude)
{
  double Q[9], dQ[9], C[9], CD[9];
  strainTensor(G, Q);
  cofactor(F, C);
  cofactor(D, CD);
  double qDotDq = 0, normDq = 0, normDqMagnitude = 0, qMagnitude = 0, traceMagnitude = 0, dTraceMagnitude = 0;
  for (int c = 0; c < 3; ++c)
    for (int r = 0; r < 3; ++r) {
      double value = 0, changeScale = 0, qScale = fabs(G[3 * c + r]) + fabs(G[3 * r + c]);
      for (int k = 0; k < 3; ++k) {
        value += F[3 * r + k] * D[3 * c + k] + D[3 * r + k] * F[3 * c + k] + D[3 * r + k] * D[3 * c + k];
        changeScale += fabs(F[3 * r + k] * D[3 * c + k]) + fabs(D[3 * r + k] * F[3 * c + k]) +
          fabs(D[3 * r + k] * D[3 * c + k]);
        qScale += fabs(G[3 * r + k] * G[3 * c + k]);
      }
      const int i = 3 * c + r;
      dQ[i] = value;
      qDotDq += Q[i] * value;
      normDq += value * value;
      normDqMagnitude += changeScale * changeScale;
      qMagnitude += qScale * changeScale;
      if (c == r) {
        traceMagnitude += qScale;
        dTraceMagnitude += changeScale;
      }
    }
  const double trace = Q[0] + Q[4] + Q[8], dTrace = dQ[0] + dQ[4] + dQ[8];
  const double strainChange = lambda * trace * dTrace / 4 + lambda * dTrace * dTrace / 8 +
    mu * qDotDq / 2 + mu * normDq / 4;
  double deltaJ = determinant(D), deltaMagnitude = determinantMagnitude(D);
  for (int i = 0; i < 9; ++i) {
    deltaJ += C[i] * D[i] + CD[i] * F[i];
    deltaMagnitude += fabs(C[i] * D[i]) + fabs(CD[i] * F[i]);
  }
  const double offset = determinantOffset(G), J0 = 1 + offset, J1 = J0 + deltaJ;
  const double absoluteDelta = (J0 < 0) == (J1 < 0) ? copysign(1.0, J0) * deltaJ : fabs(J1) - fabs(J0);
  const double compressionOffset = J0 >= 0 ? offset : -2 - offset;
  const double c0 = fmin(compressionOffset, 0.0), c1 = fmin(compressionOffset + absoluteDelta, 0.0);
  const double dc = c0 < 0 && c1 < 0 ? absoluteDelta : c1 - c0;
  const double compressionChange = -compression * dc * (c0 * c0 + c0 * c1 + c1 * c1) / 2592;
  const double volumeChange = volume * (offset * deltaJ + 0.5 * deltaJ * deltaJ);
  const double offsetMagnitude = fabs(G[0]) + fabs(G[4]) + fabs(G[8]) + determinantMagnitude(G) +
    fabs(G[0] * G[4]) + fabs(G[0] * G[8]) + fabs(G[4] * G[8]) +
    fabs(G[3] * G[1]) + fabs(G[6] * G[2]) + fabs(G[7] * G[5]);
  // Retain the uncancelled arithmetic scale: skew directions can cancel the
  // linear strain increment, and a rotated base can cancel determinant terms.
  magnitude = fabs(lambda) * (traceMagnitude * dTraceMagnitude / 4 + dTraceMagnitude * dTraceMagnitude / 8) +
    mu * (qMagnitude / 2 + normDqMagnitude / 4) +
    compression * deltaMagnitude * (c0 * c0 + fabs(c0 * c1) + c1 * c1) / 2592 +
    fabs(volume) * (offsetMagnitude * deltaMagnitude + 0.5 * deltaMagnitude * deltaMagnitude);
  return strainChange + compressionChange + volumeChange;
}

template<int NumVertices>
__global__ void energyKernel(const QuadratureData<NumVertices> *tets, int count, const double *x, double *sum)
{
  const int e = blockIdx.x * blockDim.x + threadIdx.x;
  double energy = 0.0;
  if (e < count) {
    const auto &tet = tets[e];
    double G[9];
    deformationGradient(tet, x, false, G);
    energy = tet.volume * (tet.stvkVolume ? stvkVolumeEnergy(G, tet.mu, tet.lambda, tet.compression, tet.volumePenalty) : stableNeoEnergy(G, tet.mu, tet.lambda));
  }
  __shared__ double values[kBlockSize];
  values[threadIdx.x] = energy;
  __syncthreads();
  for (int stride = blockDim.x / 2; stride; stride /= 2) {
    if (threadIdx.x < stride)
      values[threadIdx.x] += values[threadIdx.x + stride];
    __syncthreads();
  }
  if (threadIdx.x == 0)
    atomicAdd(sum, values[0]);
}

template<int NumVertices>
__global__ void energyChangeKernel(const QuadratureData<NumVertices> *points, int count,
  const double *base, const double *trial, double *sum)
{
  const int e = blockIdx.x * blockDim.x + threadIdx.x;
  double change = 0.0, magnitude = 0.0;
  if (e < count) {
    const auto &point = points[e];
    double G[9], F[9], D[9], P[9], CF[9], CD[9], stressMagnitude[9];
    deformationGradient(point, base, false, G);
    // Form the represented trial displacement first; subtracting two computed
    // deformation gradients would lose the small line-search direction.
    double localDifference[3 * NumVertices];
    for (int a = 0; a < NumVertices; ++a)
      for (int i = 0; i < 3; ++i) {
        const int global = 3 * point.vertices[a] + i;
        localDifference[3 * a + i] = trial[global] - base[global];
      }
    for (int c = 0; c < 3; ++c) {
      for (int i = 0; i < 3; ++i) {
        double value = 0.0;
        for (int a = 1; a < NumVertices; ++a)
          value += (localDifference[3 * a + i] - localDifference[i]) * point.shapeGradient[3 * a + c];
        D[3 * c + i] = value;
        F[3 * c + i] = G[3 * c + i] + (c == i ? 1.0 : 0.0);
      }
    }
    if (point.stvkVolume) {
      change = point.volume * stvkVolumeEnergyChange(G, F, D, point.mu, point.lambda, point.compression, point.volumePenalty, magnitude);
      magnitude *= point.volume;
    }
    else {
      const double deltaJ0 = stableNeoStress(G, point.mu, point.lambda, P, CF, stressMagnitude);
      cofactor(D, CD);
      double linear = 0.0, j1 = 0.0, j2 = 0.0, squaredNorm = 0.0;
      double j1Magnitude = 0.0, j2Magnitude = 0.0;
      for (int i = 0; i < 9; ++i) {
        linear += P[i] * D[i];
        j1 += CF[i] * D[i];
        j2 += CD[i] * F[i];
        squaredNorm += D[i] * D[i];
        magnitude += stressMagnitude[i] * fabs(D[i]);
        j1Magnitude += fabs(CF[i] * D[i]);
        j2Magnitude += fabs(CD[i] * F[i]);
      }
      const double j3 = determinant(D), deltaJ = j1 + j2 + j3;
      const double factor = point.lambda * deltaJ0 - point.mu;
      const double quadratic = 0.5 * point.mu * squaredNorm;
      const double higher = factor * (j2 + j3);
      const double volumeChange = 0.5 * point.lambda * deltaJ * deltaJ;
      change = point.volume * (linear + quadratic + higher + volumeChange);
      const double j3Magnitude = determinantMagnitude(D);
      const double deltaMagnitude = j1Magnitude + j2Magnitude + j3Magnitude;
      magnitude = point.volume * (magnitude + quadratic + fabs(factor) * (j2Magnitude + j3Magnitude) + 0.5 * point.lambda * deltaMagnitude * deltaMagnitude);
    }
  }
  __shared__ double changes[kBlockSize], magnitudes[kBlockSize];
  changes[threadIdx.x] = change;
  magnitudes[threadIdx.x] = magnitude;
  __syncthreads();
  for (int stride = blockDim.x / 2; stride; stride /= 2) {
    if (threadIdx.x < stride) {
      changes[threadIdx.x] += changes[threadIdx.x + stride];
      magnitudes[threadIdx.x] += magnitudes[threadIdx.x + stride];
    }
    __syncthreads();
  }
  if (threadIdx.x == 0) {
    atomicAdd(sum, changes[0]);
    atomicAdd(sum + 1, magnitudes[0]);
  }
}

template<int NumVertices>
__global__ void gradientKernel(const QuadratureData<NumVertices> *tets, int count, const double *x, double *gradient)
{
  const int e = blockIdx.x * blockDim.x + threadIdx.x;
  if (e >= count)
    return;
  const auto &tet = tets[e];
  double G[9], C[9], P[9];
  deformationGradient(tet, x, false, G);
  if (tet.stvkVolume)
    stvkVolumeStress(G, tet.mu, tet.lambda, tet.compression, tet.volumePenalty, P, C);
  else
    stableNeoStress(G, tet.mu, tet.lambda, P, C);
  for (int a = 0; a < NumVertices; ++a) {
    for (int i = 0; i < 3; ++i) {
      double value = 0.0;
      for (int c = 0; c < 3; ++c)
        value += P[3 * c + i] * tet.shapeGradient[3 * a + c];
      atomicAdd(gradient + 3 * tet.vertices[a] + i, tet.volume * value);
    }
  }
}

__device__ int permutationSign(int i, int j)
{
  return (j == (i + 1) % 3) ? 1 : -1;
}

template<int NumQuadraturePoints, int NumVertices>
__global__ void hessianKernel(const QuadratureData<NumVertices> *tets, int count, const double *x,
  bool projectMaterialTangent, double *hessians)
{
  const int e = blockIdx.x * blockDim.x + threadIdx.x;
  if (e >= count)
    return;
  const auto &tet = tets[e];
  double F[9], C[9];
  deformationGradient(tet, x, true, F);
  cofactor(F, C);
  Matrix9 tangent;
  double volumeScale = 0.0;
  if (tet.stvkVolume) {
    double G[9], Q[9], S[9], FFt[9];
    deformationGradient(tet, x, false, G);
    strainTensor(G, Q);
    const double trace = Q[0] + Q[4] + Q[8];
    for (int c = 0; c < 3; ++c)
      for (int r = 0; r < 3; ++r) {
        S[3 * c + r] = tet.mu * Q[3 * c + r] + (c == r ? 0.5 * tet.lambda * trace : 0.0);
        double value = 0;
        for (int k = 0; k < 3; ++k)
          value += F[3 * k + r] * F[3 * k + c];
        FFt[3 * c + r] = value;
      }
    const double offset = determinantOffset(G), J = 1 + offset;
    const double compressionStrain = fmin(J >= 0 ? offset : -2 - offset, 0.0);
    const double first = tet.compression > 0 && J == 0 ? nan("") :
                                                         -tet.compression * compressionStrain * compressionStrain * copysign(1.0, J) / 864;
    const double second = -tet.compression * compressionStrain / 432;
    for (int b = 0; b < 9; ++b)
      for (int a = 0; a < 9; ++a) {
        const int i = a % 3, j = a / 3, k = b % 3, l = b / 3;
        double value = (i == k ? S[3 * j + l] : 0) + tet.lambda * F[a] * F[b] +
          tet.mu * (F[3 * l + i] * F[3 * j + k] + (j == l ? FFt[3 * k + i] : 0)) + second * C[a] * C[b];
        if (i != k && j != l)
          value += first * permutationSign(i, k) * permutationSign(j, l) * F[3 * (3 - j - l) + (3 - i - k)];
        tangent(a, b) = value;
      }
    volumeScale = tet.volumePenalty * offset;
  }
  else {
    const double scale = tet.lambda * (determinant(F) - 1.0) - tet.mu;
    for (int b = 0; b < 9; ++b)
      for (int a = 0; a < 9; ++a) {
        double value = tet.lambda * C[a] * C[b] + (a == b ? tet.mu : 0.0);
        const int i = a % 3, j = a / 3, k = b % 3, l = b / 3;
        if (i != k && j != l)
          value += scale * permutationSign(i, k) * permutationSign(j, l) * F[3 * (3 - j - l) + (3 - i - k)];
        tangent(a, b) = value;
      }
  }
  // CPU StVK-vol projects its invariant-material child, then adds the raw
  // determinant-penalty child. Projecting their sum changes that material.
  if (projectMaterialTangent)
    uipc::backend::cuda::make_spd<9>(tangent);
  if (tet.stvkVolume) {
    for (int b = 0; b < 9; ++b)
      for (int a = 0; a < 9; ++a) {
        double value = tet.volumePenalty * C[a] * C[b];
        const int i = a % 3, j = a / 3, k = b % 3, l = b / 3;
        if (i != k && j != l)
          value += volumeScale * permutationSign(i, k) * permutationSign(j, l) * F[3 * (3 - j - l) + (3 - i - k)];
        tangent(a, b) += value;
      }
  }
  constexpr int NumDOFs = 3 * NumVertices;
  double *hessian = hessians + static_cast<std::size_t>(e / NumQuadraturePoints) * NumDOFs * NumDOFs;
  for (int b = 0; b < NumDOFs; ++b) {
    for (int a = 0; a < NumDOFs; ++a) {
      double value = 0.0;
      for (int d = 0; d < 3; ++d)
        for (int c = 0; c < 3; ++c)
          value += tet.shapeGradient[3 * (a / 3) + c] * tangent(3 * c + a % 3, 3 * d + b % 3) * tet.shapeGradient[3 * (b / 3) + d];
      if constexpr (NumQuadraturePoints == 1)
        hessian[NumDOFs * b + a] = tet.volume * value;
      else
        atomicAdd(hessian + NumDOFs * b + a, tet.volume * value);
    }
  }
}

template<int NumQuadraturePoints, int NumVertices>
__global__ void hessianVectorKernel(const QuadratureData<NumVertices> *tets, int count, const double *hessians,
  const double *v, double *out)
{
  const int e = blockIdx.x * blockDim.x + threadIdx.x;
  if (e >= count)
    return;
  const auto &tet = tets[e * NumQuadraturePoints];
  constexpr int NumDOFs = 3 * NumVertices;
  const double *hessian = hessians + static_cast<std::size_t>(e) * NumDOFs * NumDOFs;
  double local[NumDOFs];
  for (int a = 0; a < NumDOFs; ++a)
    local[a] = v[3 * tet.vertices[a / 3] + a % 3];
  for (int a = 0; a < NumDOFs; ++a) {
    double value = 0.0;
    for (int b = 0; b < NumDOFs; ++b)
      value += hessian[NumDOFs * b + a] * local[b];
    atomicAdd(out + 3 * tet.vertices[a / 3] + a % 3, value);
  }
}

template<int NumQuadraturePoints, int NumVertices>
__global__ void diagonalKernel(const QuadratureData<NumVertices> *tets, int count, const double *hessians, double *diagonal)
{
  const int e = blockIdx.x * blockDim.x + threadIdx.x;
  if (e >= count)
    return;
  constexpr int NumDOFs = 3 * NumVertices;
  const double *hessian = hessians + static_cast<std::size_t>(e) * NumDOFs * NumDOFs;
  for (int a = 0; a < NumDOFs; ++a)
    atomicAdd(diagonal + 3 * tets[e * NumQuadraturePoints].vertices[a / 3] + a % 3, hessian[(NumDOFs + 1) * a]);
}

struct Polynomial
{
  double c[4];
  __device__ double eval(double alpha) const { return ((c[3] * alpha + c[2]) * alpha + c[1]) * alpha + c[0]; }
};

__device__ double interiorClamp(double alpha)
{
  return alpha >= 1.0 ? 1.0 : fmin(alpha, fmax(1e-12, 0.99 * alpha));
}

__device__ double feasibleStep(const double F[9], const double dF[9])
{
  Polynomial poly{};
  // Multilinearity expands det(F + alpha*dF) with no sampling cancellation.
  for (int mask = 0; mask < 8; ++mask) {
    double mixed[9];
    int degree = 0;
    for (int col = 0; col < 3; ++col) {
      const bool direction = mask & (1 << col);
      degree += direction;
      for (int row = 0; row < 3; ++row)
        mixed[3 * col + row] = direction ? dF[3 * col + row] : F[3 * col + row];
    }
    poly.c[degree] += determinant(mixed);
  }
  if (!isfinite(poly.c[0]))
    return 0.0;
  if (poly.c[0] <= kDetEpsilon)
    return 1e-12;  // Same recovery bound as CPU material max step.
  poly.c[0] -= kDetEpsilon;
  for (int i = 1; i < 4; ++i)
    if (!isfinite(poly.c[i]))
      return 0.0;

  double cuts[4] = { 0.0, 1.0, 0.0, 0.0 };
  int count = 2;
  const double qa = 3.0 * poly.c[3], qb = 2.0 * poly.c[2], qc = poly.c[1];
  if (!isfinite(qa) || !isfinite(qb))
    return 0.0;
  double roots[2];
  int numRoots = 0;
  if (qa == 0.0) {
    if (qb != 0.0)
      roots[numRoots++] = -qc / qb;
  }
  else {
    const double disc = qb * qb - 4.0 * qa * qc;
    if (!isfinite(disc))
      return 0.0;
    if (disc >= 0.0) {
      const double root = sqrt(fmax(0.0, disc));
      roots[numRoots++] = (-qb - root) / (2.0 * qa);
      roots[numRoots++] = (-qb + root) / (2.0 * qa);
    }
  }
  for (int i = 0; i < numRoots; ++i)
    if (roots[i] > 0.0 && roots[i] < 1.0)
      cuts[count++] = roots[i];
  for (int i = 1; i < count; ++i) {
    const double value = cuts[i];
    int j = i;
    while (j > 0 && cuts[j - 1] > value) {
      cuts[j] = cuts[j - 1];
      --j;
    }
    cuts[j] = value;
  }
  int uniqueCount = 1;
  for (int i = 1; i < count; ++i)
    if (fabs(cuts[i] - cuts[uniqueCount - 1]) > kBoundaryTolerance)
      cuts[uniqueCount++] = cuts[i];
  for (int i = 1; i < uniqueCount; ++i) {
    double left = cuts[i - 1], right = cuts[i];
    const double value = poly.eval(right);
    if (value > kBoundaryTolerance)
      continue;
    if (fabs(value) <= kBoundaryTolerance)
      return interiorClamp(right);
    for (int iteration = 0; iteration < 64; ++iteration) {
      const double middle = 0.5 * (left + right);
      if (poly.eval(middle) > 0.0)
        left = middle;
      else
        right = middle;
    }
    return interiorClamp(left);
  }
  return 1.0;
}

__global__ void directionNonzeroKernel(int count, const double *dx, double *nonzero)
{
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  __shared__ double values[kBlockSize];
  // Squared norms can underflow for nonzero represented motion. This query
  // only distinguishes a stationary state from a direction, so count entries.
  values[threadIdx.x] = i < count && dx[i] != 0.0 ? 1.0 : 0.0;
  __syncthreads();
  for (int stride = blockDim.x / 2; stride; stride /= 2) {
    if (threadIdx.x < stride)
      values[threadIdx.x] += values[threadIdx.x + stride];
    __syncthreads();
  }
  if (threadIdx.x == 0)
    atomicAdd(nonzero, values[0]);
}

template<int NumVertices>
__global__ void maxStepKernel(const QuadratureData<NumVertices> *tets, int count, const double *x, const double *dx,
  const double *directionNonzero, double *alpha)
{
  const int e = blockIdx.x * blockDim.x + threadIdx.x;
  if (e >= count)
    return;
  double F[9], dF[9];
  deformationGradient(tets[e], x, true, F);
  deformationGradient(tets[e], dx, false, dF);
  // A zero direction is also the solver's initial-state feasibility query.
  // It must not report an inverted stationary state as a valid equilibrium.
  const double initialDeterminant = determinant(F);
  const double bound = *directionNonzero == 0.0 ?
    (isfinite(initialDeterminant) && initialDeterminant > kDetEpsilon ? 1.0 : 0.0) :
    feasibleStep(F, dF);
  // All bounds are nonnegative, so the IEEE representation preserves order.
  atomicMin(reinterpret_cast<unsigned long long *>(alpha), static_cast<unsigned long long>(__double_as_longlong(bound)));
}

__global__ void hexFastMaxStepKernel(const HexJacobian::Element *elements, int count,
  const double *x, const double *dx, double determinantFloor, const double *limit,
  unsigned char *needsSubdivision, double *alpha)
{
  const int e = blockIdx.x * blockDim.x + threadIdx.x;
  if (e >= count)
    return;
  HexJacobian::Motion motion;
  needsSubdivision[e] = 0;
  if (!HexJacobian::loadMotion(elements[e], x, dx, motion)) {
    atomicMin(reinterpret_cast<unsigned long long *>(alpha), 0ULL);
    return;
  }
  needsSubdivision[e] = !HexJacobian::certifyEdgeHull(elements[e], motion, *limit, determinantFloor);
}

__global__ void hexMaxStepKernel(const HexJacobian::Element *elements, int count,
  const double *x, const double *dx, double determinantFloor, const double *inputLimit,
  const unsigned char *needsSubdivision, double *alpha)
{
  const int e = blockIdx.x * blockDim.x + threadIdx.x;
  if (e >= count || !needsSubdivision[e])
    return;
  const double limit = *inputLimit;
  if (limit == 0)
    return;
  HexJacobian::Motion motion;
  double bound = 0;
  if (HexJacobian::loadMotion(elements[e], x, dx, motion)) {
    if (HexJacobian::certify(elements[e], motion, limit, determinantFloor))
      bound = limit;
    else if (HexJacobian::certify(elements[e], motion, 0, determinantFloor)) {
      // Only a complete original-prefix certificate moves the lower bound.
      // Failed or exhausted queries are not evidence of positivity. Returning
      // a smaller prefix of a certified interval preserves that certificate.
      double lower = 0, upper = limit;
      for (int iteration = 0; iteration < 16; ++iteration) {
        const double middle = lower + 0.5 * (upper - lower);
        if (HexJacobian::certify(elements[e], motion, middle, determinantFloor))
          lower = middle;
        else
          upper = middle;
      }
      bound = 0.99 * lower;
    }
  }
  atomicMin(reinterpret_cast<unsigned long long *>(alpha), static_cast<unsigned long long>(__double_as_longlong(bound)));
}

std::vector<HexJacobian::Element> prepareHexElements(const SimulationMesh &mesh)
{
  std::vector<HexJacobian::Element> elements(mesh.getNumElements());
  const int bitToVertex[8]{ 0, 1, 3, 2, 4, 5, 7, 6 };
  for (int e = 0; e < mesh.getNumElements(); ++e) {
    auto &element = elements[e];
    for (int v = 0; v < 8; ++v) {
      element.vertices[v] = mesh.getVertexIndex(e, v);
      mesh.getVertex(element.vertices[v], element.rest[v]);
    }
    Eigen::Matrix3d center = Eigen::Matrix3d::Zero();
    for (int axis = 0; axis < 3; ++axis)
      for (int bits = 0; bits < 8; ++bits) {
        if (bits & (1 << axis))
          continue;
        for (int row = 0; row < 3; ++row)
          center(row, axis) += 0.25 * (element.rest[bitToVertex[bits | (1 << axis)]][row] - element.rest[bitToVertex[bits]][row]);
      }
    const double determinant = center.determinant();
    if (!std::isfinite(determinant) || determinant == 0)
      throw std::invalid_argument("CUDA cubic FEM requires a consistently oriented nondegenerate rest cell.");
    element.orientation = determinant < 0 ? -1 : 1;
  }
  return elements;
}

template<int NumVertices>
std::vector<QuadratureData<NumVertices>> prepareQuadrature(const SimulationMesh &mesh, DeformationModelElasticMaterial materialType)
{
  constexpr int NumQuadraturePoints = NumVertices == 4 ? 1 : 8;
  if (mesh.getNumElements() > std::numeric_limits<int>::max() / NumQuadraturePoints)
    throw std::length_error("CUDA deformation energy quadrature count exceeds integer range.");
  std::vector<QuadratureData<NumVertices>> points(static_cast<std::size_t>(mesh.getNumElements()) * NumQuadraturePoints);
  const int corner[8][3] = { { 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 }, { 0, 1, 0 },
    { 0, 0, 1 }, { 1, 0, 1 }, { 1, 1, 1 }, { 0, 1, 1 } };
  const double offset = 0.5 / std::sqrt(3.0);
  const double gauss[2] = { 0.5 - offset, 0.5 + offset };
  for (int e = 0; e < mesh.getNumElements(); ++e) {
    Eigen::Matrix<double, 3, NumVertices> rest;
    int vertices[NumVertices];
    for (int a = 0; a < NumVertices; ++a) {
      vertices[a] = mesh.getVertexIndex(e, a);
      if (vertices[a] < 0 || vertices[a] >= mesh.getNumVertices())
        throw std::invalid_argument("CUDA deformation energy encountered an invalid element vertex index.");
      mesh.getVertex(vertices[a], rest.col(a).data());
    }
    if (!rest.allFinite())
      throw std::invalid_argument("CUDA deformation energy requires finite rest positions.");
    const Eigen::Vector3d origin = rest.col(0);
    rest.colwise() -= origin;
    const auto *material = mesh.getElementNumMaterials(e) > 0 ?
      dynamic_cast<const SimulationMeshENuMaterial *>(mesh.getElementMaterial(e, 0)) :
      nullptr;
    if (!material || !std::isfinite(material->getE()) || material->getE() <= 0.0 ||
      !std::isfinite(material->getNu()) || material->getNu() <= -1.0 || material->getNu() >= 0.5)
      throw std::invalid_argument("CUDA volumetric FEM requires E > 0 and -1 < nu < 0.5.");
    for (int q = 0; q < NumQuadraturePoints; ++q) {
      auto &point = points[e * NumQuadraturePoints + q];
      Eigen::Matrix<double, 3, NumVertices> dN;
      if constexpr (NumVertices == 4) {
        dN.col(0).setConstant(-1.0);
        dN.template rightCols<3>().setIdentity();
      }
      else {
        const double coordinate[3] = { gauss[q / 4], gauss[(q / 2) % 2], gauss[q % 2] };
        for (int a = 0; a < NumVertices; ++a) {
          for (int axis = 0; axis < 3; ++axis) {
            double value = corner[a][axis] ? 1.0 : -1.0;
            for (int other = 0; other < 3; ++other)
              if (other != axis)
                value *= corner[a][other] ? coordinate[other] : 1.0 - coordinate[other];
            dN(axis, a) = value;
          }
        }
      }
      const Eigen::Matrix3d Dm = rest * dN.transpose();
      const Eigen::FullPivLU<Eigen::Matrix3d> decomposition(Dm);
      if (!Dm.allFinite() || !decomposition.isInvertible())
        throw std::invalid_argument("CUDA deformation energy requires nondegenerate rest element Jacobians.");
      const Eigen::Matrix3d inverse = decomposition.inverse();
      point.volume = std::abs(Dm.determinant()) / (NumVertices == 4 ? 6.0 : 8.0);
      if (!(point.volume > 0.0) || !std::isfinite(point.volume) || !inverse.allFinite())
        throw std::invalid_argument("CUDA deformation energy encountered unrepresentable rest geometry.");
      Eigen::Matrix<double, 3, NumVertices> shape = inverse.transpose() * dN;
      // Enforce partition of unity in the same local coordinate frame used by
      // deformationGradient, including on distorted trilinear elements.
      shape.col(0) = -shape.template rightCols<NumVertices - 1>().rowwise().sum();
      for (int a = 0; a < NumVertices; ++a) {
        point.vertices[a] = vertices[a];
        for (int axis = 0; axis < 3; ++axis)
          point.shapeGradient[3 * a + axis] = shape(axis, a);
      }
      point.mu = material->getMuLame();
      point.stvkVolume = materialType == DeformationModelElasticMaterial::STVK_VOL;
      point.lambda = material->getLambdaLame() + (point.stvkVolume ? 0 : point.mu);
      point.volumePenalty = point.stvkVolume ? material->getCompressionRatio() : 0;
      point.compression = point.volumePenalty > 0 ? point.volumePenalty * material->getE() / (1 - 2 * material->getNu()) : 0;
      if (!std::isfinite(point.mu) || !std::isfinite(point.lambda) ||
        !std::isfinite(point.volumePenalty) || !std::isfinite(point.compression))
        throw std::invalid_argument("CUDA volumetric FEM parameters must be finite and representable.");
    }
  }
  return points;
}
}  // namespace

class CudaDeformationModelEnergy::Impl
{
public:
  int numDOFs = 0;
  int numElements = 0;
  bool cubic = false;
  bool enableMaterialMaxStep = true;
  CubicValidity cubicValidity = CubicValidity::Gauss;
  bool hessianPrepared = false;
  DeviceBuffer<QuadratureData<4>> tetQuadrature;
  DeviceBuffer<QuadratureData<8>> cubicQuadrature;
  DeviceBuffer<HexJacobian::Element> hexElements;
  mutable DeviceBuffer<unsigned char> hexNeedsSubdivision;
  mutable DeviceBuffer<double> hexLimit{ 1 };
  DeviceBuffer<double> hessians;
  mutable DeviceBuffer<double> scalar{ 1 };
  mutable DeviceBuffer<double> directionNonzero{ 1 };
  DeviceBuffer<double> changeScalars{ 2 };

  template<class Function>
  void visitQuadrature(Function &&function) const
  {
    if (cubic)
      function(cubicQuadrature, std::integral_constant<int, 8>{});
    else
      function(tetQuadrature, std::integral_constant<int, 1>{});
  }
};

CudaDeformationModelEnergy::CudaDeformationModelEnergy(const SimulationMesh &mesh, bool enableMaterialMaxStep,
  CubicValidity cubicValidity):
  CudaDeformationModelEnergy(mesh, DeformationModelElasticMaterial::STABLE_NEO, enableMaterialMaxStep, cubicValidity)
{
}

CudaDeformationModelEnergy::CudaDeformationModelEnergy(const SimulationMesh &mesh,
  DeformationModelElasticMaterial material, bool enableMaterialMaxStep, CubicValidity cubicValidity):
  data(std::make_unique<Impl>())
{
  if (material != DeformationModelElasticMaterial::STABLE_NEO && material != DeformationModelElasticMaterial::STVK_VOL)
    throw std::invalid_argument("CUDA volumetric FEM supports only Stable Neo-Hookean and StVK-vol materials.");
  if (cubicValidity != CubicValidity::Gauss && cubicValidity != CubicValidity::WholeCell)
    throw std::invalid_argument("CUDA cubic validity must select Gauss sampling or whole-cell certification.");
  const bool tetrahedral = mesh.getElementType() == SimulationMeshType::TET && mesh.getNumElementVertices() == 4;
  data->cubic = mesh.getElementType() == SimulationMeshType::CUBIC && mesh.getNumElementVertices() == 8;
  if (!tetrahedral && !data->cubic)
    throw std::invalid_argument("CUDA deformation energy requires tetrahedral or cubic elements.");
  if (mesh.getNumVertices() > std::numeric_limits<int>::max() / 3)
    throw std::length_error("CUDA deformation energy DOF count exceeds integer range.");
  data->numDOFs = 3 * mesh.getNumVertices();
  data->numElements = mesh.getNumElements();
  data->enableMaterialMaxStep = enableMaterialMaxStep;
  data->cubicValidity = cubicValidity;
  if (data->cubic) {
    const auto points = prepareQuadrature<8>(mesh, material);
    data->cubicQuadrature.upload(points.data(), points.size());
    if (cubicValidity == CubicValidity::WholeCell) {
      const auto elements = prepareHexElements(mesh);
      data->hexElements.upload(elements.data(), elements.size());
      data->hexNeedsSubdivision.resize(elements.size());
    }
    if (cubicValidity == CubicValidity::WholeCell && data->numElements) {
      const double one = 1;
      data->scalar.upload(&one, 1);
      data->hexLimit.upload(&one, 1);
      hexFastMaxStepKernel<<<(data->numElements + kBlockSize - 1) / kBlockSize, kBlockSize>>>(
        data->hexElements.data(), data->numElements, nullptr, nullptr, 0, data->hexLimit.data(), data->hexNeedsSubdivision.data(), data->scalar.data());
      hexMaxStepKernel<<<(data->numElements + 31) / 32, 32>>>(
        data->hexElements.data(), data->numElements, nullptr, nullptr, 0, data->hexLimit.data(), data->hexNeedsSubdivision.data(), data->scalar.data());
      checkCuda(cudaGetLastError(), "CUDA cubic rest-cell certificate");
      double certified;
      data->scalar.download(&certified, 1);
      if (certified != 1)
        throw std::invalid_argument("CUDA cubic FEM could not certify positive Jacobians throughout every rest cell.");
    }
  }
  else {
    const auto points = prepareQuadrature<4>(mesh, material);
    data->tetQuadrature.upload(points.data(), points.size());
  }
  const std::size_t elementDOFs = data->cubic ? 24 : 12;
  data->hessians.resize(static_cast<std::size_t>(data->numElements) * elementDOFs * elementDOFs);
}

CudaDeformationModelEnergy::~CudaDeformationModelEnergy() = default;

int CudaDeformationModelEnergy::getNumDOFs() const
{
  return data->numDOFs;
}

double CudaDeformationModelEnergy::func(const double *deviceX)
{
  data->scalar.zero();
  if (data->numElements) {
    data->visitQuadrature([&](const auto &points, auto) {
      const int count = static_cast<int>(points.size());
      energyKernel<<<(count + kBlockSize - 1) / kBlockSize, kBlockSize>>>(
        points.data(), count, deviceX, data->scalar.data());
    });
    checkCuda(cudaGetLastError(), "CUDA FEM energy kernel");
  }
  double value;
  data->scalar.download(&value, 1);
  return value;
}

CudaDeformationModelEnergy::EnergyChange CudaDeformationModelEnergy::energyChange(
  const double *deviceBase, const double *deviceTrial)
{
  data->changeScalars.zero();
  if (data->numElements) {
    data->visitQuadrature([&](const auto &points, auto) {
      const int count = static_cast<int>(points.size());
      energyChangeKernel<<<(count + kBlockSize - 1) / kBlockSize, kBlockSize>>>(
        points.data(), count, deviceBase, deviceTrial, data->changeScalars.data());
    });
    checkCuda(cudaGetLastError(), "CUDA FEM paired energy-change kernel");
  }
  double values[2];
  data->changeScalars.download(values, 2);
  return { values[0], values[1] };
}

void CudaDeformationModelEnergy::addGradient(const double *deviceX, double *deviceGradient)
{
  if (!data->numElements)
    return;
  data->visitQuadrature([&](const auto &points, auto) {
    const int count = static_cast<int>(points.size());
    gradientKernel<<<(count + kBlockSize - 1) / kBlockSize, kBlockSize>>>(
      points.data(), count, deviceX, deviceGradient);
  });
  checkCuda(cudaGetLastError(), "CUDA FEM gradient kernel");
}

void CudaDeformationModelEnergy::prepareHessian(const double *deviceX)
{
  prepareHessian(deviceX, true);
}

void CudaDeformationModelEnergy::prepareExactHessian(const double *deviceX)
{
  prepareHessian(deviceX, false);
}

void CudaDeformationModelEnergy::prepareHessian(const double *deviceX, bool projectMaterialTangent)
{
  if (data->numElements) {
    if (data->cubic)
      data->hessians.zero();
    data->visitQuadrature([&](const auto &points, auto quadratureCount) {
      const int count = static_cast<int>(points.size());
      hessianKernel<decltype(quadratureCount)::value><<<(count + kBlockSize - 1) / kBlockSize, kBlockSize>>>(
        points.data(), count, deviceX, projectMaterialTangent, data->hessians.data());
    });
    checkCuda(cudaGetLastError(), "CUDA FEM Hessian kernel");
  }
  data->hessianPrepared = true;
}

void CudaDeformationModelEnergy::addHessianVector(const double *deviceV, double *deviceOut) const
{
  if (!data->hessianPrepared)
    throw std::logic_error("CUDA FEM Hessian must be prepared before applying it.");
  if (!data->numElements)
    return;
  data->visitQuadrature([&](const auto &points, auto quadratureCount) {
    hessianVectorKernel<decltype(quadratureCount)::value><<<(data->numElements + kBlockSize - 1) / kBlockSize, kBlockSize>>>(
      points.data(), data->numElements, data->hessians.data(), deviceV, deviceOut);
  });
  checkCuda(cudaGetLastError(), "CUDA FEM Hessian-vector kernel");
}

void CudaDeformationModelEnergy::addDiagonal(double *deviceDiagonal) const
{
  if (!data->hessianPrepared)
    throw std::logic_error("CUDA FEM Hessian must be prepared before querying its diagonal.");
  if (!data->numElements)
    return;
  data->visitQuadrature([&](const auto &points, auto quadratureCount) {
    diagonalKernel<decltype(quadratureCount)::value><<<(data->numElements + kBlockSize - 1) / kBlockSize, kBlockSize>>>(
      points.data(), data->numElements, data->hessians.data(), deviceDiagonal);
  });
  checkCuda(cudaGetLastError(), "CUDA FEM diagonal kernel");
}

double CudaDeformationModelEnergy::computeMaxStepSize(const double *deviceX, const double *deviceDx) const
{
  if (!data->enableMaterialMaxStep || !data->numElements)
    return 1.0;
  const double one = 1.0;
  data->scalar.upload(&one, 1);
  data->directionNonzero.zero();
  directionNonzeroKernel<<<(data->numDOFs + kBlockSize - 1) / kBlockSize, kBlockSize>>>(
    data->numDOFs, deviceDx, data->directionNonzero.data());
  checkCuda(cudaGetLastError(), "CUDA FEM nonzero-direction kernel");
  data->visitQuadrature([&](const auto &points, auto) {
    const int count = static_cast<int>(points.size());
    maxStepKernel<<<(count + kBlockSize - 1) / kBlockSize, kBlockSize>>>(
      points.data(), count, deviceX, deviceDx, data->directionNonzero.data(), data->scalar.data());
  });
  checkCuda(cudaGetLastError(), "CUDA FEM material max-step kernel");
  if (data->cubic && data->cubicValidity == CubicValidity::WholeCell) {
    checkCuda(cudaMemcpyAsync(data->hexLimit.data(), data->scalar.data(), sizeof(double), cudaMemcpyDeviceToDevice), "CUDA cubic step limit");
    hexFastMaxStepKernel<<<(data->numElements + kBlockSize - 1) / kBlockSize, kBlockSize>>>(
      data->hexElements.data(), data->numElements, deviceX, deviceDx, kDetEpsilon, data->hexLimit.data(),
      data->hexNeedsSubdivision.data(), data->scalar.data());
    hexMaxStepKernel<<<(data->numElements + 31) / 32, 32>>>(
      data->hexElements.data(), data->numElements, deviceX, deviceDx, kDetEpsilon, data->hexLimit.data(),
      data->hexNeedsSubdivision.data(), data->scalar.data());
    checkCuda(cudaGetLastError(), "CUDA cubic whole-cell max-step certificate");
  }
  double alpha;
  data->scalar.download(&alpha, 1);
  return alpha;
}
}  // namespace pgo::SolidDeformationModel
