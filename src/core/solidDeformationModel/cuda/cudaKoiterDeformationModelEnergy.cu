#include "cudaKoiterDeformationModelEnergy.h"

#include "cudaDeviceBuffer.h"
#include "simulationMesh.h"
#include <libuipc/math/make_spd.h>

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/QR>
#include <cuda_runtime.h>

#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace pgo::SolidDeformationModel
{
namespace
{
using CudaUtilities::DeviceBuffer;
using CudaUtilities::checkCuda;
using V3 = Eigen::Vector3d;
using V4 = Eigen::Vector4d;
using V18 = Eigen::Matrix<double, 18, 1>;
using M3 = Eigen::Matrix3d;
using M4 = Eigen::Matrix4d;
using M18 = Eigen::Matrix<double, 18, 18>;
using J3 = Eigen::Matrix<double, 3, 18>;
using J4 = Eigen::Matrix<double, 4, 18>;
constexpr int kBlockSize = 64;

struct Patch
{
  int vertices[6];
  double rest[18];
  double restB[4];
  double membrane[16];
  double bending[16];
  bool zeroRestCurvature;
};

__host__ __device__ int oppositeVertex(int edge)
{
  return 3 + (edge + 1) % 3;
}

__host__ __device__ V4 curvatureForm(const V3 &entries)
{
  V4 b;
  b << entries[0] + entries[1], entries[0], entries[0], entries[0] + entries[2];
  return b;
}

__host__ __device__ V3 normal(const V3 x[6], int a, int b, int c)
{
  return (x[b] - x[a]).cross(x[c] - x[a]);
}

__host__ __device__ V4 secondForm(const Patch &patch, const V3 x[6])
{
  const V3 central = normal(x, 0, 1, 2);
  V3 entries;
  for (int i = 0; i < 3; ++i) {
    const int p = (i + 1) % 3, q = (i + 2) % 3, o = oppositeVertex(i);
    const V3 adjacent = patch.vertices[o] >= 0 ? normal(x, o, q, p) : V3::Zero().eval();
    entries[i] = (x[p] + x[q] - 2 * x[i]).dot(adjacent) / (central + adjacent).norm();
  }
  return curvatureForm(entries);
}

__device__ void gather(const Patch &patch, const double *u, V3 x[6], V3 displacement[6])
{
  for (int a = 0; a < 6; ++a) {
    displacement[a].setZero();
    if (patch.vertices[a] >= 0)
      for (int c = 0; c < 3; ++c)
        displacement[a][c] = u[3 * patch.vertices[a] + c] - u[3 * patch.vertices[0] + c];
    x[a] = Eigen::Map<const V3>(patch.rest + 3 * a) + displacement[a];
  }
}

__device__ double absoluteDot(const V3 &a, const V3 &b)
{
  return a.cwiseAbs().dot(b.cwiseAbs());
}

__device__ V3 crossMagnitude(const V3 &a, const V3 &b)
{
  V3 result;
  for (int i = 0; i < 3; ++i) {
    const int j = (i + 1) % 3, k = (i + 2) % 3;
    result[i] = fabs(a[j] * b[k]) + fabs(a[k] * b[j]);
  }
  return result;
}

__device__ V3 normalChange(const V3 x[6], const V3 d[6], int a, int b, int c, V3 &scale)
{
  const V3 e1 = x[b] - x[a], e2 = x[c] - x[a];
  const V3 d1 = d[b] - d[a], d2 = d[c] - d[a];
  scale = crossMagnitude(e1, d2) + crossMagnitude(d1, e2) + crossMagnitude(d1, d2);
  return e1.cross(d2) + d1.cross(e2) + d1.cross(d2);
}

// Exact polynomial metric difference and rationalized normal-length difference.
// Keeping the small changes explicit avoids subtracting two strained energies
// or two nearly equal fundamental forms during Newton's line search.
__device__ void formChange(const Patch &patch, const V3 x[6], const V3 d[6],
  V4 &da, V4 &db, V4 &scaleA, V4 &scaleB)
{
  const V3 e[2] = { x[1] - x[0], x[2] - x[0] };
  const V3 de[2] = { d[1] - d[0], d[2] - d[0] };
  for (int c = 0; c < 2; ++c)
    for (int r = 0; r < 2; ++r) {
      da[2 * c + r] = e[r].dot(de[c]) + de[r].dot(e[c]) + de[r].dot(de[c]);
      scaleA[2 * c + r] = absoluteDot(e[r], de[c]) + absoluteDot(de[r], e[c]) + absoluteDot(de[r], de[c]);
    }
  const V3 central = normal(x, 0, 1, 2);
  V3 centralScale;
  const V3 centralChange = normalChange(x, d, 0, 1, 2, centralScale);
  V3 entries, scales;
  for (int i = 0; i < 3; ++i) {
    const int p = (i + 1) % 3, q = (i + 2) % 3, o = oppositeVertex(i);
    // Missing neighbors contribute identically zero, including their derivatives.
    if (patch.vertices[o] < 0) {
      entries[i] = scales[i] = 0;
      continue;
    }
    const V3 adjacent = normal(x, o, q, p);
    V3 adjacentScale;
    const V3 adjacentChange = normalChange(x, d, o, q, p, adjacentScale);
    const V3 m = central + adjacent, dm = centralChange + adjacentChange;
    const V3 dmScale = centralScale + adjacentScale;
    const V3 v = x[p] + x[q] - 2 * x[i], dv = d[p] + d[q] - 2 * d[i];
    const double s = m.norm(), trialS = (m + dm).norm();
    const double numerator = v.dot(adjacent);
    const double numeratorChange = v.dot(adjacentChange) + dv.dot(adjacent) + dv.dot(adjacentChange);
    const double lengthChange = (2 * m.dot(dm) + dm.squaredNorm()) / (trialS + s);
    const double curvature = numerator / s;
    entries[i] = (numeratorChange - curvature * lengthChange) / trialS;
    const double lengthScale = (2 * absoluteDot(m, dmScale) + dmScale.squaredNorm()) / (trialS + s);
    const double numeratorScale = absoluteDot(v, adjacentScale) + absoluteDot(dv, adjacent) + absoluteDot(dv, adjacentScale);
    scales[i] = (numeratorScale + absoluteDot(v, adjacent) / s * lengthScale) / trialS;
  }
  db = curvatureForm(entries);
  scaleB = curvatureForm(scales);
}

__device__ void strainForms(const Patch &patch, const V3 displacement[6], V4 &a, V4 &b)
{
  V3 rest[6];
  for (int i = 0; i < 6; ++i)
    rest[i] = Eigen::Map<const V3>(patch.rest + 3 * i);
  V4 scaleA, scaleB;
  formChange(patch, rest, displacement, a, b, scaleA, scaleB);
  if (patch.zeroRestCurvature)
    b += Eigen::Map<const V4>(patch.restB);
}

__device__ M3 crossMatrix(const V3 &v)
{
  M3 result;
  result << 0, -v[2], v[1], v[2], 0, -v[0], -v[1], v[0], 0;
  return result;
}

__device__ J3 normalJacobian(const V3 x[6], int a, int b, int c)
{
  J3 result = J3::Zero();
  result.block<3, 3>(0, 3 * a) += crossMatrix(x[c] - x[b]);
  result.block<3, 3>(0, 3 * b) += crossMatrix(x[a] - x[c]);
  result.block<3, 3>(0, 3 * c) += crossMatrix(x[b] - x[a]);
  return result;
}

// Hessian of w.dot((xb-xa).cross(xc-xa)); its only nonzero blocks
// are constant skew matrices. This contracts the three normal Hessians
// directly, without storing three additional 18x18 matrices per thread.
__device__ void addNormalHessian(M18 &H, const V3 &w, int a, int b, int c)
{
  const M3 W = crossMatrix(w);
  H.block<3, 3>(3 * a, 3 * b) -= W;
  H.block<3, 3>(3 * a, 3 * c) += W;
  H.block<3, 3>(3 * b, 3 * c) -= W;
  H.block<3, 3>(3 * b, 3 * a) += W;
  H.block<3, 3>(3 * c, 3 * a) -= W;
  H.block<3, 3>(3 * c, 3 * b) += W;
}

template<bool Hessian>
__device__ void formDerivatives(const Patch &patch, const V3 x[6],
  const V4 &gradientA, const V4 &gradientB, J4 &Ja, J4 &Jb, M18 &geometricHessian)
{
  Ja.setZero();
  Jb.setZero();
  if constexpr (Hessian)
    geometricHessian.setZero();
  const V3 e1 = x[1] - x[0], e2 = x[2] - x[0];
  const int s1[3] = { -1, 1, 0 }, s2[3] = { -1, 0, 1 };
  for (int a = 0; a < 3; ++a) {
    Ja.block<1, 3>(0, 3 * a) = 2 * s1[a] * e1.transpose();
    Ja.block<1, 3>(1, 3 * a) = (s1[a] * e2 + s2[a] * e1).transpose();
    Ja.block<1, 3>(2, 3 * a) = Ja.block<1, 3>(1, 3 * a);
    Ja.block<1, 3>(3, 3 * a) = 2 * s2[a] * e2.transpose();
    if constexpr (Hessian)
      for (int b = 0; b < 3; ++b) {
        const double value = 2 * gradientA[0] * s1[a] * s1[b] + 2 * gradientA[3] * s2[a] * s2[b] +
          (gradientA[1] + gradientA[2]) * (s1[a] * s2[b] + s2[a] * s1[b]);
        for (int c = 0; c < 3; ++c)
          geometricHessian(3 * a + c, 3 * b + c) += value;
      }
  }
  const V3 central = e1.cross(e2);
  const J3 centralJ = normalJacobian(x, 0, 1, 2);
  for (int i = 0; i < 3; ++i) {
    const int p = (i + 1) % 3, q = (i + 2) % 3, o = oppositeVertex(i);
    if (patch.vertices[o] < 0)
      continue;
    const V3 adjacent = normal(x, o, q, p), m = central + adjacent;
    const V3 v = x[p] + x[q] - 2 * x[i];
    const J3 adjacentJ = normalJacobian(x, o, q, p), mJ = centralJ + adjacentJ;
    J3 vJ = J3::Zero();
    for (int c = 0; c < 3; ++c) {
      vJ(c, 3 * p + c) = 1;
      vJ(c, 3 * q + c) = 1;
      vJ(c, 3 * i + c) = -2;
    }
    const double s = m.norm(), inverse = 1 / s, inverse3 = inverse * inverse * inverse;
    const double numerator = v.dot(adjacent);
    const V18 numeratorJ = vJ.transpose() * adjacent + adjacentJ.transpose() * v;
    const V18 normJ = mJ.transpose() * m;
    const V18 curvatureJ = inverse * numeratorJ - numerator * inverse3 * normJ;
    if (i == 0)
      for (int row = 0; row < 4; ++row)
        Jb.row(row) += curvatureJ.transpose();
    else
      Jb.row(i == 1 ? 0 : 3) += curvatureJ.transpose();

    if constexpr (Hessian) {
      const double weight = i == 0 ? gradientB.sum() : gradientB[i == 1 ? 0 : 3];
      for (int col = 0; col < 18; ++col)
        for (int row = 0; row < 18; ++row) {
          const double numeratorH = vJ.col(row).dot(adjacentJ.col(col)) + adjacentJ.col(row).dot(vJ.col(col));
          const double value = inverse * numeratorH - inverse3 * (numeratorJ[row] * normJ[col] + normJ[row] * numeratorJ[col]) +
            numerator * (3 * inverse3 * inverse * inverse * normJ[row] * normJ[col] - inverse3 * mJ.col(row).dot(mJ.col(col)));
          geometricHessian(row, col) += weight * value;
        }
      addNormalHessian(geometricHessian, weight * (inverse * v - numerator * inverse3 * m), o, q, p);
      addNormalHessian(geometricHessian, -weight * numerator * inverse3 * m, 0, 1, 2);
    }
  }
}

__global__ void energyKernel(const Patch *patches, int count, const double *x, double *result)
{
  const int e = blockIdx.x * blockDim.x + threadIdx.x;
  if (e >= count)
    return;
  const Patch &patch = patches[e];
  V3 positions[6], displacement[6];
  gather(patch, x, positions, displacement);
  V4 a, b;
  strainForms(patch, displacement, a, b);
  const double energy = 0.5 * (a.dot(Eigen::Map<const M4>(patch.membrane) * a) + b.dot(Eigen::Map<const M4>(patch.bending) * b));
  atomicAdd(result, energy);
}

__global__ void energyChangeKernel(const Patch *patches, int count,
  const double *base, const double *trial, double *result)
{
  const int e = blockIdx.x * blockDim.x + threadIdx.x;
  if (e >= count)
    return;
  const Patch &patch = patches[e];
  V3 positions[6], displacement[6], increment[6];
  gather(patch, base, positions, displacement);
  for (int a = 0; a < 6; ++a) {
    increment[a].setZero();
    if (patch.vertices[a] >= 0)
      for (int c = 0; c < 3; ++c)
        increment[a][c] = (trial[3 * patch.vertices[a] + c] - base[3 * patch.vertices[a] + c]) -
          (trial[3 * patch.vertices[0] + c] - base[3 * patch.vertices[0] + c]);
  }
  V4 a, b, da, db, scaleA, scaleB;
  strainForms(patch, displacement, a, b);
  formChange(patch, positions, increment, da, db, scaleA, scaleB);
  const M4 membrane = Eigen::Map<const M4>(patch.membrane), bending = Eigen::Map<const M4>(patch.bending);
  const double change = a.dot(membrane * da) + 0.5 * da.dot(membrane * da) + b.dot(bending * db) + 0.5 * db.dot(bending * db);
  const double scale = (a.cwiseAbs() + 0.5 * scaleA).dot(membrane.cwiseAbs() * scaleA) +
    (b.cwiseAbs() + 0.5 * scaleB).dot(bending.cwiseAbs() * scaleB);
  atomicAdd(result, change);
  atomicAdd(result + 1, scale);
}

__global__ void gradientKernel(const Patch *patches, int count, const double *x, double *result)
{
  const int e = blockIdx.x * blockDim.x + threadIdx.x;
  if (e >= count)
    return;
  const Patch &patch = patches[e];
  V3 positions[6], displacement[6];
  gather(patch, x, positions, displacement);
  V4 a, b;
  strainForms(patch, displacement, a, b);
  const V4 gradientA = Eigen::Map<const M4>(patch.membrane) * a;
  const V4 gradientB = Eigen::Map<const M4>(patch.bending) * b;
  J4 Ja, Jb;
  M18 unused;
  formDerivatives<false>(patch, positions, gradientA, gradientB, Ja, Jb, unused);
  const V18 gradient = Ja.transpose() * gradientA + Jb.transpose() * gradientB;
  for (int a = 0; a < 6; ++a)
    if (patch.vertices[a] >= 0)
      for (int c = 0; c < 3; ++c)
        atomicAdd(result + 3 * patch.vertices[a] + c, gradient[3 * a + c]);
}

__global__ void hessianKernel(const Patch *patches, int count, const double *x, bool project, double *result)
{
  const int e = blockIdx.x * blockDim.x + threadIdx.x;
  if (e >= count)
    return;
  const Patch &patch = patches[e];
  V3 positions[6], displacement[6];
  gather(patch, x, positions, displacement);
  V4 a, b;
  strainForms(patch, displacement, a, b);
  const M4 membrane = Eigen::Map<const M4>(patch.membrane), bending = Eigen::Map<const M4>(patch.bending);
  J4 Ja, Jb;
  M18 H;
  formDerivatives<true>(patch, positions, membrane * a, bending * b, Ja, Jb, H);
  H += Ja.transpose() * membrane * Ja + Jb.transpose() * bending * Jb;
  if (project)
    uipc::backend::cuda::make_spd<18>(H);
  Eigen::Map<M18>(result + static_cast<std::size_t>(e) * 324) = H;
}

__global__ void hessianVectorKernel(const Patch *patches, int count, const double *hessians,
  const double *v, double *result)
{
  const int e = blockIdx.x * blockDim.x + threadIdx.x;
  if (e >= count)
    return;
  const Patch &patch = patches[e];
  V18 local = V18::Zero();
  for (int a = 0; a < 6; ++a)
    if (patch.vertices[a] >= 0)
      for (int c = 0; c < 3; ++c)
        local[3 * a + c] = v[3 * patch.vertices[a] + c];
  const V18 product = Eigen::Map<const M18>(hessians + static_cast<std::size_t>(e) * 324) * local;
  for (int a = 0; a < 6; ++a)
    if (patch.vertices[a] >= 0)
      for (int c = 0; c < 3; ++c)
        atomicAdd(result + 3 * patch.vertices[a] + c, product[3 * a + c]);
}

__global__ void diagonalKernel(const Patch *patches, int count, const double *hessians, double *result)
{
  const int e = blockIdx.x * blockDim.x + threadIdx.x;
  if (e >= count)
    return;
  const Patch &patch = patches[e];
  const Eigen::Map<const M18> H(hessians + static_cast<std::size_t>(e) * 324);
  for (int a = 0; a < 6; ++a)
    if (patch.vertices[a] >= 0)
      for (int c = 0; c < 3; ++c) {
        double value = 0;
        // A closed patch can reference the same neighboring vertex more than
        // once. Its global diagonal includes those off-diagonal local entries.
        for (int b = 0; b < 6; ++b)
          if (patch.vertices[a] == patch.vertices[b])
            value += H(3 * a + c, 3 * b + c);
        atomicAdd(result + 3 * patch.vertices[a] + c, value);
      }
}

M4 materialTensor(const Eigen::Matrix2d &inverseRest, double E, double nu, double thicknessFactor)
{
  const double lambda = E * nu / ((1 + nu) * (1 - 2 * nu)), mu = E / (2 * (1 + nu));
  M4 hessian;
  hessian << lambda + 2 * mu, 0, 0, lambda,
    0, 0, 2 * mu, 0,
    0, 2 * mu, 0, 0,
    lambda, 0, 0, lambda + 2 * mu;
  M4 mapping = M4::Zero();
  mapping.topLeftCorner<2, 2>() = inverseRest;
  mapping.bottomRightCorner<2, 2>() = inverseRest;
  return thicknessFactor * mapping.transpose() * hessian * mapping;
}
}  // namespace

class CudaKoiterDeformationModelEnergy::Impl
{
public:
  int numDOFs = 0;
  int numElements = 0;
  bool hessianPrepared = false;
  DeviceBuffer<Patch> patches;
  DeviceBuffer<double> hessians;
  DeviceBuffer<double> scalars{ 2 };
};

CudaKoiterDeformationModelEnergy::CudaKoiterDeformationModelEnergy(const SimulationMesh &mesh,
  EigenSupport::ConstRefVecXd parameters, bool zeroRestCurvature):
  data(std::make_unique<Impl>())
{
  if (mesh.getElementType() != SimulationMeshType::SHELL || mesh.getNumElementVertices() != 6)
    throw std::invalid_argument("CUDA Koiter energy requires six-vertex shell patches.");
  if (mesh.getNumVertices() > std::numeric_limits<int>::max() / 3)
    throw std::length_error("CUDA Koiter DOF count exceeds integer range.");
  data->numDOFs = 3 * mesh.getNumVertices();
  data->numElements = mesh.getNumElements();
  if (parameters.size() != 5 * static_cast<Eigen::Index>(data->numElements) || !parameters.allFinite())
    throw std::invalid_argument("CUDA Koiter requires five finite elastic parameters per element.");
  std::vector<Patch> patches(data->numElements);
  for (int e = 0; e < data->numElements; ++e) {
    Patch &patch = patches[e];
    V3 x[6];
    patch.zeroRestCurvature = zeroRestCurvature;
    for (int a = 0; a < 6; ++a) {
      patch.vertices[a] = mesh.getVertexIndex(e, a);
      if (patch.vertices[a] >= mesh.getNumVertices() || (a < 3 && patch.vertices[a] < 0))
        throw std::invalid_argument("CUDA Koiter patch has invalid vertex indices.");
      x[a].setZero();
      if (patch.vertices[a] >= 0)
        mesh.getVertex(patch.vertices[a], x[a].data());
    }
    const V3 origin = x[0];
    for (int a = 0; a < 6; ++a) {
      if (patch.vertices[a] >= 0)
        x[a] -= origin;
      if (!x[a].allFinite())
        throw std::invalid_argument("CUDA Koiter rest positions must be finite.");
      Eigen::Map<V3>(patch.rest + 3 * a) = x[a];
    }
    const V3 e1 = x[1] - x[0], e2 = x[2] - x[0];
    const double area = 0.5 * e1.cross(e2).norm();
    Eigen::Matrix2d restA;
    restA << e1.squaredNorm(), e1.dot(e2), e1.dot(e2), e2.squaredNorm();
    const auto decomposition = restA.fullPivHouseholderQr();
    const Eigen::Matrix2d inverse = decomposition.inverse();
    const V4 restB = secondForm(patch, x);
    if (!(area > 0) || !std::isfinite(area) || !decomposition.isInvertible() || !inverse.allFinite() || !restB.allFinite())
      throw std::invalid_argument("CUDA Koiter requires nondegenerate, representable rest fundamental forms.");
    Eigen::Map<V4>(patch.restB) = restB;
    const auto p = parameters.segment<5>(5 * e);
    if (!(p[0] > 0) || !(p[2] > 0) || !(p[4] > 0) || !(p[1] > -1 && p[1] < 0.5) || !(p[3] > -1 && p[3] < 0.5))
      throw std::invalid_argument("CUDA Koiter requires positive Young moduli/thickness and Poisson ratios in (-1, 0.5).");
    const M4 membrane = materialTensor(inverse, p[0], p[1], area * p[4]);
    const M4 bending = materialTensor(inverse, p[2], p[3], area * p[4] * p[4] * p[4] / 12);
    if (!membrane.allFinite() || !bending.allFinite())
      throw std::invalid_argument("CUDA Koiter material tensors must be representable.");
    Eigen::Map<M4>(patch.membrane) = membrane;
    Eigen::Map<M4>(patch.bending) = bending;
  }
  data->patches.upload(patches.data(), patches.size());
  data->hessians.resize(static_cast<std::size_t>(data->numElements) * 324);
}

CudaKoiterDeformationModelEnergy::~CudaKoiterDeformationModelEnergy() = default;

int CudaKoiterDeformationModelEnergy::getNumDOFs() const
{
  return data->numDOFs;
}

double CudaKoiterDeformationModelEnergy::func(const double *deviceX)
{
  data->scalars.zero();
  if (data->numElements) {
    energyKernel<<<(data->numElements + kBlockSize - 1) / kBlockSize, kBlockSize>>>(
      data->patches.data(), data->numElements, deviceX, data->scalars.data());
    checkCuda(cudaGetLastError(), "CUDA Koiter energy kernel");
  }
  double values[2];
  data->scalars.download(values, 2);
  return values[0];
}

CudaKoiterDeformationModelEnergy::EnergyChange CudaKoiterDeformationModelEnergy::energyChange(
  const double *deviceBase, const double *deviceTrial)
{
  data->scalars.zero();
  if (data->numElements) {
    energyChangeKernel<<<(data->numElements + kBlockSize - 1) / kBlockSize, kBlockSize>>>(
      data->patches.data(), data->numElements, deviceBase, deviceTrial, data->scalars.data());
    checkCuda(cudaGetLastError(), "CUDA Koiter paired energy-change kernel");
  }
  double values[2];
  data->scalars.download(values, 2);
  return { values[0], values[1] };
}

void CudaKoiterDeformationModelEnergy::addGradient(const double *deviceX, double *deviceGradient)
{
  if (!data->numElements)
    return;
  gradientKernel<<<(data->numElements + kBlockSize - 1) / kBlockSize, kBlockSize>>>(
    data->patches.data(), data->numElements, deviceX, deviceGradient);
  checkCuda(cudaGetLastError(), "CUDA Koiter gradient kernel");
}

void CudaKoiterDeformationModelEnergy::prepareHessian(const double *deviceX)
{
  prepareHessian(deviceX, true);
}

void CudaKoiterDeformationModelEnergy::prepareExactHessian(const double *deviceX)
{
  prepareHessian(deviceX, false);
}

void CudaKoiterDeformationModelEnergy::prepareHessian(const double *deviceX, bool project)
{
  if (data->numElements) {
    hessianKernel<<<(data->numElements + kBlockSize - 1) / kBlockSize, kBlockSize>>>(
      data->patches.data(), data->numElements, deviceX, project, data->hessians.data());
    checkCuda(cudaGetLastError(), "CUDA Koiter Hessian kernel");
  }
  data->hessianPrepared = true;
}

void CudaKoiterDeformationModelEnergy::addHessianVector(const double *deviceV, double *deviceOut) const
{
  if (!data->hessianPrepared)
    throw std::logic_error("CUDA Koiter Hessian must be prepared before applying it.");
  if (!data->numElements)
    return;
  hessianVectorKernel<<<(data->numElements + kBlockSize - 1) / kBlockSize, kBlockSize>>>(
    data->patches.data(), data->numElements, data->hessians.data(), deviceV, deviceOut);
  checkCuda(cudaGetLastError(), "CUDA Koiter Hessian-vector kernel");
}

void CudaKoiterDeformationModelEnergy::addDiagonal(double *deviceDiagonal) const
{
  if (!data->hessianPrepared)
    throw std::logic_error("CUDA Koiter Hessian must be prepared before querying its diagonal.");
  if (!data->numElements)
    return;
  diagonalKernel<<<(data->numElements + kBlockSize - 1) / kBlockSize, kBlockSize>>>(
    data->patches.data(), data->numElements, data->hessians.data(), deviceDiagonal);
  checkCuda(cudaGetLastError(), "CUDA Koiter diagonal kernel");
}
}  // namespace pgo::SolidDeformationModel
