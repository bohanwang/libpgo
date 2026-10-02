#include "cudaQuadraticEnergy.h"
#include "cudaBlas.h"

#include <cuda_runtime.h>

#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace pgo::NonlinearOptimization
{
namespace
{
using CudaUtilities::DeviceBuffer;
using CudaUtilities::checkCuda;
constexpr int blockSize = 256;

__global__ void displacementKernel(int n, const double *x, const double *reference, double *dx)
{
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < n)
    dx[i] = x[i] - reference[i];
}

__global__ void csrProductKernel(int n, const int *rowOffsets, const int *columns,
  const double *values, const double *x, double *out, bool add)
{
  const int row = blockIdx.x * blockDim.x + threadIdx.x;
  if (row >= n)
    return;
  double sum = 0;
  for (int k = rowOffsets[row]; k < rowOffsets[row + 1]; ++k)
    sum += values[k] * x[columns[k]];
  if (add)
    out[row] += sum;
  else
    out[row] = sum;
}

__global__ void addLinearGradientKernel(int n, const double *Ax, const double *linear, double *gradient)
{
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < n)
    gradient[i] += Ax[i] + linear[i];
}

__global__ void changeScaleKernel(int n, const int *rowOffsets, const int *columns,
  const double *values, const double *baseDisplacement, const double *delta,
  const double *linear, double *out)
{
  const int row = blockIdx.x * blockDim.x + threadIdx.x;
  if (row >= n)
    return;
  double magnitude = fabs(linear[row]);
  for (int k = rowOffsets[row]; k < rowOffsets[row + 1]; ++k)
    magnitude += fabs(values[k]) * (fabs(baseDisplacement[columns[k]]) + 0.5 * fabs(delta[columns[k]]));
  out[row] = fabs(delta[row]) * magnitude;
}

__global__ void addDiagonalKernel(int n, const double *diagonal, double *out)
{
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < n)
    out[i] += diagonal[i];
}
}  // namespace

class CudaQuadraticEnergy::Impl
{
public:
  Impl(const EigenSupport::SpMatD &A, EigenSupport::ConstRefVecXd ref,
    EigenSupport::ConstRefVecXd b, double c):
    constant(c)
  {
    if (A.rows() <= 0 || A.rows() != A.cols() || A.rows() > std::numeric_limits<int>::max() ||
      ref.size() != A.rows() || b.size() != A.rows() || !ref.allFinite() || !b.allFinite() || !std::isfinite(c))
      throw std::invalid_argument("CudaQuadraticEnergy requires finite, matching nonempty square matrix/vector inputs.");
    n = static_cast<int>(A.rows());
    for (Eigen::Index i = 0; i < A.nonZeros(); ++i)
      if (!std::isfinite(A.valuePtr()[i]))
        throw std::invalid_argument("CudaQuadraticEnergy matrix contains nonfinite values.");
    // A symmetric matrix is required for the gradient and PCG contract. Small
    // roundoff asymmetry is averaged once during upload.
    const EigenSupport::SpMatD transposed = A.transpose();
    const EigenSupport::SpMatD asymmetry = A - transposed;
    if (asymmetry.norm() > 1e-12 * std::max(1.0, A.norm()))
      throw std::invalid_argument("CudaQuadraticEnergy matrix must be symmetric.");
    Eigen::SparseMatrix<double, Eigen::RowMajor, int> csr = 0.5 * (A + transposed);
    csr.makeCompressed();
    if (csr.nonZeros() > std::numeric_limits<int>::max())
      throw std::length_error("CudaQuadraticEnergy matrix exceeds CUDA CSR index capacity.");
    rowOffsets.upload(csr.outerIndexPtr(), n + 1);
    columns.upload(csr.innerIndexPtr(), csr.nonZeros());
    values.upload(csr.valuePtr(), csr.nonZeros());
    std::vector<double> diagonalHost(n);
    for (int i = 0; i < n; ++i)
      diagonalHost[i] = csr.coeff(i, i);
    diagonal.upload(diagonalHost.data(), n);
    reference.upload(ref.data(), n);
    linear.upload(b.data(), n);
    dx.resize(n);
    product.resize(n);
    delta.resize(n);
    changeGradient.resize(n);
    changeScale.resize(n);
  }

  void multiply(const double *x, double *out, bool add) const
  {
    csrProductKernel<<<(n + blockSize - 1) / blockSize, blockSize>>>(
      n, rowOffsets.data(), columns.data(), values.data(), x, out, add);
    checkCuda(cudaGetLastError(), "CUDA quadratic sparse product");
  }
  void displacement(const double *x)
  {
    displacementKernel<<<(n + blockSize - 1) / blockSize, blockSize>>>(n, x, reference.data(), dx.data());
    checkCuda(cudaGetLastError(), "CUDA quadratic displacement");
  }

  int n = 0;
  double constant;
  DeviceBuffer<int> rowOffsets, columns;
  DeviceBuffer<double> values, diagonal, reference, linear, dx, product, delta, changeGradient, changeScale;
  detail::BlasHandle blas;
};

CudaQuadraticEnergy::CudaQuadraticEnergy(const EigenSupport::SpMatD &A,
  EigenSupport::ConstRefVecXd reference, EigenSupport::ConstRefVecXd linear, double constant):
  impl_(std::make_unique<Impl>(A, reference, linear, constant))
{
}
CudaQuadraticEnergy::~CudaQuadraticEnergy() = default;
int CudaQuadraticEnergy::getNumDOFs() const
{
  return impl_->n;
}

double CudaQuadraticEnergy::func(const double *deviceX)
{
  auto &d = *impl_;
  d.displacement(deviceX);
  d.multiply(d.dx.data(), d.product.data(), false);
  return 0.5 * d.blas.dot(d.n, d.dx.data(), d.product.data()) +
    d.blas.dot(d.n, d.dx.data(), d.linear.data()) + d.constant;
}

CudaPotentialEnergy::EnergyChange CudaQuadraticEnergy::energyChange(const double *deviceBase, const double *deviceTrial)
{
  auto &d = *impl_;
  // Evaluate the exact quadratic increment without subtracting two energies
  // containing the same potentially large constant or kinetic offset.
  d.displacement(deviceBase);
  displacementKernel<<<(d.n + blockSize - 1) / blockSize, blockSize>>>(d.n, deviceTrial, deviceBase, d.delta.data());
  checkCuda(cudaGetLastError(), "CUDA quadratic represented increment");
  d.multiply(d.dx.data(), d.product.data(), false);
  d.changeGradient.zero();
  addLinearGradientKernel<<<(d.n + blockSize - 1) / blockSize, blockSize>>>(d.n, d.product.data(), d.linear.data(), d.changeGradient.data());
  checkCuda(cudaGetLastError(), "CUDA quadratic base gradient");
  const double first = d.blas.dot(d.n, d.delta.data(), d.changeGradient.data());
  d.multiply(d.delta.data(), d.product.data(), false);
  const double second = 0.5 * d.blas.dot(d.n, d.delta.data(), d.product.data());
  changeScaleKernel<<<(d.n + blockSize - 1) / blockSize, blockSize>>>(d.n, d.rowOffsets.data(), d.columns.data(),
    d.values.data(), d.dx.data(), d.delta.data(), d.linear.data(), d.changeScale.data());
  checkCuda(cudaGetLastError(), "CUDA quadratic increment roundoff scale");
  double scale = 0;
  detail::checkBlas(cublasDasum(d.blas, d.n, d.changeScale.data(), 1, &scale), "CUDA quadratic increment scale reduction");
  return { first + second, scale };
}

void CudaQuadraticEnergy::addGradient(const double *deviceX, double *deviceGradient)
{
  auto &d = *impl_;
  d.displacement(deviceX);
  d.multiply(d.dx.data(), d.product.data(), false);
  addLinearGradientKernel<<<(d.n + blockSize - 1) / blockSize, blockSize>>>(
    d.n, d.product.data(), d.linear.data(), deviceGradient);
  checkCuda(cudaGetLastError(), "CUDA quadratic gradient");
}

void CudaQuadraticEnergy::prepareHessian(const double *)
{
}
void CudaQuadraticEnergy::addHessianVector(const double *deviceV, double *deviceOut) const
{
  impl_->multiply(deviceV, deviceOut, true);
}
void CudaQuadraticEnergy::addDiagonal(double *deviceDiagonal) const
{
  const auto &d = *impl_;
  addDiagonalKernel<<<(d.n + blockSize - 1) / blockSize, blockSize>>>(d.n, d.diagonal.data(), deviceDiagonal);
  checkCuda(cudaGetLastError(), "CUDA quadratic diagonal");
}

}  // namespace pgo::NonlinearOptimization
