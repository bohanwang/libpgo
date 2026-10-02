#include "cudaFrozenHessianEnergy.h"
#include "cudaBlas.h"

#include <cmath>
#include <stdexcept>

namespace pgo::NonlinearOptimization
{
namespace
{
using CudaUtilities::DeviceBuffer;
using CudaUtilities::checkCuda;
constexpr int blockSize = 256;

__global__ void differenceKernel(int n, const double *x, const double *reference, double *dx)
{
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < n)
    dx[i] = x[i] - reference[i];
}

__global__ void scaledAddKernel(int n, double scale, const double *x, double *out)
{
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < n)
    out[i] += scale * x[i];
}
__global__ void productMagnitudeKernel(int n, const double *delta, const double *product, double *out, double scale, bool add)
{
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < n) {
    const double magnitude = scale * fabs(delta[i] * product[i]);
    out[i] = (add ? out[i] : 0.0) + magnitude;
  }
}

}  // namespace

class CudaFrozenHessianEnergy::Impl
{
public:
  explicit Impl(std::shared_ptr<CudaPotentialEnergy> model_):
    model(std::move(model_))
  {
    if (!model || model->getNumDOFs() <= 0)
      throw std::invalid_argument("Frozen CUDA Hessian requires a nonempty model.");
    n = model->getNumDOFs();
    reference.resize(n);
    difference.resize(n);
    product.resize(n);
    diagonal.resize(n);
    delta.resize(n);
    changeScale.resize(n);
  }

  void requireReady() const
  {
    if (!ready)
      throw std::logic_error("Set the frozen Hessian reference before evaluation.");
  }

  void productAt(const double *x)
  {
    requireReady();
    differenceKernel<<<(n + blockSize - 1) / blockSize, blockSize>>>(n, x, reference.data(), difference.data());
    checkCuda(cudaGetLastError(), "CUDA damping displacement");
    product.zero();
    model->addHessianVector(difference.data(), product.data());
  }

  void addProduct(const double *v, double *out) const
  {
    requireReady();
    product.zero();
    model->addHessianVector(v, product.data());
    scaledAddKernel<<<(n + blockSize - 1) / blockSize, blockSize>>>(n, scale, product.data(), out);
    checkCuda(cudaGetLastError(), "CUDA damping Hessian-vector product");
  }

  std::shared_ptr<CudaPotentialEnergy> model;
  int n = 0;
  double scale = 0;
  bool ready = false;
  DeviceBuffer<double> reference, difference, diagonal, delta, changeScale;
  mutable DeviceBuffer<double> product;
  detail::BlasHandle blas;
};

CudaFrozenHessianEnergy::CudaFrozenHessianEnergy(std::shared_ptr<CudaPotentialEnergy> model):
  impl_(std::make_unique<Impl>(std::move(model)))
{
}
CudaFrozenHessianEnergy::~CudaFrozenHessianEnergy() = default;

void CudaFrozenHessianEnergy::setReference(EigenSupport::ConstRefVecXd reference, double scale)
{
  auto &data = *impl_;
  if (reference.size() != data.n || !reference.allFinite() || !std::isfinite(scale) || scale < 0)
    throw std::invalid_argument("Frozen CUDA Hessian reference and nonnegative scale must be finite.");
  data.ready = false;
  data.reference.upload(reference.data(), data.n);
  data.model->prepareHessian(data.reference.data());
  data.diagonal.zero();
  data.model->addDiagonal(data.diagonal.data());
  data.scale = scale;
  data.ready = true;
}

int CudaFrozenHessianEnergy::getNumDOFs() const
{
  return impl_->n;
}

double CudaFrozenHessianEnergy::func(const double *deviceX)
{
  auto &data = *impl_;
  data.productAt(deviceX);
  return 0.5 * data.scale * data.blas.dot(data.n, data.difference.data(), data.product.data());
}

CudaPotentialEnergy::EnergyChange CudaFrozenHessianEnergy::energyChange(const double *deviceBase, const double *deviceTrial)
{
  auto &data = *impl_;
  data.productAt(deviceBase);
  differenceKernel<<<(data.n + blockSize - 1) / blockSize, blockSize>>>(data.n, deviceTrial, deviceBase, data.delta.data());
  checkCuda(cudaGetLastError(), "CUDA damping represented increment");
  const double first = data.blas.dot(data.n, data.delta.data(), data.product.data());
  productMagnitudeKernel<<<(data.n + blockSize - 1) / blockSize, blockSize>>>(data.n, data.delta.data(), data.product.data(), data.changeScale.data(), 1.0, false);
  checkCuda(cudaGetLastError(), "CUDA damping linear increment scale");
  data.product.zero();
  data.model->addHessianVector(data.delta.data(), data.product.data());
  const double second = 0.5 * data.blas.dot(data.n, data.delta.data(), data.product.data());
  productMagnitudeKernel<<<(data.n + blockSize - 1) / blockSize, blockSize>>>(data.n, data.delta.data(), data.product.data(), data.changeScale.data(), 0.5, true);
  checkCuda(cudaGetLastError(), "CUDA damping quadratic increment scale");
  double scale = 0;
  detail::checkBlas(cublasDasum(data.blas, data.n, data.changeScale.data(), 1, &scale), "CUDA damping increment scale reduction");
  return { data.scale * (first + second), data.scale * scale };
}

void CudaFrozenHessianEnergy::addGradient(const double *deviceX, double *deviceGradient)
{
  auto &data = *impl_;
  data.productAt(deviceX);
  scaledAddKernel<<<(data.n + blockSize - 1) / blockSize, blockSize>>>(data.n, data.scale, data.product.data(), deviceGradient);
  checkCuda(cudaGetLastError(), "CUDA damping gradient");
}

void CudaFrozenHessianEnergy::prepareHessian(const double *)
{
  impl_->requireReady();
}
void CudaFrozenHessianEnergy::addHessianVector(const double *deviceV, double *deviceOut) const
{
  impl_->addProduct(deviceV, deviceOut);
}

void CudaFrozenHessianEnergy::addDiagonal(double *deviceDiagonal) const
{
  const auto &data = *impl_;
  data.requireReady();
  scaledAddKernel<<<(data.n + blockSize - 1) / blockSize, blockSize>>>(data.n, data.scale, data.diagonal.data(), deviceDiagonal);
  checkCuda(cudaGetLastError(), "CUDA damping Hessian diagonal");
}

}  // namespace pgo::NonlinearOptimization
