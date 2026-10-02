#include "cudaDeviceBuffer.h"
#include "cudaFrozenHessianEnergy.h"
#include "cudaQuadraticEnergy.h"

#include <gtest/gtest.h>

#include <cmath>
#include <memory>

namespace
{
namespace ES = pgo::EigenSupport;
namespace NO = pgo::NonlinearOptimization;
using pgo::CudaUtilities::DeviceBuffer;
using pgo::CudaUtilities::checkCuda;

class CudaEnergyChangeTest : public testing::Test
{
protected:
  void SetUp() override
  {
    int count = 0;
    if (cudaGetDeviceCount(&count) != cudaSuccess || count == 0)
      GTEST_SKIP() << "No CUDA device is accessible.";
    checkCuda(cudaSetDevice(0), "CUDA energy-change test device selection");
  }
};

TEST_F(CudaEnergyChangeTest, SparseQuadraticResolvesChangesBelowAbsoluteEnergySpacing)
{
  ES::MXd matrix(4, 4);
  matrix << 5, 1, 0, -0.5,
    1, 4, 0.5, 0,
    0, 0.5, 3, 1,
    -0.5, 0, 1, 2;
  const ES::VXd reference = (ES::VXd(4) << 100, -200, 300, 400).finished();
  const ES::VXd linear = (ES::VXd(4) << -2, 3, 1, -4).finished();
  const ES::VXd base = reference + (ES::VXd(4) << 0.5, 1, -1, 0.5).finished();
  const ES::VXd trial = base + (ES::VXd(4) << 1e-7, 0.5e-7, -2e-7, 3e-7).finished();
  NO::CudaQuadraticEnergy energy(matrix.sparseView(), reference, linear, 1e20);
  DeviceBuffer<double> deviceBase, deviceTrial;
  deviceBase.upload(base.data(), base.size());
  deviceTrial.upload(trial.data(), trial.size());
  EXPECT_DOUBLE_EQ(energy.func(deviceBase.data()), energy.func(deviceTrial.data()));
  long double expected = 0;
  for (int i = 0; i < 4; ++i) {
    const long double delta = static_cast<long double>(trial[i]) - base[i];
    long double gradient = linear[i], hDelta = 0;
    for (int j = 0; j < 4; ++j) {
      gradient += static_cast<long double>(matrix(i, j)) * (static_cast<long double>(base[j]) - reference[j]);
      hDelta += static_cast<long double>(matrix(i, j)) * (static_cast<long double>(trial[j]) - base[j]);
    }
    expected += delta * (gradient + 0.5L * hDelta);
  }
  const auto change = energy.energyChange(deviceBase.data(), deviceTrial.data());
  EXPECT_NEAR(change.value, static_cast<double>(expected), 3e-21);
  EXPECT_GE(change.roundoffScale, std::abs(change.value));
  EXPECT_LT(change.roundoffScale, 1e-3);
  EXPECT_DOUBLE_EQ(energy.energyChange(deviceBase.data(), deviceBase.data()).value, 0);
}

TEST_F(CudaEnergyChangeTest, FrozenHessianResolvesPureSecondOrderIncrement)
{
  ES::SpMatD matrix(2, 2);
  matrix.insert(0, 0) = 2;
  matrix.insert(1, 1) = 3;
  auto model = std::make_shared<NO::CudaQuadraticEnergy>(matrix, ES::VXd::Zero(2), ES::VXd::Zero(2));
  NO::CudaFrozenHessianEnergy energy(model);
  const ES::VXd reference = (ES::VXd(2) << 100, -100).finished();
  energy.setReference(reference, 0.4);
  const ES::VXd base = reference + (ES::VXd(2) << 3, 2).finished();
  const double h = std::ldexp(1.0, -30);
  const ES::VXd trial = base + (ES::VXd(2) << h, -h).finished();
  DeviceBuffer<double> deviceBase, deviceTrial;
  deviceBase.upload(base.data(), base.size());
  deviceTrial.upload(trial.data(), trial.size());
  EXPECT_DOUBLE_EQ(energy.func(deviceBase.data()), energy.func(deviceTrial.data()));
  // The first-order terms are +6h and -6h. The exact remaining increment is
  // scale * (2+3)/2 * h^2, far below the spacing of the absolute energy.
  const auto change = energy.energyChange(deviceBase.data(), deviceTrial.data());
  EXPECT_NEAR(change.value, h * h, 1e-30);
  EXPECT_GE(change.roundoffScale, std::abs(change.value));
  EXPECT_LT(change.roundoffScale, 1e-6);
  EXPECT_DOUBLE_EQ(energy.energyChange(deviceBase.data(), deviceBase.data()).value, 0);
}
}  // namespace
