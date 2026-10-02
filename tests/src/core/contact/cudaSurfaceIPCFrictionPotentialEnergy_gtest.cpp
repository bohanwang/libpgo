#include "ipc/cuda/cudaSurfaceIPCFrictionPotentialEnergy.h"
#include "ipc/core/surfaceIPCFriction.h"
#include "cudaDeviceBuffer.h"

#include <gtest/gtest.h>
#include <cmath>
#include <limits>

namespace
{
namespace ES = pgo::EigenSupport;
using pgo::Contact::CIPC::CudaSurfaceIPCFrictionPotentialEnergy;
using pgo::Contact::CIPC::SurfaceIPCFriction;
using pgo::Contact::CIPC::SurfaceIPCCore;
using pgo::CudaUtilities::DeviceBuffer;

class CudaSurfaceIPCFrictionPotentialEnergyGTest : public ::testing::Test
{
protected:
  void SetUp() override
  {
    int count = 0;
    if (cudaGetDeviceCount(&count) != cudaSuccess || count == 0) {
      GTEST_SKIP() << "A CUDA device is required.";
    }
  }
};
struct Mesh
{
  ES::MXd vertices;
  ES::MXi faces;
  ES::SpMatD map;
  std::vector<uint8_t> mask;
  ES::VXd rest() const
  {
    ES::VXd result(3 * vertices.rows());
    for (int i = 0; i < vertices.rows(); ++i)
      result.segment<3>(3 * i) = vertices.row(i).transpose();
    return result;
  }
};
Mesh makeMesh(bool obstacle, bool coupled)
{
  Mesh m;
  m.vertices.resize(6, 3);
  m.vertices << 0, 0, 0, 1, 0, 0, 0, 1, 0,
    .12, .08, .045, .86, .13, .055, .09, .82, .05;
  m.faces.resize(2, 3);
  m.faces << 0, 1, 2, 3, 4, 5;
  m.mask.assign(6, 1);
  if (obstacle)
    std::fill(m.mask.begin(), m.mask.begin() + 3, 0);
  const int n = obstacle ? 9 : 18;
  m.map.resize(18, n);
  std::vector<ES::TripletD> entries;
  for (int row = obstacle ? 9 : 0; row < 18; ++row) {
    const int col = obstacle ? row - 9 : row;
    entries.emplace_back(row, col, coupled ? .7 : 1.0);
    if (coupled)
      entries.emplace_back(row, (col + 4) % n, .3);
  }
  m.map.setFromTriplets(entries.begin(), entries.end());
  return m;
}
SurfaceIPCCore::Parameters parameters(double epsEE = 0.0)
{
  SurfaceIPCCore::Parameters p;
  p.dhat = .1;
  p.kappa = 2.3;
  p.eps_ee = epsEE;
  p.frictionCoeff = .43;
  p.frictionEpsV = .02;
  return p;
}
void updateCpu(SurfaceIPCFriction &friction, const Mesh &mesh, const SurfaceIPCCore::Parameters &p,
  const ES::VXd &reference, const ES::VXd &lagged, double timestep)
{
  SurfaceIPCCore core(p);
  core.setMesh(mesh.vertices, mesh.faces, mesh.mask);
  const ES::VXd worldReference = mesh.rest() + mesh.map * reference;
  const ES::VXd worldLagged = mesh.rest() + mesh.map * lagged;
  core.prepareForSurfacePositions(worldLagged);
  friction.update(worldReference, worldLagged, core.getPTPairs(), core.getEEPairs(),
    p.dhat, p.kappa, p.eps_ee, p.frictionCoeff, p.frictionEpsV, timestep);
}
void nearVector(const ES::VXd &actual, const ES::VXd &expected, double tol = 2e-8)
{
  EXPECT_LE((actual - expected).norm(), tol * std::max(1.0, expected.norm()));
}
void parity(bool obstacle, bool coupled, double epsEE)
{
  auto mesh = makeMesh(obstacle, coupled);
  auto p = parameters(epsEE);
  CudaSurfaceIPCFrictionPotentialEnergy gpu(mesh.vertices, mesh.faces, mesh.map, mesh.mask, p);
  SurfaceIPCFriction cpu;
  const int n = mesh.map.cols();
  const double h = .05;
  ES::VXd reference(n), lagged(n), direction(n), offset(n), actual(n), surfaceGradient(18);
  for (int i = 0; i < n; ++i) {
    reference[i] = .0003 * std::sin(.4 * i);
    lagged[i] = .0002 * std::cos(.6 * i);
    direction[i] = std::cos(.7 * i + .2);
    offset[i] = .03 * std::sin(.3 * i);
  }
  DeviceBuffer<double> dref, dlagged, dx, dv, out;
  dref.upload(reference.data(), n);
  dlagged.upload(lagged.data(), n);
  dv.upload(direction.data(), n);
  gpu.updateFriction(dref.data(), dlagged.data(), h);
  updateCpu(cpu, mesh, p, reference, lagged, h);
  ASSERT_GT(gpu.getNumFrictionPairs(), 0u);
  EXPECT_EQ(gpu.getNumFrictionPairs(), cpu.numPairs());
  for (double slipScale : { 0.0, .0002, .003 }) {
    const ES::VXd x = reference + slipScale * direction;
    const ES::VXd world = mesh.rest() + mesh.map * x;
    dx.upload(x.data(), n);
    EXPECT_NEAR(gpu.func(dx.data()), cpu.computeEnergy(world), 2e-10 * std::max(1.0, cpu.computeEnergy(world)));
    surfaceGradient.setZero();
    cpu.addGradient(world, surfaceGradient);
    const ES::VXd expectedGradient = mesh.map.transpose() * surfaceGradient;
    out.upload(offset.data(), n);
    gpu.addGradient(dx.data(), out.data());
    out.download(actual.data(), n);
    nearVector(actual, offset + expectedGradient);
    ES::SpMatD surfaceHessian;
    cpu.computeHessian(world, surfaceHessian);
    const ES::SpMatD hessian = mesh.map.transpose() * surfaceHessian * mesh.map;
    gpu.prepareHessian(dx.data());
    out.upload(offset.data(), n);
    gpu.addHessianVector(dv.data(), out.data());
    out.download(actual.data(), n);
    nearVector(actual, offset + hessian * direction);
    out.upload(offset.data(), n);
    gpu.addDiagonal(out.data());
    out.download(actual.data(), n);
    nearVector(actual, offset + ES::VXd(hessian.diagonal()));
    EXPECT_GE(direction.dot(hessian * direction), -1e-9);
    if (slipScale == 0.0)
      EXPECT_GT((hessian * direction).norm(), 1.0);
    const double step = 1e-7;
    ES::VXd trial = x + step * direction;
    dx.upload(trial.data(), n);
    const double plus = gpu.func(dx.data());
    trial = x - step * direction;
    dx.upload(trial.data(), n);
    const double minus = gpu.func(dx.data());
    EXPECT_NEAR((plus - minus) / (2 * step), expectedGradient.dot(direction), 2e-5 * std::max(1.0, expectedGradient.norm()));
    // Energy and gradient evaluations must retain the last prepared Hessian.
    gpu.addGradient(dx.data(), out.data());
    out.zero();
    gpu.addHessianVector(dv.data(), out.data());
    out.download(actual.data(), n);
    nearVector(actual, hessian * direction);
  }
  // Re-lagging changes closest geometry and normal force, keeping x_n fixed.
  lagged += .0004 * direction;
  dlagged.upload(lagged.data(), n);
  gpu.updateFriction(dref.data(), dlagged.data(), h);
  updateCpu(cpu, mesh, p, reference, lagged, h);
  dx.upload(reference.data(), n);
  EXPECT_NEAR(gpu.func(dx.data()), cpu.computeEnergy(mesh.rest() + mesh.map * reference), 1e-10);
  out.zero();
  gpu.addGradient(dx.data(), out.data());
  out.download(actual.data(), n);
  EXPECT_EQ(actual.norm(), 0.0);
  EXPECT_THROW(gpu.addDiagonal(out.data()), std::logic_error);
}
TEST_F(CudaSurfaceIPCFrictionPotentialEnergyGTest, MatchesSelfContactAcrossSlipRegimes)
{
  parity(false, false, 0.0);
}
TEST_F(CudaSurfaceIPCFrictionPotentialEnergyGTest, MatchesEmbeddedObstacleAndMollifier)
{
  parity(true, true, .01);
}
TEST_F(CudaSurfaceIPCFrictionPotentialEnergyGTest, MatchesGeneralEmbeddingAndMollifier)
{
  parity(false, true, .01);
}

TEST_F(CudaSurfaceIPCFrictionPotentialEnergyGTest, StableChangesPreserveTinyWorkAndZeroMotion)
{
  auto mesh = makeMesh(true, true);
  auto p = parameters();
  CudaSurfaceIPCFrictionPotentialEnergy gpu(mesh.vertices, mesh.faces, mesh.map, mesh.mask, p);
  const int n = mesh.map.cols();
  ES::VXd reference = ES::VXd::Constant(n, -10.0), lagged = ES::VXd::Zero(n);
  // Large reference displacement tests cancellation without translating lagged geometry.
  DeviceBuffer<double> dr, dl, dx, dt, dg;
  dr.upload(reference.data(), n);
  dl.upload(lagged.data(), n);
  gpu.updateFriction(dr.data(), dl.data(), .05);
  dg.resize(n);
  for (double magnitude : { .0001, .001, .004 }) {
    ES::VXd x = reference, direction(n);
    for (int i = 0; i < n; ++i) {
      direction[i] = std::cos(.4 * i);
      x[i] += magnitude * direction[i];
    }
    dx.upload(x.data(), n);
    const auto zero = gpu.energyChange(dx.data(), dx.data());
    EXPECT_EQ(zero.value, 0.0);
    EXPECT_EQ(zero.roundoffScale, 0.0);
    const ES::VXd trial = x + 1e-11 * direction;
    dt.upload(trial.data(), n);
    dg.zero();
    gpu.addGradient(dx.data(), dg.data());
    ES::VXd gradient(n);
    dg.download(gradient.data(), n);
    const auto change = gpu.energyChange(dx.data(), dt.data());
    EXPECT_NEAR(change.value, gradient.dot(trial - x), 2e-7 * std::max(1e-16, std::abs(change.value)));
    EXPECT_GE(change.roundoffScale, std::abs(change.value));
    const auto reverse = gpu.energyChange(dt.data(), dx.data());
    EXPECT_NEAR(change.value, -reverse.value, 1e-8 * std::max(1e-16, std::abs(change.value)));
  }
}
TEST_F(CudaSurfaceIPCFrictionPotentialEnergyGTest, PrescribedPlacementUpdatesWeightsAndInvalidatesLagging)
{
  auto mesh = makeMesh(true, true);
  auto p = parameters();
  CudaSurfaceIPCFrictionPotentialEnergy gpu(mesh.vertices, mesh.faces, mesh.map, mesh.mask, p);
  ES::VXd reference = ES::VXd::Zero(mesh.map.cols()), lagged = reference, trial = reference;
  trial[0] = .001;
  DeviceBuffer<double> dr, dl, dx;
  dr.upload(reference.data(), reference.size());
  dl.upload(lagged.data(), lagged.size());
  dx.upload(trial.data(), trial.size());
  gpu.updateFriction(dr.data(), dl.data(), .05);
  ES::MXd obstacle = mesh.vertices.topRows(3);
  obstacle.col(0) *= 1.1;
  obstacle.col(2).array() += .003;
  gpu.updateExternalSurface(0, obstacle);
  EXPECT_EQ(gpu.getNumFrictionPairs(), 0u);
  EXPECT_THROW(gpu.func(dx.data()), std::logic_error);
  mesh.vertices.topRows(3) = obstacle;
  SurfaceIPCFriction cpu;
  updateCpu(cpu, mesh, p, reference, lagged, .05);
  gpu.updateFriction(dr.data(), dl.data(), .05);
  EXPECT_NEAR(gpu.func(dx.data()), cpu.computeEnergy(mesh.rest() + mesh.map * trial), 1e-10);
  EXPECT_THROW(gpu.updateExternalSurface(3, obstacle), std::invalid_argument);
}
TEST_F(CudaSurfaceIPCFrictionPotentialEnergyGTest, ValidatesLifecycleAndParametersAndAllowsZeroCoefficient)
{
  auto mesh = makeMesh(true, false);
  auto p = parameters();
  p.frictionCoeff = 0.0;
  CudaSurfaceIPCFrictionPotentialEnergy gpu(mesh.vertices, mesh.faces, mesh.map, mesh.mask, p);
  ES::VXd zero = ES::VXd::Zero(mesh.map.cols());
  DeviceBuffer<double> dx;
  dx.upload(zero.data(), zero.size());
  EXPECT_THROW(gpu.func(dx.data()), std::logic_error);
  EXPECT_THROW(gpu.updateFriction(dx.data(), dx.data(), 0.0), std::invalid_argument);
  gpu.updateFriction(dx.data(), dx.data(), .05);
  EXPECT_EQ(gpu.func(dx.data()), 0.0);
  EXPECT_EQ(gpu.getNumFrictionPairs(), 0u);
  gpu.clear();
  EXPECT_THROW(gpu.func(dx.data()), std::logic_error);
  zero[0] = std::numeric_limits<double>::quiet_NaN();
  dx.upload(zero.data(), zero.size());
  EXPECT_THROW(gpu.updateFriction(dx.data(), dx.data(), .05), std::invalid_argument);
}
TEST_F(CudaSurfaceIPCFrictionPotentialEnergyGTest, EnergyChangesCrossSmoothingBoundaryWithoutInvalidatingHessian)
{
  auto mesh = makeMesh(true, false);
  mesh.vertices << -5, -5, 0, 5, -5, 0, 0, 5, 0,
    -.1, 0, .05, .1, 0, .05, 0, .2, .05;
  auto p = parameters();
  const double timestep = .05, epsilon = p.frictionEpsV * timestep;
  CudaSurfaceIPCFrictionPotentialEnergy gpu(mesh.vertices, mesh.faces, mesh.map, mesh.mask, p);
  SurfaceIPCFriction cpu;
  ES::VXd zero = ES::VXd::Zero(9), base = zero, trial = zero, direction = zero;
  for (int vertex = 0; vertex < 3; ++vertex) {
    base[3 * vertex] = epsilon * (1.0 - 1e-6);
    trial[3 * vertex] = epsilon * (1.0 + 1e-6);
    direction[3 * vertex + 1] = 1.0;
  }
  DeviceBuffer<double> dz, db, dt, dv, out;
  dz.upload(zero.data(), 9);
  db.upload(base.data(), 9);
  dt.upload(trial.data(), 9);
  dv.upload(direction.data(), 9);
  out.resize(9);
  gpu.updateFriction(dz.data(), dz.data(), timestep);
  updateCpu(cpu, mesh, p, zero, zero, timestep);
  ASSERT_EQ(gpu.getNumFrictionPairs(), 3u);
  gpu.prepareHessian(db.data());
  out.zero();
  gpu.addHessianVector(dv.data(), out.data());
  ES::VXd before(9), after(9);
  out.download(before.data(), 9);
  const double expected = cpu.computeEnergy(mesh.rest() + mesh.map * trial) - cpu.computeEnergy(mesh.rest() + mesh.map * base);
  const auto forward = gpu.energyChange(db.data(), dt.data());
  const auto reverse = gpu.energyChange(dt.data(), db.data());
  EXPECT_NEAR(forward.value, expected, 2e-9 * std::abs(expected));
  EXPECT_NEAR(forward.value, -reverse.value, 1e-12 * std::abs(expected));
  const auto unchanged = gpu.energyChange(db.data(), db.data());
  EXPECT_EQ(unchanged.value, 0.0);
  EXPECT_EQ(unchanged.roundoffScale, 0.0);
  out.zero();
  gpu.addHessianVector(dv.data(), out.data());
  out.download(after.data(), 9);
  nearVector(after, before, 1e-13);
}
}  // namespace
