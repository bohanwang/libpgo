#pragma once

// Test-only measurements: compare the same discrete elastic model, with exact
// and CPU-matching projected operators kept separate in the report.
#include "cudaDeviceBuffer.h"
#include "cudaPotentialEnergy.h"
#include "deformationModelEnergy.h"
#include "deformationModelManager.h"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <string>
#include <vector>

namespace pgo::Tests::ElasticEquivalence
{
namespace ES = EigenSupport;
using Json = nlohmann::json;
using CudaUtilities::DeviceBuffer;

inline Json error(ES::ConstRefVecXd actual, ES::ConstRefVecXd expected,
  double absoluteTolerance, double relativeTolerance)
{
  const bool finite = actual.allFinite() && expected.allFinite();
  const double absolute = finite ? (actual - expected).cwiseAbs().maxCoeff() : 0;
  const double referenceInfinity = finite ? expected.cwiseAbs().maxCoeff() : 0;
  const double referenceL2 = finite ? expected.norm() : 0;
  const double differenceL2 = finite ? (actual - expected).norm() : 0;
  const double threshold = absoluteTolerance + relativeTolerance * referenceInfinity;
  return { { "finite", finite }, { "passed", finite && absolute <= threshold },
    { "maximum_absolute_error", finite ? Json(absolute) : Json(nullptr) },
    { "relative_infinity_error", finite && referenceInfinity > 0 ? Json(absolute / referenceInfinity) : Json(nullptr) },
    { "relative_l2_error", finite && referenceL2 > 0 ? Json(differenceL2 / referenceL2) : Json(nullptr) },
    { "reference_infinity_norm", referenceInfinity }, { "reference_l2_norm", referenceL2 },
    { "absolute_tolerance", absoluteTolerance }, { "relative_tolerance", relativeTolerance },
    { "acceptance_threshold", threshold },
    { "criterion", "max_abs_error <= absolute_tolerance + relative_tolerance * reference_infinity_norm" } };
}

inline Json timed(const std::function<void()> &operation, const std::function<void()> &synchronize)
{
  using Clock = std::chrono::steady_clock;
  constexpr int samples = 5, maximumCalls = 4096;
  constexpr double targetSeconds = 0.1;
  auto batch = [&](int calls) {
    synchronize();
    const auto start = Clock::now();
    for (int i = 0; i < calls; ++i)
      operation();
    synchronize();
    return std::chrono::duration<double>(Clock::now() - start).count();
  };
  operation();
  operation();
  synchronize();
  int calls = 1;
  double pilot = batch(calls);
  while (pilot < targetSeconds && calls < maximumCalls) {
    calls *= 2;
    pilot = batch(calls);
  }
  std::vector<double> elapsed, perCall;
  for (int i = 0; i < samples; ++i) {
    elapsed.push_back(batch(calls));
    perCall.push_back(elapsed.back() / calls);
  }
  std::sort(perCall.begin(), perCall.end());
  return { { "samples", samples }, { "calls_per_sample", calls },
    { "warmup_calls", 2 }, { "target_seconds_per_sample", targetSeconds },
    { "maximum_calls_per_sample", maximumCalls },
    { "call_cap_reached", calls == maximumCalls }, { "calibration_batch_seconds", pilot },
    { "sample_seconds", elapsed }, { "median_seconds_per_call", perCall[samples / 2] },
    { "minimum_seconds_per_call", perCall.front() }, { "maximum_seconds_per_call", perCall.back() } };
}

inline Json compare(SolidDeformationModel::DeformationModelEnergy &cpu,
  SolidDeformationModel::DeformationModelManager &manager,
  NonlinearOptimization::CudaPotentialEnergy &gpu, ES::ConstRefVecXd state,
  bool fullColumns, bool includeTimings = false, int probeCount = 8)
{
  const int n = state.size();
  if (n == 0 || cpu.getNumDOFs() != n || gpu.getNumDOFs() != n || !state.allFinite())
    throw std::invalid_argument("Elastic equivalence requires matching, nonempty finite states.");
  DeviceBuffer<double> deviceState, deviceDirection(n), deviceResult(n);
  deviceState.upload(state.data(), n);
  ES::VXd expected(n), actual(n), direction(n);
  const double cpuEnergy = cpu.func(state), gpuEnergy = gpu.func(deviceState.data());
  ES::VXd cpuScalar(1), gpuScalar(1);
  cpuScalar[0] = cpuEnergy;
  gpuScalar[0] = gpuEnergy;
  Json report = {
    { "num_dofs", n }, { "comparison", "double-precision numerical equivalence, not bitwise identity" },
    { "energy", error(gpuScalar, cpuScalar, 1e-10, 5e-10) },
    { "cpu_energy", cpuEnergy }, { "cuda_energy", gpuEnergy },
    { "timings_enabled", includeTimings },
    { "hessian_coverage", fullColumns ? "all coordinate columns reconstructed through CUDA HVPs" : "deterministic matrix-free HVP probes; not a full matrix comparison" }
  };
  cpu.gradient(state, expected);
  deviceResult.zero();
  gpu.addGradient(deviceState.data(), deviceResult.data());
  deviceResult.download(actual.data(), n);
  report["gradient"] = error(actual, expected, 1e-9, 1e-9);
  bool passed = report["energy"]["passed"].get<bool>() && report["gradient"]["passed"].get<bool>();
  auto synchronize = [] {
    CudaUtilities::checkCuda(cudaDeviceSynchronize(), "Synchronize elastic measurement");
  };
  auto cpuSynchronize = [] {
  };
  if (includeTimings) {
    volatile double energySink = 0;
    report["timings"]["cpu_energy"] = timed([&] { energySink = cpu.func(state); }, cpuSynchronize);
    report["timings"]["cuda_energy"] = timed([&] { energySink = gpu.func(deviceState.data()); }, synchronize);
    report["timings"]["cpu_gradient"] = timed([&] { cpu.gradient(state, expected); }, cpuSynchronize);
    report["timings"]["cuda_gradient_with_output_zero"] = timed([&] {
      deviceResult.zero();
      gpu.addGradient(deviceState.data(), deviceResult.data());
    },
      synchronize);
    report["timings"]["host_to_device_one_vector"] = timed([&] { deviceState.upload(state.data(), n); }, synchronize);
    report["timings"]["device_to_host_one_vector"] = timed([&] { deviceResult.download(actual.data(), n); }, synchronize);
    report["timing_scope"] = "Resident input/output buffers. GPU regions synchronize before and after every measured batch; gradient/HVP include output clearing. Upload/download are separate one-vector measurements. CPU Hessian assembly and GPU cached preparation are reported separately from applications.";
    report["transfer_bytes_per_call"] = sizeof(double) * n;
    (void)energySink;
  }
  ES::SpMatD hessian;
  cpu.createHessian(hessian);
  for (bool projected : { false, true }) {
    const std::string label = projected ? "cpu_matching_projected" : "exact_unprojected";
    manager.setEnforceSPD(projected ? 1 : 0);
    cpu.hessian(state, hessian);
    if (projected)
      gpu.prepareHessian(deviceState.data());
    else
      gpu.prepareExactHessian(deviceState.data());
    Json result = { { "cpu_sparse_nonzeros", hessian.nonZeros() },
      { "directions", Json::array() }, { "direction_count", fullColumns ? n : probeCount } };
    double maxAbsolute = 0, maxRelativeInfinity = 0, maxRelativeL2 = 0;
    double sumSquaredError = 0, sumSquaredReference = 0, maxReference = 0;
    bool operatorPassed = true;
    for (int probe = 0; probe < (fullColumns ? n : probeCount); ++probe) {
      if (fullColumns) {
        direction.setZero();
        direction[probe] = 1;
      }
      else {
        for (int i = 0; i < n; ++i)
          direction[i] = std::sin(0.17 * (i + 1) * (probe + 1)) + 0.3 * std::cos(0.13 * i + 0.7 * probe);
        direction.normalize();
      }
      expected = hessian * direction;
      deviceDirection.upload(direction.data(), n);
      deviceResult.zero();
      gpu.addHessianVector(deviceDirection.data(), deviceResult.data());
      deviceResult.download(actual.data(), n);
      Json measurement = error(actual, expected, 1e-8, projected ? 5e-8 : 3e-9);
      measurement["direction_index"] = probe;
      operatorPassed = operatorPassed && measurement["passed"].get<bool>();
      if (measurement["finite"].get<bool>()) {
        maxAbsolute = std::max(maxAbsolute, measurement["maximum_absolute_error"].get<double>());
        if (!measurement["relative_infinity_error"].is_null())
          maxRelativeInfinity = std::max(maxRelativeInfinity, measurement["relative_infinity_error"].get<double>());
        if (!measurement["relative_l2_error"].is_null())
          maxRelativeL2 = std::max(maxRelativeL2, measurement["relative_l2_error"].get<double>());
        sumSquaredError += (actual - expected).squaredNorm();
        sumSquaredReference += expected.squaredNorm();
        maxReference = std::max(maxReference, expected.cwiseAbs().maxCoeff());
      }
      result["directions"].push_back(measurement);
    }
    result["passed"] = operatorPassed;
    result["maximum_absolute_error"] = maxAbsolute;
    result["maximum_direction_relative_infinity_error"] = maxRelativeInfinity;
    result["maximum_direction_relative_l2_error"] = maxRelativeL2;
    if (fullColumns) {
      result["full_matrix_relative_frobenius_error"] = sumSquaredReference > 0 ? Json(std::sqrt(sumSquaredError / sumSquaredReference)) : Json(nullptr);
      result["full_matrix_maximum_entry_relative_error"] = maxReference > 0 ? Json(maxAbsolute / maxReference) : Json(nullptr);
    }
    expected = hessian.diagonal();
    deviceResult.zero();
    gpu.addDiagonal(deviceResult.data());
    deviceResult.download(actual.data(), n);
    result["diagonal"] = error(actual, expected, 1e-8, projected ? 5e-8 : 3e-9);
    passed = passed && operatorPassed && result["diagonal"]["passed"].get<bool>();
    if (includeTimings) {
      result["timings"]["cpu_sparse_hessian_assembly"] = timed([&] { cpu.hessian(state, hessian); }, cpuSynchronize);
      result["timings"]["cuda_cached_hessian_preparation"] = timed([&] {
        if (projected)
          gpu.prepareHessian(deviceState.data());
        else
          gpu.prepareExactHessian(deviceState.data());
      },
        synchronize);
      result["timings"]["cpu_sparse_matrix_vector"] = timed([&] { expected.noalias() = hessian * direction; }, cpuSynchronize);
      result["timings"]["cuda_matrix_free_hessian_vector_with_output_zero"] = timed([&] {
        deviceResult.zero();
        gpu.addHessianVector(deviceDirection.data(), deviceResult.data());
      },
        synchronize);
    }
    report[label] = std::move(result);
  }
  manager.setEnforceSPD(1);
  report["passed"] = passed;
  return report;
}
}  // namespace pgo::Tests::ElasticEquivalence
