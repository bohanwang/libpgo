#pragma once

#include <memory>
#include <string>

namespace pgo::Simulation
{
class ImplicitBackwardEulerTimeIntegrator;
}

namespace pgo::NonlinearOptimization
{
class PotentialEnergies;
}

namespace pgo::RunIPCSim
{
struct IpcSimulationContext;
struct RunIPCSimRuntimeConfig;
struct RunIPCSimSession;

enum class ContactBackendKind
{
  Ipc,
  LegacyPenalty,
};

class RunIPCSimContactBackend
{
public:
  virtual ~RunIPCSimContactBackend() = default;

  virtual ContactBackendKind kind() const = 0;
  virtual std::string description() const = 0;

  virtual void initializeAfterRestart(const RunIPCSimRuntimeConfig &runtimeConfig,
    IpcSimulationContext &context, RunIPCSimSession &session) = 0;

  virtual void beginFrame(int frame, const RunIPCSimRuntimeConfig &runtimeConfig,
    IpcSimulationContext &context, RunIPCSimSession &session) = 0;

  virtual void addForces(int frame, const RunIPCSimRuntimeConfig &runtimeConfig,
    IpcSimulationContext &context, RunIPCSimSession &session) = 0;

  virtual void afterStep(int frame, const RunIPCSimRuntimeConfig &runtimeConfig,
    IpcSimulationContext &context, RunIPCSimSession &session) = 0;

  virtual void addStaticEnergies(const RunIPCSimRuntimeConfig &runtimeConfig,
    IpcSimulationContext &context, NonlinearOptimization::PotentialEnergies &energyAll) = 0;

  virtual void logSummary(const IpcSimulationContext &context,
    const RunIPCSimSession &session) const = 0;
};
}  // namespace pgo::RunIPCSim
