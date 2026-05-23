#pragma once

#include <memory>

namespace pgo
{
class ConfigFileJSON;

namespace Contact::CIPC
{
class EmbeddedSurfaceIPCPotentialEnergy;
}

namespace RunSim
{
class ScopedRunSimCliLogRedirect;
}

namespace Simulation
{
class ImplicitBackwardEulerTimeIntegrator;
}

namespace SolidDeformationModel
{
class DeformationModelEnergy;
}
}  // namespace pgo

namespace pgo::RunIPCSim
{
struct RunIPCSimOptions;
struct RunIPCSimRuntimeConfig;
class RunIPCSimOutput;

class RunIPCSimRunScope
{
public:
  RunIPCSimRunScope(const pgo::ConfigFileJSON &config,
    const RunIPCSimRuntimeConfig &runtimeConfig,
    const RunIPCSimOptions &options,
    const RunIPCSimOutput &output);
  ~RunIPCSimRunScope();

  RunIPCSimRunScope(const RunIPCSimRunScope &) = delete;
  RunIPCSimRunScope &operator=(const RunIPCSimRunScope &) = delete;

  void logProfileSummaryIfEnabled() const;

private:
  bool profilingEnabled_ = false;
  std::unique_ptr<RunSim::ScopedRunSimCliLogRedirect> logRedirect_;
};

void logRunIPCSimMaxStepSummary(
  const std::shared_ptr<pgo::SolidDeformationModel::DeformationModelEnergy> &elasticEnergy,
  const std::shared_ptr<pgo::Contact::CIPC::EmbeddedSurfaceIPCPotentialEnergy> &collisionHandler,
  const std::shared_ptr<pgo::Simulation::ImplicitBackwardEulerTimeIntegrator> &integrator);
}  // namespace pgo::RunIPCSim
