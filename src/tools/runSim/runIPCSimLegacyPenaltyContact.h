#pragma once

#include "EigenSupport.h"
#include "runIPCSimContactBackend.h"
#include "triMeshGeo.h"

#include <memory>
#include <vector>

namespace pgo
{
class ConfigFileJSON;
}

namespace pgo::RunIPCSim
{
struct LegacyPenaltyContactConfig
{
  double stiffness = 0.0;
  int samples = 0;
  double frictionCoeff = 0.0;
  double velocityEps = 0.0;
};

LegacyPenaltyContactConfig parseLegacyPenaltyContactConfig(const pgo::ConfigFileJSON &config);

std::shared_ptr<RunIPCSimContactBackend> makeLegacyPenaltyContactBackend(
  const pgo::ConfigFileJSON &config,
  const LegacyPenaltyContactConfig &contactConfig,
  const pgo::Mesh::TriMeshGeo &surfaceMesh,
  const std::vector<int> &embeddingVertexIndices,
  const std::vector<double> &embeddingWeights,
  int simulationDofCount,
  double scale);
}  // namespace pgo::RunIPCSim
