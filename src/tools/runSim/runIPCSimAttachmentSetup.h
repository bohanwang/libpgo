#pragma once

#include "EigenSupport.h"

#include <memory>
#include <string>
#include <vector>

namespace pgo
{
class ConfigFileJSON;
namespace ConstraintPotentialEnergies
{
class MultipleVertexPulling;
}
}

namespace pgo::RunIPCSim
{
void buildPullingConstraints(const pgo::ConfigFileJSON &jconfig, const std::vector<std::string> &fixedVertexFilenames,
  const EigenSupport::VXd &simulationRestPosition, const EigenSupport::SpMatD &K,
  std::vector<std::shared_ptr<ConstraintPotentialEnergies::MultipleVertexPulling>> &pullingEnergies,
  std::vector<EigenSupport::VXd> &pullingTargets, std::vector<EigenSupport::VXd> &pullingTargetRests);
}  // namespace pgo::RunIPCSim
