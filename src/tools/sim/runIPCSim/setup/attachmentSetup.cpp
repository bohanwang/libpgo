#include "setup/attachmentSetup.h"

#include "basicIO.h"
#include "configFileJSON.h"
#include "multiVertexPullingSoftConstraints.h"

#include <algorithm>
#include <array>
#include <stdexcept>

namespace pgo::RunIPCSim
{
namespace ES = pgo::EigenSupport;

void buildPullingConstraints(const pgo::ConfigFileJSON &jconfig, const std::vector<std::string> &fixedVertexFilenames,
  const ES::VXd &simulationRestPosition, const ES::SpMatD &K,
  std::vector<std::shared_ptr<ConstraintPotentialEnergies::MultipleVertexPulling>> &pullingEnergies,
  std::vector<ES::VXd> &pullingTargets, std::vector<ES::VXd> &pullingTargetRests)
{
  int fixedVertexFileIndex = 0;
  for (const auto &fv : jconfig.handle()["fixed-vertices"]) {
    const std::string &filename = fixedVertexFilenames.at(fixedVertexFileIndex++);
    const std::array<double, 3> movement = fv["movement"].get<std::array<double, 3>>();
    const double attachmentCoeff = fv["coeff"].get<double>();

    std::vector<int> fixedVertices;
    if (BasicIO::read1DText(filename.c_str(), std::back_inserter(fixedVertices)) != 0) {
      throw std::runtime_error("Failed to read fixed vertex file: " + filename);
    }
    std::sort(fixedVertices.begin(), fixedVertices.end());

    ES::VXd tgtVertexPositions(fixedVertices.size() * 3);
    ES::VXd tgtVertexRests(fixedVertices.size() * 3);
    const ES::V3d movementVec(movement[0], movement[1], movement[2]);
    for (int vi = 0; vi < static_cast<int>(fixedVertices.size()); ++vi) {
      tgtVertexPositions.segment<3>(vi * 3) =
        simulationRestPosition.segment<3>(fixedVertices[vi] * 3) + movementVec;
      tgtVertexRests.segment<3>(vi * 3) =
        simulationRestPosition.segment<3>(fixedVertices[vi] * 3);
    }

    auto pullingEnergy = std::make_shared<ConstraintPotentialEnergies::MultipleVertexPulling>(
      K, simulationRestPosition.data(), static_cast<int>(fixedVertices.size()), fixedVertices.data(), tgtVertexPositions.data(), nullptr, 1);
    pullingEnergy->setCoeff(attachmentCoeff);
    pullingEnergies.push_back(pullingEnergy);
    pullingTargets.push_back(tgtVertexPositions);
    pullingTargetRests.push_back(tgtVertexRests);
  }
}
}  // namespace pgo::RunIPCSim
