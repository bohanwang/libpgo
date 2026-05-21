#include "runIPCSimSetup.h"

#include "barycentricCoordinates.h"
#include "configFileJSON.h"
#include "deformationModelEnergy.h"
#include "generateMassMatrix.h"
#include "runIPCSimAttachmentSetup.h"
#include "runIPCSimLegacyPenaltyContact.h"
#include "runIPCSimSetupCommon.h"
#include "runSimFEMSetup.h"
#include "runSimVolumeMeshIO.h"
#include "volumetricMesh.h"

#include <iostream>
#include <memory>

namespace pgo::RunIPCSim
{
namespace ES = pgo::EigenSupport;

IpcSimulationContext buildVolumeLegacyPenaltySimulation(const pgo::ConfigFileJSON &jconfig)
{
  validateZeroInitialDisplacement(jconfig);
  if (!jconfig.exist("tet-mesh") && !jconfig.exist("cubic-mesh"))
    throwConfigError("runIPCSim --legacy only supports volume legacy configs with `tet-mesh` or `cubic-mesh`; shell legacy has been removed.");
  if (jconfig.exist("tet-mesh") && jconfig.exist("cubic-mesh"))
    throwConfigError("runIPCSim --legacy expects exactly one of `tet-mesh` or `cubic-mesh`.");
  if (!jconfig.exist("fixed-vertices"))
    throwConfigError("Missing required field `fixed-vertices`.");

  const double scale = jconfig.getDouble("scale", 1);
  const SolidDeformationModel::DeformationModelElasticMaterial elasticMat = parseVolumeElasticMaterial(jconfig);
  const bool enableMaterialMaxStep = parseEnableMaterialMaxStep(jconfig);
  const RunSim::ResolvedRunSimPaths resolvedPaths = RunSim::resolveRunSimPaths(jconfig);
  std::unique_ptr<VolumetricMeshes::VolumetricMesh> volumetricMesh =
    RunSim::loadValidatedVolumeMesh(RunSim::parseVolumeMeshInputConfig(jconfig), scale);

  pgo::Mesh::TriMeshGeo surfaceMesh;
  ES::VXd surfaceRestPositions;
  loadSurfaceMeshAndRestPositions(resolvedPaths.surfaceMeshFilename, scale, surfaceMesh, surfaceRestPositions);

  pgo::InterpolationCoordinates::BarycentricCoordinates bc(
    surfaceMesh.numVertices(), surfaceRestPositions.data(), volumetricMesh.get());
  ES::SpMatD W = bc.generateInterpolationMatrix();
  if (W.rows() != surfaceMesh.numVertices() * 3)
    throwConfigError("runIPCSim --legacy volume setup produced an embedding matrix with unexpected row count.");
  if (W.cols() != volumetricMesh->getNumVertices() * 3)
    throwConfigError("runIPCSim --legacy volume setup produced an embedding matrix with unexpected column count.");
  if (W.nonZeros() <= 0)
    throwConfigError("runIPCSim --legacy volume setup produced an empty embedding matrix.");

  ES::SpMatD M;
  VolumetricMeshes::GenerateMassMatrix::computeMassMatrix(volumetricMesh.get(), M, true);

  RunSim::InitializedVolumetricSimulation initialized =
    RunSim::initializeVolumetricSimulation(*volumetricMesh, elasticMat,
      SolidDeformationModel::DeformationModelPlasticMaterial::VOLUMETRIC_DOF6,
      enableMaterialMaxStep);
  if (W.cols() != initialized.restPosition.size())
    throwConfigError("runIPCSim --legacy volume setup produced an embedding matrix incompatible with simulation DOFs.");

  ES::VXd zero = ES::VXd::Zero(initialized.restPosition.size());
  ES::SpMatD K;
  initialized.elasticEnergy->createHessian(K);
  initialized.elasticEnergy->hessian(zero, K);

  std::vector<std::shared_ptr<ConstraintPotentialEnergies::MultipleVertexPulling>> pullingEnergies;
  std::vector<ES::VXd> pullingTargets;
  std::vector<ES::VXd> pullingTargetRests;
  buildPullingConstraints(jconfig, resolvedPaths.fixedVertexFilenames, initialized.restPosition, K,
    pullingEnergies, pullingTargets, pullingTargetRests);

  LegacyPenaltyContactConfig legacyContactConfig = parseLegacyPenaltyContactConfig(jconfig);
  std::cout << "runIPCSim legacy volume penalty contact parameters: "
            << "contact-stiffness=" << legacyContactConfig.stiffness << ", "
            << "contact-samples=" << legacyContactConfig.samples << ", "
            << "contact-friction-coeff=" << legacyContactConfig.frictionCoeff << ", "
            << "contact-vel-eps=" << legacyContactConfig.velocityEps << std::endl;

  IpcSimulationContext context;
  context.M = std::move(M);
  context.simulationRestPosition = std::move(initialized.restPosition);
  context.surfaceRestPositions = std::move(surfaceRestPositions);
  context.plasticParams = std::move(initialized.plasticity);
  context.surfaceFromSimulationDispMap = std::move(W);
  context.simulationMeshOwner = initialized.simMesh;
  context.deformationModelManagerOwner = initialized.dmm;
  context.deformationModelAssemblerOwner = initialized.assembler;
  context.elasticEnergy = initialized.elasticEnergy;
  context.pullingEnergies = std::move(pullingEnergies);
  context.pullingTargets = std::move(pullingTargets);
  context.pullingTargetRests = std::move(pullingTargetRests);
  context.surfaceMesh = std::move(surfaceMesh);
  context.contactBackend = makeLegacyPenaltyContactBackend(jconfig, legacyContactConfig,
    context.surfaceMesh, bc.getEmbeddingVertexIndices(), bc.getEmbeddingWeights(),
    static_cast<int>(context.simulationRestPosition.size()), scale);
  return context;
}
}  // namespace pgo::RunIPCSim
