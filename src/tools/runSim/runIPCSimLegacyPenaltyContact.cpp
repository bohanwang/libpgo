#include "runIPCSimLegacyPenaltyContact.h"

#include "configFileJSON.h"
#include "implicitBackwardEulerTimeIntegrator.h"
#include "potentialEnergies.h"
#include "pointPenetrationEnergy.h"
#include "pointTrianglePairCouplingEnergyWithCollision.h"
#include "runIPCSimConfig.h"
#include "runIPCSimLogging.h"
#include "runIPCSimSession.h"
#include "runIPCSimSetup.h"
#include "triangleMeshExternalContactHandler.h"
#include "triangleMeshSelfContactHandler.h"

#include <array>
#include <stdexcept>
#include <string>

namespace pgo::RunIPCSim
{
namespace
{
namespace ES = pgo::EigenSupport;

struct LegacyKinematicObject
{
  pgo::Mesh::TriMeshGeo mesh;
  ES::V3d movement = ES::V3d::Zero();
};

std::vector<LegacyKinematicObject> loadLegacyKinematicObjects(const pgo::ConfigFileJSON &config, double scale)
{
  std::vector<LegacyKinematicObject> objects;
  if (!config.exist("external-objects"))
    return objects;

  const auto &externalObjects = config.handle()["external-objects"];
  if (!externalObjects.is_array())
    throw std::invalid_argument("runIPCSim --legacy expects `external-objects` to be an array.");

  objects.reserve(externalObjects.size());
  for (const auto &objectJson : externalObjects) {
    if (!objectJson.is_object())
      throw std::invalid_argument("runIPCSim --legacy expects each `external-objects[]` entry to be an object.");
    if (!objectJson.contains("filename"))
      throw std::invalid_argument("Missing required field `external-objects[].filename`.");
    if (!objectJson.contains("movement"))
      throw std::invalid_argument("Missing required field `external-objects[].movement`.");

    LegacyKinematicObject object;
    const std::string filename = config.resolvePath(objectJson["filename"].get<std::string>());
    if (!object.mesh.load(filename))
      throw std::runtime_error("Failed to load external object mesh: " + filename);
    for (int vi = 0; vi < object.mesh.numVertices(); ++vi)
      object.mesh.pos(vi) *= scale;

    const std::array<double, 3> movement = objectJson["movement"].get<std::array<double, 3>>();
    object.movement = ES::V3d(movement[0], movement[1], movement[2]);
    objects.push_back(std::move(object));
  }

  return objects;
}

std::vector<pgo::Mesh::TriMeshRef> makeObjectRefs(std::vector<LegacyKinematicObject> &objects)
{
  std::vector<pgo::Mesh::TriMeshRef> refs;
  refs.reserve(objects.size());
  for (auto &object : objects)
    refs.emplace_back(object.mesh);
  return refs;
}

class LegacyPenaltyContactBackend final : public RunIPCSimContactBackend
{
public:
  LegacyPenaltyContactBackend(const LegacyPenaltyContactConfig &config,
    const pgo::Mesh::TriMeshGeo &surfaceMesh,
    std::vector<LegacyKinematicObject> objects,
    const std::vector<int> &embeddingVertexIndices,
    const std::vector<double> &embeddingWeights,
    int simulationDofCount):
    config_(config), objects_(std::move(objects))
  {
    if (config_.stiffness <= 0.0)
      return;

    std::vector<pgo::Mesh::TriMeshRef> objectRefs = makeObjectRefs(objects_);
    if (!objectRefs.empty()) {
      externalContactHandler_ = std::make_shared<Contact::TriangleMeshExternalContactHandler>(
        surfaceMesh.positions(), surfaceMesh.triangles(), simulationDofCount,
        objectRefs, config_.samples, &embeddingVertexIndices, &embeddingWeights);
    }

    selfContactHandler_ = std::make_shared<Contact::TriangleMeshSelfContactHandler>(
      surfaceMesh.positions(), surfaceMesh.triangles(), simulationDofCount,
      config_.samples, &embeddingVertexIndices, &embeddingWeights);
  }

  ~LegacyPenaltyContactBackend() override
  {
    releaseActiveBuffers();
  }

  ContactBackendKind kind() const override { return ContactBackendKind::LegacyPenalty; }
  std::string description() const override { return "legacy-penalty"; }

  void initializeAfterRestart(const RunIPCSimRuntimeConfig &runtimeConfig,
    IpcSimulationContext &context, RunIPCSimSession &session) override
  {
    ES::mv(context.surfaceFromSimulationDispMap, session.u, session.usurf);
    const double denom = runtimeConfig.numSimSteps > 1 ? static_cast<double>(runtimeConfig.numSimSteps - 1) : 1.0;
    const double advanceFrames = static_cast<double>(session.frameStart + 1);
    for (std::size_t oi = 0; oi < objects_.size(); ++oi) {
      const ES::V3d movement = objects_[oi].movement / denom * advanceFrames;
      for (int vi = 0; vi < objects_[oi].mesh.numVertices(); ++vi)
        objects_[oi].mesh.pos(vi) += movement;
      updateExternalSurface(oi);
    }
  }

  void beginFrame(int, const RunIPCSimRuntimeConfig &, IpcSimulationContext &, RunIPCSimSession &) override
  {
    releaseActiveBuffers();
  }

  void addForces(int, const RunIPCSimRuntimeConfig &runtimeConfig,
    IpcSimulationContext &context, RunIPCSimSession &session) override
  {
    addExternalContactForce(runtimeConfig, context, session);
    addSelfContactForce(runtimeConfig, context, session);
  }

  void afterStep(int, const RunIPCSimRuntimeConfig &runtimeConfig,
    IpcSimulationContext &context, RunIPCSimSession &session) override
  {
    releaseActiveBuffers();
    ES::mv(context.surfaceFromSimulationDispMap, session.u, session.usurf);

    const double denom = runtimeConfig.numSimSteps > 1 ? static_cast<double>(runtimeConfig.numSimSteps - 1) : 1.0;
    for (std::size_t oi = 0; oi < objects_.size(); ++oi) {
      const ES::V3d movement = objects_[oi].movement / denom;
      for (int vi = 0; vi < objects_[oi].mesh.numVertices(); ++vi)
        objects_[oi].mesh.pos(vi) += movement;
      updateExternalSurface(oi);
    }
  }

  void addStaticEnergies(const RunIPCSimRuntimeConfig &,
    IpcSimulationContext &, NonlinearOptimization::PotentialEnergies &) override
  {
    // Old runSim static mode did not add legacy penalty contact energies.
  }

  void logSummary(const IpcSimulationContext &context, const RunIPCSimSession &session) const override
  {
    logRunIPCSimMaxStepSummary(context.elasticEnergy, nullptr, session.integrator);
  }

private:
  void releaseActiveBuffers()
  {
    if (activeExternalEnergy_ && activeExternalBuffer_) {
      activeExternalEnergy_->freeBuffer(activeExternalBuffer_);
      activeExternalEnergy_->setBuffer(nullptr);
    }
    activeExternalBuffer_ = nullptr;
    activeExternalEnergy_.reset();

    if (activeSelfEnergy_ && activeSelfBuffer_) {
      activeSelfEnergy_->freeBuffer(activeSelfBuffer_);
      activeSelfEnergy_->setBuffer(nullptr);
    }
    activeSelfBuffer_ = nullptr;
    activeSelfEnergy_.reset();
  }

  void updateExternalSurface(std::size_t objectIndex)
  {
    if (!externalContactHandler_)
      return;
    externalContactHandler_->updateExternalSurface(
      static_cast<int>(objectIndex), pgo::Mesh::TriMeshRef(objects_[objectIndex].mesh));
  }

  void addExternalContactForce(const RunIPCSimRuntimeConfig &runtimeConfig,
    IpcSimulationContext &context, RunIPCSimSession &session)
  {
    if (!externalContactHandler_)
      return;

    externalContactHandler_->execute(session.usurf.data());
    if (externalContactHandler_->getNumCollidingSamples() <= 0)
      return;

    activeExternalEnergy_ = externalContactHandler_->buildContactEnergy();
    activeExternalBuffer_ = activeExternalEnergy_->allocateBuffer();

    auto posFunc = [&context](const ES::V3d &u, ES::V3d &p, int dofStart) {
      p = u + context.simulationRestPosition.segment<3>(dofStart);
    };

    auto lastPosFunc = [&context, &session](const ES::V3d &, ES::V3d &p, int dofStart) {
      p = session.u.segment<3>(dofStart) + context.simulationRestPosition.segment<3>(dofStart);
    };

    activeExternalEnergy_->setComputePosFunction(posFunc);
    activeExternalEnergy_->setComputeLastPosFunction(lastPosFunc);
    activeExternalEnergy_->setBuffer(activeExternalBuffer_);
    activeExternalEnergy_->setCoeff(config_.stiffness);
    activeExternalEnergy_->setFrictionCoeff(config_.frictionCoeff);
    activeExternalEnergy_->setTimestep(runtimeConfig.timestep);
    activeExternalEnergy_->setVelEps(config_.velocityEps);
    session.integrator->addGeneralImplicitForceModel(activeExternalEnergy_, 0, 0);
  }

  void addSelfContactForce(const RunIPCSimRuntimeConfig &runtimeConfig,
    IpcSimulationContext &context, RunIPCSimSession &session)
  {
    if (!selfContactHandler_)
      return;

    selfContactHandler_->execute(session.usurf.data());
    if (selfContactHandler_->getCollidingTrianglePair().empty())
      return;

    selfContactHandler_->handleContactDCD(0, 100);
    activeSelfEnergy_ = selfContactHandler_->buildContactEnergy();
    activeSelfEnergy_->setToPosFunction([&context](const ES::V3d &x, ES::V3d &p, int offset) {
      p = x + context.simulationRestPosition.segment<3>(offset);
    });
    activeSelfEnergy_->setToLastPosFunction([&context, &session](const ES::V3d &, ES::V3d &p, int offset) {
      p = context.simulationRestPosition.segment<3>(offset) + session.u.segment<3>(offset);
    });

    activeSelfBuffer_ = activeSelfEnergy_->allocateBuffer();
    activeSelfEnergy_->setBuffer(activeSelfBuffer_);
    activeSelfEnergy_->setCoeff(config_.stiffness);
    activeSelfEnergy_->computeClosestPosition(session.u.data());
    activeSelfEnergy_->setFrictionCoeff(config_.frictionCoeff);
    activeSelfEnergy_->setTimestep(runtimeConfig.timestep);
    activeSelfEnergy_->setVelEps(config_.velocityEps);
    session.integrator->addGeneralImplicitForceModel(activeSelfEnergy_, 0, 0);
  }

  LegacyPenaltyContactConfig config_;
  std::vector<LegacyKinematicObject> objects_;
  std::shared_ptr<Contact::TriangleMeshExternalContactHandler> externalContactHandler_;
  std::shared_ptr<Contact::TriangleMeshSelfContactHandler> selfContactHandler_;
  std::shared_ptr<Contact::PointPenetrationEnergy> activeExternalEnergy_;
  Contact::PointPenetrationEnergyBuffer *activeExternalBuffer_ = nullptr;
  std::shared_ptr<Contact::PointTrianglePairCouplingEnergyWithCollision> activeSelfEnergy_;
  Contact::PointTrianglePairCouplingEnergyWithCollisionBuffer *activeSelfBuffer_ = nullptr;
};
}  // namespace

LegacyPenaltyContactConfig parseLegacyPenaltyContactConfig(const pgo::ConfigFileJSON &config)
{
  const bool hasContactSample = config.exist("contact-sample");
  const bool hasContactSamples = config.exist("contact-samples");
  if (hasContactSample && hasContactSamples)
    throw std::invalid_argument("runIPCSim --legacy accepts either `contact-sample` or `contact-samples`, not both.");

  LegacyPenaltyContactConfig parsed;
  parsed.stiffness = config.getDouble("contact-stiffness", 1);
  parsed.samples = hasContactSamples ? config.getInt("contact-samples", 1) : config.getInt("contact-sample", 1);
  parsed.frictionCoeff = config.getDouble("contact-friction-coeff", 1);
  parsed.velocityEps = config.getDouble("contact-vel-eps", 1);

  if (parsed.stiffness < 0.0)
    throw std::invalid_argument("runIPCSim --legacy requires non-negative `contact-stiffness`.");
  if (parsed.samples <= 0)
    throw std::invalid_argument("runIPCSim --legacy requires positive contact sample count.");
  if (parsed.velocityEps <= 0.0)
    throw std::invalid_argument("runIPCSim --legacy requires positive `contact-vel-eps`.");

  return parsed;
}

std::shared_ptr<RunIPCSimContactBackend> makeLegacyPenaltyContactBackend(
  const pgo::ConfigFileJSON &config,
  const LegacyPenaltyContactConfig &contactConfig,
  const pgo::Mesh::TriMeshGeo &surfaceMesh,
  const std::vector<int> &embeddingVertexIndices,
  const std::vector<double> &embeddingWeights,
  int simulationDofCount,
  double scale)
{
  return std::make_shared<LegacyPenaltyContactBackend>(
    contactConfig, surfaceMesh, loadLegacyKinematicObjects(config, scale),
    embeddingVertexIndices, embeddingWeights, simulationDofCount);
}
}  // namespace pgo::RunIPCSim
