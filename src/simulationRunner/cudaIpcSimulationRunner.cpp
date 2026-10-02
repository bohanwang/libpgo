#include "cudaIpcSimulationRunner.h"

#include "configFileJSON.h"
#include "cudaDeformationModelEnergy.h"
#include "cudaKoiterDeformationModelEnergy.h"
#include "cudaDeviceBuffer.h"
#include "cudaFrozenHessianEnergy.h"
#include "cudaEmbeddedSurfaceIPCPotentialEnergy.h"
#include "cudaSurfaceIPCFrictionPotentialEnergy.h"
#include "cudaInit.h"
#include "cudaNewtonSolver.h"
#include "cudaQuadraticEnergy.h"
#include "deformationModelEnergy.h"
#include "deformationModelManager.h"
#include "EigenSupport.h"
#include "initPredicates.h"
#include "multiVertexPullingSoftConstraints.h"
#include "pgoLogging.h"
#include "runIPCSimSetup.h"
#include "runSimCliLogging.h"
#include "runSimTiming.h"

#include <fmt/format.h>
#include <nlohmann/json.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace
{
namespace ES = pgo::EigenSupport;
namespace NO = pgo::NonlinearOptimization;
using Context = pgo::RunIPCSim::IpcSimulationContext;

struct NewtonResult
{
  NO::CudaNewtonSolver::Diagnostics diagnostics;
  std::size_t frictionPairs = 0;
  int status = NO::CudaNewtonSolver::SOLVE_NOT_CONVERGED;
};

struct SubstepResult
{
  std::vector<NewtonResult> solves;
  int maxNewtonIterations = 0;
  ES::MXd state;
  NO::CudaNewtonSolver::Diagnostics diagnostics;
  std::size_t activePairs = 0;
  int status = NO::CudaNewtonSolver::SOLVE_NOT_CONVERGED;
};

// Counts cover every solve; terminal objective diagnostics refer to the last
// solve. The maximum gradient ensures that every lagging iteration satisfies
// the configured convergence tolerance.
void appendDiagnostics(NO::CudaNewtonSolver::Diagnostics &total,
  const NO::CudaNewtonSolver::Diagnostics &current)
{
  total.newtonIterations += current.newtonIterations;
  total.pcgIterations += current.pcgIterations;
  total.regularizationAttempts += current.regularizationAttempts;
  total.roundoffAcceptedSteps += current.roundoffAcceptedSteps;
  total.exactHessianAttempts += current.exactHessianAttempts;
  total.exactHessianDirections += current.exactHessianDirections;
  total.exactHessianFallbacks += current.exactHessianFallbacks;
  total.inputChordRejections += current.inputChordRejections;
  total.trialChordRejections += current.trialChordRejections;
  total.stagnationTriggers += current.stagnationTriggers;
  total.residualCorrectionAttempts += current.residualCorrectionAttempts;
  total.residualCorrectionCandidates += current.residualCorrectionCandidates;
  total.residualCorrectionSteps += current.residualCorrectionSteps;
  total.maxRegularization = std::max(total.maxRegularization, current.maxRegularization);
  total.finalGradientInfinityNorm = std::max(total.finalGradientInfinityNorm, current.finalGradientInfinityNorm);
  total.finalEnergy = current.finalEnergy;
  total.lastLinearResidualNorm = current.lastLinearResidualNorm;
  total.lastFeasibleAlpha = current.lastFeasibleAlpha;
  total.lastAcceptedAlpha = current.lastAcceptedAlpha;
  total.stopReason = current.stopReason;
}

// The contact model certifies the entire original chord with interval bounds.
// A safe prefix is insufficient: reconstructing new endpoints from a prefix
// would lose the original interpolation uncertainty near a grazing contact.
bool certifyCudaChord(const std::array<NO::CudaPotentialEnergy *, 2> &energies, ES::ConstRefVecXd from,
  ES::ConstRefVecXd to, pgo::CudaUtilities::DeviceBuffer<double> &deviceFrom,
  pgo::CudaUtilities::DeviceBuffer<double> &deviceDelta)
{
  const ES::VXd delta = to - from;
  deviceFrom.upload(from.data(), from.size());
  deviceDelta.upload(delta.data(), delta.size());
  double bound = 1;
  for (const auto *energy : energies) {
    const double candidate = energy->computeMaxStepSize(deviceFrom.data(), deviceDelta.data());
    if (!std::isfinite(candidate) || candidate < 0 || candidate > 1)
      return false;
    bound = std::min(bound, candidate);
  }
  return bound == 1.0;
}

void writeStateFile(const std::filesystem::path &filename, const ES::MXd &state)
{
  if (ES::writeMatrix(filename.string().c_str(), state) != 0)
    throw std::runtime_error("Failed to write CUDA IPC state: " + filename.string());
}

void writeSurfaceFile(const Context &context, ES::ConstRefVecXd u,
  const std::filesystem::path &filename, double scale)
{
  const ES::VXd positions = context.surfaceRestPositions + context.surfaceFromSimulationDispMap * u;
  pgo::Mesh::TriMeshGeo surface = context.surfaceMesh;
  for (int i = 0; i < surface.numVertices(); ++i)
    surface.pos(i) = positions.segment<3>(3 * i) / scale;
  if (!surface.save(filename.string()))
    throw std::runtime_error("Failed to write CUDA IPC surface: " + filename.string());
}

void saveCudaState(const Context &context, ES::ConstRefVecXd u,
  ES::ConstRefVecXd velocity, ES::ConstRefVecXd acceleration,
  const std::filesystem::path &output, int frame, int dumpInterval, double scale)
{
  ES::MXd state(u.size(), 3);
  state.col(0) = u;
  state.col(1) = velocity;
  state.col(2) = acceleration;
  const auto stateFile = output / fmt::format("deform{:04d}.u", frame);
  writeStateFile(stateFile, state);

  if (frame % dumpInterval == 0) {
    const auto meshFile = output / fmt::format("ret{:04d}.obj", frame / dumpInterval);
    writeSurfaceFile(context, u, meshFile, scale);
  }
}

// Pulling energies are quadratic. Extract their coefficients when targets
// change, then evaluate them entirely on the device throughout the solve.
void addPullingQuadratic(const Context &context, ES::ConstRefVecXd reference,
  ES::SpMatD &matrix, ES::VXd &linear, double &constant)
{
  auto add = [&](const auto &energy) {
    ES::SpMatD hessian;
    energy->createHessian(hessian);
    energy->hessian(reference, hessian);
    matrix += hessian;
    ES::VXd gradient(reference.size());
    energy->gradient(reference, gradient);
    linear += gradient;
    constant += energy->func(reference);
  };
  for (const auto &energy : context.pullingEnergies)
    add(energy);
  for (const auto &pulling : context.trajectoryPullings)
    add(pulling.energy);
}
}  // namespace

int pgo::SimulationRunner::runCudaIPCSimulation(const ConfigFileJSON &config, bool enableCliLog,
  const IPCSimulationRunOptions &options)
{
  std::unique_ptr<RunSim::ScopedRunSimCliLogRedirect> logRedirect;
  try {
    RunSim::IPCSimulationTiming timing(options);
    const std::string materialName = config.getString("elastic-material", 1);
    const bool shell = !config.exist("tet-mesh") && !config.exist("cubic-mesh");
    const std::string cubicValidityName = config.getValue<std::string>("cuda-cubic-validity", 0, "gauss");
    if (cubicValidityName != "gauss" && cubicValidityName != "whole-cell")
      throw std::invalid_argument("`cuda-cubic-validity` must be `gauss` or `whole-cell`.");
    if (config.exist("cuda-cubic-validity") && !config.exist("cubic-mesh"))
      throw std::invalid_argument("`cuda-cubic-validity` requires a cubic mesh.");
    using CubicValidity = SolidDeformationModel::CudaDeformationModelEnergy::CubicValidity;
    const auto cubicValidity = cubicValidityName == "gauss" ? CubicValidity::Gauss : CubicValidity::WholeCell;
    if ((shell && materialName != "koiter-stvk") ||
      (!shell && materialName != "stable-neo" && materialName != "stvk-vol"))
      throw std::invalid_argument("CUDA IPC supports `koiter-stvk` shells and `stable-neo` or `stvk-vol` tetrahedral/cubic solids.");
    const double frictionCoefficient = config.getValue<double>("contact-friction-coeff", 0, 0.0);
    if (!std::isfinite(frictionCoefficient) || frictionCoefficient < 0)
      throw std::invalid_argument("CUDA IPC requires a finite, nonnegative friction coefficient.");
    const bool frictionEnabled = frictionCoefficient > 0;
    const int frictionIterations = config.getValue<int>("ipc-friction-iterations", 0, 1);
    if (frictionIterations < 1)
      throw std::invalid_argument("`ipc-friction-iterations` must be at least one.");
    if (config.getValue<bool>("deterministic-assembly", 0, false))
      throw std::invalid_argument("CUDA IPC deterministic assembly is not implemented yet.");
    if (config.exist("linear-solver"))
      throw std::invalid_argument("CUDA IPC uses GPU PCG; remove the CPU `linear-solver` selection.");

    const std::string simType = config.getString("sim-type", 1);
    if (simType != "static" && simType != "dynamic")
      throw std::invalid_argument("CUDA IPC requires `sim-type` equal to `static` or `dynamic`.");
    const bool dynamic = simType == "dynamic";
    const bool restart = config.getValue<bool>("restart-from-u", 0, false);
    if (restart && !dynamic)
      throw std::invalid_argument("`restart-from-u` requires a dynamic simulation.");
    const double timestep = (dynamic || frictionEnabled) ? config.getDouble("timestep", 1) : 1.0;
    const int numSteps = dynamic ? config.getInt("num-timestep", 1) : 1;
    const int dumpInterval = dynamic ? config.getInt("dump-interval", 1) : 1;
    const int maxIterations = config.getInt("solver-max-iter", 1);
    const double tolerance = config.getDouble("solver-eps", 1);
    if (!std::isfinite(timestep) || timestep <= 0 || numSteps < 0 || dumpInterval < 1 || maxIterations < 1 || !std::isfinite(tolerance) || tolerance <= 0)
      throw std::invalid_argument("Invalid CUDA IPC timestep, iteration count, dump interval, or solver tolerance.");
    const auto damping = dynamic ? config.getValue<std::array<double, 2>>("damping-params", 1) : std::array<double, 2>{ 0, 0 };
    if (!std::isfinite(damping[0]) || damping[0] < 0 || !std::isfinite(damping[1]) || damping[1] < 0)
      throw std::invalid_argument("CUDA IPC requires finite, nonnegative Rayleigh damping coefficients.");

    const double scale = config.getDouble("scale", 1);
    const std::filesystem::path output = config.getResolvedPath("output", 1);
    std::filesystem::create_directories(output);
    if (enableCliLog)
      logRedirect = std::make_unique<RunSim::ScopedRunSimCliLogRedirect>((output / "runIPCSim.log").string());
    Logging::setLevel(RunSim::resolveConfiguredLogLevel(config));
    const int device = config.getValue<int>("cuda-device", 0, 0);
    if (CudaUtilities::initDevice(device) != 0)
      throw std::runtime_error("Unable to initialize the requested CUDA device.");
    Mesh::initPredicates();

    auto context = shell ? RunIPCSim::buildShellIpcSimulation(config) : RunIPCSim::buildVolumeIpcSimulation(config);
    const int n3 = static_cast<int>(context.simulationRestPosition.size());
    auto makeElastic = [&]() -> std::shared_ptr<NO::CudaPotentialEnergy> {
      if (shell)
        return std::make_shared<SolidDeformationModel::CudaKoiterDeformationModelEnergy>(
          *context.simulationMeshOwner, context.elasticParameters, context.zeroShellRestCurvature);
      const auto material = materialName == "stable-neo" ? SolidDeformationModel::DeformationModelElasticMaterial::STABLE_NEO :
                                                           SolidDeformationModel::DeformationModelElasticMaterial::STVK_VOL;
      return std::make_shared<SolidDeformationModel::CudaDeformationModelEnergy>(
        *context.simulationMeshOwner, material, context.elasticEnergy->isMaterialMaxStepEnabled(), cubicValidity);
    };
    auto elastic = makeElastic();
    {
      const nlohmann::json validity = {
        { "schema_version", 1 },
        { "element_type", shell ? "shell" : config.exist("cubic-mesh") ? "cubic" :
                                                                         "tet" },
        { "cubic_validity", cubicValidityName },
        { "material_max_step_enabled", context.elasticEnergy->isMaterialMaxStepEnabled() },
        { "determinant_floor", 1e-8 }
      };
      const auto validityPath = output / "material-validity.json";
      if (restart && config.exist("cubic-mesh") && !std::filesystem::exists(validityPath) &&
        !config.exist("cuda-cubic-validity"))
        throw std::invalid_argument("A cubic restart without material-validity.json requires an explicit `cuda-cubic-validity` choice.");
      if (restart && std::filesystem::exists(validityPath)) {
        std::ifstream previousStream(validityPath);
        const auto previous = nlohmann::json::parse(previousStream);
        if (previous != validity)
          throw std::invalid_argument("A CUDA restart must retain its recorded material validity criterion.");
      }
      std::ofstream stream(validityPath);
      stream << std::setw(2) << validity << '\n';
      stream.flush();
      if (!stream)
        throw std::runtime_error("Failed to write material validity metadata.");
    }
    std::shared_ptr<NO::CudaFrozenHessianEnergy> stiffnessDamping;
    if (dynamic && damping[1] > 0) {
      // Independent Hessian storage preserves K(u_n) while the elastic model
      // prepares K(u) at successive Newton iterates, matching CPU updateD().
      auto dampingModel = makeElastic();
      stiffnessDamping = std::make_shared<NO::CudaFrozenHessianEnergy>(std::move(dampingModel));
    }
    const auto parameters = context.collisionHandler->getIPCParameters();
    auto barrierParameters = parameters;
    barrierParameters.frictionCoeff = 0;  // Friction has its own frozen lagging state.
    using ContactEnergy = Contact::CIPC::CudaEmbeddedSurfaceIPCPotentialEnergy;
    auto makeContact = [&]() {
      const auto &surface = context.collisionSurface;
      return std::make_shared<ContactEnergy>(surface.vertices, surface.triangles,
        surface.displacementMap, surface.vertexIsDeformable, barrierParameters);
    };
    auto contact = makeContact();
    std::shared_ptr<Contact::CIPC::CudaSurfaceIPCFrictionPotentialEnergy> friction;
    if (frictionEnabled) {
      const auto &surface = context.collisionSurface;
      friction = std::make_shared<Contact::CIPC::CudaSurfaceIPCFrictionPotentialEnergy>(
        surface.vertices, surface.triangles, surface.displacementMap, surface.vertexIsDeformable, parameters);
    }

    NO::CudaNewtonSolver::SolverParam solverParameters;
    // Dynamic steps follow their input-to-output physical chord. A static
    // solve follows its recorded Newton path toward a full-load equilibrium.
    solverParameters.checkInputChord = dynamic;
    solverParameters.checkTrialChord = !dynamic;
    const int maxRefinements = config.getValue<int>("cuda-max-timestep-refinements", 0, 8);
    if (maxRefinements < 0 || maxRefinements > 12)
      throw std::invalid_argument("`cuda-max-timestep-refinements` must be between 0 and 12.");
    const int solverVerbose = config.getValue<int>("cuda-solver-verbose", 0, 0);
    if (solverVerbose < 0 || solverVerbose > 3)
      throw std::invalid_argument("`cuda-solver-verbose` must be between 0 and 3.");
    if (config.exist("cuda-pcg-max-iter"))
      solverParameters.maxPCGIterations = config.getInt("cuda-pcg-max-iter", 1);
    if (config.exist("cuda-pcg-tolerance"))
      solverParameters.pcgRelativeTolerance = config.getDouble("cuda-pcg-tolerance", 1);
    if (config.exist("cuda-pcg-batch-size"))
      solverParameters.pcgReductionBatchSize = config.getInt("cuda-pcg-batch-size", 1);
    if (config.exist("cuda-exact-hessian-start-iteration"))
      solverParameters.exactHessianStartIteration = config.getInt("cuda-exact-hessian-start-iteration", 1);
    if (config.exist("cuda-stagnation-window"))
      solverParameters.stagnationWindow = config.getInt("cuda-stagnation-window", 1);
    if (config.exist("cuda-stagnation-tolerance"))
      solverParameters.stagnationRelativeTolerance = config.getDouble("cuda-stagnation-tolerance", 1);
    if (solverParameters.maxPCGIterations < 0 || !std::isfinite(solverParameters.pcgRelativeTolerance) ||
      solverParameters.pcgRelativeTolerance <= 0 || solverParameters.pcgRelativeTolerance >= 1)
      throw std::invalid_argument("CUDA PCG requires a nonnegative iteration limit and a relative tolerance in (0, 1).");
    const auto gravity = config.getValue<std::array<double, 3>>("g", 1);
    const ES::V3d acceleration = ES::Mp<const ES::V3d>(gravity.data());
    if (!acceleration.allFinite())
      throw std::invalid_argument("Gravity must be finite.");
    ES::VXd gravityAll(n3);
    for (int i = 0; i < n3 / 3; ++i)
      gravityAll.segment<3>(3 * i) = acceleration;
    const ES::VXd externalForce = context.M * gravityAll;
    ES::VXd u = context.initialDisplacement;
    ES::VXd velocity = ES::VXd::Zero(n3), accelerationState = ES::VXd::Zero(n3);
    if (dynamic) {
      const auto initialVelocity = config.getValue<std::array<double, 3>>("init-vel", 1);
      const ES::V3d v = ES::Mp<const ES::V3d>(initialVelocity.data());
      if (!v.allFinite())
        throw std::invalid_argument("Initial velocity must be finite.");
      for (int i = 0; i < n3 / 3; ++i)
        velocity.segment<3>(3 * i) = v;
    }

    int firstFrame = 0;
    if (restart) {
      for (int frame = numSteps - 1; frame >= 0; --frame) {
        const auto filename = output / fmt::format("deform{:04d}.u", frame);
        if (!std::filesystem::exists(filename))
          continue;
        ES::MXd state;
        if (ES::readMatrix(filename.string().c_str(), state) != 0 || state.rows() != n3 || state.cols() < 3 || !state.leftCols(3).allFinite())
          throw std::runtime_error("Invalid CUDA IPC restart state: " + filename.string());
        u = state.col(0);
        velocity = state.col(1);
        accelerationState = state.col(2);
        firstFrame = frame + 1;
        break;
      }
    }

    // Like runIPCSim.log, this file records the current invocation, including
    // the original frame indices when continuing from a checkpoint.
    std::ofstream solverLog(output / "newton-solves.csv");
    if (!solverLog)
      throw std::runtime_error("Unable to create CUDA Newton diagnostics.");
    solverLog << std::setprecision(17)
              << "frame,newton_iterations,pcg_iterations,active_pairs,gradient_inf,tolerance,status,stop,energy,linear_residual,last_feasible_alpha,last_accepted_alpha,regularizations,roundoff_steps,exact_hessian_attempts,exact_hessian_directions,exact_hessian_fallbacks,substeps,refinement_retries,max_substep_newton_iterations,rejected_newton_solves,input_chord_rejections,macroframe_chord_certified,residual_correction_attempts,residual_correction_candidates,residual_correction_steps,friction_solves,trial_chord_rejections,stagnation_triggers\n";
    const std::string substepHeader = "frame,substep,substeps,fraction,physical_time,state_file,newton_iterations,pcg_iterations,gradient_inf,tolerance,status,stop,active_pairs";
    const auto substepPath = output / "newton-substeps.csv";
    const bool appendSubsteps = restart && firstFrame > 0 && std::filesystem::exists(substepPath);
    if (appendSubsteps) {
      std::ifstream previousLog(substepPath);
      std::string line;
      if (!std::getline(previousLog, line) || line != substepHeader)
        throw std::runtime_error("Invalid accepted-substep manifest at restart.");
      int previousFrame = -1, previousStep = -1, previousCount = 0;
      while (std::getline(previousLog, line)) {
        std::istringstream row(line);
        int recordedFrame, recordedStep, recordedCount;
        char separator1, separator2;
        if (!(row >> recordedFrame >> separator1 >> recordedStep >> separator2 >> recordedCount) ||
          separator1 != ',' || separator2 != ',' || recordedFrame < 0 || recordedCount < 1 ||
          recordedStep < 0 || recordedStep >= recordedCount ||
          (previousFrame == -1 && recordedStep != 0) ||
          (previousFrame >= 0 && !(recordedFrame == previousFrame && recordedCount == previousCount && recordedStep == previousStep + 1) &&
            !(recordedFrame == previousFrame + 1 && previousStep + 1 == previousCount && recordedStep == 0)))
          throw std::runtime_error("Incomplete or unordered accepted-substep manifest at restart.");
        previousFrame = recordedFrame;
        previousStep = recordedStep;
        previousCount = recordedCount;
      }
      if (previousFrame != firstFrame - 1 || previousStep + 1 != previousCount)
        throw std::runtime_error("The accepted-substep manifest does not end at the restart state.");
    }
    std::ofstream attemptLog(output / "newton-attempts.csv");
    std::ofstream substepLog(substepPath, appendSubsteps ? std::ios::app : std::ios::out);
    if (!attemptLog || !substepLog)
      throw std::runtime_error("Unable to create CUDA timestep diagnostics.");
    attemptLog << std::setprecision(17)
               << "frame,refinement,substep,planned_substeps,accepted,newton_iterations,pcg_iterations,gradient_inf,tolerance,status,stop,rejection,input_chord_rejections,residual_correction_attempts,residual_correction_candidates,residual_correction_steps,friction_iteration,friction_pairs,trial_chord_rejections,stagnation_triggers,stagnation_step_ratio,stagnation_norm2_relative_change,stagnation_infinity_relative_change\n";
    substepLog << std::setprecision(17);
    if (!appendSubsteps)
      substepLog << substepHeader << '\n';
    nlohmann::json staticNewtonPath = {
      { "schema_version", 1 }, { "path_type", "static_newton_iterates" },
      { "parameter", "accepted_newton_iteration" }, { "physical_time", false },
      { "num_dofs", n3 }, { "complete", false },
      { "states", nlohmann::json::array() }, { "solves", nlohmann::json::array() }
    };
    auto writeStaticNewtonPath = [&]() {
      std::ofstream stream(output / "static-newton-path.json");
      stream << std::setw(2) << staticNewtonPath << '\n';
      stream.flush();
      if (!stream)
        throw std::runtime_error("Failed to write the static Newton path manifest.");
    };
    if (!dynamic) {
      std::filesystem::create_directories(output / "newton-path");
      writeStaticNewtonPath();
    }
    CudaUtilities::DeviceBuffer<double> chordFrom(n3), chordDelta(n3), frictionReference, frictionLagged;
    std::cout << "simulation-backend=cuda; " << (shell ? "shell" : config.exist("cubic-mesh") ? "cubic" :
                                                                                                "tet")
              << ' ' << materialName << "; IPC friction=" << frictionCoefficient
              << "; lagging-iterations=" << (frictionEnabled ? frictionIterations : 1) << "; DOFs=" << n3
              << "; material-max-step=" << context.elasticEnergy->isMaterialMaxStepEnabled()
              << "; cubic-validity=" << cubicValidityName << '\n';
    const auto start = std::chrono::steady_clock::now();
    std::size_t peakContacts = 0;
    int totalNewton = 0, totalPCG = 0, nonconverged = 0;
    // A completed restart still needs the saved external-object placement for
    // final state validation, even though no timestep will be executed.
    if (dynamic && firstFrame == numSteps && firstFrame > 0) {
      RunIPCSim::placeIPCKinematicObjects(*context.collisionHandler,
        context.kinematicExternalObjects, firstFrame - 1, numSteps);
    }
    int executedFrames = 0;
    for (int frame = firstFrame; frame < numSteps &&
         (options.maxFrames < 0 || executedFrames < options.maxFrames);
         ++frame) {
      RunSim::IPCSimulationFrameTimer frameTiming(timing, frame);
      const double objectRatio = numSteps > 1 ? static_cast<double>(frame) / (numSteps - 1) : 0.0;
      bool moved = false;
      for (const auto &object : context.kinematicExternalObjects) {
        if (!object.moves() || !dynamic)
          continue;
        const ES::MXd positions = object.positionsAt(objectRatio);
        context.collisionSurface.vertices.middleRows(object.firstSurfaceVertex, positions.rows()) = positions;
        context.collisionHandler->updateExternalSurface(object.firstSurfaceVertex, positions);
        if (friction)
          friction->updateExternalSurface(object.firstSurfaceVertex, positions);
        moved = true;
      }
      if (moved)
        contact = makeContact();
      if (frame == firstFrame || moved)
        context.collisionHandler->validateCollisionFreeState(u);

      const ES::VXd frameStart = u;
      const ES::VXd frameVelocity = velocity;
      const ES::VXd frameAcceleration = accelerationState;
      std::vector<SubstepResult> accepted;
      bool macroframeChordCertified = false;
      int refinement = 0, rejectedSolves = 0;
      for (; refinement <= (dynamic ? maxRefinements : 0); ++refinement) {
        const int count = 1 << refinement;
        const double h = timestep / count;
        u = frameStart;
        velocity = frameVelocity;
        accelerationState = frameAcceleration;
        std::vector<SubstepResult> trial;
        std::string rejection;
        for (int step = 0; step < count; ++step) {
          // Interpolate prescribed pulling targets between macroframes. The
          // first frame retains the original runner's initial target time.
          const double frameCoordinate = std::max(0.0, frame - 1.0 + (step + 1.0) / count);
          const double targetRatio = numSteps > 1 ? frameCoordinate / (numSteps - 1) : 1.0;
          for (std::size_t i = 0; i < context.pullingEnergies.size(); ++i) {
            const ES::VXd target = dynamic ? ES::VXd((1.0 - targetRatio) * context.pullingTargetRests[i] + targetRatio * context.pullingTargets[i]) : context.pullingTargets[i];
            context.pullingEnergies[i]->setTargetPos(target.data());
          }
          for (auto &pulling : context.trajectoryPullings) {
            const double time = dynamic ? frameCoordinate * timestep : pulling.track.times.back();
            const bool active = !dynamic || pulling.track.isActiveAt(time);
            pulling.energy->setCoeff(active ? pulling.coeff : 0.0);
            pulling.released = !active;
            if (active) {
              const ES::VXd target = pulling.targetAt(time);
              pulling.energy->setTargetPos(target.data());
            }
          }

          const ES::VXd previous = u;
          const ES::VXd previousVelocity = velocity;
          ES::SpMatD matrix(n3, n3);
          ES::VXd linear = -externalForce;
          if (dynamic) {
            matrix = context.M * (1.0 / (h * h) + damping[0] / h);
            linear -= (context.M * velocity) / h;
          }
          double constant = 0;
          addPullingQuadratic(context, previous, matrix, linear, constant);
          matrix.prune(0.0);
          matrix.makeCompressed();
          auto quadratic = std::make_shared<NO::CudaQuadraticEnergy>(matrix, previous, linear, constant);
          std::vector<std::shared_ptr<NO::CudaPotentialEnergy>> energies{ elastic, contact, quadratic };
          if (stiffnessDamping) {
            stiffnessDamping->setReference(previous, damping[1] / h);
            energies.push_back(stiffnessDamping);
          }
          if (friction) {
            energies.push_back(friction);
            frictionReference.upload(previous.data(), n3);
          }
          SubstepResult result;
          // Every lagging iteration uses the original x_n, v_n, frozen damping
          // and incremental inertia. Only the final solution advances time.
          for (int iteration = 0; iteration < (frictionEnabled ? frictionIterations : 1); ++iteration) {
            NewtonResult solveResult;
            if (friction) {
              frictionLagged.upload(u.data(), n3);
              friction->updateFriction(frictionReference.data(), frictionLagged.data(), h);
              solveResult.frictionPairs = friction->getNumFrictionPairs();
            }
            NO::CudaNewtonSolver solver(energies, solverParameters);
            const std::size_t firstNewtonState = staticNewtonPath["states"].size();
            if (!dynamic) {
              solver.setAcceptedIterateObserver([&](int newtonIteration, const double *deviceState,
                                                  double energy, double gradientInfinity) {
                frameTiming.output([&]() {
                  if (!std::isfinite(energy) || !std::isfinite(gradientInfinity))
                    throw std::runtime_error("Nonfinite static Newton path diagnostics.");
                  ES::MXd state = ES::MXd::Zero(n3, 3);
                  CudaUtilities::checkCuda(cudaMemcpy(state.col(0).data(), deviceState,
                                             sizeof(double) * n3, cudaMemcpyDeviceToHost),
                    "Download static Newton iterate");
                  if (!state.allFinite())
                    throw std::runtime_error("Nonfinite static Newton path state.");
                  const auto file = std::filesystem::path("newton-path") /
                    fmt::format("solve{:04d}-iterate{:06d}.u", iteration, newtonIteration);
                  writeStateFile(output / file, state);
                  staticNewtonPath["states"].push_back({ { "solve", iteration }, { "iteration", newtonIteration },
                    { "state_file", file.generic_string() }, { "energy", energy },
                    { "gradient_inf", gradientInfinity } });
                });
              });
            }
            if (solverVerbose)
              std::cout << "CUDA solve frame " << frame << ", refinement=" << refinement
                        << ", substep=" << step + 1 << '/' << count
                        << ", friction-iteration=" << iteration << ", friction-pairs=" << solveResult.frictionPairs << std::endl;
            solveResult.status = solver.solve(u.data(), maxIterations, tolerance, solverVerbose);
            frameTiming.recordSolve(solveResult.status);
            solveResult.diagnostics = solver.getDiagnostics();
            if (!dynamic) {
              frameTiming.output([&]() {
                const auto &d = solveResult.diagnostics;
                const auto stateCount = staticNewtonPath["states"].size();
                staticNewtonPath["solves"].push_back({ { "solve", iteration }, { "status", solveResult.status },
                  { "stop", d.stopReason }, { "tolerance", tolerance },
                  { "gradient_inf", d.finalGradientInfinityNorm },
                  { "newton_iterations", d.newtonIterations },
                  { "first_state", firstNewtonState }, { "last_state", static_cast<long long>(stateCount) - 1 } });
                writeStaticNewtonPath();
                if (solveResult.status == NO::CudaNewtonSolver::SOLVE_CONVERGED &&
                  stateCount - firstNewtonState != static_cast<std::size_t>(d.newtonIterations) + 1)
                  throw std::runtime_error("Incomplete accepted-iterate coverage in the static Newton path.");
              });
            }
            result.status = solveResult.status;
            result.maxNewtonIterations = std::max(result.maxNewtonIterations, solveResult.diagnostics.newtonIterations);
            appendDiagnostics(result.diagnostics, solveResult.diagnostics);
            result.solves.push_back(std::move(solveResult));
            if (result.status != NO::CudaNewtonSolver::SOLVE_CONVERGED)
              break;
          }
          result.activePairs = contact->getNumActivePairs();
          if (dynamic) {
            velocity = (u - previous) / h;
            accelerationState = (velocity - previousVelocity) / h;
          }
          result.state.resize(n3, 3);
          result.state.col(0) = u;
          result.state.col(1) = velocity;
          result.state.col(2) = accelerationState;
          trial.push_back(std::move(result));
          if (trial.back().status != NO::CudaNewtonSolver::SOLVE_CONVERGED) {
            rejection = "newton_not_converged";
            break;
          }
          if (!trial.back().state.allFinite()) {
            rejection = "nonfinite_state";
            break;
          }
          if (dynamic && !certifyCudaChord({ contact.get(), elastic.get() }, previous, u, chordFrom, chordDelta)) {
            rejection = "uncertified_physical_substep";
            break;
          }
        }
        // The recorded physical substeps define the trajectory. A direct
        // shortcut between macroframe endpoints can cross a surface even when
        // that actual path is clear. Retain its certificate as a diagnostic;
        // consumers must follow the complete recorded path when it is false.
        if (rejection.empty())
          macroframeChordCertified = (dynamic && count == 1) ||
            certifyCudaChord({ contact.get(), elastic.get() }, frameStart, u, chordFrom, chordDelta);
        frameTiming.output([&]() {
          for (std::size_t step = 0; step < trial.size(); ++step) {
            const auto &result = trial[step];
            for (std::size_t iteration = 0; iteration < result.solves.size(); ++iteration) {
              const auto &solveResult = result.solves[iteration];
              const auto &d = solveResult.diagnostics;
              attemptLog << frame << ',' << refinement << ',' << step << ',' << count << ',' << rejection.empty()
                         << ',' << d.newtonIterations << ',' << d.pcgIterations << ',' << d.finalGradientInfinityNorm
                         << ',' << tolerance << ',' << solveResult.status << ',' << d.stopReason << ',' << rejection
                         << ',' << d.inputChordRejections << ',' << d.residualCorrectionAttempts
                         << ',' << d.residualCorrectionCandidates << ',' << d.residualCorrectionSteps
                         << ',' << iteration << ',' << solveResult.frictionPairs
                         << ',' << d.trialChordRejections << ',' << d.stagnationTriggers
                         << ',' << d.lastStagnationStepRatio << ',' << d.lastStagnationNorm2RelativeChange
                         << ',' << d.lastStagnationInfinityRelativeChange << '\n';
            }
          }
          attemptLog.flush();
          if (!attemptLog)
            throw std::runtime_error("Failed to write CUDA timestep attempt diagnostics.");
        });
        if (rejection.empty()) {
          accepted = std::move(trial);
          break;
        }
        for (const auto &result : trial)
          rejectedSolves += static_cast<int>(result.solves.size());
        if (solverVerbose >= 3) {
          frameTiming.output([&]() {
            const auto directory = output / "rejected" / fmt::format("frame{:04d}-refinement{:02d}", frame, refinement);
            std::filesystem::create_directories(directory);
            ES::MXd initial(n3, 3);
            initial.col(0) = frameStart;
            initial.col(1) = frameVelocity;
            initial.col(2) = frameAcceleration;
            writeStateFile(directory / "start.u", initial);
            for (std::size_t step = 0; step < trial.size(); ++step)
              writeStateFile(directory / fmt::format("step{:04d}.u", step), trial[step].state);
          });
        }
        std::cout << "CUDA frame " << frame << " retry: " << rejection
                  << ", attempted-substeps=" << count << std::endl;
      }
      if (accepted.empty())
        throw std::runtime_error(fmt::format("CUDA frame {} could not converge with certified physical motion after {} refinements; no state was committed.", frame, dynamic ? maxRefinements : 0));

      NO::CudaNewtonSolver::Diagnostics diagnostics;
      int maxSubstepIterations = 0, frictionSolves = 0;
      for (const auto &result : accepted) {
        appendDiagnostics(diagnostics, result.diagnostics);
        maxSubstepIterations = std::max(maxSubstepIterations, result.maxNewtonIterations);
        if (frictionEnabled)
          frictionSolves += static_cast<int>(result.solves.size());
        peakContacts = std::max(peakContacts, result.activePairs);
      }
      const auto activePairs = accepted.back().activePairs;
      frameTiming.output([&]() {
        solverLog << frame << ',' << diagnostics.newtonIterations << ',' << diagnostics.pcgIterations
                  << ',' << activePairs << ',' << diagnostics.finalGradientInfinityNorm
                  << ',' << tolerance << ',' << NO::CudaNewtonSolver::SOLVE_CONVERGED << ',' << diagnostics.stopReason
                  << ',' << diagnostics.finalEnergy << ',' << diagnostics.lastLinearResidualNorm
                  << ',' << diagnostics.lastFeasibleAlpha << ',' << diagnostics.lastAcceptedAlpha
                  << ',' << diagnostics.regularizationAttempts << ',' << diagnostics.roundoffAcceptedSteps
                  << ',' << diagnostics.exactHessianAttempts << ',' << diagnostics.exactHessianDirections
                  << ',' << diagnostics.exactHessianFallbacks << ',' << accepted.size() << ',' << refinement
                  << ',' << maxSubstepIterations << ',' << rejectedSolves << ',' << diagnostics.inputChordRejections
                  << ',' << macroframeChordCertified << ',' << diagnostics.residualCorrectionAttempts
                  << ',' << diagnostics.residualCorrectionCandidates << ',' << diagnostics.residualCorrectionSteps
                  << ',' << frictionSolves << ',' << diagnostics.trialChordRejections
                  << ',' << diagnostics.stagnationTriggers << std::endl;
        if (!solverLog)
          throw std::runtime_error("Failed to write CUDA Newton diagnostics.");
        totalNewton += diagnostics.newtonIterations;
        totalPCG += diagnostics.pcgIterations;
        saveCudaState(context, u, velocity, accelerationState, output, frame, dumpInterval, scale);
        std::filesystem::create_directories(output / "substeps");
        for (std::size_t step = 0; step < accepted.size(); ++step) {
          const auto &result = accepted[step];
          const auto &d = result.diagnostics;
          std::filesystem::path stateFile = fmt::format("deform{:04d}.u", frame);
          if (step + 1 < accepted.size()) {
            stateFile = std::filesystem::path("substeps") / fmt::format("frame{:04d}-step{:04d}.u", frame, step);
            std::filesystem::create_directories((output / stateFile).parent_path());
            writeStateFile(output / stateFile, result.state);
          }
          const auto surfaceFile = output / "substeps" / fmt::format("frame{:04d}-step{:04d}.obj", frame, step);
          writeSurfaceFile(context, result.state.col(0), surfaceFile, scale);
          const double fraction = static_cast<double>(step + 1) / accepted.size();
          substepLog << frame << ',' << step << ',' << accepted.size() << ',' << fraction
                     << ',' << (dynamic ? (frame + fraction) * timestep : 0.0) << ',' << stateFile.generic_string()
                     << ',' << d.newtonIterations << ',' << d.pcgIterations << ',' << d.finalGradientInfinityNorm
                     << ',' << tolerance << ',' << result.status << ',' << d.stopReason << ',' << result.activePairs << '\n';
        }
        substepLog.flush();
        if (!substepLog)
          throw std::runtime_error("Failed to write CUDA accepted-substep diagnostics.");
      });
      std::cout << "CUDA frame " << frame << ": Newton=" << diagnostics.newtonIterations
                << ", PCG=" << diagnostics.pcgIterations << ", contacts=" << activePairs
                << ", |g|inf=" << diagnostics.finalGradientInfinityNorm
                << ", roundoff-steps=" << diagnostics.roundoffAcceptedSteps
                << ", exact-directions=" << diagnostics.exactHessianDirections
                << ", exact-fallbacks=" << diagnostics.exactHessianFallbacks
                << ", substeps=" << accepted.size() << ", refinement-retries=" << refinement
                << ", macroframe-chord-certified=" << macroframeChordCertified
                << ", friction-solves=" << frictionSolves
                << ", residual-correction-steps=" << diagnostics.residualCorrectionSteps
                << ", input-chord-rejections=" << diagnostics.inputChordRejections
                << ", trial-chord-rejections=" << diagnostics.trialChordRejections
                << ", stagnation-triggers=" << diagnostics.stagnationTriggers
                << ", stop=" << diagnostics.stopReason << std::endl;
      if (!dynamic) {
        context.collisionHandler->validateCollisionFreeState(u);
        frameTiming.output([&]() {
          staticNewtonPath["complete"] = true;
          writeStaticNewtonPath();
        });
      }
      frameTiming.finish(true, true);
      ++executedFrames;
    }
    // Independent CPU geometry validation at the output boundary catches errors
    // in the GPU path without introducing a CPU contact evaluation in Newton.
    if (dynamic)
      context.collisionHandler->validateCollisionFreeState(u);
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    std::cout << "CUDA IPC complete: frames=" << executedFrames << ", peak-contacts=" << peakContacts
              << ", Newton=" << totalNewton << ", PCG=" << totalPCG << ", nonconverged=" << nonconverged
              << ", elapsed=" << seconds << " s" << std::endl;
    return 0;
  }
  catch (const std::exception &error) {
    SPDLOG_LOGGER_ERROR(Logging::lgr(), "{}", error.what());
    if (logRedirect) {
      logRedirect.reset();
      std::cerr << error.what() << '\n';
    }
    return 1;
  }
}
