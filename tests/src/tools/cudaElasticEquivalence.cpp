// Optional test/measurement executable; does not run a simulation or write
// into its configured output. All geometry and material setup uses the same
// repository setup path as the CPU and CUDA IPC runners.
#include "cudaElasticEquivalence.h"
#include "cudaDeformationModelEnergy.h"
#include "cudaKoiterDeformationModelEnergy.h"
#include "configFileJSON.h"
#include "deformationModelAssembler.h"
#include "initPredicates.h"
#include "pgoLogging.h"
#include "runIPCSimSetup.h"
#include "simulationMesh.h"
#include "triMeshGeo.h"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;
namespace ES = pgo::EigenSupport;
namespace SDM = pgo::SolidDeformationModel;
using Json = nlohmann::json;

int main(int argc, char **argv)
{
  try {
    if (argc < 2)
      throw std::invalid_argument("Usage: cudaElasticEquivalence CONFIG --report JSON [--state STATE.u ...] [--skip-initial] [--timings] [--full-columns] [--device N] [--cpu-assembly default|deterministic]");
    const fs::path configFile = fs::absolute(argv[1]);
    fs::path reportFile;
    std::vector<fs::path> states;
    bool timings = false, fullColumns = false, initial = true;
    std::string cpuAssembly = "deterministic";
    int device = 0;
    for (int i = 2; i < argc; ++i) {
      const std::string option = argv[i];
      if (option == "--report" && i + 1 < argc)
        reportFile = fs::absolute(argv[++i]);
      else if (option == "--state" && i + 1 < argc)
        states.push_back(fs::absolute(argv[++i]));
      else if (option == "--device" && i + 1 < argc)
        device = std::stoi(argv[++i]);
      else if (option == "--cpu-assembly" && i + 1 < argc)
        cpuAssembly = argv[++i];
      else if (option == "--timings")
        timings = true;
      else if (option == "--full-columns")
        fullColumns = true;
      else if (option == "--skip-initial")
        initial = false;
      else
        throw std::invalid_argument("Unknown or incomplete argument: " + option);
    }
    if (reportFile.empty() || (!initial && states.empty()))
      throw std::invalid_argument("A report path and at least one state are required.");
    if (cpuAssembly != "default" && cpuAssembly != "deterministic")
      throw std::invalid_argument("CPU assembly must be default or deterministic.");
    pgo::Logging::init();
    pgo::CudaUtilities::checkCuda(cudaSetDevice(device), "Select elastic measurement GPU");
    pgo::CudaUtilities::checkCuda(cudaFree(nullptr), "Initialize elastic measurement GPU");
    cudaDeviceProp properties;
    pgo::CudaUtilities::checkCuda(cudaGetDeviceProperties(&properties, device), "Read elastic measurement GPU");
    pgo::Mesh::initPredicates();
    pgo::ConfigFileJSON config(configFile.string().c_str());
    const bool shell = !config.exist("tet-mesh") && !config.exist("cubic-mesh");
    const std::string material = config.getString("elastic-material", 1);
    if ((shell && material != "koiter-stvk") || (!shell && material != "stable-neo" && material != "stvk-vol"))
      throw std::invalid_argument("Only stable-neo/StVK-vol volumes and Koiter-StVK shells are supported.");
    auto context = shell ? pgo::RunIPCSim::buildShellIpcSimulation(config) : pgo::RunIPCSim::buildVolumeIpcSimulation(config);
    context.deformationModelAssemblerOwner->setDeterministic(cpuAssembly == "deterministic");
    std::unique_ptr<pgo::NonlinearOptimization::CudaPotentialEnergy> gpu;
    if (shell)
      gpu = std::make_unique<SDM::CudaKoiterDeformationModelEnergy>(*context.simulationMeshOwner,
        context.elasticParameters, context.zeroShellRestCurvature);
    else
      gpu = std::make_unique<SDM::CudaDeformationModelEnergy>(*context.simulationMeshOwner,
        material == "stable-neo" ? SDM::DeformationModelElasticMaterial::STABLE_NEO : SDM::DeformationModelElasticMaterial::STVK_VOL,
        context.elasticEnergy->isMaterialMaxStepEnabled(), SDM::CudaDeformationModelEnergy::CubicValidity::Gauss);
    if (fullColumns && gpu->getNumDOFs() > 256)
      throw std::invalid_argument("Full Hessian reconstruction is limited to small fixtures (at most 256 DOFs).");
    Json report = {
      { "schema_version", 1 }, { "config", configFile.string() }, { "material", material },
      { "element_type", shell ? "koiter_shell" : config.exist("cubic-mesh") ? "cubic" :
                                                                              "tet" },
      { "num_dofs", gpu->getNumDOFs() }, { "num_elements", context.simulationMeshOwner->getNumElements() },
      { "gpu_name", properties.name }, { "logical_cuda_device", device },
      { "precision", "IEEE double" }, { "scope", "elastic model only; excludes contact, inertia, damping, external forces, and time integration" },
      { "cpu_assembly", cpuAssembly },
      { "cpu_assembly_description", cpuAssembly == "deterministic" ? "fixed-order assembly of per-element contributions" : "repository default parallel atomic assembly; deterministic assembly disabled" },
      { "timings_enabled", timings },
      { "timing_scope", timings ? "operation measurements; require an externally reserved isolated resource slot for performance claims" : "disabled" },
      { "projection_scope", shell ? "full 18x18 Koiter patch Hessian" : material == "stable-neo" ? "9x9 tangent at each material quadrature point" :
                                                                                                   "StVK child 9x9 tangent projection plus raw volumetric penalty, matching CPU; not necessarily PSD" },
      { "states", Json::array() }
    };
    bool passed = true;
    auto compareState = [&](const std::string &kind, const fs::path &source, ES::ConstRefVecXd displacement) {
      auto measurement = pgo::Tests::ElasticEquivalence::compare(*context.elasticEnergy,
        *context.deformationModelManagerOwner, *gpu, displacement, fullColumns, timings);
      measurement["state_kind"] = kind;
      if (!source.empty())
        measurement["state_file"] = source.string();
      measurement["displacement_infinity_norm"] = displacement.cwiseAbs().maxCoeff();
      passed = passed && measurement.at("passed").get<bool>();
      std::cout << kind << ": elastic equivalence " << (measurement.at("passed").get<bool>() ? "PASS" : "FAIL") << '\n';
      report["states"].push_back(std::move(measurement));
      report["passed"] = passed;
      if (!reportFile.parent_path().empty())
        fs::create_directories(reportFile.parent_path());
      std::ofstream stream(reportFile);
      stream << report.dump(2) << '\n';
      if (!stream)
        throw std::runtime_error("Cannot write elastic equivalence report.");
    };
    if (initial)
      compareState("configured_initial_displacement", {}, context.initialDisplacement);
    for (const auto &stateFile : states) {
      ES::MXd state;
      if (ES::readMatrix(stateFile.string().c_str(), state) != 0 || state.rows() != gpu->getNumDOFs() ||
        state.cols() != 3 || !state.allFinite())
        throw std::invalid_argument("Invalid saved u/v/a state: " + stateFile.string());
      compareState("saved_displacement", stateFile, state.col(0));
    }
    return passed ? 0 : 2;
  }
  catch (const std::exception &error) {
    std::cerr << "Elastic equivalence failed: " << error.what() << '\n';
    return 1;
  }
}
