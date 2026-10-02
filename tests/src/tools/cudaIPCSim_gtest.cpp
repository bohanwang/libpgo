#include <gtest/gtest.h>

#include "configFileJSON.h"
#include "deformationModelEnergy.h"
#include "EigenSupport.h"
#include "pgoLogging.h"
#include "runIPCSimSetup.h"
#include "simulationRunner.h"

#include <cuda_runtime_api.h>
#include <nlohmann/json.hpp>

#include <chrono>
#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <limits>
#include <map>
#include <sstream>
#include <vector>

namespace
{
namespace fs = std::filesystem;
namespace ES = pgo::EigenSupport;

std::vector<std::map<std::string, std::string>> readCSV(const fs::path &path)
{
  std::ifstream stream(path);
  std::string line, value;
  std::getline(stream, line);
  std::istringstream header(line);
  std::vector<std::string> names;
  while (std::getline(header, value, ','))
    names.push_back(value);
  std::vector<std::map<std::string, std::string>> rows;
  while (std::getline(stream, line)) {
    std::istringstream record(line);
    std::map<std::string, std::string> row;
    for (const auto &name : names) {
      std::getline(record, value, ',');
      row[name] = value;
    }
    rows.push_back(std::move(row));
  }
  return rows;
}

std::string shellQuote(const fs::path &path)
{
  std::string result = "'";
  for (char c : path.string())
    result += c == '\'' ? "'\\''" : std::string(1, c);
  return result + "'";
}

fs::path substepSurfaceFile(const fs::path &output, const std::map<std::string, std::string> &row)
{
  std::ostringstream name;
  name << "frame" << std::setfill('0') << std::setw(4) << std::stoi(row.at("frame"))
       << "-step" << std::setw(4) << std::stoi(row.at("substep")) << ".obj";
  return output / "substeps" / name.str();
}

std::string readFile(const fs::path &path)
{
  std::ifstream stream(path, std::ios::binary);
  return { std::istreambuf_iterator<char>(stream), std::istreambuf_iterator<char>() };
}

class CudaIPCRunnerTest : public ::testing::Test
{
protected:
  void SetUp() override
  {
    int count = 0;
    if (cudaGetDeviceCount(&count) != cudaSuccess || count == 0)
      GTEST_SKIP() << "A CUDA GPU is required for the IPC integration tests.";
    pgo::Logging::init();
    directory = fs::temp_directory_path() / ("pgo-cuda-ipc-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    ASSERT_TRUE(fs::create_directory(directory));
    std::ofstream plane(directory / "plane.obj");
    plane << "v -1 0 -1\nv 1 0 -1\nv 1 0 1\nv -1 0 1\nf 1 3 2\nf 1 4 3\n";
  }

  void TearDown() override
  {
    if (!directory.empty()) {
      std::error_code error;
      fs::remove_all(directory, error);
    }
  }

  nlohmann::json configuration(const std::string &backend, int steps = 3)
  {
    nlohmann::json config = {
      { "tet-mesh", std::string(PGO_CUDA_TEST_BOX_DIR) + "/box-tet.veg" },
      { "surface-mesh", std::string(PGO_CUDA_TEST_BOX_DIR) + "/box.obj" },
      { "fixed-vertices", nlohmann::json::array() },
      { "ipc-external-objects", nlohmann::json::array({ { { "filename", (directory / "plane.obj").string() }, { "scale", 1.0 }, { "initial-translation", { 0, 0, 0 } }, { "movement", { 0, 0, 0 } } } }) },
      { "g", { 0, -9.81, 0 } }, { "init-vel", { 0, 0, 0 } },
      { "init-disp", { 0, -0.377, 0 } }, { "scale", 1.0 },
      { "timestep", 0.001 }, { "num-timestep", steps }, { "damping-params", { 0, 0 } },
      { "sim-type", "dynamic" }, { "solver-eps", 1e-6 }, { "solver-max-iter", 200 },
      { "elastic-material", "stable-neo" }, { "dump-interval", 1 },
      { "output", (directory / backend).string() },
      { "ipc-dhat", 0.002 }, { "ipc-kappa", 3000.0 },
      { "enable-material-max-step", true }, { "contact-model", "ipc" },
      { "simulation-backend", backend }, { "loglevel", "warn" }
    };
    if (backend == "cpu")
      config["linear-solver"] = "eigen";
    else if (backend == "cuda") {
      config["cuda-pcg-batch-size"] = 8;
      config["cuda-stagnation-window"] = 4;
    }
    return config;
  }

  fs::path writeConfiguration(const nlohmann::json &config, const std::string &name)
  {
    const auto filename = directory / name;
    std::ofstream(filename) << config.dump(2);
    return filename;
  }

  void checkPhysicalSurfaceExports(const fs::path &configFile, const nlohmann::json &config)
  {
    const fs::path output = config.at("output").get<std::string>();
    const auto rows = readCSV(output / "newton-substeps.csv");
    ASSERT_FALSE(rows.empty());
    pgo::ConfigFileJSON parsed(configFile.string().c_str());
    const bool shell = !config.contains("tet-mesh") && !config.contains("cubic-mesh");
    const auto context = shell ? pgo::RunIPCSim::buildShellIpcSimulation(parsed) : pgo::RunIPCSim::buildVolumeIpcSimulation(parsed);
    const double scale = config.at("scale").get<double>();
    const double h = config.at("timestep").get<double>();
    double previousTime = -1;
    for (const auto &row : rows) {
      SCOPED_TRACE("frame " + row.at("frame") + ", substep " + row.at("substep"));
      ES::MXd state;
      ASSERT_EQ(ES::readMatrix((output / row.at("state_file")).string().c_str(), state), 0);
      ASSERT_EQ(state.cols(), 3);
      ASSERT_TRUE(state.allFinite());
      pgo::Mesh::TriMeshGeo surface;
      ASSERT_TRUE(surface.load(substepSurfaceFile(output, row).string()));
      ASSERT_EQ(surface.numVertices(), context.surfaceMesh.numVertices());
      ASSERT_EQ(surface.numTriangles(), context.surfaceMesh.numTriangles());
      const ES::VXd expected = (context.surfaceRestPositions + context.surfaceFromSimulationDispMap * state.col(0)) / scale;
      for (int vertex = 0; vertex < surface.numVertices(); ++vertex)
        EXPECT_LT((surface.pos(vertex) - expected.segment<3>(3 * vertex)).cwiseAbs().maxCoeff(), 1e-13);
      const double fraction = std::stod(row.at("fraction"));
      EXPECT_DOUBLE_EQ(fraction, (std::stoi(row.at("substep")) + 1.0) / std::stoi(row.at("substeps")));
      const double time = std::stod(row.at("physical_time"));
      EXPECT_DOUBLE_EQ(time, config.at("sim-type") == "static" ? 0.0 : (std::stoi(row.at("frame")) + fraction) * h);
      EXPECT_GT(time, previousTime);
      previousTime = time;
    }
  }

  fs::path directory;
};

TEST_F(CudaIPCRunnerTest, ContactTrajectoryMatchesCPUAndStaysSeparated)
{
  auto gpuConfig = configuration("cuda");
  const auto gpuFile = writeConfiguration(gpuConfig, "cuda.json");
  const auto cpuFile = writeConfiguration(configuration("cpu"), "cpu.json");
  ASSERT_EQ(pgo::SimulationRunner::runSimulationFromConfig(gpuFile), 0);
  ASSERT_EQ(pgo::SimulationRunner::runSimulationFromConfig(cpuFile), 0);

  ES::MXd gpu, cpu;
  ASSERT_EQ(ES::readMatrix((directory / "cuda/deform0002.u").string().c_str(), gpu), 0);
  ASSERT_EQ(ES::readMatrix((directory / "cpu/deform0002.u").string().c_str(), cpu), 0);
  ASSERT_EQ(gpu.rows(), cpu.rows());
  ASSERT_EQ(gpu.cols(), 3);
  ASSERT_TRUE(gpu.allFinite());
  EXPECT_LT((gpu.col(0) - cpu.col(0)).cwiseAbs().maxCoeff(), 2e-7);
  EXPECT_LT((gpu.col(1) - cpu.col(1)).cwiseAbs().maxCoeff(), 2e-4);

  pgo::ConfigFileJSON parsed(gpuFile.string().c_str());
  auto context = pgo::RunIPCSim::buildVolumeIpcSimulation(parsed);
  EXPECT_NO_THROW(context.collisionHandler->validateCollisionFreeState(gpu.col(0)));
  auto minimumGap = [&](ES::ConstRefVecXd displacement) {
    const ES::VXd surface = context.surfaceRestPositions + context.surfaceFromSimulationDispMap * displacement;
    double gap = 1.0;
    for (int i = 1; i < surface.size(); i += 3)
      gap = std::min(gap, surface[i]);
    return gap;
  };
  EXPECT_GT(minimumGap(gpu.col(0)), 0.0);
  // The first frame exercises active contact. Later frames may leave the
  // barrier region as the elastic box rebounds from the plane.
  ES::MXd first;
  ASSERT_EQ(ES::readMatrix((directory / "cuda/deform0000.u").string().c_str(), first), 0);
  EXPECT_NO_THROW(context.collisionHandler->validateCollisionFreeState(first.col(0)));
  EXPECT_GT(minimumGap(first.col(0)), 0.0);
  EXPECT_LT(minimumGap(first.col(0)), 0.002);
  EXPECT_TRUE(fs::exists(directory / "cuda/ret0002.obj"));
}

TEST_F(CudaIPCRunnerTest, RestartMatchesContinuousRun)
{
  auto config = configuration("cuda", 2);
  auto file = writeConfiguration(config, "restart.json");
  ASSERT_EQ(pgo::SimulationRunner::runSimulationFromConfig(file), 0);
  const auto originalRows = readCSV(directory / "cuda/newton-substeps.csv");
  ASSERT_FALSE(originalRows.empty());
  std::map<fs::path, std::string> originalAssets;
  for (const auto &row : originalRows) {
    for (const fs::path &asset : { directory / "cuda" / row.at("state_file"), substepSurfaceFile(directory / "cuda", row) }) {
      ASSERT_TRUE(fs::exists(asset));
      originalAssets.emplace(asset, readFile(asset));
    }
  }
  config["num-timestep"] = 4;
  config["restart-from-u"] = true;
  file = writeConfiguration(config, "restart.json");
  ASSERT_EQ(pgo::SimulationRunner::runSimulationFromConfig(file), 0);
  const auto resumedRows = readCSV(directory / "cuda/newton-substeps.csv");
  ASSERT_GE(resumedRows.size(), originalRows.size());
  EXPECT_TRUE(std::equal(originalRows.begin(), originalRows.end(), resumedRows.begin()));
  for (const auto &[asset, contents] : originalAssets) {
    ASSERT_TRUE(fs::exists(asset));
    EXPECT_EQ(readFile(asset), contents);
  }
  checkPhysicalSurfaceExports(file, config);
  // A restart must retain the accepted physical path of the earlier frames,
  // including its manifest, so a full independent replay remains possible.
  const std::string validationCommand = shellQuote(PGO_TEST_VALIDATE_IPC_TRAJECTORY_BIN) + " " + shellQuote(file) +
    " --swept --report " + shellQuote(directory / "restart-collision-validation.json") +
    " > " + shellQuote(directory / "restart-validation.log") + " 2>&1";
  ASSERT_EQ(std::system(validationCommand.c_str()), 0);
  config["restart-from-u"] = false;
  config["output"] = (directory / "continuous").string();
  file = writeConfiguration(config, "continuous.json");
  ASSERT_EQ(pgo::SimulationRunner::runSimulationFromConfig(file), 0);
  ES::MXd resumed, continuous;
  ASSERT_EQ(ES::readMatrix((directory / "cuda/deform0003.u").string().c_str(), resumed), 0);
  ASSERT_EQ(ES::readMatrix((directory / "continuous/deform0003.u").string().c_str(), continuous), 0);
  EXPECT_LT((resumed.col(0) - continuous.col(0)).cwiseAbs().maxCoeff(), 1e-9);
  EXPECT_LT((resumed.col(1) - continuous.col(1)).cwiseAbs().maxCoeff(), 1e-6);
}

TEST_F(CudaIPCRunnerTest, BoundedInvocationPreservesMovingObjectsAndPullingSchedule)
{
  using pgo::SimulationRunner::IPCSimulationRunOptions;
  const auto cloth = directory / "bounded-cloth.obj";
  std::ofstream(cloth) << "v -.1 .019 -.1\nv .1 .019 -.1\nv .1 .019 .1\nv -.1 .019 .1\nf 1 3 2\nf 1 4 3\n";
  // A single large triangle avoids an almost-parallel floor diagonal and
  // isolates the original frame-dependent pulling and obstacle schedules.
  std::ofstream(directory / "plane.obj") << "v -10 0 -10\nv 10 0 -10\nv 0 0 10\nf 1 3 2\n";
  std::ofstream(directory / "bounded-fixed.txt") << "0\n";
  for (const std::string backend : { "cpu", "cuda" }) {
    SCOPED_TRACE(backend);
    auto config = configuration(backend, 5);
    config.erase("tet-mesh");
    config["surface-mesh"] = cloth.string();
    config["elastic-material"] = "koiter-stvk";
    config["shell-youngs-modulus"] = 1000.0;
    config["shell-bending-youngs-modulus"] = 1000.0;
    config["shell-areal-density"] = 1.0;
    config["init-disp"] = { 0, 0, 0 };
    config["init-vel"] = { .2, 0, 0 };
    config["timestep"] = .002;
    config["ipc-dhat"] = .02;
    config["ipc-kappa"] = 1.0;
    config["contact-friction-coeff"] = .1;
    config["ipc-friction-iterations"] = 2;
    config["fixed-vertices"] = nlohmann::json::array({ { { "filename", (directory / "bounded-fixed.txt").string() },
      { "movement", { .002, 0, 0 } }, { "coeff", 10.0 } } });
    config["ipc-external-objects"][0]["movement"] = { 0, .002, 0 };
    const auto fullOutput = directory / (backend + "-full");
    config["output"] = fullOutput.string();
    auto file = writeConfiguration(config, backend + "-bounded.json");
    ASSERT_EQ(pgo::SimulationRunner::runIPCSimulationFromConfig(file), 0);

    const auto boundedOutput = directory / (backend + "-bounded");
    config["output"] = boundedOutput.string();
    file = writeConfiguration(config, backend + "-bounded.json");
    IPCSimulationRunOptions options;
    options.maxFrames = 2;
    options.timingOutput = directory / (backend + "-prefix-timing.csv");
    ASSERT_EQ(pgo::SimulationRunner::runIPCSimulationFromConfig(file, false, options), 0);
    EXPECT_TRUE(fs::exists(boundedOutput / "deform0001.u"));
    EXPECT_FALSE(fs::exists(boundedOutput / "deform0002.u"));
    auto checkTiming = [&](const fs::path &path, int first, int count) {
      const auto rows = readCSV(path);
      ASSERT_EQ(rows.size(), static_cast<std::size_t>(count));
      for (int i = 0; i < count; ++i) {
        const auto &row = rows[i];
        EXPECT_EQ(std::stoi(row.at("frame")), first + i);
        EXPECT_EQ(row.at("terminal_status"), "0");
        EXPECT_EQ(row.at("all_solver_status_converged"), "1");
        EXPECT_EQ(row.at("accepted_solvers_converged"), "1");
        EXPECT_EQ(row.at("committed"), "1");
        EXPECT_EQ(std::stoi(row.at("newton_solves")), 2);
        const double compute = std::stod(row.at("compute_seconds"));
        const double output = std::stod(row.at("output_seconds"));
        const double total = std::stod(row.at("total_seconds"));
        EXPECT_GE(compute, 0);
        EXPECT_GT(output, 0);
        EXPECT_NEAR(compute + output, total, 1e-12);
      }
    };
    checkTiming(options.timingOutput, 0, 2);
    config["restart-from-u"] = true;
    file = writeConfiguration(config, backend + "-bounded.json");
    options.timingOutput = directory / (backend + "-middle-timing.csv");
    ASSERT_EQ(pgo::SimulationRunner::runIPCSimulationFromConfig(file, false, options), 0);
    EXPECT_TRUE(fs::exists(boundedOutput / "deform0003.u"));
    EXPECT_FALSE(fs::exists(boundedOutput / "deform0004.u"));
    checkTiming(options.timingOutput, 2, 2);
    options.timingOutput = directory / (backend + "-last-timing.csv");
    ASSERT_EQ(pgo::SimulationRunner::runIPCSimulationFromConfig(file, backend == "cuda", options), 0);
    checkTiming(options.timingOutput, 4, 1);
    if (backend == "cuda") {
      EXPECT_NE(readFile(boundedOutput / "runIPCSim.log").find("CUDA IPC complete: frames=1,"), std::string::npos);
    }
    for (int frame = 0; frame < 5; ++frame) {
      ES::MXd full, bounded;
      const auto filename = "deform000" + std::to_string(frame) + ".u";
      ASSERT_EQ(ES::readMatrix((fullOutput / filename).string().c_str(), full), 0);
      ASSERT_EQ(ES::readMatrix((boundedOutput / filename).string().c_str(), bounded), 0);
      EXPECT_LT((full.col(0) - bounded.col(0)).cwiseAbs().maxCoeff(), 1e-9);
      EXPECT_LT((full.col(1) - bounded.col(1)).cwiseAbs().maxCoeff(), 1e-6);
      EXPECT_LT((full.col(2) - bounded.col(2)).cwiseAbs().maxCoeff(), 1e-3);
    }
    // Shortening the configured schedule is deliberately different: both the
    // pulling target and plane already reach their final placement at frame1.
    config["restart-from-u"] = false;
    config["num-timestep"] = 2;
    const auto shortenedOutput = directory / (backend + "-shortened");
    config["output"] = shortenedOutput.string();
    file = writeConfiguration(config, backend + "-shortened.json");
    ASSERT_EQ(pgo::SimulationRunner::runIPCSimulationFromConfig(file), 0);
    ES::MXd original, shortened;
    ASSERT_EQ(ES::readMatrix((fullOutput / "deform0001.u").string().c_str(), original), 0);
    ASSERT_EQ(ES::readMatrix((shortenedOutput / "deform0001.u").string().c_str(), shortened), 0);
    EXPECT_GT((original.col(0) - shortened.col(0)).cwiseAbs().maxCoeff(), 1e-8);
  }
}

TEST_F(CudaIPCRunnerTest, TimingReportsUnconvergedCPUAndUncommittedCUDA)
{
  for (const std::string backend : { "cpu", "cuda" }) {
    SCOPED_TRACE(backend);
    auto config = configuration(backend, 3);
    config["solver-max-iter"] = 1;
    config["solver-eps"] = 1e-30;
    if (backend == "cuda")
      config["cuda-max-timestep-refinements"] = 0;
    const auto file = writeConfiguration(config, backend + "-failure-timing.json");
    pgo::SimulationRunner::IPCSimulationRunOptions options;
    options.maxFrames = 1;
    options.timingOutput = directory / (backend + "-failure-timing.csv");
    const int code = pgo::SimulationRunner::runIPCSimulationFromConfig(file, false, options);
    EXPECT_EQ(code, backend == "cpu" ? 0 : 1);
    const auto rows = readCSV(options.timingOutput);
    ASSERT_EQ(rows.size(), 1);
    EXPECT_EQ(rows[0].at("frame"), "0");
    EXPECT_EQ(rows[0].at("terminal_status"), "1");
    EXPECT_EQ(rows[0].at("all_solver_status_converged"), "0");
    EXPECT_EQ(rows[0].at("accepted_solvers_converged"), "0");
    EXPECT_EQ(rows[0].at("committed"), backend == "cpu" ? "1" : "0");
    EXPECT_EQ(rows[0].at("newton_solves"), "1");
    EXPECT_TRUE(std::isfinite(std::stod(rows[0].at("total_seconds"))));
  }
}

TEST_F(CudaIPCRunnerTest, BoundedStaticTimingIncludesZeroIterationSolve)
{
  for (const std::string backend : { "cpu", "cuda" }) {
    SCOPED_TRACE(backend);
    auto config = configuration(backend, 5);
    config["sim-type"] = "static";
    config["g"] = { 0, 0, 0 };
    config["init-disp"] = { 0, 0, 0 };
    const auto file = writeConfiguration(config, backend + "-static-timing.json");
    pgo::SimulationRunner::IPCSimulationRunOptions options;
    options.maxFrames = 3;
    options.timingOutput = directory / (backend + "-static-timing.csv");
    ASSERT_EQ(pgo::SimulationRunner::runIPCSimulationFromConfig(file, false, options), 0);
    const auto rows = readCSV(options.timingOutput);
    ASSERT_EQ(rows.size(), 1);
    EXPECT_EQ(rows[0].at("frame"), "0");
    EXPECT_EQ(rows[0].at("terminal_status"), "0");
    EXPECT_EQ(rows[0].at("all_solver_status_converged"), "1");
    EXPECT_EQ(rows[0].at("accepted_solvers_converged"), "1");
    EXPECT_EQ(rows[0].at("committed"), "1");
    EXPECT_EQ(rows[0].at("newton_solves"), "1");
    EXPECT_GT(std::stod(rows[0].at("output_seconds")), 0);
    EXPECT_FALSE(fs::exists(directory / backend / "deform0001.u"));
    if (backend == "cuda") {
      const auto path = nlohmann::json::parse(readFile(directory / backend / "static-newton-path.json"));
      EXPECT_TRUE(path.at("complete"));
      EXPECT_EQ(path.at("solves").at(0).at("newton_iterations"), 0);
      EXPECT_EQ(path.at("states").size(), 1);
    }
  }
}

TEST_F(CudaIPCRunnerTest, RayleighDampingMatchesIndependentCPUReference)
{
  // A small tensile state with shear isolates K(u_previous) damping. Build a
  // CPU reference directly from the backward-Euler residual and sparse LDLT:
  // the production CPU energy line search can stall on scalar cancellation
  // before its residual converges. The original near-parallel contact fixture
  // is retained separately below with strict accepted-substep checks.
  auto gpuConfig = configuration("cuda", 4);
  gpuConfig["damping-params"] = { 0.01, 0.01 };
  gpuConfig["init-disp"] = { 0, 0, 0 };
  gpuConfig["restart-from-u"] = true;
  const auto gpuFile = writeConfiguration(gpuConfig, "cuda-damping.json");
  pgo::ConfigFileJSON parsed(gpuFile.string().c_str());
  auto context = pgo::RunIPCSim::buildVolumeIpcSimulation(parsed);
  const int n = static_cast<int>(context.simulationRestPosition.size());
  ES::MXd initial = ES::MXd::Zero(n, 3);
  ES::M3d strain = ES::M3d::Zero();
  strain.diagonal() << 0.02, 0.04, 0.03;
  strain(0, 1) = 0.025;
  ES::V3d center = ES::V3d::Zero();
  for (int i = 0; i < n / 3; ++i)
    center += context.simulationRestPosition.segment<3>(3 * i);
  center /= n / 3;
  for (int i = 0; i < n / 3; ++i) {
    const ES::V3d relative = context.simulationRestPosition.segment<3>(3 * i) - center;
    initial.col(0).segment<3>(3 * i) = strain * relative;
    initial.col(1).segment<3>(3 * i) = 0.1 * relative;
  }
  fs::create_directories(directory / "cuda");
  ASSERT_EQ(ES::writeMatrix((directory / "cuda/deform0000.u").string().c_str(), initial), 0);
  ASSERT_EQ(pgo::SimulationRunner::runSimulationFromConfig(gpuFile), 0);
  ES::MXd gpu, cpu;
  ES::MXd gpuPrevious = initial, cpuPrevious = initial;
  ES::VXd gravity = ES::VXd::Zero(n);
  for (int i = 0; i < n / 3; ++i)
    gravity[3 * i + 1] = -9.81;
  const double h = 0.001, alpha = 0.01, beta = 0.01;
  auto residual = [&](const ES::MXd &state, const ES::MXd &previous) {
    ES::SpMatD stiffness;
    context.elasticEnergy->createHessian(stiffness);
    context.elasticEnergy->hessian(previous.col(0), stiffness);
    ES::VXd gradient(n);
    context.elasticEnergy->gradient(state.col(0), gradient);
    const ES::VXd delta = state.col(0) - previous.col(0);
    gradient += context.M * ((1 / (h * h) + alpha / h) * delta - previous.col(1) / h - gravity);
    gradient += (beta / h) * (stiffness * delta);
    return gradient.cwiseAbs().maxCoeff();
  };
  ES::SpMatD initialStiffness;
  context.elasticEnergy->createHessian(initialStiffness);
  context.elasticEnergy->hessian(initial.col(0), initialStiffness);
  for (int frame = 1; frame < 4; ++frame) {
    ES::SpMatD previousStiffness;
    context.elasticEnergy->createHessian(previousStiffness);
    context.elasticEnergy->hessian(cpuPrevious.col(0), previousStiffness);
    if (frame > 1) {
      EXPECT_GT((previousStiffness - initialStiffness).norm(), 1e-8 * initialStiffness.norm());
    }
    const ES::SpMatD fixedHessian = (1 / (h * h) + alpha / h) * context.M + (beta / h) * previousStiffness;
    const ES::VXd linear = -context.M * (cpuPrevious.col(1) / h + gravity);
    ES::VXd reference = cpuPrevious.col(0), gradient(n);
    for (int iteration = 0; iteration < 200; ++iteration) {
      context.elasticEnergy->gradient(reference, gradient);
      gradient += fixedHessian * (reference - cpuPrevious.col(0)) + linear;
      if (gradient.cwiseAbs().maxCoeff() < 1e-6)
        break;
      ES::SpMatD hessian;
      context.elasticEnergy->createHessian(hessian);
      context.elasticEnergy->hessian(reference, hessian);
      hessian += fixedHessian;
      ES::SymSolver factorization;
      factorization.compute(hessian);
      ASSERT_EQ(factorization.info(), Eigen::Success);
      const ES::VXd correction = factorization.solve(-gradient);
      ASSERT_TRUE(correction.allFinite());
      const double bound = context.elasticEnergy->computeMaxStepSize(reference, correction);
      ASSERT_GT(bound, 0);
      ASSERT_LE(bound, 1);
      reference += bound * correction;
    }
    cpu.resize(n, 3);
    cpu.col(0) = reference;
    cpu.col(1) = (reference - cpuPrevious.col(0)) / h;
    cpu.col(2) = (cpu.col(1) - cpuPrevious.col(1)) / h;
    const std::string name = "deform000" + std::to_string(frame) + ".u";
    ASSERT_EQ(ES::readMatrix((directory / "cuda" / name).string().c_str(), gpu), 0);
    ASSERT_TRUE(gpu.allFinite());
    const ES::VXd surface = context.surfaceRestPositions + context.surfaceFromSimulationDispMap * reference;
    for (int i = 1; i < surface.size(); i += 3)
      ASSERT_GT(surface[i], 0.1);  // The reference contains no plane contact.
    EXPECT_LT((gpu.col(0) - cpu.col(0)).cwiseAbs().maxCoeff(), 2e-7);
    EXPECT_LT((gpu.col(1) - cpu.col(1)).cwiseAbs().maxCoeff(), 2e-4);
    EXPECT_LT(residual(gpu, gpuPrevious), 1e-6);
    EXPECT_LT(residual(cpu, cpuPrevious), 1e-6);
    gpuPrevious = gpu;
    cpuPrevious = cpu;
  }
  const auto rows = readCSV(directory / "cuda/newton-solves.csv");
  ASSERT_EQ(rows.size(), 3);
  for (const auto &row : rows) {
    EXPECT_EQ(row.at("status"), "0");
    EXPECT_LT(std::stod(row.at("gradient_inf")), 1e-6);
    EXPECT_EQ(row.at("substeps"), "1");
  }
}

TEST_F(CudaIPCRunnerTest, NearParallelDampedContactCommitsOnlyConvergedCertifiedSubsteps)
{
  auto config = configuration("cuda", 1);
  config["damping-params"] = { 0.01, 0.01 };
  const auto file = writeConfiguration(config, "damped-near-parallel.json");
  ASSERT_EQ(pgo::SimulationRunner::runSimulationFromConfig(file), 0);
  const auto macros = readCSV(directory / "cuda/newton-solves.csv");
  ASSERT_EQ(macros.size(), 1);
  EXPECT_EQ(macros.front().at("status"), "0");
  // Coordinate recovery can make the original full timestep converge. The
  // same strict residual and geometry requirements apply to one or many steps.
  EXPECT_GE(std::stoi(macros.front().at("substeps")), 1);
  const auto substeps = readCSV(directory / "cuda/newton-substeps.csv");
  ASSERT_EQ(substeps.size(), std::stoul(macros.front().at("substeps")));
  for (const auto &step : substeps) {
    EXPECT_EQ(step.at("status"), "0");
    EXPECT_LT(std::stod(step.at("gradient_inf")), 1e-6);
    EXPECT_EQ(std::stod(step.at("tolerance")), 1e-6);
    EXPECT_LE(std::stoi(step.at("newton_iterations")), 200);
    EXPECT_TRUE(fs::exists(directory / "cuda" / step.at("state_file")));
  }
  checkPhysicalSurfaceExports(file, config);
  const auto report = directory / "collision-validation.json";
  const std::string command = shellQuote(PGO_TEST_VALIDATE_IPC_TRAJECTORY_BIN) + " " + shellQuote(file) +
    " --swept --report " + shellQuote(report) + " > " + shellQuote(directory / "validation.log") + " 2>&1";
  ASSERT_EQ(std::system(command.c_str()), 0);
  std::ifstream reportStream(report);
  const auto validation = nlohmann::json::parse(reportStream);
  EXPECT_TRUE(validation.at("passed").get<bool>());
}

TEST_F(CudaIPCRunnerTest, StVKVolumeRunsNonrigidDampedStep)
{
  auto config = configuration("cuda", 2);
  config["elastic-material"] = "stvk-vol";
  config["damping-params"] = { 0.01, 0.01 };
  config["init-disp"] = { 0, 0, 0 };
  config["restart-from-u"] = true;
  const auto file = writeConfiguration(config, "stvk-damping.json");
  pgo::ConfigFileJSON parsed(file.string().c_str());
  auto context = pgo::RunIPCSim::buildVolumeIpcSimulation(parsed);
  const int n = static_cast<int>(context.simulationRestPosition.size());
  ES::MXd initial = ES::MXd::Zero(n, 3);
  const ES::V3d strain(0.02, -0.01, 0.01);
  for (int vertex = 0; vertex < n / 3; ++vertex)
    initial.col(0).segment<3>(3 * vertex) = strain.cwiseProduct(context.simulationRestPosition.segment<3>(3 * vertex));
  const double initialElasticEnergy = context.elasticEnergy->func(initial.col(0));
  ASSERT_GT(initialElasticEnergy, 0);
  fs::create_directories(directory / "cuda");
  ASSERT_EQ(ES::writeMatrix((directory / "cuda/deform0000.u").string().c_str(), initial), 0);
  ASSERT_EQ(pgo::SimulationRunner::runSimulationFromConfig(file), 0);
  ES::MXd finalState;
  ASSERT_EQ(ES::readMatrix((directory / "cuda/deform0001.u").string().c_str(), finalState), 0);
  ASSERT_TRUE(finalState.allFinite());
  EXPECT_GT((finalState.col(0) - initial.col(0)).norm(), 1e-8);
  EXPECT_LT(context.elasticEnergy->func(finalState.col(0)), initialElasticEnergy);
  const auto rows = readCSV(directory / "cuda/newton-solves.csv");
  ASSERT_EQ(rows.size(), 1);
  EXPECT_EQ(rows.front().at("status"), "0");
  EXPECT_LT(std::stod(rows.front().at("gradient_inf")), 1e-6);
  checkPhysicalSurfaceExports(file, config);
}

TEST_F(CudaIPCRunnerTest, KoiterShellUsesInitialSurfaceAndDamping)
{
  // A curved two-triangle hinge exercises membrane and bending parameters,
  // boundary stencils, nonuniform initial displacement and the independent
  // frozen stiffness cache used by Rayleigh damping.
  const std::vector<ES::V3d> rest = {
    { -0.1, 0.3, -0.1 }, { 0.1, 0.3, -0.1 },
    { -0.1, 0.3, 0.1 }, { 0.1, 0.31, 0.1 }
  };
  std::vector<ES::V3d> initial = rest;
  for (auto &position : initial) {
    position[0] *= 1.025;
    position[2] *= 0.99;
  }
  initial[3][1] += 0.004;
  auto writeSurface = [&](const fs::path &path, const std::vector<ES::V3d> &positions) {
    std::ofstream stream(path);
    stream << std::setprecision(17);
    for (const auto &position : positions)
      stream << "v " << position.transpose() << '\n';
    stream << "f 1 3 2\nf 2 3 4\n";
  };
  writeSurface(directory / "shell-rest.obj", rest);
  writeSurface(directory / "shell-initial.obj", initial);
  auto config = configuration("cuda", 3);
  config.erase("tet-mesh");
  config["surface-mesh"] = (directory / "shell-rest.obj").string();
  config["initial-surface-mesh"] = (directory / "shell-initial.obj").string();
  config["elastic-material"] = "koiter-stvk";
  config["shell-youngs-modulus"] = 10000.0;
  config["shell-bending-youngs-modulus"] = 30000.0;
  config["shell-thickness"] = 0.01;
  config["shell-areal-density"] = 2.0;
  config["shell-zero-rest-curvature"] = true;
  config["damping-params"] = { 0.01, 0.003 };
  config["g"] = { 0, 0, 0 };
  config["init-disp"] = { 0, 0.01, 0 };
  config["scale"] = 1.25;
  config["timestep"] = 0.002;
  const auto file = writeConfiguration(config, "koiter-damping.json");
  pgo::ConfigFileJSON parsed(file.string().c_str());
  auto context = pgo::RunIPCSim::buildShellIpcSimulation(parsed);
  ASSERT_TRUE(context.zeroShellRestCurvature);
  for (int vertex = 0; vertex < static_cast<int>(rest.size()); ++vertex) {
    const ES::V3d expected = 1.25 * (initial[vertex] - rest[vertex]) + ES::V3d(0, 0.01, 0);
    EXPECT_LT((context.initialDisplacement.segment<3>(3 * vertex) - expected).norm(), 1e-15);
  }
  const double initialElasticEnergy = context.elasticEnergy->func(context.initialDisplacement);
  ASSERT_GT(initialElasticEnergy, 0);
  ASSERT_EQ(pgo::SimulationRunner::runSimulationFromConfig(file), 0);
  ES::MXd finalState;
  ASSERT_EQ(ES::readMatrix((directory / "cuda/deform0002.u").string().c_str(), finalState), 0);
  ASSERT_TRUE(finalState.allFinite());
  EXPECT_GT((finalState.col(0) - context.initialDisplacement).norm(), 1e-8);
  EXPECT_LT(context.elasticEnergy->func(finalState.col(0)), initialElasticEnergy);
  const auto rows = readCSV(directory / "cuda/newton-solves.csv");
  ASSERT_EQ(rows.size(), 3);
  for (const auto &row : rows) {
    EXPECT_EQ(row.at("status"), "0");
    EXPECT_LT(std::stod(row.at("gradient_inf")), 1e-6);
  }
  const auto physicalRows = readCSV(directory / "cuda/newton-substeps.csv");
  ASSERT_GE(physicalRows.size(), rows.size());
  for (const auto &row : physicalRows) {
    EXPECT_EQ(row.at("status"), "0");
    EXPECT_DOUBLE_EQ(std::stod(row.at("tolerance")), 1e-6);
    EXPECT_LT(std::stod(row.at("gradient_inf")), 1e-6);
    EXPECT_LE(std::stoi(row.at("newton_iterations")), 200);
  }
  checkPhysicalSurfaceExports(file, config);
  const auto report = directory / "shell-collision-validation.json";
  const std::string command = shellQuote(PGO_TEST_VALIDATE_IPC_TRAJECTORY_BIN) + " " + shellQuote(file) +
    " --swept --report " + shellQuote(report) + " > " + shellQuote(directory / "shell-validation.log") + " 2>&1";
  ASSERT_EQ(std::system(command.c_str()), 0);
  std::ifstream reportStream(report);
  EXPECT_TRUE(nlohmann::json::parse(reportStream).at("passed").get<bool>());

  // A static solve records its accepted optimization path, independently of
  // physical timesteps. Intermediate knots may have a large residual; only
  // the terminal solve state must meet the original tolerance.
  config["sim-type"] = "static";
  config["output"] = (directory / "static-shell").string();
  const auto staticFile = writeConfiguration(config, "static-koiter.json");
  ASSERT_EQ(pgo::SimulationRunner::runSimulationFromConfig(staticFile), 0);
  std::ifstream pathStream(directory / "static-shell/static-newton-path.json");
  const auto path = nlohmann::json::parse(pathStream);
  EXPECT_TRUE(path.at("complete").get<bool>());
  EXPECT_FALSE(path.at("physical_time").get<bool>());
  EXPECT_EQ(path.at("path_type"), "static_newton_iterates");
  const auto &solves = path.at("solves");
  ASSERT_EQ(solves.size(), 1);
  EXPECT_EQ(solves[0].at("status"), 0);
  EXPECT_LT(solves[0].at("gradient_inf").get<double>(), 1e-6);
  const auto &knots = path.at("states");
  ASSERT_EQ(knots.size(), solves[0].at("newton_iterations").get<int>() + 1);
  ASSERT_GT(knots.size(), 1);
  EXPECT_EQ(solves[0].at("first_state"), 0);
  EXPECT_EQ(solves[0].at("last_state"), knots.size() - 1);
  ES::MXd lastKnot;
  for (std::size_t index = 0; index < knots.size(); ++index) {
    EXPECT_EQ(knots[index].at("iteration"), index);
    EXPECT_EQ(knots[index].at("solve"), 0);
    EXPECT_TRUE(std::isfinite(knots[index].at("energy").get<double>()));
    EXPECT_TRUE(std::isfinite(knots[index].at("gradient_inf").get<double>()));
    ASSERT_EQ(ES::readMatrix((directory / "static-shell" / knots[index].at("state_file").get<std::string>()).string().c_str(), lastKnot), 0);
    ASSERT_TRUE(lastKnot.allFinite());
    if (index == 0) {
      EXPECT_LT((lastKnot.col(0) - context.initialDisplacement).norm(), 1e-15);
    }
  }
  EXPECT_GT(knots[0].at("gradient_inf").get<double>(), 1e-6);
  ES::MXd staticState;
  ASSERT_EQ(ES::readMatrix((directory / "static-shell/deform0000.u").string().c_str(), staticState), 0);
  EXPECT_EQ(lastKnot, staticState);
  checkPhysicalSurfaceExports(staticFile, config);
  const std::string staticCommand = shellQuote(PGO_TEST_VALIDATE_IPC_TRAJECTORY_BIN) + " " + shellQuote(staticFile) +
    " --swept --report " + shellQuote(directory / "static-shell-validation.json") +
    " > " + shellQuote(directory / "static-shell-validation.log") + " 2>&1";
  ASSERT_EQ(std::system(staticCommand.c_str()), 0);
}

TEST_F(CudaIPCRunnerTest, QuantizedDampingCheckpointConvergesAtOriginalTolerance)
{
  // Preserve the original frame 172->173 force history while making this a
  // one-step restart test: the original linear track has 1000 intervals and
  // movement -.1. A 173-interval track with movement -.1*173/1000 gives identical
  // targets throughout the tested interval. The full example keeps 1001 frames.
  constexpr int frame = 173;
  constexpr double originalMovement = -0.1;
  const double shortenedMovement = originalMovement * frame / 1000;
  for (int endpoint : { frame - 1, frame }) {
    const double originalTarget = originalMovement * endpoint / 1000;
    const double shortenedTarget = shortenedMovement * endpoint / frame;
    EXPECT_NEAR(shortenedTarget, originalTarget, 4 * std::numeric_limits<double>::epsilon() * std::abs(originalTarget));
  }
  auto config = configuration("cuda", frame + 1);
  const auto assets = fs::path(PGO_CUDA_TEST_BOX_DIR).parent_path() / "box-contact";
  config["fixed-vertices"] = nlohmann::json::array({ { { "filename", (assets / "top-fixed.txt").string() },
    { "movement", { 0, shortenedMovement, 0 } }, { "coeff", 100000 } } });
  config["ipc-external-objects"][0]["filename"] = (assets / "plane.obj").string();
  config["g"] = { 0, 0, 0 };
  config["init-disp"] = { 0, 0, 0 };
  config["damping-params"] = { 0.01, 0.01 };
  config["solver-max-iter"] = 300;
  config["restart-from-u"] = true;
  fs::create_directories(directory / "cuda");
  fs::copy_file(fs::path(PGO_CUDA_TEST_FIXTURE_DIR) / "damping-contact-0172.u", directory / "cuda/deform0172.u");
  const auto file = writeConfiguration(config, "damping-checkpoint.json");
  ASSERT_EQ(pgo::SimulationRunner::runSimulationFromConfig(file), 0);
  const auto rows = readCSV(directory / "cuda/newton-solves.csv");
  ASSERT_EQ(rows.size(), 1);
  const auto &macro = rows.front();
  EXPECT_EQ(std::stoi(macro.at("frame")), frame);
  EXPECT_EQ(macro.at("status"), "0");
  EXPECT_DOUBLE_EQ(std::stod(macro.at("tolerance")), 1e-6);
  EXPECT_TRUE(std::isfinite(std::stod(macro.at("gradient_inf"))));
  EXPECT_LT(std::stod(macro.at("gradient_inf")), 1e-6);
  const int count = std::stoi(macro.at("substeps"));
  ASSERT_GE(count, 1);
  ASSERT_LE(count, 256);
  ASSERT_EQ(count & (count - 1), 0);
  const int refinement = std::stoi(macro.at("refinement_retries"));
  ASSERT_GE(refinement, 0);
  ASSERT_LE(refinement, 8);
  EXPECT_EQ(count, 1 << refinement);

  // Rounded residual recovery may succeed at the original timestep or after
  // adaptive refinement. Only the committed chain contributes to its counts;
  // converged solves in a rolled-back attempt are not accepted physical steps.
  const auto attempts = readCSV(directory / "cuda/newton-attempts.csv");
  std::map<int, std::map<std::string, std::string>> committed;
  int rejectedSolves = 0;
  for (const auto &attempt : attempts) {
    EXPECT_EQ(std::stoi(attempt.at("frame")), frame);
    if (attempt.at("accepted") == "0") {
      EXPECT_LT(std::stoi(attempt.at("refinement")), refinement);
      ++rejectedSolves;
      continue;
    }
    ASSERT_EQ(attempt.at("accepted"), "1");
    EXPECT_EQ(std::stoi(attempt.at("refinement")), refinement);
    EXPECT_EQ(std::stoi(attempt.at("planned_substeps")), count);
    EXPECT_EQ(attempt.at("friction_iteration"), "0");
    EXPECT_TRUE(attempt.at("rejection").empty());
    const int step = std::stoi(attempt.at("substep"));
    ASSERT_GE(step, 0);
    ASSERT_LT(step, count);
    ASSERT_TRUE(committed.emplace(step, attempt).second);
  }
  ASSERT_EQ(committed.size(), static_cast<std::size_t>(count));
  EXPECT_EQ(rejectedSolves, std::stoi(macro.at("rejected_newton_solves")));
  const auto substeps = readCSV(directory / "cuda/newton-substeps.csv");
  ASSERT_EQ(substeps.size(), static_cast<std::size_t>(count));
  ES::MXd previous;
  ASSERT_EQ(ES::readMatrix((fs::path(PGO_CUDA_TEST_FIXTURE_DIR) / "damping-contact-0172.u").string().c_str(), previous), 0);
  ASSERT_EQ(previous.cols(), 3);
  ASSERT_TRUE(previous.allFinite());
  const double h = 0.001 / count;
  int recoverySteps = 0, newtonIterations = 0, pcgIterations = 0, maxNewtonIterations = 0;
  double maxGradient = 0;
  for (int step = 0; step < count; ++step) {
    SCOPED_TRACE("accepted substep " + std::to_string(step));
    const auto &row = substeps[step];
    const auto &attempt = committed.at(step);
    EXPECT_EQ(std::stoi(row.at("frame")), frame);
    EXPECT_EQ(std::stoi(row.at("substep")), step);
    EXPECT_EQ(std::stoi(row.at("substeps")), count);
    for (const char *field : { "newton_iterations", "pcg_iterations", "gradient_inf", "tolerance", "status", "stop" })
      EXPECT_EQ(row.at(field), attempt.at(field));
    EXPECT_EQ(row.at("status"), "0");
    EXPECT_EQ(row.at("stop"), "gradient_tolerance");
    EXPECT_DOUBLE_EQ(std::stod(row.at("tolerance")), 1e-6);
    const double gradient = std::stod(row.at("gradient_inf"));
    EXPECT_TRUE(std::isfinite(gradient));
    EXPECT_GE(gradient, 0);
    EXPECT_LT(gradient, 1e-6);
    const int iterations = std::stoi(row.at("newton_iterations"));
    EXPECT_GE(iterations, 0);
    EXPECT_LE(iterations, 300);
    const int corrections = std::stoi(attempt.at("residual_correction_steps"));
    EXPECT_GE(corrections, 0);
    EXPECT_LE(corrections, iterations);
    EXPECT_GE(std::stoi(row.at("pcg_iterations")), 0);
    recoverySteps += corrections;
    newtonIterations += iterations;
    pcgIterations += std::stoi(row.at("pcg_iterations"));
    maxNewtonIterations = std::max(maxNewtonIterations, iterations);
    maxGradient = std::max(maxGradient, gradient);

    ES::MXd state;
    ASSERT_EQ(ES::readMatrix((directory / "cuda" / row.at("state_file")).string().c_str(), state), 0);
    ASSERT_EQ(state.rows(), previous.rows());
    ASSERT_EQ(state.cols(), 3);
    ASSERT_TRUE(state.allFinite());
    const ES::VXd expectedVelocity = (state.col(0) - previous.col(0)) / h;
    const ES::VXd expectedAcceleration = (state.col(1) - previous.col(1)) / h;
    ASSERT_TRUE(expectedVelocity.allFinite());
    ASSERT_TRUE(expectedAcceleration.allFinite());
    const double roundoff = 16 * std::numeric_limits<double>::epsilon();
    EXPECT_LE((state.col(1) - expectedVelocity).cwiseAbs().maxCoeff(),
      roundoff * std::max(1.0, expectedVelocity.cwiseAbs().maxCoeff()));
    EXPECT_LE((state.col(2) - expectedAcceleration).cwiseAbs().maxCoeff(),
      roundoff * std::max(1.0, expectedAcceleration.cwiseAbs().maxCoeff()));
    previous = state;
  }
  EXPECT_GT(recoverySteps, 0);
  EXPECT_EQ(recoverySteps, std::stoi(macro.at("residual_correction_steps")));
  EXPECT_EQ(newtonIterations, std::stoi(macro.at("newton_iterations")));
  EXPECT_EQ(pcgIterations, std::stoi(macro.at("pcg_iterations")));
  EXPECT_EQ(maxNewtonIterations, std::stoi(macro.at("max_substep_newton_iterations")));
  EXPECT_DOUBLE_EQ(maxGradient, std::stod(macro.at("gradient_inf")));
  EXPECT_EQ(substeps.back().at("state_file"), "deform0173.u");
  ES::MXd corrected;
  ASSERT_EQ(ES::readMatrix((directory / "cuda/deform0173.u").string().c_str(), corrected), 0);
  ASSERT_TRUE(corrected.allFinite());
  EXPECT_EQ(corrected, previous);
  if (count == 1) {
    // This reference fixture belongs to the original full-h objective.
    // Refined steps have different inertia and must be checked as a chain.
    ES::MXd failed;
    ASSERT_EQ(ES::readMatrix((fs::path(PGO_CUDA_TEST_FIXTURE_DIR) / "damping-contact-0173-failed.u").string().c_str(), failed), 0);
    ASSERT_EQ(failed.rows(), corrected.rows());
    ASSERT_EQ(failed.cols(), 3);
    EXPECT_LT((corrected.col(0) - failed.col(0)).cwiseAbs().maxCoeff(), 2e-7);
  }
  checkPhysicalSurfaceExports(file, config);
  const std::string command = shellQuote(PGO_TEST_VALIDATE_IPC_TRAJECTORY_BIN) + " " + shellQuote(file) +
    " --swept --first-frame 173 --last-frame 173 --report " + shellQuote(directory / "collision-validation.json") +
    " > " + shellQuote(directory / "validation.log") + " 2>&1";
  ASSERT_EQ(std::system(command.c_str()), 0);
  std::ifstream reportStream(directory / "collision-validation.json");
  const auto validation = nlohmann::json::parse(reportStream);
  EXPECT_TRUE(validation.at("passed").get<bool>());
  EXPECT_EQ(validation.at("first_frame"), frame);
  EXPECT_EQ(validation.at("last_frame"), frame);
  EXPECT_TRUE(validation.at("physical_substeps").at("passed").get<bool>());
  EXPECT_EQ(validation.at("physical_substeps").at("checked"), count);
}

TEST_F(CudaIPCRunnerTest, FrictionLaggingCommitsOnceAndMatchesCPUAndRestart)
{
  // Keep all floor edges outside the barrier region. Nearly parallel
  // diagonals with eps_ee=0 give an ill-conditioned contact fixture in which
  // the existing CPU runner itself does not reach its gradient tolerance.
  // A single covering triangle isolates the friction update lifecycle.
  std::ofstream(directory / "friction-floor.obj") << "v -1 0 -1\nv 1 0 -1\nv 0 0 1\nf 1 3 2\n";
  std::ofstream(directory / "sliding-shell.obj") << "v -.1 .019 -.1\nv .1 .019 -.1\nv .1 .019 .1\nv -.1 .019 .1\nf 1 3 2\nf 1 4 3\n";
  auto config = configuration("cuda", 40);
  config.erase("tet-mesh");
  config["surface-mesh"] = (directory / "sliding-shell.obj").string();
  config["ipc-external-objects"][0]["filename"] = (directory / "friction-floor.obj").string();
  config["elastic-material"] = "koiter-stvk";
  config["shell-youngs-modulus"] = 1000.0;
  config["shell-bending-youngs-modulus"] = 1000.0;
  config["shell-areal-density"] = 1.0;
  config["init-disp"] = { 0, 0, 0 };
  config["init-vel"] = { 0.2, 0, 0 };
  config["timestep"] = 0.002;
  config["ipc-dhat"] = 0.02;
  config["ipc-kappa"] = 1.0;
  config["solver-eps"] = 1e-8;
  config["solver-max-iter"] = 100;
  config["contact-friction-coeff"] = 0.5;
  config["ipc-friction-iterations"] = 2;
  const auto file = writeConfiguration(config, "friction.json");
  ASSERT_EQ(pgo::SimulationRunner::runSimulationFromConfig(file), 0);
  const auto macros = readCSV(directory / "cuda/newton-solves.csv");
  const auto substeps = readCSV(directory / "cuda/newton-substeps.csv");
  const auto attempts = readCSV(directory / "cuda/newton-attempts.csv");
  ASSERT_EQ(macros.size(), 40);
  ASSERT_EQ(substeps.size(), 40);
  ASSERT_EQ(attempts.size(), 80);
  ES::MXd previous = ES::MXd::Zero(12, 3), finalState;
  for (int vertex = 0; vertex < 4; ++vertex)
    previous(3 * vertex, 1) = 0.2;
  for (int frame = 0; frame < 40; ++frame) {
    const auto &macro = macros[frame];
    EXPECT_EQ(macro.at("status"), "0");
    EXPECT_EQ(macro.at("substeps"), "1");
    EXPECT_EQ(macro.at("friction_solves"), "2");
    EXPECT_LT(std::stod(macro.at("gradient_inf")), 1e-8);
    for (int iteration = 0; iteration < 2; ++iteration) {
      const auto &attempt = attempts[2 * frame + iteration];
      EXPECT_EQ(std::stoi(attempt.at("frame")), frame);
      EXPECT_EQ(std::stoi(attempt.at("friction_iteration")), iteration);
      EXPECT_EQ(attempt.at("accepted"), "1");
      EXPECT_EQ(attempt.at("status"), "0");
      EXPECT_GT(std::stoul(attempt.at("friction_pairs")), 0);
      EXPECT_LT(std::stod(attempt.at("gradient_inf")), 1e-8);
      EXPECT_LE(std::stoi(attempt.at("newton_iterations")), 100);
    }
    ASSERT_EQ(ES::readMatrix((directory / "cuda" / substeps[frame].at("state_file")).string().c_str(), finalState), 0);
    ASSERT_TRUE(finalState.allFinite());
    // Both lagging solves share x_n and v_n. Physical time advances once,
    // so exported velocity and acceleration must use the full physical step.
    EXPECT_LT((finalState.col(1) - (finalState.col(0) - previous.col(0)) / 0.002).norm(), 1e-12);
    EXPECT_LT((finalState.col(2) - (finalState.col(1) - previous.col(1)) / 0.002).norm(), 1e-10);
    previous = finalState;
  }
  checkPhysicalSurfaceExports(file, config);
  const std::string validationCommand = shellQuote(PGO_TEST_VALIDATE_IPC_TRAJECTORY_BIN) + " " + shellQuote(file) +
    " --swept --report " + shellQuote(directory / "friction-collision-validation.json") +
    " > " + shellQuote(directory / "friction-validation.log") + " 2>&1";
  ASSERT_EQ(std::system(validationCommand.c_str()), 0);

  // Use CPU gradients and sparse LDLT as an independent reference. The CPU
  // runner's absolute-energy line search can stop on roundoff before 1e-8;
  // require the actual residual to converge at every lagging iteration.
  pgo::ConfigFileJSON parsed(file.string().c_str());
  auto cpuContext = pgo::RunIPCSim::buildShellIpcSimulation(parsed);
  cpuContext.collisionHandler->setDeterministic(true);
  ES::MXd cpuState = ES::MXd::Zero(12, 3);
  for (int vertex = 0; vertex < 4; ++vertex)
    cpuState(3 * vertex, 1) = 0.2;
  ES::VXd gravity = ES::VXd::Zero(12);
  for (int vertex = 0; vertex < 4; ++vertex)
    gravity[3 * vertex + 1] = -9.81;
  constexpr double h = 0.002;
  const ES::SpMatD inertia = cpuContext.M * (1.0 / (h * h));
  int convergedCPUOuterSolves = 0;
  for (int frame = 0; frame < 40; ++frame) {
    const ES::MXd previousCPU = cpuState;
    const ES::VXd linear = -cpuContext.M * (previousCPU.col(1) / h + gravity);
    ES::VXd reference = previousCPU.col(0);
    auto residual = [&](const ES::VXd &at) {
      ES::VXd gradient(12), contactGradient(12);
      cpuContext.elasticEnergy->gradient(at, gradient);
      cpuContext.collisionHandler->gradient(at, contactGradient);
      gradient += contactGradient + inertia * (at - previousCPU.col(0)) + linear;
      return gradient;
    };
    for (int outer = 0; outer < 2; ++outer) {
      SCOPED_TRACE("CPU reference frame " + std::to_string(frame) + ", outer " + std::to_string(outer));
      cpuContext.collisionHandler->updateFriction(previousCPU.col(0), reference, h);
      for (int iteration = 0; iteration < 100; ++iteration) {
        const ES::VXd gradient = residual(reference);
        ASSERT_TRUE(gradient.allFinite());
        if (gradient.cwiseAbs().maxCoeff() < 1e-8)
          break;
        ES::SpMatD hessian, contactHessian;
        cpuContext.elasticEnergy->createHessian(hessian);
        cpuContext.elasticEnergy->hessian(reference, hessian);
        cpuContext.collisionHandler->hessianDirect(reference, contactHessian);
        hessian += contactHessian + inertia;
        ES::SymSolver factorization;
        factorization.compute(hessian);
        ASSERT_EQ(factorization.info(), Eigen::Success);
        const ES::VXd direction = factorization.solve(-gradient);
        ASSERT_TRUE(direction.allFinite());
        double alpha = std::min(cpuContext.elasticEnergy->computeMaxStepSize(reference, direction),
          cpuContext.collisionHandler->computeMaxStepSize(reference, direction));
        ASSERT_GT(alpha, 0);
        bool accepted = false;
        for (int backtrack = 0; backtrack < 40; ++backtrack, alpha *= 0.5) {
          const ES::VXd trial = reference + alpha * direction;
          const ES::VXd trialGradient = residual(trial);
          if (trialGradient.allFinite() && trialGradient.norm() < gradient.norm()) {
            reference = trial;
            accepted = true;
            break;
          }
        }
        ASSERT_TRUE(accepted);
      }
      ASSERT_LT(residual(reference).cwiseAbs().maxCoeff(), 1e-8);
      ++convergedCPUOuterSolves;
    }
    EXPECT_NO_THROW(cpuContext.collisionHandler->validateCollisionFreeState(reference));
    EXPECT_DOUBLE_EQ(cpuContext.collisionHandler->computeMaxStepSize(previousCPU.col(0), reference - previousCPU.col(0)), 1.0);
    cpuState.col(0) = reference;
    cpuState.col(1) = (reference - previousCPU.col(0)) / h;
    cpuState.col(2) = (cpuState.col(1) - previousCPU.col(1)) / h;
  }
  EXPECT_EQ(convergedCPUOuterSolves, 80);
  EXPECT_LT((finalState.col(0) - cpuState.col(0)).norm(), 2e-7);
  EXPECT_LT((finalState.col(1) - cpuState.col(1)).norm(), 2e-4);
  double meanSlidingVelocity = 0;
  for (int vertex = 0; vertex < 4; ++vertex)
    meanSlidingVelocity += finalState(3 * vertex, 1) / 4;
  EXPECT_LT(std::abs(meanSlidingVelocity), 0.04);

  // Retain the first half of the accepted physical path and restart from its
  // final state. The friction reference must be reconstructed from that state.
  const auto restarted = directory / "restarted-friction";
  fs::create_directories(restarted / "substeps");
  std::ifstream originalManifest(directory / "cuda/newton-substeps.csv");
  std::ofstream restartedManifest(restarted / "newton-substeps.csv");
  std::string line;
  std::getline(originalManifest, line);
  restartedManifest << line << '\n';
  for (int frame = 0; frame < 20; ++frame) {
    std::getline(originalManifest, line);
    restartedManifest << line << '\n';
    const auto &row = substeps[frame];
    const fs::path state = row.at("state_file");
    fs::copy_file(directory / "cuda" / state, restarted / state);
    fs::copy_file(substepSurfaceFile(directory / "cuda", row), substepSurfaceFile(restarted, row));
  }
  restartedManifest.close();
  // Macrostate files and substep files can have distinct names.
  fs::copy_file(directory / "cuda/deform0019.u", restarted / "deform0019.u", fs::copy_options::overwrite_existing);
  config["output"] = restarted.string();
  config["restart-from-u"] = true;
  const auto restartFile = writeConfiguration(config, "restarted-friction.json");
  ASSERT_EQ(pgo::SimulationRunner::runSimulationFromConfig(restartFile), 0);
  ES::MXd resumed;
  ASSERT_EQ(ES::readMatrix((restarted / "deform0039.u").string().c_str(), resumed), 0);
  EXPECT_LT((resumed.col(0) - finalState.col(0)).norm(), 1e-8);
  EXPECT_LT((resumed.col(1) - finalState.col(1)).norm(), 1e-6);
  EXPECT_LT((resumed.col(2) - finalState.col(2)).norm(), 1e-4);
  checkPhysicalSurfaceExports(restartFile, config);
  const auto resumedSubsteps = readCSV(restarted / "newton-substeps.csv");
  ASSERT_EQ(resumedSubsteps.size(), 40);
  EXPECT_TRUE(std::equal(substeps.begin(), substeps.begin() + 20, resumedSubsteps.begin()));
  const std::string restartValidationCommand = shellQuote(PGO_TEST_VALIDATE_IPC_TRAJECTORY_BIN) + " " + shellQuote(restartFile) +
    " --swept --report " + shellQuote(directory / "friction-restart-validation.json") +
    " > " + shellQuote(directory / "friction-restart-validation.log") + " 2>&1";
  ASSERT_EQ(std::system(restartValidationCommand.c_str()), 0);
}

TEST_F(CudaIPCRunnerTest, CubicValidityIsRecordedAndRetainedOnRestart)
{
  auto config = configuration("cuda", 1);
  config.erase("tet-mesh");
  config["cubic-mesh"] = PGO_CUDA_TEST_CUBIC_BOX_VEG;
  const auto file = writeConfiguration(config, "cubic.json");
  ASSERT_EQ(pgo::SimulationRunner::runSimulationFromConfig(file), 0);
  const auto metadataFile = directory / "cuda/material-validity.json";
  const auto originalMetadata = readFile(metadataFile);
  const auto metadata = nlohmann::json::parse(originalMetadata);
  EXPECT_EQ(metadata.at("cubic_validity"), "gauss");
  EXPECT_EQ(metadata.at("element_type"), "cubic");
  EXPECT_EQ(metadata.at("material_max_step_enabled"), true);
  config["restart-from-u"] = true;
  fs::remove(metadataFile);
  EXPECT_NE(pgo::SimulationRunner::runSimulationFromConfig(writeConfiguration(config, "legacy-ambiguous-validity.json")), 0);
  EXPECT_FALSE(fs::exists(metadataFile));
  std::ofstream(metadataFile) << originalMetadata;
  config["cuda-cubic-validity"] = "whole-cell";
  EXPECT_NE(pgo::SimulationRunner::runSimulationFromConfig(writeConfiguration(config, "changed-validity.json")), 0);
  EXPECT_EQ(readFile(metadataFile), originalMetadata);
  config["restart-from-u"] = false;
  config["output"] = (directory / "whole-cell").string();
  ASSERT_EQ(pgo::SimulationRunner::runSimulationFromConfig(writeConfiguration(config, "whole-cell.json")), 0);
  EXPECT_EQ(nlohmann::json::parse(readFile(directory / "whole-cell/material-validity.json")).at("cubic_validity"), "whole-cell");
  ES::MXd gauss, whole;
  ASSERT_EQ(ES::readMatrix((directory / "cuda/deform0000.u").string().c_str(), gauss), 0);
  ASSERT_EQ(ES::readMatrix((directory / "whole-cell/deform0000.u").string().c_str(), whole), 0);
  EXPECT_LT((gauss - whole).cwiseAbs().maxCoeff(), 1e-8);
}

TEST_F(CudaIPCRunnerTest, UnsupportedOptionsFailExplicitly)
{
  auto config = configuration("cuda");
  config["contact-friction-coeff"] = -0.2;
  EXPECT_NE(pgo::SimulationRunner::runSimulationFromConfig(writeConfiguration(config, "friction.json")), 0);
  config["contact-friction-coeff"] = 0;
  config["elastic-material"] = "koiter-stvk";
  EXPECT_NE(pgo::SimulationRunner::runSimulationFromConfig(writeConfiguration(config, "material.json")), 0);
  config["elastic-material"] = "stable-neo";
  config["cuda-cubic-validity"] = "disabled";
  EXPECT_NE(pgo::SimulationRunner::runSimulationFromConfig(writeConfiguration(config, "invalid-validity.json")), 0);
  config["cuda-cubic-validity"] = "gauss";
  EXPECT_NE(pgo::SimulationRunner::runSimulationFromConfig(writeConfiguration(config, "tet-validity.json")), 0);
  config.erase("cuda-cubic-validity");
  config["cuda-stagnation-tolerance"] = 0.0;
  EXPECT_NE(pgo::SimulationRunner::runSimulationFromConfig(writeConfiguration(config, "invalid-stagnation-tolerance.json")), 0);
  config.erase("cuda-stagnation-tolerance");
  config["simulation-backend"] = "unavailable";
  EXPECT_NE(pgo::SimulationRunner::runSimulationFromConfig(writeConfiguration(config, "backend.json")), 0);
}
}  // namespace
