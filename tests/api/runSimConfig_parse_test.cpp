#include "configFileJSON.h"
#include "runSimCore.h"
#include "pgoLogging.h"

#include <gtest/gtest.h>
#include <nlohmann/json.hpp>

#include <filesystem>
#include <fstream>
#include <random>
#include <regex>
#include <stdexcept>
#include <string>
#include <vector>

namespace test {

namespace {

class ScopedTempDir {
public:
    ScopedTempDir() {
        std::random_device rd;
        path_ = std::filesystem::temp_directory_path() /
                std::filesystem::path("libpgo-runsim-config-test-" + std::to_string(rd()) + "-" +
                                      std::to_string(rd()));
        std::filesystem::create_directories(path_);
    }

    ~ScopedTempDir() {
        std::error_code ec;
        std::filesystem::remove_all(path_, ec);
    }

    const std::filesystem::path& path() const { return path_; }

private:
    std::filesystem::path path_;
};

bool writeTextFile(const std::filesystem::path& filePath, const std::string& content) {
    std::filesystem::create_directories(filePath.parent_path());
    std::ofstream os(filePath);
    if (!os) {
        return false;
    }

    os << content;
    return static_cast<bool>(os);
}

std::string readTextFile(const std::filesystem::path& filePath) {
    std::ifstream is(filePath);
    if (!is) {
        throw std::runtime_error("Failed to open text file for reading.");
    }

    return std::string((std::istreambuf_iterator<char>(is)), std::istreambuf_iterator<char>());
}

nlohmann::json makeBaseConfig() {
    return {
        {"surface-mesh", "surface.obj"},
        {"fixed-vertices", nlohmann::json::array()},
        {"external-objects", nlohmann::json::array()},
        {"g", {0.0, 0.0, 0.0}},
        {"init-vel", {0.0, 0.0, 0.0}},
        {"init-disp", {0.0, 0.0, 0.0}},
        {"scale", 1.0},
        {"timestep", 0.01},
        {"num-timestep", 1},
        {"damping-params", {0.0, 0.0}},
        {"sim-type", "static"},
        {"contact-stiffness", 0.0},
        {"contact-sample", 1},
        {"contact-friction-coeff", 0.0},
        {"contact-vel-eps", 1e-8},
        {"solver-eps", 1e-8},
        {"solver-max-iter", 2},
        {"elastic-material", "stable-neo"},
        {"dump-interval", 1},
        {"output", "result.obj"},
    };
}

pgo::api::RunSimConfig parseConfigFile(const std::filesystem::path& configFile) {
    pgo::ConfigFileJSON jconfig;
    if (!jconfig.open(configFile.string().c_str())) {
        throw std::runtime_error("Failed to open config json file.");
    }

    return pgo::api::parseRunSimConfig(jconfig, configFile.string());
}

std::string normalizedAbsolutePath(const std::filesystem::path& input) {
    return std::filesystem::absolute(input).lexically_normal().string();
}

std::filesystem::path repoRootPath() {
    return std::filesystem::absolute(std::filesystem::path(__FILE__)).parent_path().parent_path().parent_path();
}

bool hasFeasibleAlphaUpperBoundBelowOne(const std::string& logContents) {
    const std::regex alphaPattern(
        R"(IPC feasible alpha upper bound: ([+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?))");

    for (std::sregex_iterator it(logContents.begin(), logContents.end(), alphaPattern), end; it != end; ++it) {
        const double alphaUpper = std::stod((*it)[1].str());
        if (alphaUpper < 1.0 - 1e-12) {
            return true;
        }
    }

    return false;
}

std::vector<int> extractIntSequence(const std::string& logContents, const std::regex& pattern) {
    std::vector<int> sequence;

    for (std::sregex_iterator it(logContents.begin(), logContents.end(), pattern), end; it != end; ++it) {
        sequence.push_back(std::stoi((*it)[1].str()));
    }

    return sequence;
}

std::vector<std::string> extractSolverStatusSequence(const std::string& logContents) {
    const std::regex statusPattern(R"(Frame [0-9]+ solver status: ([^\n\r]+)\.)");
    std::vector<std::string> statuses;

    for (std::sregex_iterator it(logContents.begin(), logContents.end(), statusPattern), end; it != end; ++it) {
        statuses.push_back((*it)[1].str());
    }

    return statuses;
}

}  // namespace

TEST(RunSimConfigParseTest, TetMeshConfigSelectsTetInputType) {
    ScopedTempDir tempDir;

    const std::filesystem::path configDir  = tempDir.path() / "project" / "config";
    const std::filesystem::path configFile = configDir / "sim.json";

    nlohmann::json j = makeBaseConfig();
    j["tet-mesh"]   = "../assets/tet/simple.veg";

    ASSERT_TRUE(writeTextFile(configFile, j.dump(2)));

    const pgo::api::RunSimConfig config = parseConfigFile(configFile);

    EXPECT_EQ(config.mesh.volumetricMeshType, pgo::api::VolumetricMeshInputType::TET);
    EXPECT_EQ(config.mesh.tetMeshFilename, normalizedAbsolutePath(configDir / "../assets/tet/simple.veg"));
    EXPECT_TRUE(config.mesh.cubicMeshFilename.empty());
    EXPECT_EQ(config.mesh.surfaceMeshFilename, normalizedAbsolutePath(configDir / "surface.obj"));
}

TEST(RunSimConfigParseTest, CubicMeshConfigSelectsCubicInputType) {
    ScopedTempDir tempDir;

    const std::filesystem::path configDir  = tempDir.path() / "project" / "config";
    const std::filesystem::path configFile = configDir / "sim.json";

    nlohmann::json j = makeBaseConfig();
    j["cubic-mesh"] = "../assets/cubic/simple.veg";

    ASSERT_TRUE(writeTextFile(configFile, j.dump(2)));

    const pgo::api::RunSimConfig config = parseConfigFile(configFile);

    EXPECT_EQ(config.mesh.volumetricMeshType, pgo::api::VolumetricMeshInputType::CUBIC);
    EXPECT_TRUE(config.mesh.tetMeshFilename.empty());
    EXPECT_EQ(config.mesh.cubicMeshFilename, normalizedAbsolutePath(configDir / "../assets/cubic/simple.veg"));
    EXPECT_EQ(config.mesh.surfaceMeshFilename, normalizedAbsolutePath(configDir / "surface.obj"));
}

TEST(RunSimConfigParseTest, RejectsConfigWithBothTetAndCubicMeshKeys) {
    ScopedTempDir tempDir;

    const std::filesystem::path configFile = tempDir.path() / "sim.json";

    nlohmann::json j = makeBaseConfig();
    j["tet-mesh"]   = "tet.veg";
    j["cubic-mesh"] = "cubic.veg";

    ASSERT_TRUE(writeTextFile(configFile, j.dump(2)));

    EXPECT_THROW((void)parseConfigFile(configFile), std::runtime_error);
}

TEST(RunSimConfigParseTest, RejectsConfigWithMissingVolumetricMeshKey) {
    ScopedTempDir tempDir;

    const std::filesystem::path configFile = tempDir.path() / "sim.json";

    nlohmann::json j = makeBaseConfig();

    ASSERT_TRUE(writeTextFile(configFile, j.dump(2)));

    EXPECT_THROW((void)parseConfigFile(configFile), std::runtime_error);
}

TEST(RunSimConfigParseTest, RejectsVolumetricMeshTypeHintMismatch) {
    ScopedTempDir tempDir;

    const std::filesystem::path configFile = tempDir.path() / "sim.json";

    nlohmann::json j             = makeBaseConfig();
    j["cubic-mesh"]             = "cubic.veg";
    j["volumetric-mesh-type"]   = "tet";

    ASSERT_TRUE(writeTextFile(configFile, j.dump(2)));

    EXPECT_THROW((void)parseConfigFile(configFile), std::runtime_error);
}

TEST(RunSimConfigParseTest, AcceptsVolumetricMeshTypeHintWhenConsistent) {
    ScopedTempDir tempDir;

    const std::filesystem::path configFile = tempDir.path() / "sim.json";

    nlohmann::json j             = makeBaseConfig();
    j["cubic-mesh"]             = "cubic.veg";
    j["volumetric-mesh-type"]   = "cubic";

    ASSERT_TRUE(writeTextFile(configFile, j.dump(2)));

    const pgo::api::RunSimConfig config = parseConfigFile(configFile);
    EXPECT_EQ(config.mesh.volumetricMeshType, pgo::api::VolumetricMeshInputType::CUBIC);
}

TEST(RunSimConfigParseTest, DefaultsToPenaltyContactModelAndDefaultIpcParameters) {
    ScopedTempDir tempDir;

    const std::filesystem::path configFile = tempDir.path() / "sim.json";

    nlohmann::json j = makeBaseConfig();
    j["cubic-mesh"] = "cubic.veg";

    ASSERT_TRUE(writeTextFile(configFile, j.dump(2)));

    const pgo::api::RunSimConfig config = parseConfigFile(configFile);

    EXPECT_EQ(config.contact.contactModel, "penalty");
    EXPECT_DOUBLE_EQ(config.contact.externalIpcDhat, 1e-3);
    EXPECT_DOUBLE_EQ(config.contact.selfIpcDhat, 1e-3);
    EXPECT_DOUBLE_EQ(config.contact.externalIpcKappa, 1e4);
    EXPECT_DOUBLE_EQ(config.contact.selfIpcKappa, 1e4);
    EXPECT_DOUBLE_EQ(config.contact.ipcAlphaSafety, 0.9);
    EXPECT_TRUE(config.contact.ipcEnableFeasibleLineSearch);
}

TEST(RunSimConfigParseTest, ParsesIpcBarrierContactParameters) {
    ScopedTempDir tempDir;

    const std::filesystem::path configFile = tempDir.path() / "sim.json";

    nlohmann::json j        = makeBaseConfig();
    j["cubic-mesh"]         = "cubic.veg";
    j["contact-model"]      = "ipc-barrier";
    j["external-ipc-dhat"]  = 0.02;
    j["self-ipc-dhat"]      = 0.02;
    j["external-ipc-kappa"] = 2.5e5;
    j["self-ipc-kappa"]     = 2.5e5;
    j["ipc-alpha-safety"]   = 0.8;
    j["ipc-enable-feasible-line-search"] = false;

    ASSERT_TRUE(writeTextFile(configFile, j.dump(2)));

    const pgo::api::RunSimConfig config = parseConfigFile(configFile);

    EXPECT_EQ(config.contact.contactModel, "ipc-barrier");
    EXPECT_DOUBLE_EQ(config.contact.externalIpcDhat, 0.02);
    EXPECT_DOUBLE_EQ(config.contact.selfIpcDhat, 0.02);
    EXPECT_DOUBLE_EQ(config.contact.externalIpcKappa, 2.5e5);
    EXPECT_DOUBLE_EQ(config.contact.selfIpcKappa, 2.5e5);
    EXPECT_DOUBLE_EQ(config.contact.ipcAlphaSafety, 0.8);
    EXPECT_FALSE(config.contact.ipcEnableFeasibleLineSearch);
}

TEST(RunSimConfigParseTest, ParsesDecoupledExternalAndSelfIpcBarrierParameters) {
    ScopedTempDir tempDir;

    const std::filesystem::path configFile = tempDir.path() / "sim.json";

    nlohmann::json j             = makeBaseConfig();
    j["cubic-mesh"]              = "cubic.veg";
    j["contact-model"]           = "ipc-barrier";
    j["external-ipc-dhat"]       = 0.03;
    j["self-ipc-dhat"]           = 0.04;
    j["external-ipc-kappa"]      = 3.5e5;
    j["self-ipc-kappa"]          = 4.5e5;

    ASSERT_TRUE(writeTextFile(configFile, j.dump(2)));

    const pgo::api::RunSimConfig config = parseConfigFile(configFile);

    EXPECT_DOUBLE_EQ(config.contact.externalIpcDhat, 0.03);
    EXPECT_DOUBLE_EQ(config.contact.selfIpcDhat, 0.04);
    EXPECT_DOUBLE_EQ(config.contact.externalIpcKappa, 3.5e5);
    EXPECT_DOUBLE_EQ(config.contact.selfIpcKappa, 4.5e5);
}

TEST(RunSimConfigParseTest, RejectsLegacySharedIpcBarrierParameters) {
    ScopedTempDir tempDir;

    const std::filesystem::path configFile = tempDir.path() / "sim.json";

    nlohmann::json j      = makeBaseConfig();
    j["cubic-mesh"]       = "cubic.veg";
    j["contact-model"]    = "ipc-barrier";
    j["ipc-dhat"]         = 0.02;
    j["ipc-kappa"]        = 2.5e5;

    ASSERT_TRUE(writeTextFile(configFile, j.dump(2)));

    EXPECT_THROW((void)parseConfigFile(configFile), std::runtime_error);
}

TEST(RunSimConfigParseTest, IpcBarrierFeasibleLineSearchDefaultsToEnabled) {
    ScopedTempDir tempDir;

    const std::filesystem::path configFile = tempDir.path() / "sim.json";

    nlohmann::json j      = makeBaseConfig();
    j["cubic-mesh"]       = "cubic.veg";
    j["contact-model"]    = "ipc-barrier";
    j["external-ipc-dhat"]  = 0.02;
    j["self-ipc-dhat"]      = 0.02;
    j["external-ipc-kappa"] = 2.5e5;
    j["self-ipc-kappa"]     = 2.5e5;
    j["ipc-alpha-safety"] = 0.8;

    ASSERT_TRUE(writeTextFile(configFile, j.dump(2)));

    const pgo::api::RunSimConfig config = parseConfigFile(configFile);

    EXPECT_TRUE(config.contact.ipcEnableFeasibleLineSearch);
}

TEST(RunSimConfigParseTest, RejectsUnsupportedContactModel) {
    ScopedTempDir tempDir;

    const std::filesystem::path configFile = tempDir.path() / "sim.json";

    nlohmann::json j        = makeBaseConfig();
    j["cubic-mesh"]         = "cubic.veg";
    j["contact-model"]      = "unsupported";

    ASSERT_TRUE(writeTextFile(configFile, j.dump(2)));

    EXPECT_THROW((void)parseConfigFile(configFile), std::runtime_error);
}

TEST(RunSimConfigParseTest, RejectsNonPositiveExternalIpcDhat) {
    ScopedTempDir tempDir;

    const std::filesystem::path configFile = tempDir.path() / "sim.json";

    nlohmann::json j        = makeBaseConfig();
    j["cubic-mesh"]         = "cubic.veg";
    j["contact-model"]      = "ipc-barrier";
    j["external-ipc-dhat"]  = 0.0;

    ASSERT_TRUE(writeTextFile(configFile, j.dump(2)));

    EXPECT_THROW((void)parseConfigFile(configFile), std::runtime_error);
}

TEST(RunSimConfigParseTest, RejectsNonPositiveSelfIpcKappa) {
    ScopedTempDir tempDir;

    const std::filesystem::path configFile = tempDir.path() / "sim.json";

    nlohmann::json j        = makeBaseConfig();
    j["cubic-mesh"]         = "cubic.veg";
    j["contact-model"]      = "ipc-barrier";
    j["self-ipc-kappa"]     = 0.0;

    ASSERT_TRUE(writeTextFile(configFile, j.dump(2)));

    EXPECT_THROW((void)parseConfigFile(configFile), std::runtime_error);
}

TEST(RunSimConfigParseTest, RejectsOutOfRangeIpcAlphaSafety) {
    ScopedTempDir tempDir;

    const std::filesystem::path configFile = tempDir.path() / "sim-a.json";
    const std::filesystem::path configFile2 = tempDir.path() / "sim-b.json";

    nlohmann::json j        = makeBaseConfig();
    j["cubic-mesh"]         = "cubic.veg";
    j["contact-model"]      = "ipc-barrier";
    j["ipc-alpha-safety"]   = 0.0;

    ASSERT_TRUE(writeTextFile(configFile, j.dump(2)));

    nlohmann::json j2      = j;
    j2["ipc-alpha-safety"] = 1.1;
    ASSERT_TRUE(writeTextFile(configFile2, j2.dump(2)));

    EXPECT_THROW((void)parseConfigFile(configFile), std::runtime_error);
    EXPECT_THROW((void)parseConfigFile(configFile2), std::runtime_error);
}

TEST(RunSimConfigParseTest, RejectsFrictionWithIpcBarrierContactModel) {
    ScopedTempDir tempDir;

    const std::filesystem::path configFile = tempDir.path() / "sim.json";

    nlohmann::json j                 = makeBaseConfig();
    j["cubic-mesh"]                  = "cubic.veg";
    j["contact-model"]               = "ipc-barrier";
    j["contact-friction-coeff"]      = 0.3;

    ASSERT_TRUE(writeTextFile(configFile, j.dump(2)));

    EXPECT_THROW((void)parseConfigFile(configFile), std::runtime_error);
}

TEST(RunSimConfigParseTest, RunSimFromConfigCubicDynamicSmokeTest) {
    ScopedTempDir tempDir;
    const std::filesystem::path repoRoot = repoRootPath();
    const std::filesystem::path cubicMeshFile = repoRoot / "examples" / "cubic-box" / "cubic-box.veg";
    const std::filesystem::path surfaceMeshFile = repoRoot / "examples" / "cubic-box" / "cubic-box.obj";
    const std::filesystem::path outputDir = tempDir.path() / "output";

    ASSERT_TRUE(std::filesystem::exists(cubicMeshFile));
    ASSERT_TRUE(std::filesystem::exists(surfaceMeshFile));

    pgo::Logging::init();

    pgo::api::RunSimConfig config;
    config.mesh.cubicMeshFilename = cubicMeshFile.string();
    config.mesh.surfaceMeshFilename = surfaceMeshFile.string();
    config.mesh.volumetricMeshType = pgo::api::VolumetricMeshInputType::CUBIC;
    config.scene.extAcc = pgo::EigenSupport::V3d::Zero();
    config.scene.initialVel = pgo::EigenSupport::V3d::Zero();
    config.scene.initialDisp = pgo::EigenSupport::V3d::Zero();
    config.mesh.scale = 1.0;
    config.simulation.timestep = 0.01;
    config.contact.contactStiffness = 0.0;
    config.contact.contactSamples = 1;
    config.contact.contactFrictionCoeff = 0.0;
    config.contact.contactVelEps = 1e-8;
    config.solver.solverEps = 1e-8;
    config.solver.solverMaxIter = 4;
    config.solver.dampingParams = {0.0, 0.0};
    config.solver.elasticMaterial = pgo::SolidDeformationModel::DeformationModelElasticMaterial::STABLE_NEO;
    config.simulation.numSimSteps = 2;
    config.simulation.dumpInterval = 1;
    config.simulation.simType = "dynamic";
    config.runtime.outputFolder = outputDir.string();
    config.runtime.deterministicMode = true;

    EXPECT_EQ(pgo::api::runSimFromConfig(config), 0);
    EXPECT_TRUE(std::filesystem::exists(outputDir / "deform0001.u"));
    EXPECT_TRUE(std::filesystem::exists(outputDir / "ret0001.obj"));

    std::error_code ec;
    std::filesystem::remove("fv.obj", ec);
}

TEST(RunSimConfigParseTest, RunSimFromConfigCubicDynamicIpcSmokeTest) {
    ScopedTempDir tempDir;
    const std::filesystem::path repoRoot = repoRootPath();
    const std::filesystem::path cubicMeshFile = repoRoot / "examples" / "cubic-box" / "cubic-box.veg";
    const std::filesystem::path surfaceMeshFile = repoRoot / "examples" / "cubic-box" / "cubic-box.obj";
    const std::filesystem::path externalObjectFile = repoRoot / "examples" / "bottom.obj";
    const std::filesystem::path outputDir = tempDir.path() / "output";

    ASSERT_TRUE(std::filesystem::exists(cubicMeshFile));
    ASSERT_TRUE(std::filesystem::exists(surfaceMeshFile));
    ASSERT_TRUE(std::filesystem::exists(externalObjectFile));

    pgo::Logging::init();

    pgo::api::RunSimConfig config;
    config.mesh.cubicMeshFilename = cubicMeshFile.string();
    config.mesh.surfaceMeshFilename = surfaceMeshFile.string();
    config.mesh.volumetricMeshType = pgo::api::VolumetricMeshInputType::CUBIC;
    config.scene.externalObjects.push_back({externalObjectFile.string(), {0.0, 0.0, 0.0}});
    config.scene.extAcc = pgo::EigenSupport::V3d::Zero();
    config.scene.initialVel = pgo::EigenSupport::V3d::Zero();
    config.scene.initialDisp = pgo::EigenSupport::V3d::Zero();
    config.mesh.scale = 1.0;
    config.simulation.timestep = 0.01;
    config.contact.contactStiffness = 0.0;
    config.contact.contactSamples = 1;
    config.contact.enableSelfContact = false;
    config.contact.contactFrictionCoeff = 0.0;
    config.contact.contactVelEps = 1e-8;
    config.contact.contactModel = "ipc-barrier";
    config.contact.externalIpcDhat = 0.01;
    config.contact.selfIpcDhat = 0.01;
    config.contact.externalIpcKappa = 1e5;
    config.contact.selfIpcKappa = 1e5;
    config.contact.ipcAlphaSafety = 0.9;
    config.solver.solverEps = 1e-8;
    config.solver.solverMaxIter = 4;
    config.solver.dampingParams = {0.0, 0.0};
    config.solver.elasticMaterial = pgo::SolidDeformationModel::DeformationModelElasticMaterial::STABLE_NEO;
    config.simulation.numSimSteps = 2;
    config.simulation.dumpInterval = 1;
    config.simulation.simType = "dynamic";
    config.runtime.outputFolder = outputDir.string();
    config.runtime.deterministicMode = true;

    EXPECT_EQ(pgo::api::runSimFromConfig(config), 0);
    EXPECT_TRUE(std::filesystem::exists(outputDir / "deform0001.u"));
    EXPECT_TRUE(std::filesystem::exists(outputDir / "ret0001.obj"));

    std::error_code ec;
    std::filesystem::remove("fv.obj", ec);
}

TEST(RunSimConfigParseTest, RunSimFromConfigCubicDynamicIpcNearContactActivatesExternalBarrier) {
    ScopedTempDir tempDir;
    const std::filesystem::path repoRoot = repoRootPath();
    const std::filesystem::path cubicMeshFile = repoRoot / "examples" / "cubic-box" / "cubic-box.veg";
    const std::filesystem::path surfaceMeshFile = repoRoot / "examples" / "cubic-box" / "cubic-box.obj";
    const std::filesystem::path externalObjectFile = repoRoot / "examples" / "bottom.obj";
    const std::filesystem::path outputDir = tempDir.path() / "output";
    const std::filesystem::path logFile = tempDir.path() / "runsim-ipc.log";
    const std::string logFilename = logFile.string();

    ASSERT_TRUE(std::filesystem::exists(cubicMeshFile));
    ASSERT_TRUE(std::filesystem::exists(surfaceMeshFile));
    ASSERT_TRUE(std::filesystem::exists(externalObjectFile));

    pgo::Logging::init(logFilename.c_str());

    pgo::api::RunSimConfig config;
    config.mesh.cubicMeshFilename = cubicMeshFile.string();
    config.mesh.surfaceMeshFilename = surfaceMeshFile.string();
    config.mesh.volumetricMeshType = pgo::api::VolumetricMeshInputType::CUBIC;
    config.scene.externalObjects.push_back({externalObjectFile.string(), {0.0, 0.0, 0.0}});
    config.scene.extAcc = pgo::EigenSupport::V3d::Zero();
    config.scene.initialVel = pgo::EigenSupport::V3d::Zero();
    config.scene.initialDisp = pgo::EigenSupport::V3d(0.0, -0.2, 0.0);
    config.mesh.scale = 1.0;
    config.simulation.timestep = 0.01;
    config.contact.contactStiffness = 0.0;
    config.contact.contactSamples = 1;
    config.contact.enableSelfContact = false;
    config.contact.contactFrictionCoeff = 0.0;
    config.contact.contactVelEps = 1e-8;
    config.contact.contactModel = "ipc-barrier";
    config.contact.externalIpcDhat = 0.01;
    config.contact.selfIpcDhat = 0.01;
    config.contact.externalIpcKappa = 1e5;
    config.contact.selfIpcKappa = 1e5;
    config.contact.ipcAlphaSafety = 0.9;
    config.solver.solverEps = 1e-8;
    config.solver.solverMaxIter = 4;
    config.solver.dampingParams = {0.0, 0.0};
    config.solver.elasticMaterial = pgo::SolidDeformationModel::DeformationModelElasticMaterial::STABLE_NEO;
    config.simulation.numSimSteps = 2;
    config.simulation.dumpInterval = 1;
    config.simulation.simType = "dynamic";
    config.runtime.outputFolder = outputDir.string();
    config.runtime.deterministicMode = true;

    EXPECT_EQ(pgo::api::runSimFromConfig(config), 0);
    EXPECT_TRUE(std::filesystem::exists(outputDir / "deform0001.u"));
    EXPECT_TRUE(std::filesystem::exists(outputDir / "ret0001.obj"));
    ASSERT_TRUE(std::filesystem::exists(logFile));

    const std::string logContents = readTextFile(logFile);
    const std::regex activeSamplesPattern(R"(# external active samples: ([1-9][0-9]*))");

    EXPECT_TRUE(std::regex_search(logContents, activeSamplesPattern));
    EXPECT_NE(logContents.find("IPC feasible alpha callback active."), std::string::npos);
    EXPECT_EQ(logContents.find("IPC barrier energy not yet implemented for external contact."), std::string::npos);

    std::error_code ec;
    std::filesystem::remove("fv.obj", ec);
}

TEST(RunSimConfigParseTest, RunSimFromConfigCubicDynamicIpcSelfBarrierSmokeTest) {
    ScopedTempDir tempDir;
    const std::filesystem::path repoRoot = repoRootPath();
    const std::filesystem::path exampleConfigFile =
        repoRoot / "examples" / "pulled-cubic-box-self-ipc" / "pulled-cubic-box-self-ipc.json";
    const std::filesystem::path outputDir = tempDir.path() / "output";
    const std::filesystem::path logFile = tempDir.path() / "runsim-self-ipc.log";
    const std::string logFilename = logFile.string();

    ASSERT_TRUE(std::filesystem::exists(exampleConfigFile));

    pgo::Logging::init(logFilename.c_str());

    pgo::api::RunSimConfig config = parseConfigFile(exampleConfigFile);
    config.runtime.outputFolder = outputDir.string();
    config.runtime.deterministicMode = true;
    config.simulation.numSimSteps = 5;
    config.simulation.dumpInterval = 1;
    config.contact.externalIpcDhat = 0.1;
    config.contact.selfIpcDhat = 0.1;
    config.contact.externalIpcKappa = 10.0;
    config.contact.selfIpcKappa = 10.0;
    ASSERT_GE(config.scene.fixedVertices.size(), 2u);
    config.scene.fixedVertices[1].movement = {0.0, -1.2, 0.0};

    EXPECT_EQ(pgo::api::runSimFromConfig(config), 0);
    EXPECT_TRUE(std::filesystem::exists(outputDir / "deform0001.u"));
    EXPECT_TRUE(std::filesystem::exists(outputDir / "ret0001.obj"));
    EXPECT_TRUE(std::filesystem::exists(outputDir / "deform0004.u"));
    EXPECT_TRUE(std::filesystem::exists(outputDir / "ret0004.obj"));
    ASSERT_TRUE(std::filesystem::exists(logFile));

    const std::string logContents = readTextFile(logFile);
    const std::regex activePairsPattern(R"(# self active pairs: ([1-9][0-9]*))");

    EXPECT_TRUE(std::regex_search(logContents, activePairsPattern));
    EXPECT_TRUE(hasFeasibleAlphaUpperBoundBelowOne(logContents));
    EXPECT_EQ(logContents.find("IPC barrier energy not yet implemented for self contact."), std::string::npos);

    std::error_code ec;
    std::filesystem::remove("fv.obj", ec);
}

TEST(RunSimConfigParseTest, RunSimFromConfigCubicDynamicIpcMergedBarrierSmokeTest) {
    ScopedTempDir tempDir;
    const std::filesystem::path repoRoot = repoRootPath();
    const std::filesystem::path exampleConfigFile =
        repoRoot / "examples" / "pulled-cubic-box-self-ipc" / "pulled-cubic-box-self-ipc.json";
    const std::filesystem::path externalObjectFile = repoRoot / "examples" / "bottom.obj";
    const std::filesystem::path outputDir = tempDir.path() / "output";
    const std::filesystem::path logFile = tempDir.path() / "runsim-merged-ipc.log";
    const std::string logFilename = logFile.string();

    ASSERT_TRUE(std::filesystem::exists(exampleConfigFile));
    ASSERT_TRUE(std::filesystem::exists(externalObjectFile));

    pgo::Logging::init(logFilename.c_str());

    pgo::api::RunSimConfig config = parseConfigFile(exampleConfigFile);
    config.scene.externalObjects.push_back({externalObjectFile.string(), {0.0, 0.0, 0.0}});
    config.scene.initialDisp = pgo::EigenSupport::V3d(0.0, -0.2, 0.0);
    config.runtime.outputFolder = outputDir.string();
    config.runtime.deterministicMode = true;
    config.simulation.numSimSteps = 8;
    config.simulation.dumpInterval = 1;
    config.contact.externalIpcDhat = 0.1;
    config.contact.selfIpcDhat = 0.1;
    config.contact.externalIpcKappa = 10.0;
    config.contact.selfIpcKappa = 10.0;
    config.contact.ipcAlphaSafety = 0.9;
    config.contact.ipcEnableFeasibleLineSearch = true;
    ASSERT_GE(config.scene.fixedVertices.size(), 2u);
    config.scene.fixedVertices[1].movement = {0.0, -1.2, 0.0};

    EXPECT_EQ(pgo::api::runSimFromConfig(config), 0);
    EXPECT_TRUE(std::filesystem::exists(outputDir / "deform0001.u"));
    EXPECT_TRUE(std::filesystem::exists(outputDir / "ret0001.obj"));
    EXPECT_TRUE(std::filesystem::exists(outputDir / "deform0007.u"));
    EXPECT_TRUE(std::filesystem::exists(outputDir / "ret0007.obj"));
    ASSERT_TRUE(std::filesystem::exists(logFile));

    const std::string logContents = readTextFile(logFile);
    const std::regex activeSamplesPattern(R"(# external active samples: ([1-9][0-9]*))");
    const std::regex activePairsPattern(R"(# self active pairs: ([1-9][0-9]*))");

    EXPECT_TRUE(std::regex_search(logContents, activeSamplesPattern));
    EXPECT_TRUE(std::regex_search(logContents, activePairsPattern));
    EXPECT_NE(logContents.find("IPC feasible alpha merged callback active."), std::string::npos);
    EXPECT_TRUE(hasFeasibleAlphaUpperBoundBelowOne(logContents));
    EXPECT_EQ(logContents.find("terminated with solver status"), std::string::npos);
    EXPECT_EQ(logContents.find("nan"), std::string::npos);
    EXPECT_EQ(logContents.find("NaN"), std::string::npos);
    EXPECT_EQ(logContents.find("IPC barrier energy not yet implemented for external contact."), std::string::npos);
    EXPECT_EQ(logContents.find("IPC barrier energy not yet implemented for self contact."), std::string::npos);

    std::error_code ec;
    std::filesystem::remove("fv.obj", ec);
}

TEST(RunSimConfigParseTest, RunSimFromConfigCubicDynamicIpcMergedDeterministicSmokeTest) {
    ScopedTempDir tempDir;
    const std::filesystem::path repoRoot = repoRootPath();
    const std::filesystem::path exampleConfigFile =
        repoRoot / "examples" / "pulled-cubic-box-self-ipc" / "pulled-cubic-box-self-ipc.json";
    const std::filesystem::path externalObjectFile = repoRoot / "examples" / "bottom.obj";
    const std::filesystem::path outputDirA = tempDir.path() / "output-a";
    const std::filesystem::path outputDirB = tempDir.path() / "output-b";
    const std::filesystem::path logFileA = tempDir.path() / "runsim-merged-ipc-a.log";
    const std::filesystem::path logFileB = tempDir.path() / "runsim-merged-ipc-b.log";

    ASSERT_TRUE(std::filesystem::exists(exampleConfigFile));
    ASSERT_TRUE(std::filesystem::exists(externalObjectFile));

    pgo::api::RunSimConfig config = parseConfigFile(exampleConfigFile);
    config.scene.externalObjects.push_back({externalObjectFile.string(), {0.0, 0.0, 0.0}});
    config.scene.initialDisp = pgo::EigenSupport::V3d(0.0, -0.2, 0.0);
    config.runtime.deterministicMode = true;
    config.simulation.numSimSteps = 8;
    config.simulation.dumpInterval = 1;
    config.contact.externalIpcDhat = 0.1;
    config.contact.selfIpcDhat = 0.1;
    config.contact.externalIpcKappa = 10.0;
    config.contact.selfIpcKappa = 10.0;
    config.contact.ipcAlphaSafety = 0.9;
    config.contact.ipcEnableFeasibleLineSearch = true;
    ASSERT_GE(config.scene.fixedVertices.size(), 2u);
    config.scene.fixedVertices[1].movement = {0.0, -1.2, 0.0};

    pgo::Logging::init(logFileA.string().c_str());
    config.runtime.outputFolder = outputDirA.string();
    EXPECT_EQ(pgo::api::runSimFromConfig(config), 0);
    ASSERT_TRUE(std::filesystem::exists(logFileA));

    pgo::Logging::init(logFileB.string().c_str());
    config.runtime.outputFolder = outputDirB.string();
    EXPECT_EQ(pgo::api::runSimFromConfig(config), 0);
    ASSERT_TRUE(std::filesystem::exists(logFileB));

    const std::string logA = readTextFile(logFileA);
    const std::string logB = readTextFile(logFileB);
    const std::regex externalSeqPattern(R"(# external active samples: ([0-9]+))");
    const std::regex selfSeqPattern(R"(# self active pairs: ([0-9]+))");

    const std::vector<int> externalSeqA = extractIntSequence(logA, externalSeqPattern);
    const std::vector<int> externalSeqB = extractIntSequence(logB, externalSeqPattern);
    const std::vector<int> selfSeqA = extractIntSequence(logA, selfSeqPattern);
    const std::vector<int> selfSeqB = extractIntSequence(logB, selfSeqPattern);
    const std::vector<std::string> statusSeqA = extractSolverStatusSequence(logA);
    const std::vector<std::string> statusSeqB = extractSolverStatusSequence(logB);

    ASSERT_FALSE(externalSeqA.empty());
    ASSERT_FALSE(selfSeqA.empty());
    ASSERT_FALSE(statusSeqA.empty());
    EXPECT_EQ(externalSeqA, externalSeqB);
    EXPECT_EQ(selfSeqA, selfSeqB);
    EXPECT_EQ(statusSeqA, statusSeqB);

    std::error_code ec;
    std::filesystem::remove("fv.obj", ec);
}

}  // namespace test
