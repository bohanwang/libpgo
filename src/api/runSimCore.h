#pragma once

#include "EigenSupport.h"
#include "configFileJSON.h"
#include "deformationModelManager.h"

#include <array>
#include <string>
#include <vector>

namespace pgo {
namespace api {

struct FixedVertexAttachment {
    std::string           filename;
    std::array<double, 3> movement{};
    double                coeff = 0.0;
};

struct ExternalObject {
    std::string           filename;
    std::array<double, 3> movement{};
};

enum class VolumetricMeshInputType {
    TET,
    CUBIC,
};

struct RunSimMeshConfig {
    std::string             tetMeshFilename;
    std::string             cubicMeshFilename;
    VolumetricMeshInputType volumetricMeshType = VolumetricMeshInputType::TET;
    std::string             surfaceMeshFilename;
    double                  scale = 1.0;
};

struct RunSimSceneConfig {
    std::vector<FixedVertexAttachment> fixedVertices;
    std::vector<ExternalObject>        externalObjects;
    EigenSupport::V3d                  extAcc;
    EigenSupport::V3d                  initialVel;
    EigenSupport::V3d                  initialDisp;
};

struct RunSimSimulationConfig {
    double      timestep     = 1.0;
    int         numSimSteps  = 0;
    int         dumpInterval = 1;
    std::string simType;
};

struct RunSimContactConfig {
    std::string contactModel         = "penalty";
    double      contactStiffness     = 0.0;
    int         contactSamples       = 1;
    bool        enableSelfContact    = true;
    double      contactFrictionCoeff = 0.0;
    double      contactVelEps        = 0.0;
    double      externalIpcDhat      = 1e-3;
    double      selfIpcDhat          = 1e-3;
    double      externalIpcKappa     = 1e4;
    double      selfIpcKappa         = 1e4;
    double      ipcAlphaSafety       = 0.9;
    bool        ipcEnableFeasibleLineSearch = true;
};

struct RunSimSolverConfig {
    double                                                 solverEps     = 0.0;
    int                                                    solverMaxIter = 0;
    std::array<double, 2>                                  dampingParams{};
    SolidDeformationModel::DeformationModelElasticMaterial elasticMaterial{};
};

struct RunSimRuntimeConfig {
    std::string outputFolder;
    bool        deterministicMode = false;
    bool        restartEnabled    = false;
    bool        restartPause      = false;
};

struct RunSimConfig {
    RunSimMeshConfig       mesh;
    RunSimSceneConfig      scene;
    RunSimSimulationConfig simulation;
    RunSimContactConfig    contact;
    RunSimSolverConfig     solver;
    RunSimRuntimeConfig    runtime;
};

struct RunSimState {
    bool              hasRestart = false;
    int               frameStart = 0;
    EigenSupport::VXd u;
    EigenSupport::VXd uvel;
    EigenSupport::VXd uacc;
};

RunSimConfig parseRunSimConfig(const ConfigFileJSON& jconfig, const std::string& configFilePath);

int runSimFromConfig(const RunSimConfig& config);

}  // namespace api
}  // namespace pgo
