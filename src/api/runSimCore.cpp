#include "runSimCore.h"

#include "fileService.h"
#include "basicIO.h"
#include "generateTetMeshMatrix.h"
#include "cubicMesh.h"
#include "tetMesh.h"
#include "triMeshGeo.h"
#include "geometryQuery.h"
#include "simulationMesh.h"
#include "deformationModelManager.h"
#include "tetMeshDeformationModel.h"
#include "plasticModel3DDeformationGradient.h"
#include "deformationModelAssembler.h"
#include "deformationModelEnergy.h"
#include "multiVertexPullingSoftConstraints.h"
#include "implicitBackwardEulerTimeIntegrator.h"
#include "TRBDF2TimeIntegrator.h"
#include "generateMassMatrix.h"
#include "barycentricCoordinates.h"
#include "triangleMeshExternalContactHandler.h"
#include "pointPenetrationBarrierEnergy.h"
#include "pointPenetrationEnergy.h"
#include "triangleMeshSelfContactHandler.h"
#include "pointTrianglePairBarrierEnergy.h"
#include "pointTrianglePairCouplingEnergyWithCollision.h"
#include "linearPotentialEnergy.h"
#include "NewtonRaphsonSolver.h"
#include "pgoLogging.h"

#include <tbb/global_control.h>

#include <algorithm>
#include <fmt/format.h>
#include <iostream>
#include <iterator>
#include <memory>
#include <stdexcept>
#include <thread>
#include <vector>

namespace pgo {
namespace api {

namespace {

namespace ES = EigenSupport;

bool isDynamicContactEnabled(const RunSimContactConfig& contact) {
    return (contact.contactModel == "penalty" && contact.contactStiffness > 0.0) ||
           (contact.contactModel == "ipc-barrier" &&
            (contact.externalIpcKappa > 0.0 || contact.selfIpcKappa > 0.0));
}

Contact::PosFunction makeCurrentPositionFunction(const ES::VXd& restPosition) {
    return [&restPosition](const ES::V3d& u, ES::V3d& p, int dofStart) { p = u + restPosition.segment<3>(dofStart); };
}

Contact::PosFunction makeLastPositionFunction(const ES::VXd& restPosition, const ES::VXd& u) {
    return [&restPosition, &u](const ES::V3d& x, ES::V3d& p, int dofStart) {
        (void)x;
        p = u.segment<3>(dofStart) + restPosition.segment<3>(dofStart);
    };
}

std::string describeSolverStatus(Simulation::TimeIntegratorSolverOption option, int status) {
    if (option == Simulation::TimeIntegratorSolverOption::SO_NEWTON) {
        return fmt::format("{} ({})", status, NonlinearOptimization::NewtonRaphsonSolver::statusToString(status));
    }

    return fmt::format("{}", status);
}

}  // namespace

RunSimConfig parseRunSimConfig(const ConfigFileJSON& jconfig, const std::string& configFilePath) {
    using namespace EigenSupport;
    namespace ES = EigenSupport;

    if (configFilePath.empty()) {
        throw std::runtime_error("Config file path is empty.");
    }
    const FileService::ConfigPathResolver configPathResolver(configFilePath);

    RunSimConfig config;

    const bool hasTetMesh   = jconfig.exist("tet-mesh");
    const bool hasCubicMesh = jconfig.exist("cubic-mesh");

    if (hasTetMesh && hasCubicMesh) {
        throw std::runtime_error("Config cannot contain both 'tet-mesh' and 'cubic-mesh'.");
    }

    if (!hasTetMesh && !hasCubicMesh) {
        throw std::runtime_error("Config must contain either 'tet-mesh' or 'cubic-mesh'.");
    }

    if (hasTetMesh) {
        config.mesh.tetMeshFilename = configPathResolver.resolve(jconfig.getString("tet-mesh", 1));
        if (config.mesh.tetMeshFilename.empty()) {
            throw std::runtime_error("'tet-mesh' resolves to an empty path.");
        }
        config.mesh.volumetricMeshType = VolumetricMeshInputType::TET;
    } else {
        config.mesh.cubicMeshFilename = configPathResolver.resolve(jconfig.getString("cubic-mesh", 1));
        if (config.mesh.cubicMeshFilename.empty()) {
            throw std::runtime_error("'cubic-mesh' resolves to an empty path.");
        }
        config.mesh.volumetricMeshType = VolumetricMeshInputType::CUBIC;
    }

    if (jconfig.exist("volumetric-mesh-type")) {
        const std::string hintedType = jconfig.getString("volumetric-mesh-type", 1);

        VolumetricMeshInputType hintedMeshType;
        if (hintedType == "tet") {
            hintedMeshType = VolumetricMeshInputType::TET;
        } else if (hintedType == "cubic") {
            hintedMeshType = VolumetricMeshInputType::CUBIC;
        } else {
            throw std::runtime_error(fmt::format("Unsupported volumetric-mesh-type: {}", hintedType));
        }

        if (hintedMeshType != config.mesh.volumetricMeshType) {
            throw std::runtime_error(
                "volumetric-mesh-type does not match the volumetric mesh filename key.");
        }
    }

    config.mesh.surfaceMeshFilename = configPathResolver.resolve(jconfig.getString("surface-mesh", 1));
    auto readVec3              = [&](const char* key) {
        auto values = jconfig.getValue<std::array<double, 3>>(key, 1);
        return ES::V3d(values[0], values[1], values[2]);
    };

    config.scene.extAcc               = readVec3("g");
    config.scene.initialVel           = readVec3("init-vel");
    config.scene.initialDisp          = readVec3("init-disp");
    config.mesh.scale                 = jconfig.getDouble("scale", 1);
    config.simulation.timestep        = jconfig.getDouble("timestep", 1);
    config.contact.contactStiffness   = jconfig.getDouble("contact-stiffness", 1);
    config.contact.contactSamples     = jconfig.getInt("contact-sample", 1);
    config.contact.enableSelfContact  = jconfig.handle().value("enable-self-contact", true);
    config.contact.contactFrictionCoeff = jconfig.getDouble("contact-friction-coeff", 1);
    config.contact.contactVelEps      = jconfig.getDouble("contact-vel-eps", 1);
    config.contact.contactModel       = jconfig.handle().value("contact-model", std::string("penalty"));
    if (config.contact.contactModel != "penalty" && config.contact.contactModel != "ipc-barrier") {
        throw std::runtime_error(fmt::format(
            "Unsupported contact-model: '{}'. Expected 'penalty' or 'ipc-barrier'.", config.contact.contactModel));
    }
    if (config.contact.contactModel == "ipc-barrier") {
        if (jconfig.exist("ipc-dhat") || jconfig.exist("ipc-kappa")) {
            throw std::runtime_error(
                "ipc-dhat/ipc-kappa are no longer supported. "
                "Use external-ipc-dhat, self-ipc-dhat, external-ipc-kappa, and self-ipc-kappa.");
        }

        config.contact.externalIpcDhat  = jconfig.handle().value("external-ipc-dhat", 1e-3);
        config.contact.selfIpcDhat      = jconfig.handle().value("self-ipc-dhat", 1e-3);
        config.contact.externalIpcKappa = jconfig.handle().value("external-ipc-kappa", 1e4);
        config.contact.selfIpcKappa     = jconfig.handle().value("self-ipc-kappa", 1e4);
        config.contact.ipcAlphaSafety = jconfig.handle().value("ipc-alpha-safety", 0.9);
        config.contact.ipcEnableFeasibleLineSearch =
            jconfig.handle().value("ipc-enable-feasible-line-search", true);

        if (config.contact.externalIpcDhat <= 0.0) {
            throw std::runtime_error("external-ipc-dhat must be positive.");
        }
        if (config.contact.selfIpcDhat <= 0.0) {
            throw std::runtime_error("self-ipc-dhat must be positive.");
        }
        if (config.contact.externalIpcKappa <= 0.0) {
            throw std::runtime_error("external-ipc-kappa must be positive.");
        }
        if (config.contact.selfIpcKappa <= 0.0) {
            throw std::runtime_error("self-ipc-kappa must be positive.");
        }
        if (config.contact.ipcAlphaSafety <= 0.0 || config.contact.ipcAlphaSafety > 1.0) {
            throw std::runtime_error("ipc-alpha-safety must be in (0, 1].");
        }
        if (config.contact.contactFrictionCoeff > 0.0) {
            throw std::runtime_error(
                "contact-model 'ipc-barrier' does not support friction yet. "
                "Set contact-friction-coeff to 0.");
        }
    }
    config.solver.solverEps            = jconfig.getDouble("solver-eps", 1);
    config.solver.solverMaxIter        = jconfig.getInt("solver-max-iter", 1);
    config.solver.dampingParams        = jconfig.getValue<std::array<double, 2>>("damping-params", 1);

    std::string material = jconfig.getString("elastic-material");
    if (material == "stable-neo") {
        config.solver.elasticMaterial = SolidDeformationModel::DeformationModelElasticMaterial::STABLE_NEO;
    } else if (material == "stvk-vol") {
        config.solver.elasticMaterial = SolidDeformationModel::DeformationModelElasticMaterial::STVK_VOL;
    } else {
        throw std::runtime_error(fmt::format("Unsupported elastic material: {}", material));
    }

    config.simulation.numSimSteps  = jconfig.getInt("num-timestep", 1);
    config.simulation.dumpInterval = jconfig.getInt("dump-interval", 1);
    config.simulation.simType      = jconfig.getString("sim-type");
    config.runtime.outputFolder    = configPathResolver.resolve(jconfig.getString("output", 1));
    config.runtime.deterministicMode = jconfig.handle().value("deterministic", false);

    if (jconfig.exist("fixed-vertices")) {
        for (const auto& fv : jconfig.handle()["fixed-vertices"]) {
            FixedVertexAttachment attachment;
            attachment.filename = configPathResolver.resolve(fv["filename"].get<std::string>());
            attachment.movement = fv["movement"].get<std::array<double, 3>>();
            attachment.coeff    = fv["coeff"].get<double>();
            config.scene.fixedVertices.push_back(attachment);
        }
    }

    if (jconfig.exist("external-objects")) {
        for (const auto& obj : jconfig.handle()["external-objects"]) {
            ExternalObject ext;
            ext.filename = configPathResolver.resolve(obj["filename"].get<std::string>());
            ext.movement = obj["movement"].get<std::array<double, 3>>();
            config.scene.externalObjects.push_back(ext);
        }
    }

    return config;
}

int runSimFromConfig(const RunSimConfig& config) {
    using namespace EigenSupport;
    namespace ES = EigenSupport;

    const int defaultParallelism =
        std::min(64, static_cast<int>(std::max(1u, std::thread::hardware_concurrency())));
    const auto& mesh       = config.mesh;
    const auto& scene      = config.scene;
    const auto& simulation = config.simulation;
    const auto& contact    = config.contact;
    const auto& solver     = config.solver;
    const auto& runtime    = config.runtime;
    const double externalIpcDhat  = contact.externalIpcDhat;
    const double selfIpcDhat      = contact.selfIpcDhat;
    const double externalIpcKappa = contact.externalIpcKappa;
    const double selfIpcKappa     = contact.selfIpcKappa;

    const int effectiveParallelism = runtime.deterministicMode ? 1 : defaultParallelism;

    tbb::global_control threadLimit(tbb::global_control::max_allowed_parallelism, effectiveParallelism);
    tbb::global_control stackSizeLimit(tbb::global_control::thread_stack_size, 16 * 1024 * 1024);

    if (runtime.deterministicMode) {
        SPDLOG_LOGGER_INFO(Logging::lgr(), "Deterministic mode enabled; forcing single-threaded execution.");
    }

    std::unique_ptr<VolumetricMeshes::VolumetricMesh>        volumetricMesh;
    std::shared_ptr<SolidDeformationModel::SimulationMesh>   simMesh;
    if (mesh.volumetricMeshType == VolumetricMeshInputType::TET) {
        auto tetMesh = std::make_unique<VolumetricMeshes::TetMesh>(mesh.tetMeshFilename.c_str());
        for (int vi = 0; vi < tetMesh->getNumVertices(); vi++) {
            tetMesh->setVertex(vi, tetMesh->getVertex(vi) * mesh.scale);
        }

        simMesh        = SolidDeformationModel::SimulationMesh::createFromTetMesh(*tetMesh);
        volumetricMesh = std::move(tetMesh);
    } else {
        auto cubicMesh = std::make_unique<VolumetricMeshes::CubicMesh>(mesh.cubicMeshFilename.c_str());
        for (int vi = 0; vi < cubicMesh->getNumVertices(); vi++) {
            cubicMesh->setVertex(vi, cubicMesh->getVertex(vi) * mesh.scale);
        }

        simMesh        = SolidDeformationModel::SimulationMesh::createFromCubicMesh(*cubicMesh);
        volumetricMesh = std::move(cubicMesh);
    }

    Mesh::TriMeshGeo surfaceMesh;
    if (surfaceMesh.load(mesh.surfaceMeshFilename) != true)
        return 1;

    for (int vi = 0; vi < surfaceMesh.numVertices(); vi++) {
        surfaceMesh.pos(vi) *= mesh.scale;
    }

    int     surfn  = surfaceMesh.numVertices();
    int     surfn3 = surfn * 3;
    ES::VXd surfaceRestPositions(surfn3);
    for (int vi = 0; vi < surfn; vi++) {
        surfaceRestPositions.segment<3>(vi * 3) = surfaceMesh.pos(vi);
    }

    InterpolationCoordinates::BarycentricCoordinates bc(surfaceMesh.numVertices(), surfaceRestPositions.data(),
                                                        volumetricMesh.get());
    ES::SpMatD                                       W = bc.generateInterpolationMatrix();

    std::shared_ptr<SolidDeformationModel::DeformationModelManager> dmm =
        std::make_shared<SolidDeformationModel::DeformationModelManager>();

    dmm->setMesh(simMesh.get(), nullptr, nullptr);
    dmm->init(pgo::SolidDeformationModel::DeformationModelPlasticMaterial::VOLUMETRIC_DOF6, solver.elasticMaterial, 1);

    std::vector<double>                                               elementWeights(simMesh->getNumElements(), 1.0);
    std::shared_ptr<SolidDeformationModel::DeformationModelAssembler> assembler =
        std::make_shared<SolidDeformationModel::DeformationModelAssembler>(dmm, elementWeights.data());

    int n    = simMesh->getNumVertices();
    int n3   = n * 3;
    int nele = simMesh->getNumElements();

    ES::VXd plasticity(nele * 6);
    ES::M3d I = ES::M3d::Identity();
    for (int ei = 0; ei < nele; ei++) {
        const SolidDeformationModel::PlasticModel3DDeformationGradient* pm =
            dynamic_cast<const SolidDeformationModel::PlasticModel3DDeformationGradient*>(
                dmm->getDeformationModel(ei)->getPlasticModel());
        if (!pm) {
            SPDLOG_LOGGER_ERROR(Logging::lgr(), "Plastic model is not of type PlasticModel3DDeformationGradient.");
            return 1;
        }
        pm->toParam(I.data(), plasticity.data() + ei * dmm->getNumPlasticParameters());
    }

    ES::VXd restPosition(n3);
    for (int vi = 0; vi < n; vi++) {
        double p[3];
        simMesh->getVertex(vi, p);
        restPosition.segment<3>(vi * 3) = ES::V3d(p[0], p[1], p[2]);
    }

    std::shared_ptr<SolidDeformationModel::DeformationModelEnergy> elasticEnergy =
        std::make_shared<SolidDeformationModel::DeformationModelEnergy>(assembler, &restPosition, 0);
    elasticEnergy->setPlasticParams(plasticity);

    ES::VXd zero(n3);
    zero.setZero();

    ES::SpMatD K;
    elasticEnergy->createHessian(K);
    elasticEnergy->hessian(zero, K);

    std::vector<std::shared_ptr<ConstraintPotentialEnergies::MultipleVertexPulling>> pullingEnergies;
    std::vector<ES::VXd>                                                             pullingTargets, pullingTargetRests;
    Mesh::TriMeshGeo                                                                 tempMesh;
    for (const auto& fv : scene.fixedVertices) {
        std::vector<int> fixedVertices;
        if (BasicIO::read1DText(fv.filename.c_str(), std::back_inserter(fixedVertices)) != 0) {
            return 1;
        }
        std::sort(fixedVertices.begin(), fixedVertices.end());

        ES::VXd tgtVertexPositions(fixedVertices.size() * 3);
        ES::VXd tgtVertexRests(fixedVertices.size() * 3);
        for (int vi = 0; vi < (int)fixedVertices.size(); vi++) {
            ES::V3d movement(fv.movement[0], fv.movement[1], fv.movement[2]);
            tgtVertexPositions.segment<3>(vi * 3) = restPosition.segment<3>(fixedVertices[vi] * 3) + movement;
            tgtVertexRests.segment<3>(vi * 3)     = restPosition.segment<3>(fixedVertices[vi] * 3);
            tempMesh.addPos(tgtVertexPositions.segment<3>(vi * 3));
        }

        auto pullingEnergy = std::make_shared<ConstraintPotentialEnergies::MultipleVertexPulling>(
            K, restPosition.data(), (int)fixedVertices.size(), fixedVertices.data(), tgtVertexPositions.data(), nullptr,
            1);
        pullingEnergy->setCoeff(fv.coeff);
        pullingEnergies.push_back(pullingEnergy);
        pullingTargets.push_back(tgtVertexPositions);
        pullingTargetRests.push_back(tgtVertexRests);
    }

    tempMesh.save("fv.obj");

    std::vector<std::string> kinematicObjectFilenames;
    std::vector<ES::V3d>     kinematicObjectMovements;
    for (const auto& obj : scene.externalObjects) {
        kinematicObjectFilenames.push_back(obj.filename);
        kinematicObjectMovements.emplace_back(obj.movement[0], obj.movement[1], obj.movement[2]);
    }

    ES::SpMatD M;
    VolumetricMeshes::GenerateMassMatrix::computeMassMatrix(volumetricMesh.get(), M, true);

    ES::VXd g(n3);
    for (int vi = 0; vi < n; vi++) {
        g.segment<3>(vi * 3) = scene.extAcc;
    }

    ES::VXd fext(n3);
    ES::mv(M, g, fext);

    if (simulation.simType == "dynamic") {
        std::vector<Mesh::TriMeshGeo> kinematicObjects;
        for (const auto& filename : kinematicObjectFilenames) {
            kinematicObjects.emplace_back();
            if (kinematicObjects.back().load(filename) != true) {
                return 1;
            }

            for (int vi = 0; vi < kinematicObjects.back().numVertices(); vi++) {
                kinematicObjects.back().pos(vi) *= mesh.scale;
            }
        }

        std::vector<Mesh::TriMeshRef> kinematicObjectsRef;
        for (int i = 0; i < (int)kinematicObjects.size(); i++) {
            kinematicObjectsRef.emplace_back(kinematicObjects[i]);
        }

        const bool contactEnabled = isDynamicContactEnabled(contact);

        std::shared_ptr<Contact::TriangleMeshExternalContactHandler> externalContactHandler;
        std::shared_ptr<Contact::TriangleMeshSelfContactHandler>     selfCD;
        if (kinematicObjectsRef.size() && contactEnabled) {
            externalContactHandler = std::make_shared<Contact::TriangleMeshExternalContactHandler>(
                surfaceMesh.positions(), surfaceMesh.triangles(), n3, kinematicObjectsRef, contact.contactSamples,
                &bc.getEmbeddingVertexIndices(), &bc.getEmbeddingWeights());
        }

        if (contactEnabled && contact.enableSelfContact) {
            selfCD = std::make_shared<Contact::TriangleMeshSelfContactHandler>(
                surfaceMesh.positions(), surfaceMesh.triangles(), n3, contact.contactSamples,
                &bc.getEmbeddingVertexIndices(), &bc.getEmbeddingWeights());
        }

        std::shared_ptr<Simulation::ImplicitBackwardEulerTimeIntegrator> intg =
            std::make_shared<Simulation::ImplicitBackwardEulerTimeIntegrator>(M, elasticEnergy, solver.dampingParams[0],
                                                                              solver.dampingParams[1], simulation.timestep,
                                                                              solver.solverMaxIter, solver.solverEps);

#if defined(PGO_HAS_KNITRO)
        intg->setSolverOption(Simulation::TimeIntegratorSolverOption::SO_KNITRO);
        intg->setSolverConfigFile("config.opt");
#endif

        for (auto pullingEnergy : pullingEnergies)
            intg->addImplicitForceModel(pullingEnergy, 0, 0);

        intg->setExternalForce(fext.data());

        ES::VXd x = restPosition, u(n3);
        ES::VXd uvel(n3), uacc(n3), usurf(surfn3);

        u.setZero();
        uvel.setZero();
        uacc.setZero();
        usurf.setZero();

        for (int i = 0; i < n; i++) {
            uvel.segment<3>(i * 3) = scene.initialVel;
            u.segment<3>(i * 3)    = scene.initialDisp;
        }

        if (!std::filesystem::exists(runtime.outputFolder)) {
            std::filesystem::create_directories(runtime.outputFolder);
        }

        int frameStart = 0;
        if (runtime.restartEnabled) {
            for (int framei = simulation.numSimSteps - 1; framei >= 0; framei--) {
                std::string filename = fmt::format("{}/deform{:04d}.u", runtime.outputFolder, framei);
                if (!std::filesystem::exists(filename)) {
                    continue;
                }

                ES::MXd uMat(n3, 3);
                if (ES::readMatrix(filename.c_str(), uMat) == 0) {
                    frameStart     = framei;
                    u.noalias()    = uMat.col(0);
                    uvel.noalias() = uMat.col(1);
                    uacc.noalias() = uMat.col(2);
                    std::cout << "Restarting from frame " << framei << std::endl;
                    break;
                }
            }

            if (runtime.restartPause) {
                std::cout << frameStart << std::endl;
                std::cin.get();
            }
        }

        ES::mv(W, u, usurf);

        for (int framei = frameStart + 1; framei < simulation.numSimSteps; framei++) {
            intg->clearGeneralImplicitForceModel();
            intg->clearAlphaTestFunc();

            double ratio = (double)framei / (simulation.numSimSteps - 1);
            for (size_t pi = 0; pi < pullingEnergies.size(); pi++) {
                ES::VXd restTgt = pullingTargetRests[pi];
                ES::VXd curTgt  = restTgt * (1 - ratio) + pullingTargets[pi] * ratio;
                pullingEnergies[pi]->setTargetPos(curTgt.data());

                std::cout << "Frame " << framei << ", attachment " << pi << " target: " << curTgt.transpose().head(3)
                          << std::endl;
            }

            std::shared_ptr<Contact::PointPenetrationEnergy> extContactEnergy;
            Contact::PointPenetrationEnergyBuffer*           extContactBuffer = nullptr;
            std::shared_ptr<Contact::PointPenetrationBarrierEnergy> extBarrierEnergy;
            Contact::PointPenetrationBarrierEnergyBuffer*           extBarrierBuffer = nullptr;
            const bool shouldUseIpcAlphaFilter =
                contact.contactModel == "ipc-barrier" && contact.ipcEnableFeasibleLineSearch;
            bool   hasExternalAlphaUpperBound = false;
            bool   hasSelfAlphaUpperBound     = false;
            double externalDSafe              = 0.0;
            double selfDSafe                  = 0.0;
            const auto currentPositionFunc = makeCurrentPositionFunction(restPosition);
            if (externalContactHandler) {
                if (contact.contactModel == "penalty") {
                    externalContactHandler->execute(usurf.data());
                } else if (contact.contactModel == "ipc-barrier") {
                    externalContactHandler->execute(usurf.data(), externalIpcDhat);
                }

                if (externalContactHandler->getNumActiveSamples()) {
                    if (contact.contactModel == "penalty") {
                        extContactEnergy = externalContactHandler->buildContactEnergy(1);
                        extContactBuffer = extContactEnergy->allocateBuffer();
                        extContactEnergy->setComputePosFunction(currentPositionFunc);
                        extContactEnergy->setBuffer(extContactBuffer);
                        extContactEnergy->setCoeff(contact.contactStiffness);

                        extContactEnergy->setFrictionCoeff(contact.contactFrictionCoeff);
                        extContactEnergy->setComputeLastPosFunction(makeLastPositionFunction(restPosition, u));
                        extContactEnergy->setVelEps(contact.contactVelEps);
                        extContactEnergy->setTimestep(simulation.timestep);

                        intg->addGeneralImplicitForceModel(extContactEnergy, 0, 0);
                    } else if (contact.contactModel == "ipc-barrier") {
                        extBarrierEnergy = externalContactHandler->buildBarrierEnergy(externalIpcDhat);
                        extBarrierBuffer = extBarrierEnergy->allocateBuffer();
                        extBarrierEnergy->setComputePosFunction(currentPositionFunc);
                        extBarrierEnergy->setBuffer(extBarrierBuffer);
                        extBarrierEnergy->setCoeff(externalIpcKappa);

                        intg->addGeneralImplicitForceModel(extBarrierEnergy, 0, 0);

                        if (shouldUseIpcAlphaFilter) {
                            hasExternalAlphaUpperBound = true;
                            externalDSafe =
                                Contact::PointPenetrationBarrierEnergy::normalizePositiveZero(externalIpcDhat, -1.0);
                        }
                    }
                }
            }

            std::shared_ptr<Contact::PointTrianglePairCouplingEnergyWithCollision> selfContactEnergy;
            Contact::PointTrianglePairCouplingEnergyWithCollisionBuffer*           selfContactEnergyBuf = nullptr;
            std::shared_ptr<Contact::PointTrianglePairBarrierEnergy> selfBarrierEnergy;
            Contact::PointTrianglePairBarrierEnergyBuffer*           selfBarrierBuffer = nullptr;

            if (selfCD) {
                if (contact.contactModel == "penalty") {
                    selfCD->execute(usurf.data());

                    if (selfCD->getCollidingTrianglePair().size() > 0) {
                        selfCD->handleContactDCD(0, 100);

                        selfContactEnergy = selfCD->buildContactEnergy();
                        selfContactEnergy->setToPosFunction([&restPosition](const ES::V3d& x, ES::V3d& p, int offset) {
                            p = x + restPosition.segment<3>(offset);
                        });

                        selfContactEnergy->setToLastPosFunction(
                            [&restPosition, &u](const ES::V3d& x, ES::V3d& p, int offset) {
                                p = restPosition.segment<3>(offset) + u.segment<3>(offset);
                            });

                        selfContactEnergyBuf = selfContactEnergy->allocateBuffer();
                        selfContactEnergy->setBuffer(selfContactEnergyBuf);
                        selfContactEnergy->setCoeff(contact.contactStiffness);
                        selfContactEnergy->computeClosestPosition(u.data());

                        selfContactEnergy->setFrictionCoeff(contact.contactFrictionCoeff);
                        selfContactEnergy->setTimestep(simulation.timestep);
                        selfContactEnergy->setVelEps(contact.contactVelEps);

                        intg->addGeneralImplicitForceModel(selfContactEnergy, 0, 0);
                    }
                } else if (contact.contactModel == "ipc-barrier") {
                    selfCD->execute(usurf.data(), selfIpcDhat);
                    SPDLOG_LOGGER_INFO(Logging::lgr(), "# self active pairs: {}", selfCD->getNumActivePairs());

                    if (selfCD->getNumActivePairs() > 0) {
                        selfBarrierEnergy = selfCD->buildBarrierEnergy(selfIpcDhat);
                        selfBarrierBuffer = selfBarrierEnergy->allocateBuffer();
                        selfBarrierEnergy->setToPosFunction(currentPositionFunc);
                        selfBarrierEnergy->setBuffer(selfBarrierBuffer);
                        selfBarrierEnergy->setCoeff(selfIpcKappa);

                        intg->addGeneralImplicitForceModel(selfBarrierEnergy, 0, 0);

                        if (shouldUseIpcAlphaFilter) {
                            hasSelfAlphaUpperBound = true;
                            selfDSafe =
                                Contact::PointTrianglePairBarrierEnergy::normalizePositiveZero(selfIpcDhat, -1.0);
                        }
                    }
                }
            }

            if (shouldUseIpcAlphaFilter && (hasExternalAlphaUpperBound || hasSelfAlphaUpperBound)) {
                intg->setAlphaTestFunc(
                    [&u, &contact, externalContactHandler, selfCD, hasExternalAlphaUpperBound, hasSelfAlphaUpperBound,
                     externalDSafe, selfDSafe, mergedAlphaCallbackLogged = false,
                     alphaCallbackLogged = false](const ES::VXd& z, const ES::VXd& dz) mutable {
                        ES::VXd currentU = u + z;
                        double  alphaUpper = 1.0;
                        double  externalAlphaUpper = 1.0;
                        double  selfAlphaUpper     = 1.0;

                        if (!alphaCallbackLogged) {
                            SPDLOG_LOGGER_INFO(Logging::lgr(), "IPC feasible alpha callback active.");
                            alphaCallbackLogged = true;
                        }

                        if (hasExternalAlphaUpperBound) {
                            externalAlphaUpper = externalContactHandler->computeEmbeddedAlphaUpperBound(
                                currentU, dz, contact.ipcAlphaSafety, externalDSafe);
                            alphaUpper = std::min(alphaUpper, externalAlphaUpper);
                        }

                        if (hasSelfAlphaUpperBound) {
                            selfAlphaUpper = selfCD->computeEmbeddedAlphaUpperBound(
                                currentU, dz, contact.ipcAlphaSafety, selfDSafe);
                            alphaUpper = std::min(alphaUpper, selfAlphaUpper);
                        }

                        if (hasExternalAlphaUpperBound && hasSelfAlphaUpperBound && !mergedAlphaCallbackLogged) {
                            SPDLOG_LOGGER_INFO(Logging::lgr(), "IPC feasible alpha merged callback active.");
                            mergedAlphaCallbackLogged = true;
                        }

                        if (alphaUpper < 1.0 - 1e-12) {
                            SPDLOG_LOGGER_INFO(Logging::lgr(),
                                               "IPC feasible alpha upper bound: {} (external: {}, self: {})",
                                               alphaUpper, externalAlphaUpper, selfAlphaUpper);
                        }

                        return alphaUpper;
                    });
            }

            intg->setqState(u, uvel, uacc);

            intg->doTimestep(1, 2, 1);
            const int solverRet = intg->getSolverReturn();
            SPDLOG_LOGGER_INFO(Logging::lgr(), "Frame {} solver status: {}.", framei,
                               describeSolverStatus(intg->getSolverOption(), solverRet));

            if (extContactBuffer && extContactEnergy) {
                extContactEnergy->freeBuffer(extContactBuffer);
            }

            if (extBarrierBuffer && extBarrierEnergy) {
                extBarrierEnergy->freeBuffer(extBarrierBuffer);
            }

            if (selfContactEnergyBuf && selfContactEnergy) {
                selfContactEnergy->freeBuffer(selfContactEnergyBuf);
            }

            if (selfBarrierBuffer && selfBarrierEnergy) {
                selfBarrierEnergy->freeBuffer(selfBarrierBuffer);
            }

            if (NonlinearOptimization::NewtonRaphsonSolver::isHardFailure(solverRet)) {
                SPDLOG_LOGGER_ERROR(Logging::lgr(), "Frame {} terminated with solver status {}.", framei,
                                    describeSolverStatus(intg->getSolverOption(), solverRet));
                return solverRet;
            }

            const bool isNewtonSmallStepStall =
                intg->getSolverOption() == Simulation::TimeIntegratorSolverOption::SO_NEWTON &&
                solverRet == NonlinearOptimization::NewtonRaphsonSolver::SR_STALLED_SMALL_STEP;
            if (solverRet > 0 && !isNewtonSmallStepStall) {
                SPDLOG_LOGGER_WARN(Logging::lgr(), "Frame {} terminated with non-converged solver status {}.", framei,
                                   describeSolverStatus(intg->getSolverOption(), solverRet));
            }

            intg->getq(u);
            intg->getqvel(uvel);
            intg->getqacc(uacc);

            ES::mv(W, u, usurf);

            if (framei % simulation.dumpInterval == 0) {
                ES::VXd psurf = surfaceRestPositions + usurf;

                Mesh::TriMeshGeo mesh = surfaceMesh;
                for (int vi = 0; vi < mesh.numVertices(); vi++) {
                    mesh.pos(vi) = psurf.segment<3>(vi * 3) / config.mesh.scale;
                }
                mesh.save(fmt::format("{}/ret{:04d}.obj", runtime.outputFolder, framei / simulation.dumpInterval));
            }

            for (size_t eobji = 0; eobji < kinematicObjects.size(); eobji++) {
                ES::V3d movement = kinematicObjectMovements[eobji] / (simulation.numSimSteps - 1);
                for (int vi = 0; vi < kinematicObjects[eobji].numVertices(); vi++) {
                    kinematicObjects[eobji].pos(vi) += movement;
                }
                externalContactHandler->updateExternalSurface(eobji, kinematicObjectsRef[eobji]);
            }

            ES::MXd uMat(n3, 3);
            uMat.col(0) = u;
            uMat.col(1) = uvel;
            uMat.col(2) = uacc;

            ES::writeMatrix(fmt::format("{}/deform{:04d}.u", runtime.outputFolder, framei).c_str(), uMat);
        }
    } else if (simulation.simType == "static") {
        std::shared_ptr<PredefinedPotentialEnergies::LinearPotentialEnergy> externalForcesEnergy =
            std::make_shared<PredefinedPotentialEnergies::LinearPotentialEnergy>(fext);

        std::shared_ptr<NonlinearOptimization::PotentialEnergies> energyAll =
            std::make_shared<NonlinearOptimization::PotentialEnergies>(n3);
        energyAll->addPotentialEnergy(elasticEnergy);
        for (auto eng : pullingEnergies)
            energyAll->addPotentialEnergy(eng, 1.0);
        energyAll->addPotentialEnergy(externalForcesEnergy, -1.0);
        energyAll->init();

        NonlinearOptimization::NewtonRaphsonSolver::SolverParam solverParam;

        ES::VXd u(n3);
        u.setZero();

        energyAll->printEnergy(u);

        NonlinearOptimization::NewtonRaphsonSolver solver(u.data(), solverParam, energyAll, std::vector<int>(),
                                                          nullptr);
        const int solverRet = solver.solve(u.data(), config.solver.solverMaxIter, config.solver.solverEps, 2);
        if (NonlinearOptimization::NewtonRaphsonSolver::isHardFailure(solverRet)) {
            SPDLOG_LOGGER_ERROR(Logging::lgr(), "Static solve failed with solver status {}.",
                                describeSolverStatus(Simulation::TimeIntegratorSolverOption::SO_NEWTON, solverRet));
            return solverRet;
        }
        if (solverRet > 0) {
            SPDLOG_LOGGER_WARN(Logging::lgr(), "Static solve terminated with non-converged solver status {}.",
                               describeSolverStatus(Simulation::TimeIntegratorSolverOption::SO_NEWTON, solverRet));
        }

        ES::VXd x = restPosition + u;

        ES::VXd xsurf(surfn3);
        ES::mv(W, x, xsurf);

        Mesh::TriMeshGeo meshOut = surfaceMesh;
        for (int vi = 0; vi < meshOut.numVertices(); vi++) {
            meshOut.pos(vi) = xsurf.segment<3>(vi * 3);
        }
        meshOut.save(runtime.outputFolder);
    }

    return 0;
}

}  // namespace api
}  // namespace pgo
