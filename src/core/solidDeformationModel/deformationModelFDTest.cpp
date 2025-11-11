/*
author: Bohan Wang
copyright to USC,MIT,NUS
*/

#include "deformationModelFDTest.h"

#include "simulationMesh.h"
#include "deformationModel.h"
#include "deformationModelManager.h"

// #include "elementLocalDirection.h"
#include "tetMesh.h"
#include "deformationModelAssembler.h"
#include "deformationModelEnergy.h"

#include "finiteDifference.h"
#include "pgoLogging.h"
#include "EigenSupport.h"
#include "triMeshGeo.h"

#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>

#include <random>
#include <chrono>

using namespace pgo;
using namespace pgo::SolidDeformationModel;
using namespace pgo::NonlinearOptimization;

namespace ES = pgo::EigenSupport;

using hclock = std::chrono::high_resolution_clock;

inline double dur(const hclock::time_point &t1, const hclock::time_point &t2)
{
  return static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count()) / 1000.0;
}

int SolidDeformationModel::fdTestTetMesh(const char *tetMeshFilename, int numTestDOFs, int fdTestFlag)
{
  SPDLOG_LOGGER_INFO(Logging::lgr(), "Loading mesh: {} ...", tetMeshFilename);

  VolumetricMeshes::TetMesh tetMesh(tetMeshFilename);
  std::shared_ptr<SimulationMesh> mesh(loadTetMesh(&tetMesh));
  if (mesh == nullptr)
    return 1;

  std::vector<DeformationModelPlasticMaterial> testingPlasticDOFs = { DeformationModelPlasticMaterial::VOLUMETRIC_DOF6 };
  testingPlasticDOFs.push_back(DeformationModelPlasticMaterial::VOLUMETRIC_DOF0);

  std::vector<DeformationModelElasticMaterial> testMaterials;

  testMaterials.push_back(DeformationModelElasticMaterial::STABLE_NEO);
  testMaterials.push_back(DeformationModelElasticMaterial::VOLUME);
  testMaterials.push_back(DeformationModelElasticMaterial::MOONEY_RIVLIN);
  testMaterials.push_back(DeformationModelElasticMaterial::STVK);
  testMaterials.push_back(DeformationModelElasticMaterial::INV_STVK);

  // testMaterials.push_back(DeformationModelElasticMaterial::STVK_VOL);
  // testMaterials.push_back(DeformationModelElasticMaterial::VOLUME);

  std::random_device rd;
  std::mt19937 eng(rd());

  int n3 = mesh->getNumVertices() * 3;
  int nele = mesh->getNumElements();

  ES::M3d S;
  S << 1.8673, 0.6765, 1.8324,
    0.6765, 0.4383, 0.9276,
    1.8324, 0.9276, 2.4211;

  for (DeformationModelPlasticMaterial plasticMat : testingPlasticDOFs) {
    for (const auto elasticMat : testMaterials) {
      SPDLOG_LOGGER_INFO(Logging::lgr(), "*******************************************");
      SPDLOG_LOGGER_INFO(Logging::lgr(), "*******************************************");
      SPDLOG_LOGGER_INFO(Logging::lgr(), "*******************************************");
      SPDLOG_LOGGER_INFO(Logging::lgr(), "Testing PM: {}; EM: {}", (int)plasticMat, (int)elasticMat);

      if (elasticMat == DeformationModelElasticMaterial::MOONEY_RIVLIN) {
        ES::M3d Cpq;
        Cpq << 0.0, 1.0, 1.0,
          1, 1.0, 1.0,
          1.0, 1.0, 1.0;

        ES::V2d D;
        D << 1, 1;

        SimulationMeshMooneyRivlinMaterial *mat = new SimulationMeshMooneyRivlinMaterial(2, 2, Cpq.data(), D.data());
        mesh->setMaterial(-1, mat);
      }
      else {
        SimulationMeshENuMaterial *mat = new SimulationMeshENuMaterial(1000.0, 0.45);
        mesh->setMaterial(-1, mat);
      }

      std::shared_ptr<DeformationModelManager> dmm = std::make_shared<DeformationModelManager>();
      dmm->setMesh(mesh.get());
      dmm->init(plasticMat, elasticMat, 0);

      int nplastic = dmm->getNumPlasticParameters();
      int nelastic = 0;
      ES::VXd scaleParam(nplastic);

      if (nplastic == 3) {
        scaleParam.noalias() = S.diagonal();
      }
      else if (nplastic == 6) {
        double s[6] = {
          S(0, 0), S(0, 1), S(0, 2),
          S(1, 1), S(1, 2), S(2, 2)
        };
        scaleParam.noalias() = Eigen::Map<const Eigen::Matrix<double, 6, 1>>(s);
      }
      else if (nplastic == 12) {
        Eigen::SelfAdjointEigenSolver<ES::M3d> eig(S);
        ES::M3d R = eig.eigenvectors();
        ES::V3d s = eig.eigenvalues();

        scaleParam.segment<3>(0) = s;
        scaleParam.segment<9>(3) = Eigen::Map<const ES::V9d>(R.data());
      }

      FiniteDifference fd(FiniteDifference::M_FIVE_POINT, 1e-7);
      std::uniform_real_distribution<double> distrib(-0.1, 0.1);

      ES::VXd x(12 + nplastic + nelastic);
      int offset[4] = {
        0,
        12,
        12 + nplastic,
        12 + nplastic + nelastic,
      };

      double errors[3][2] = { { 0, 0 }, { 0, 0 }, { 0, 0 } };
      tbb::enumerable_thread_specific<ES::MXd> d2EdxdaTLS(ES::MXd(12, nplastic));
      tbb::enumerable_thread_specific<ES::MXd> d2EdadaTLS(ES::MXd(nplastic, nplastic));
      tbb::enumerable_thread_specific<ES::MXd> d2EdadbTLS(ES::MXd(nplastic, nelastic));
      tbb::enumerable_thread_specific<ES::MXd> d2EdbdbTLS(ES::MXd(nelastic, nelastic));
      tbb::enumerable_thread_specific<ES::MXd> d2EdxdbTLS(ES::MXd(12, nelastic));

      int numTestElements = std::min(nele, 500);

      if ((fdTestFlag & FD_TT_ELEMENT)) {
        hclock::time_point t1 = hclock::now();

        for (int ele = 0; ele < numTestElements; ele++) {
          const DeformationModel *fem = dmm->getDeformationModel(ele);
          DeformationModelCacheData *cache = fem->allocateCacheData();

          for (int j = 0; j < 4; j++) {
            ES::V3d p;
            mesh->getVertex(ele, j, p.data());
            x.segment<3>(j * 3) = p;
          }

          for (int j = 0; j < 12; j++) {
            x[j] += distrib(eng);
          }
          x.segment(12, nplastic) = scaleParam;

          FiniteDifference::EvalFunc evalFuncE = [&](const double *x_param, double *E, double *grad, double *hess) {
            fem->prepareData(x_param, scaleParam.data(), nullptr, cache);

            if (E) {
              *E = fem->computeEnergy(cache);
            }

            if (grad) {
              memset(grad, 0, sizeof(double) * 12);
              fem->compute_dE_dx(cache, grad);
            }

            if (hess) {
              fem->compute_d2E_dx2(cache, hess);
            }
          };

          double gradError, hessError;
          fd.testEnergy(evalFuncE, 12, true, true, -1, x.data(), &gradError, &hessError);

          // LGI << "Ele " << ele << "; Error 12: " << gradError << ',' << hessError;
          errors[0][0] = std::max(errors[0][0], gradError);
          errors[0][1] = std::max(errors[0][1], hessError);

          if (nplastic > 0) {
            FiniteDifference::EvalFunc evalFuncEP = [&](const double *x_param, double *E, double *grad, double *hess) {
              fem->prepareData(x_param, x_param + offset[1], nullptr, cache);

              if (E) {
                *E = fem->computeEnergy(cache);
              }

              if (grad) {
                memset(grad, 0, sizeof(double) * offset[2]);
                fem->compute_dE_dx(cache, grad);
                fem->compute_dE_da(cache, grad + offset[1]);
              }

              if (hess) {
                ES::M12d hess0;
                ES::MXd &d2E_dxda = d2EdxdaTLS.local();
                ES::MXd &d2E_da2 = d2EdadaTLS.local();

                fem->compute_d2E_dx2(cache, hess0.data());
                fem->compute_d2E_dxda(cache, d2E_dxda.data());
                fem->compute_d2E_da2(cache, d2E_da2.data());

                Eigen::Map<ES::MXd> hMap(hess, offset[2], offset[2]);
                hMap.block<12, 12>(0, 0) = hess0;
                hMap.block(0, 12, 12, nplastic) = d2E_dxda;
                hMap.block(12, 0, nplastic, 12) = d2E_dxda.transpose();
                hMap.block(12, 12, nplastic, nplastic) = d2E_da2;
              }
            };

            fd.testEnergy(evalFuncEP, offset[2], true, true, -1, x.data(), &gradError, &hessError);
            // LGI << "Ele " << ele << "; Error: " << gradError << ',' << hessError;
            errors[1][0] = std::max(errors[1][0], gradError);
            errors[1][1] = std::max(errors[1][1], hessError);
          }

          if (nelastic > 0) {
            FiniteDifference::EvalFunc evalFuncEPA = [&](const double *x_param, double *E, double *grad, double *hess) {
              fem->prepareData(x_param, x_param + offset[1], x_param + offset[2], cache);

              if (E) {
                *E = fem->computeEnergy(cache);
              }

              if (grad) {
                memset(grad, 0, sizeof(double) * offset[3]);
                fem->compute_dE_dx(cache, grad);
                fem->compute_dE_da(cache, grad + offset[1]);
                fem->compute_dE_db(cache, grad + offset[2]);
              }

              if (hess) {
                ES::M12d hess0;
                ES::MXd &d2E_dxda = d2EdxdaTLS.local();
                ES::MXd &d2E_da2 = d2EdadaTLS.local();

                fem->compute_d2E_dx2(cache, hess0.data());
                fem->compute_d2E_dxda(cache, d2E_dxda.data());
                fem->compute_d2E_da2(cache, d2E_da2.data());

                ES::MXd &d2E_dxdb = d2EdxdbTLS.local();
                ES::MXd &d2E_db2 = d2EdbdbTLS.local();
                ES::MXd &d2E_dadb = d2EdadbTLS.local();

                fem->compute_d2E_dxdb(cache, d2E_dxdb.data());
                fem->compute_d2E_dadb(cache, d2E_dadb.data());
                fem->compute_d2E_db2(cache, d2E_db2.data());

                Eigen::Map<ES::MXd> hMap(hess, offset[3], offset[3]);
                hMap.block<12, 12>(0, 0) = hess0;
                hMap.block(0, 12, 12, nplastic) = d2E_dxda;
                hMap.block(12, 0, nplastic, 12) = d2E_dxda.transpose();
                hMap.block(12, 12, nplastic, nplastic) = d2E_da2;

                hMap.block(0, 12 + nplastic, 12, nelastic) = d2E_dxdb;
                hMap.block(12 + nplastic, 0, nelastic, 12) = d2E_dxdb.transpose();
                hMap.block(12, 12 + nplastic, nplastic, nelastic) = d2E_dadb;
                hMap.block(12 + nplastic, 12, nelastic, nplastic) = d2E_dadb.transpose();

                hMap.block(12 + nplastic, 12 + nplastic, nelastic, nelastic) = d2E_db2;
              }
            };

            fd.testEnergy(evalFuncEPA, offset[3], true, true, -1, x.data(), &gradError, &hessError);
            // LGI << "Ele " << ele << "; Error: " << gradError << ',' << hessError;
            errors[2][0] = std::max(errors[2][0], gradError);
            errors[2][1] = std::max(errors[2][1], hessError);
          }

          if (ele % 10 == 0)
            std::cout << ele << ' ' << std::flush;

          fem->freeCacheData(cache);
        }

        hclock::time_point t2 = hclock::now();

        std::cout << "\n\n"
                  << std::endl;
        SPDLOG_LOGGER_INFO(Logging::lgr(), "Error: {}, {}", errors[0][0], errors[0][1]);
        SPDLOG_LOGGER_INFO(Logging::lgr(), "Error: {}, {}", errors[1][0], errors[1][1]);
        SPDLOG_LOGGER_INFO(Logging::lgr(), "Error: {}, {}", errors[2][0], errors[2][1]);
        SPDLOG_LOGGER_INFO(Logging::lgr(), "Time(per tet) : {} ms", dur(t1, t2) / (double)numTestElements);
      }  // if element test

      ES::VXd tetMeshRestPositions(n3);
      for (int i = 0; i < mesh->getNumVertices(); i++) {
        ES::V3d p;
        mesh->getVertex(i, p.data());
        tetMeshRestPositions.segment<3>(i * 3) = p;
      }
      double gradError = 0, hessError = 0;

      if (fdTestFlag & FD_TT_ASSEM_E) {
        SPDLOG_LOGGER_INFO(Logging::lgr(), "Assembler E test...");

        ES::VXd x = tetMeshRestPositions;
        for (int i = 0; i < n3; i++) {
          x[i] = tetMeshRestPositions[i] + distrib(eng);
        }

        ES::VXd scales(nele * nplastic);
        for (int hi = 0; hi < nele; hi++) {
          scales.segment(hi * nplastic, nplastic) = scaleParam;
        }

        ES::VXd elasticParams;

        std::shared_ptr<DeformationModelAssembler> forceModelAssembler;
        forceModelAssembler = std::make_shared<DeformationModelAssembler>(dmm, nullptr);

        std::shared_ptr<DeformationModelEnergy> energy = std::make_shared<DeformationModelEnergy>(forceModelAssembler);
        energy->setPlasticParams(scales);

        fd.testEnergy(energy, true, false, -1.0, x.data(), -1, &gradError, nullptr);
        fd.testEnergy(energy, false, true, -1.0, x.data(), numTestDOFs, nullptr, &hessError);

        SPDLOG_LOGGER_INFO(Logging::lgr(), "Assembler error: {}, {}", gradError, hessError);

        for (int hi = 0; hi < nele; hi++) {
          scales.segment(hi * nplastic, nplastic) = scaleParam + ES::VXd::Constant(nplastic, 1.0 * hi / nele);
        }

        energy->setPlasticParams(scales);

        fd.testEnergy(energy, true, false, -1.0, x.data(), -1, &gradError, nullptr);
        fd.testEnergy(energy, false, true, -1.0, x.data(), numTestDOFs, nullptr, &hessError);

        SPDLOG_LOGGER_INFO(Logging::lgr(), "Assembler error: {}, {}", gradError, hessError);
      }

      // if ((fdTestFlag & FD_TT_ASSEM_EFP) && nplastic > 0) {
      //   LGI << "Assembler EFP test...";
      //   std::shared_ptr<DeformationModelAssemblerElasticFullPlastic> forceModelAssemblerEFP = std::make_shared<DeformationModelAssemblerElasticFullPlastic>(dmm, nullptr, 0);

      //  ES::VXd x(forceModelAssemblerEFP->getNumDOFs());
      //  for (int i = 0; i < n3; i++) {
      //    x[i] = tetMeshRestPositions[i] + distrib(eng);
      //  }

      //  for (int hi = 0; hi < nele; hi++) {
      //    x.segment(n3 + hi * nplastic, nplastic) = scaleParam;
      //  }

      //  std::shared_ptr<DeformationModelEnergy> energy = std::make_shared<DeformationModelEnergy>(forceModelAssemblerEFP);

      //  fd.testEnergy(energy, true, false, -1.0, x.data(), -1, &gradError, nullptr);
      //  fd.testEnergy(energy, false, true, -1.0, x.data(), numDOFs, nullptr, &hessError);

      //  LGI << gradError << ',' << hessError;

      //  for (int hi = 0; hi < nele; hi++) {
      //    x.segment(n3 + hi * nplastic, nplastic) = scaleParam + ES::VXd::Constant(nplastic, 0.1 * (hi % 5));
      //    dmm->getDeformationModel(0)->getPlasticModel()->compute_paramfull(x.data() + n3 + hi * nplastic);
      //  }

      //  fd.testEnergy(energy, true, false, -1.0, x.data(), -1, &gradError, nullptr);
      //  fd.testEnergy(energy, false, true, -1.0, x.data(), numDOFs, nullptr, &hessError);

      //  LGI << gradError << ',' << hessError;
      //}

      std::cin.get();
    }
  }

  return 0;
}

int SolidDeformationModel::fdTestShellMesh(const char *surfaceMeshFilename, int numTestDOFs, int fdTestFlag)
{
  Mesh::TriMeshGeo surfaceMesh;
  if (!surfaceMesh.load(surfaceMeshFilename)) {
    SPDLOG_LOGGER_ERROR(Logging::lgr(), "Failed to load surface mesh: {}", surfaceMeshFilename);
    return 1;
  }

  SolidDeformationModel::SimulationMeshENuhMaterial mat(1e5, 0.4, 1e-3);

  std::shared_ptr<SimulationMesh> mesh(SolidDeformationModel::loadShellMesh(surfaceMesh, &mat));
  if (mesh == nullptr) {
    return 1;
  }

  std::vector<DeformationModelPlasticMaterial> testingPlasticDOFs = {
    DeformationModelPlasticMaterial::SHELL_FF_DOF1,
    DeformationModelPlasticMaterial::SHELL_FF_DOF0,
  };

  std::vector<DeformationModelElasticMaterial> testMaterials = {
    DeformationModelElasticMaterial::KOITER_STVK,
    // DeformationModelElasticMaterial::KOITER_FABRIC
  };

  std::random_device rd;
  std::mt19937 eng(rd());
  std::uniform_real_distribution<double> distrib0(1.1, 2);

  int n3 = mesh->getNumVertices() * 3;
  int nele = mesh->getNumElements();

  for (DeformationModelPlasticMaterial plasticMat : testingPlasticDOFs) {
    for (const auto elasticMat : testMaterials) {
      SPDLOG_LOGGER_INFO(Logging::lgr(), "*******************************************");
      SPDLOG_LOGGER_INFO(Logging::lgr(), "*******************************************");
      SPDLOG_LOGGER_INFO(Logging::lgr(), "*******************************************");
      SPDLOG_LOGGER_INFO(Logging::lgr(), "Testing PM: {}; EM: {}", (int)plasticMat, (int)elasticMat);

      if (elasticMat == DeformationModelElasticMaterial::MOONEY_RIVLIN) {
        ES::M3d Cpq;
        Cpq << 0.0, 1.0, 1.0,
          1, 1.0, 1.0,
          1.0, 1.0, 1.0;

        ES::V2d D;
        D << 1, 1;

        SimulationMeshMooneyRivlinhMaterial *mat = new SimulationMeshMooneyRivlinhMaterial(2, 2, Cpq.data(), D.data(), 1e-3);
        mesh->setMaterial(-1, mat);
      }
      else {
        SimulationMeshENuhMaterial *mat = new SimulationMeshENuhMaterial(1000.0, 0.45, 1e-3);
        mesh->setMaterial(-1, mat);
      }

      std::shared_ptr<DeformationModelManager> dmm = std::make_shared<DeformationModelManager>();
      dmm->setMesh(mesh.get());
      dmm->init(plasticMat, elasticMat, 0);

      int nplastic = dmm->getNumPlasticParameters();
      int nelastic = dmm->getNumElasticParameters();

      ES::VXd scaleParam(nplastic);
      if (nplastic == 1) {
        scaleParam[0] = distrib0(eng);
      }

      ES::VXd elasticParam(nelastic);
      elasticParam.setOnes();

      if (elasticMat == DeformationModelElasticMaterial::KOITER_FABRIC) {
        // elasticParam << 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1e-3;
        elasticParam << 1, 1, 1, 1, 1, 1, 1, 1000, 1000, 1000, 1, 1e-3;
      }
      else if (elasticMat == DeformationModelElasticMaterial::KOITER_STVK) {
        elasticParam << 20000, 0.45, 10000, 0.3, 1e-3;
      }

      FiniteDifference fd(FiniteDifference::M_FIVE_POINT, 1e-7);
      std::uniform_real_distribution<double> distrib(-0.1, 0.1);

      ES::VXd x(mesh->getNumElementVertices() * 3 + nplastic + nelastic);
      int offset[4] = {
        0,
        mesh->getNumElementVertices() * 3,
        mesh->getNumElementVertices() * 3 + nplastic,
        mesh->getNumElementVertices() * 3 + nplastic + nelastic,
      };

      double errors[3][2] = { { 0, 0 }, { 0, 0 }, { 0, 0 } };
      tbb::enumerable_thread_specific<ES::MXd> d2EdxdaTLS(ES::MXd::Zero(mesh->getNumElementVertices() * 3, nplastic));
      tbb::enumerable_thread_specific<ES::MXd> d2EdadaTLS(ES::MXd::Zero(nplastic, nplastic));
      tbb::enumerable_thread_specific<ES::MXd> d2EdadbTLS(ES::MXd::Zero(nplastic, nelastic));
      tbb::enumerable_thread_specific<ES::MXd> d2EdbdbTLS(ES::MXd::Zero(nelastic, nelastic));
      tbb::enumerable_thread_specific<ES::MXd> d2EdxdbTLS(ES::MXd::Zero(mesh->getNumElementVertices() * 3, nelastic));

      int numTestElements = std::min(nele, 500);

      if ((fdTestFlag & FD_TT_ELEMENT)) {
        hclock::time_point t1 = hclock::now();

        for (int ele = 0; ele < numTestElements; ele++) {
          const DeformationModel *fem = dmm->getDeformationModel(ele);
          DeformationModelCacheData *cache = fem->allocateCacheData();

          for (int j = 0; j < mesh->getNumElementVertices(); j++) {
            ES::V3d p;
            int vid = mesh->getVertexIndex(ele, j);
            if (vid >= 0) {
              mesh->getVertex(ele, j, p.data());
              x.segment<3>(j * 3) = p;
            }
            else {
              x.segment<3>(j * 3).setZero();
            }
          }

          for (int j = 0; j < mesh->getNumElementVertices() * 3; j++) {
            x[j] += distrib(eng);
          }

          x.segment(mesh->getNumElementVertices() * 3, nplastic) = scaleParam;

          FiniteDifference::EvalFunc evalFuncE = [&](const double *x_param, double *E, double *grad, double *hess) {
            fem->prepareData(x_param, scaleParam.data(), elasticParam.data(), cache);

            if (E) {
              *E = fem->computeEnergy(cache);
            }

            if (grad) {
              memset(grad, 0, sizeof(double) * mesh->getNumElementVertices() * 3);
              ;
              fem->compute_dE_dx(cache, grad);
            }

            if (hess) {
              fem->compute_d2E_dx2(cache, hess);
            }
          };

          double gradError, hessError;
          fd.testEnergy(evalFuncE, mesh->getNumElementVertices() * 3, true, true, -1, x.data(), &gradError, &hessError);

          // LGI << "Ele " << ele << "; Error 12: " << gradError << ',' << hessError;
          errors[0][0] = std::max(errors[0][0], gradError);
          errors[0][1] = std::max(errors[0][1], hessError);

          if (nplastic > 0) {
            FiniteDifference::EvalVecFunc evalFuncP = [&](const double *x_param, double *g, double *jac) {
              fem->prepareData(x_param, x_param + offset[1], elasticParam.data(), cache);

              if (g) {
                memset(g, 0, sizeof(double) * offset[1]);
                fem->compute_dE_dx(cache, g);
              }

              if (jac) {
                fem->compute_d2E_dx2(cache, jac);
                fem->compute_d2E_dxda(cache, jac + offset[1] * offset[1]);
              }
            };

            fd.testVecFunc(evalFuncP, offset[1], offset[2], -1, x.data(), &gradError);
            // LGI << "Ele " << ele << "; Error: " << gradError << ',' << hessError;
            errors[1][0] = std::max(errors[1][0], gradError);
          }

          if (nelastic > 0) {
            FiniteDifference::EvalVecFunc evalFuncA = [&](const double *x_param, double *g, double *jac) {
              fem->prepareData(x_param, scaleParam.data(), x_param + offset[1], cache);

              if (g) {
                memset(g, 0, sizeof(double) * offset[3]);
                fem->compute_dE_dx(cache, g);
              }

              if (jac) {
                fem->compute_d2E_dx2(cache, jac);
                fem->compute_d2E_dxdb(cache, jac + offset[1] * offset[1]);
              }
            };

            ES::VXd xparam(offset[1] + nelastic);
            xparam.head(offset[1]) = x.head(offset[1]);
            xparam.tail(nelastic) = elasticParam;

            fd.testVecFunc(evalFuncA, offset[1], xparam.size(), -1, xparam.data(), &gradError);
            // LGI << "Ele " << ele << "; Error: " << gradError << ',' << hessError;
            errors[2][0] = std::max(errors[2][0], gradError);
            errors[2][1] = std::max(errors[2][1], hessError);
          }

          if (ele % 10 == 0)
            std::cout << ele << ' ' << std::flush;

          fem->freeCacheData(cache);
        }

        hclock::time_point t2 = hclock::now();

        std::cout << "\n\n"
                  << std::endl;

        SPDLOG_LOGGER_INFO(Logging::lgr(), "Error: {}, {}", errors[0][0], errors[0][1]);
        SPDLOG_LOGGER_INFO(Logging::lgr(), "Error: {}, {}", errors[1][0], errors[1][1]);
        SPDLOG_LOGGER_INFO(Logging::lgr(), "Error: {}, {}", errors[2][0], errors[2][1]);
        SPDLOG_LOGGER_INFO(Logging::lgr(), "Time(per tet) : {} ms", dur(t1, t2) / (double)numTestElements);
      }  // if element test

      ES::VXd tetMeshRestPositions(n3);
      for (int i = 0; i < mesh->getNumVertices(); i++) {
        ES::V3d p;
        mesh->getVertex(i, p.data());
        tetMeshRestPositions.segment<3>(i * 3) = p;
      }
      double gradError = 0, hessError = 0;

      if (fdTestFlag & FD_TT_ASSEM_E) {
        SPDLOG_LOGGER_INFO(Logging::lgr(), "Assembler E test...");

        ES::VXd x = tetMeshRestPositions;
        for (int i = 0; i < n3; i++) {
          x[i] = tetMeshRestPositions[i] + distrib(eng);
        }

        ES::VXd scales(nele * nplastic);
        for (int hi = 0; hi < nele; hi++) {
          scales.segment(hi * nplastic, nplastic) = scaleParam;
        }

        ES::VXd elasticParams;

        std::shared_ptr<DeformationModelAssembler> forceModelAssembler;
        forceModelAssembler = std::make_shared<DeformationModelAssembler>(dmm, nullptr);

        std::shared_ptr<DeformationModelEnergy> energy = std::make_shared<DeformationModelEnergy>(forceModelAssembler);

        fd.testEnergy(energy, true, false, -1.0, x.data(), -1, &gradError, nullptr);
        fd.testEnergy(energy, false, true, -1.0, x.data(), numTestDOFs, nullptr, &hessError);

        SPDLOG_LOGGER_INFO(Logging::lgr(), "Assembler error: {}, {}", gradError, hessError);

        for (int hi = 0; hi < nele; hi++) {
          scales.segment(hi * nplastic, nplastic) = scaleParam + ES::VXd::Constant(nplastic, 1.0 * hi / nele);
        }

        energy->setPlasticParams(scales);

        fd.testEnergy(energy, true, false, -1.0, x.data(), -1, &gradError, nullptr);
        fd.testEnergy(energy, false, true, -1.0, x.data(), numTestDOFs, nullptr, &hessError);

        SPDLOG_LOGGER_INFO(Logging::lgr(), "Assembler error: {}, {}", gradError, hessError);
      }
    }
  }

  return 0;
}