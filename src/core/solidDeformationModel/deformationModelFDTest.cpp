/*
author: Bohan Wang
copyright to USC,MIT
*/

#include "deformationModelFDTest.h"
#include "simulationMesh.h"
#include "deformationModelManager.h"
#include "deformationModel.h"
// #include "elementLocalDirection.h"
#include "tetMesh.h"
#include "deformationModelAssemblerElastic.h"
#include "deformationModelAssemblerElasticFullPlastic.h"
#include "deformationModelEnergy.h"
#include "plasticModel.h"

#include "finiteDifference.h"
#include "pgoLogging.h"
#include "EigenSupport.h"
#include "matrixIO.h"
#include "objMesh.h"
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

int SolidDeformationModel::fdTestTetMesh(const char *tetMeshFilename, const char *surfaceMeshFilename,
  const char *basisFilename, const char *fiberDirectionFilename, int numDOFs, int fdTestType)
{
  LGI << "Loading mesh: " << tetMeshFilename << " ...";
  TetMesh *tetMesh = loadTetMesh(tetMeshFilename);
  SimulationMesh *mesh = loadTetMesh(tetMesh);
  if (mesh == nullptr)
    return 1;

  ObjMesh surfaceMesh(surfaceMeshFilename);
  TriMeshGeo surfaceMeshGeo;
  surfaceMesh.exportGeometry(surfaceMeshGeo);

  std::vector<DeformationModelPlasticMaterial> testingPlasticDOFs = { DeformationModelPlasticMaterial::VOLUMETRIC_DOF0 };
  if (fiberDirectionFilename && strlen(fiberDirectionFilename)) {
    testingPlasticDOFs.push_back(DeformationModelPlasticMaterial::VOLUMETRIC_DOF3);
  }
  testingPlasticDOFs.push_back(DeformationModelPlasticMaterial::VOLUMETRIC_DOF6);

  std::vector<DeformationModelElasticMaterial> testMaterials;
  testMaterials.push_back(DeformationModelElasticMaterial::STABLE_NEO);
  //testMaterials.push_back(DeformationModelElasticMaterial::STVK);
  //testMaterials.push_back(DeformationModelElasticMaterial::STVK_VOL);
  //testMaterials.push_back(DeformationModelElasticMaterial::VOLUME);

  std::random_device rd;
  std::mt19937 eng(rd());

  int n3 = mesh->getNumVertices() * 3;
  int nele = mesh->getNumElements();

  ES::M3d S;
  S << 1.8673, 0.6765, 1.8324,
    0.6765, 0.4383, 0.9276,
    1.8324, 0.9276, 2.4211;

  ES::MXd basis;
  if (basisFilename && strlen(basisFilename)) {
    std::vector<double> data;
    int mm, nn;
    int r = ReadMatrixFromDisk(basisFilename, mm, nn, data);
    if (r == 0 && mm == nele) {
      basis = Eigen::Map<ES::MXd>(data.data(), mm, nn);
    }
  }

  for (DeformationModelPlasticMaterial nPlastic : testingPlasticDOFs) {
    for (const auto mat : testMaterials) {
      LGI << "\n*******************************************\n"
          << "*******************************************\n"
          << "*******************************************\n"
          << "Testing PM: " << (int)nPlastic << ',' << (int)mat;

      DeformationModelManager *dmm = new DeformationModelManager;

      if (fiberDirectionFilename && strlen(fiberDirectionFilename)) {
        ElementLocalDirection fiberDirection;
        fiberDirection.addMuscle(surfaceMeshGeo.numVertices(), surfaceMeshGeo.positions().data(),
          surfaceMeshGeo.numTriangles(), surfaceMeshGeo.triangles().data());
        fiberDirection.addWrapperTetMesh(tetMesh);
        int ret = fiberDirection.load(fiberDirectionFilename);
        if (ret != 0)
          return 1;

        std::vector<double> elementFiberDirections, vertexFiberDirections;
        for (int ele = 0; ele < mesh->getNumElements(); ele++) {
          Vec3d dir = fiberDirection.getFiberDirection(ele);
          elementFiberDirections.push_back(dir[0]);
          elementFiberDirections.push_back(dir[1]);
          elementFiberDirections.push_back(dir[2]);
        }

        for (int vi = 0; vi < mesh->getNumVertices(); vi++) {
          Vec3d dir = fiberDirection.getVertexFiberDirection(vi);
          vertexFiberDirections.push_back(dir[0]);
          vertexFiberDirections.push_back(dir[1]);
          vertexFiberDirections.push_back(dir[2]);
        }

        dmm->setMesh(mesh, elementFiberDirections.data(), vertexFiberDirections.data());
      }
      else {
        dmm->setMesh(mesh, nullptr, nullptr);
      }

      dmm->init(nPlastic, mat, nullptr);

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

      if ((fdTestType & FD_TT_ELEMENT)) {
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
            LG_ << ele << ' ';

          fem->freeCacheData(cache);
        }

        hclock::time_point t2 = hclock::now();

        LG_ << '\n';
        LGI << '\n'
            << "Error: " << errors[0][0] << ',' << errors[0][1] << '\n'
            << "Error: " << errors[1][0] << ',' << errors[1][1] << '\n'
            << "Error: " << errors[2][0] << ',' << errors[2][1] << '\n'
            << "Time (per tet): " << dur(t1, t2) / (double)numTestElements << "ms";
      }  // if element test

      ES::VXd tetMeshRestPositions(n3);
      for (int i = 0; i < mesh->getNumVertices(); i++) {
        ES::V3d p;
        mesh->getVertex(i, p.data());
        tetMeshRestPositions.segment<3>(i * 3) = p;
      }
      double gradError = 0, hessError = 0;

      if (fdTestType & FD_TT_ASSEM_E) {
        LGI << "Assembler E test...";

        ES::VXd x = tetMeshRestPositions;
        for (int i = 0; i < n3; i++) {
          x[i] = tetMeshRestPositions[i] + distrib(eng);
        }

        int numHandles = (int)basis.cols();
        ES::VXd scales(numHandles * nplastic);

        for (int hi = 0; hi < numHandles; hi++) {
          scales.segment(hi * nplastic, nplastic) = scaleParam;
        }

        std::shared_ptr<DeformationModelAssemblerElastic> forceModelAssemblerE;
        if (basis.size()) {
          forceModelAssemblerE = std::make_shared<DeformationModelAssemblerElastic>(dmm, basis.data(), numHandles);
        }
        else {
          forceModelAssemblerE = std::make_shared<DeformationModelAssemblerElastic>(dmm, nullptr, 0);
        }
        forceModelAssemblerE->setFixedParameters(scales.data(), numHandles * nplastic);

        std::shared_ptr<DeformationModelEnergy> energy = std::make_shared<DeformationModelEnergy>(forceModelAssemblerE);

        fd.testEnergy(energy, true, false, -1.0, x.data(), -1, &gradError, nullptr);
        fd.testEnergy(energy, false, true, -1.0, x.data(), numDOFs, nullptr, &hessError);

        LGI << gradError << ',' << hessError;

        for (int hi = 0; hi < numHandles; hi++) {
          scales.segment(hi * nplastic, nplastic) = scaleParam + ES::VXd::Constant(nplastic, 0.1 * hi);
          dmm->getDeformationModel(0)->getPlasticModel()->compute_paramfull(scales.data() + hi * nplastic);
        }
        forceModelAssemblerE->setFixedParameters(scales.data(), numHandles * nplastic);

        fd.testEnergy(energy, true, false, -1.0, x.data(), -1, &gradError, nullptr);
        fd.testEnergy(energy, false, true, -1.0, x.data(), numDOFs, nullptr, &hessError);

        LGI << gradError << ',' << hessError;
      }

      if ((fdTestType & FD_TT_ASSEM_EFP) && nplastic > 0) {
        LGI << "Assembler EFP test...";
        std::shared_ptr<DeformationModelAssemblerElasticFullPlastic> forceModelAssemblerEFP = std::make_shared<DeformationModelAssemblerElasticFullPlastic>(dmm, nullptr, 0);

        ES::VXd x(forceModelAssemblerEFP->getNumDOFs());
        for (int i = 0; i < n3; i++) {
          x[i] = tetMeshRestPositions[i] + distrib(eng);
        }

        for (int hi = 0; hi < nele; hi++) {
          x.segment(n3 + hi * nplastic, nplastic) = scaleParam;
        }

        std::shared_ptr<DeformationModelEnergy> energy = std::make_shared<DeformationModelEnergy>(forceModelAssemblerEFP);

        fd.testEnergy(energy, true, false, -1.0, x.data(), -1, &gradError, nullptr);
        fd.testEnergy(energy, false, true, -1.0, x.data(), numDOFs, nullptr, &hessError);

        LGI << gradError << ',' << hessError;

        for (int hi = 0; hi < nele; hi++) {
          x.segment(n3 + hi * nplastic, nplastic) = scaleParam + ES::VXd::Constant(nplastic, 0.1 * (hi % 5));
          dmm->getDeformationModel(0)->getPlasticModel()->compute_paramfull(x.data() + n3 + hi * nplastic);
        }

        fd.testEnergy(energy, true, false, -1.0, x.data(), -1, &gradError, nullptr);
        fd.testEnergy(energy, false, true, -1.0, x.data(), numDOFs, nullptr, &hessError);

        LGI << gradError << ',' << hessError;
      }

      delete dmm;
      std::cin.get();
    }
  }

  return 0;
}

int SolidDeformationModel::fdTestBending(const char *surfaceMeshFilename, const char *basisFilename, const char *fiberDirectionFilename, int numDOFs, int testType)
{
  ObjMesh mesh(surfaceMeshFilename);
  SimulationMeshENuhMaterial mat(1e8, 0.4, 3e-4);
  SimulationMesh *sm = loadObjMesh(&mesh, &mat, 0);
  DeformationModelManager *dmm = new DeformationModelManager;
  dmm->setMesh(sm, nullptr, nullptr);
  dmm->init(DeformationModelPlasticMaterial::VOLUMETRIC_DOF0, DeformationModelElasticMaterial::CLOTH_BENDING_QUADRATIC, nullptr);

  std::default_random_engine rd;
  //std::mt19937 eng(rd());

  hclock::time_point t1 = hclock::now();

  int numTestElements = sm->getNumElements();

  FiniteDifference fd(FiniteDifference::M_FIVE_POINT, 1e-7);
  std::uniform_real_distribution<double> distrib(-0.01, 0.01);

  int nelastic = 1;
  ES::VXd x(12 + nelastic);
  int offset[3] = {
    0,
    12,
    12 + nelastic,
  };

  double errors[3][2] = { { 0, 0 }, { 0, 0 }, { 0, 0 } };
  tbb::enumerable_thread_specific<ES::MXd> d2EdbdbTLS(ES::MXd(nelastic, nelastic));
  tbb::enumerable_thread_specific<ES::MXd> d2EdxdbTLS(ES::MXd(12, nelastic));

  for (int ele = 0; ele < numTestElements; ele++) {
    const DeformationModel *fem = dmm->getDeformationModel(ele);
    DeformationModelCacheData *cache = fem->allocateCacheData();

    for (int j = 0; j < 4; j++) {
      ES::V3d p;
      sm->getVertex(ele, j, p.data());
      x.segment<3>(j * 3) = p;
    }

    for (int j = 0; j < 12; j++) {
      x[j] += distrib(rd);
    }
    x[12] = 1;

    //std::vector<Vec3d> pt;
    //std::vector<Vec3i> tri;
    //for (int i = 0; i < 4; i++) {
    //  pt.emplace_back(Vec3d(x.data() + i * 3));
    //}
    //tri.push_back(Vec3i(0, 1, 2));
    //tri.push_back(Vec3i(3, 2, 1));
    //ObjMesh(pt, tri).save("f:/abc.obj");

    double elasticCoeff[1] = { 1.0 };

    FiniteDifference::EvalFunc evalFuncE = [&](const double *x_param, double *E, double *grad, double *hess) {
      fem->prepareData(x_param, nullptr, elasticCoeff, cache);

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

    if (nelastic > 0) {
      FiniteDifference::EvalFunc evalFuncEA = [&](const double *x_param, double *E, double *grad, double *hess) {
        fem->prepareData(x_param, nullptr, x_param + offset[1], cache);

        if (E) {
          *E = fem->computeEnergy(cache);
        }

        if (grad) {
          memset(grad, 0, sizeof(double) * offset[2]);
          fem->compute_dE_dx(cache, grad);
          fem->compute_dE_db(cache, grad + offset[1]);
        }

        if (hess) {
          ES::M12d hess0;
          ES::MXd &d2E_dxdb = d2EdxdbTLS.local();
          ES::MXd &d2E_db2 = d2EdbdbTLS.local();

          fem->compute_d2E_dx2(cache, hess0.data());
          fem->compute_d2E_dxdb(cache, d2E_dxdb.data());
          fem->compute_d2E_db2(cache, d2E_db2.data());

          Eigen::Map<ES::MXd> hMap(hess, offset[2], offset[2]);
          hMap.block<12, 12>(0, 0) = hess0;
          hMap.block(0, 12, 12, nelastic) = d2E_dxdb;
          hMap.block(12, 0, nelastic, 12) = d2E_dxdb.transpose();
          hMap.block(12, 12, nelastic, nelastic) = d2E_db2;
        }
      };

      fd.testEnergy(evalFuncEA, offset[2], true, true, -1, x.data(), &gradError, &hessError);
      // LGI << "Ele " << ele << "; Error: " << gradError << ',' << hessError;
      errors[2][0] = std::max(errors[2][0], gradError);
      errors[2][1] = std::max(errors[2][1], hessError);
    }

    if (ele % 10 == 0)
      LG_ << ele << ' ';

    fem->freeCacheData(cache);
  }

  hclock::time_point t2 = hclock::now();

  LG_ << '\n';
  LGI << '\n'
      << "Error: " << errors[0][0] << ',' << errors[0][1] << '\n'
      << "Error: " << errors[1][0] << ',' << errors[1][1] << '\n'
      << "Error: " << errors[2][0] << ',' << errors[2][1] << '\n'
      << "Time (per tet): " << dur(t1, t2) / (double)numTestElements << "ms";

  return 0;
}

/*
forceModelAssemblerEMP->setFixedParameters(&actParam, 1);

x.resize(forceModelAssemblerEMP->getNumDOFs());
for (int i = 0; i < n3; i++) {
  x[i] = tetMeshRestPositions[i] + distrib(eng);
}

LGI << "Assembler EMP test...";

double gradError, hessError;
int numDOFs = 300;
std::shared_ptr<MuscleElasticEnergy> muscleEnergy = std::make_shared<MuscleElasticEnergy>(forceModelAssemblerEMP);

std::chrono::high_resolution_clock::time_point t3 = std::chrono::high_resolution_clock::now();

if (ASSM_EMP_TEST) {
  fd.testEnergy(muscleEnergy, true, false, -1.0, x.data(), -1, &gradError, nullptr);
  fd.testEnergy(muscleEnergy, false, true, -1.0, x.data(), numDOFs, nullptr, &hessError);
}

std::chrono::high_resolution_clock::time_point t4 = std::chrono::high_resolution_clock::now();
LGI << gradError << ',' << hessError << "; time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t4 - t3).count() << "ms";

for (int hi = 0; hi < numHandles; hi++) {
  x.segment(n3 + hi * nplastic, nplastic) = scaleParam + ES::VXd::Constant(nplastic, 0.1 * hi);
  muscleSimModel->getMuscleFEM(materialControlElements[hi])->getPlasticModel()->compute_paramfull(x.data() + n3 + hi * nplastic);
}

if (ASSM_EMP_TEST) {
  fd.testEnergy(muscleEnergy, true, false, -1.0, x.data(), -1, &gradError, nullptr);
  fd.testEnergy(muscleEnergy, false, true, -1.0, x.data(), numDOFs, nullptr, &hessError);
}

LGI << gradError << ',' << hessError;

*/
