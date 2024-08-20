/*
author: Bohan Wang
copyright to USC,MIT
*/

#include "deformationModelManager.h"
#include "deformationModel.h"
#include "tetMeshDeformationModel.h"
// #include "triangleMeshDeformationModel.h"
// #include "edgeDeformationModel.h"
// #include "koiterShellDeformationModel.h"
#include "simulationMesh.h"

#include "elasticModelCombinedMaterial.h"
#include "elasticModelHillTypeMaterial.h"
#include "elasticModelInvariantBasedMaterial.h"
// #include "elasticModelLambdaBasedMaterial.h"
// #include "elasticModelLowerInvariantBasedMaterial.h"
// #include "lowerInvariantBasedMaterialSpline.h"
#include "elasticModelStableNeoHookeanMaterial.h"
#include "elasticModelVolumeMaterial.h"
#include "invariantBasedMaterialStVK.h"
#include "elasticModel1DQuadratic.h"
#include "elasticModel1DCubicSpline.h"
#include "elasticModelLinearMaterial.h"
// #include "elasticModelTensileMaterial.h"
// #include "elasticModelKoiterShellMaterial.h"

// #include "plasticModel2DConstant.h"
// #include "plasticModel2D3DOF.h"
// #include "plasticModelBend1DOF.h"
// #include "plasticModelBendConstant.h"
// #include "plasticModelKoiterShell6DOF.h"
#include "plasticModel3D3DOF.h"
#include "plasticModel3D6DOF.h"
#include "plasticModel3DConstant.h"

#include "configFileJSON.h"
#include "pgoLogging.h"
#include "volumetricMesh.h"
#include "EigenSupport.h"

#include "volumetricMeshENuMaterial.h"

#include <fmt/format.h>

#include <tbb/parallel_for.h>
#include <tbb/enumerable_thread_specific.h>

#include <memory>

namespace ES = pgo::EigenSupport;

namespace pgo
{
namespace SolidDeformationModel
{
class DeformationModelManagerImpl
{
public:
  ~DeformationModelManagerImpl();

  const SimulationMesh *simulationMesh;

  std::vector<DeformationModel *> elementFEMs;

  std::vector<ElasticModelStableNeoHookeanMaterial *> stableNeoHookeanMaterials;
  std::vector<ElasticModelLinearMaterial *> linearMaterials;
  std::vector<ElasticModelHillTypeMaterial *> hillTypeMaterials;
  std::vector<ElasticModelInvariantBasedMaterial *> invariantBasedMaterials;
  std::vector<ElasticModelVolumeMaterial *> volumeMaterials;
  std::vector<ElasticModel1DQuadratic *> quadratic1D;
  std::vector<ElasticModel1DCubicSpline *> spline1D;
  std::vector<ElasticModelCombinedMaterial<2> *> combined2Materials;
  std::vector<ElasticModelCombinedMaterial<3> *> combined3Materials;

  // std::vector<ElasticModelTensileMaterial *> clothStVKMaterials;
  // std::vector<ElasticModelKoiterShellMaterial *> koiterShellMaterials;
  // std::vector<ElasticModelLowerInvariantBasedMaterial *> lowerInvariantMaterials;

  std::vector<ElasticModel *> elementMaterials;
  std::vector<InvariantBasedMaterial *> invariantModels;

  std::vector<PlasticModel3DConstant *> plasticVolConstant;
  std::vector<PlasticModel3D3DOF *> plasticVol3DOF;
  std::vector<PlasticModel3D6DOF *> plasticVol6DOF;

  // std::vector<LowerInvariantBasedMaterialSpline *> lowerInvariantSplineModels;
  // std::vector<PlasticModel2DConstant *> plasticSurfaceConstant;
  // std::vector<PlasticModel2D3DOF *> plasticSurface3DOF;
  // std::vector<PlasticModelBendConstant *> plasticBendConstant;
  // std::vector<PlasticModelBend1DOF *> plasticBend1DOF;
  // std::vector<PlasticModelKoiterShell6DOF *> plasticKoiterShell6DOF;

  ES::VXd fiberDirections;
  ES::VXd vertexFiberDirections;
  ES::M3Xd fiberAxesRest, vertexFiberAxesRest;
  ES::M3Xd fiberAxes, vertexFiberAxes;

  ES::M3d globalRotation;

  double Eact = 0.1e6;
  double Etendon = 6000;
  double Epass = 6000;
  double lo = 0.6;
  double gamma = 1.0;
  double nu = 0.49;
  double h = 3e-4;

  int enforceSPD = 1;
  int hasUserDefinedMaterial = 0;
  double compressionRatio = 0;

  int numPlasticParams;
  int nele;
  int nvtx;

  double *x = nullptr, *y1 = nullptr, *y2 = nullptr, *y3 = nullptr;
  int nSplinePt = 0;
};

DeformationModelManagerImpl::~DeformationModelManagerImpl()
{
  for (auto ptr : elementFEMs)
    delete ptr;

  for (auto ptr : stableNeoHookeanMaterials)
    if (ptr)
      delete ptr;

  for (auto ptr : linearMaterials)
    if (ptr)
      delete ptr;

  for (auto ptr : hillTypeMaterials)
    if (ptr)
      delete ptr;

  for (auto ptr : invariantBasedMaterials)
    if (ptr)
      delete ptr;

  // for (auto ptr : clothStVKMaterials)
  //   if (ptr)
  //     delete ptr;

  // for (auto ptr : koiterShellMaterials)
  //   if (ptr)
  //     delete ptr;

  for (auto ptr : volumeMaterials)
    if (ptr)
      delete ptr;

  for (auto ptr : quadratic1D)
    if (ptr)
      delete ptr;

  for (auto ptr : spline1D)
    if (ptr)
      delete ptr;

  for (auto ptr : combined2Materials)
    if (ptr)
      delete ptr;

  for (auto ptr : plasticVol3DOF)
    if (ptr)
      delete ptr;

  for (auto ptr : plasticVol6DOF)
    if (ptr)
      delete ptr;

  for (auto ptr : plasticVolConstant)
    if (ptr)
      delete ptr;

  // for (auto ptr : plasticSurface3DOF)
  //   if (ptr)
  //     delete ptr;

  // for (auto ptr : plasticKoiterShell6DOF)
  //   if (ptr)
  //     delete ptr;

  // for (auto ptr : plasticSurfaceConstant)
  //   if (ptr)
  //     delete ptr;
  // for (auto ptr : plasticBendConstant)
  //   if (ptr)
  //     delete ptr;

  // for (auto ptr : plasticBend1DOF)
  //   if (ptr)
  //     delete ptr;

  // for (auto ptr : lowerInvariantMaterials)
  //   if (ptr)
  //     delete ptr;

  // for (auto ptr : lowerInvariantSplineModels)
  //   if (ptr)
  //     delete ptr;
}

int numPlasticDOFsAry[] = {
  0, 3, 6, 0, 3, 0, 1, 6
};

}  // namespace SolidDeformationModel
}  // namespace pgo

using namespace pgo::SolidDeformationModel;

DeformationModelManager::DeformationModelManager()
{
  data = new DeformationModelManagerImpl;
}

DeformationModelManager::~DeformationModelManager()
{
  delete data;
}

void DeformationModelManager::setMesh(const SimulationMesh *simulationMesh, const double *elementFiberDirections, const double *vertexFiberDirections)
{
  data->simulationMesh = simulationMesh;
  data->nele = data->simulationMesh->getNumElements();
  data->nvtx = data->simulationMesh->getNumVertices();

  if (elementFiberDirections) {
    data->fiberDirections = Eigen::Map<const ES::VXd>(elementFiberDirections, data->nele * 3);
  }
  else {
    data->fiberDirections = ES::VXd::Zero(data->nele * 3);
    for (int i = 0; i < data->nele; i++) {
      data->fiberDirections.segment<3>(i * 3) = ES::V3d::UnitY();
    }
  }

  if (vertexFiberDirections) {
    data->vertexFiberDirections = Eigen::Map<const ES::VXd>(vertexFiberDirections, data->nvtx * 3);
  }
  else {
    data->vertexFiberDirections = ES::VXd::Zero(data->nvtx * 3);
    for (int i = 0; i < data->nvtx; i++) {
      data->vertexFiberDirections.segment<3>(i * 3) = ES::V3d::UnitY();
    }
  }
}

void DeformationModelManager::setSplineFunction(int numPts, double *x, double *y1, double *y2, double *y3)
{
  data->x = x;
  data->y1 = y1;
  data->y2 = y2;
  data->y3 = y3;
  data->nSplinePt = numPts;
}

void DeformationModelManager::init(DeformationModelPlasticMaterial plasticModelType, DeformationModelElasticMaterial elasticMaterialType, const char *muscleParameterConfigFilename /* = nullptr */)
{
  SPDLOG_LOGGER_INFO(pgo::Logging::lgr(), "Computing the type of each element...");

  data->numPlasticParams = numPlasticDOFsAry[(int)plasticModelType];
  data->globalRotation = ES::M3d::Identity();

  SPDLOG_LOGGER_INFO(pgo::Logging::lgr(), "Computing the fiber alignment transformation...");

  std::vector<std::vector<int>> vertexNearbyElements(data->nvtx);
  std::vector<std::vector<int>> elementNearbyElements(data->nele);

  for (int i = 0; i < data->nele; i++) {
    for (int j = 0; j < data->simulationMesh->getNumElementVertices(); j++) {
      vertexNearbyElements[data->simulationMesh->getVertexIndex(i, j)].push_back(i);
    }
  }

  for (int i = 0; i < data->nele; i++) {
    for (int j = 0; j < data->simulationMesh->getNumElementVertices(); j++) {
      const auto &eles = vertexNearbyElements[data->simulationMesh->getVertexIndex(i, j)];
      elementNearbyElements[i].insert(elementNearbyElements[i].end(), eles.begin(), eles.end());
    }

    // sort and remove duplications
    std::sort(elementNearbyElements[i].begin(), elementNearbyElements[i].end());
    auto itt = std::unique(elementNearbyElements[i].begin(), elementNearbyElements[i].end());
    elementNearbyElements[i].erase(itt, elementNearbyElements[i].end());

    // remove itself
    itt = std::lower_bound(elementNearbyElements[i].begin(), elementNearbyElements[i].end(), i);
    PGO_ALOG(itt != elementNearbyElements[i].end() && *itt == i);
    elementNearbyElements[i].erase(itt);
  }

  data->fiberAxesRest.resize(3, data->nele * 3);
  ES::V3d guideVector(0, 0, 1);
  std::vector<int> axisComputed(data->nele, 0);
  while (1) {
    int allComputed = 1;
    for (int f : axisComputed)
      allComputed &= f;

    if (allComputed)
      break;

    for (int ele = 0; ele < data->nele; ele++) {
      if (axisComputed[ele] == 1)
        continue;

      ES::V3d hintY = guideVector;
      ES::V3d X = data->fiberDirections.segment<3>(ele * 3);
      X.normalize();

      // if it is the degenerated case
      if (fabs(X.dot(hintY)) > 1 - 1e-6) {
        // we first search nearby elements
        bool nearbyFound = false;
        for (int elej : elementNearbyElements[ele]) {
          if (axisComputed[elej] == 0) {
            continue;
          }

          // fetch the computed direction
          ES::M3d R = data->fiberAxesRest.block<3, 3>(0, elej * 3);
          hintY = R.row(1);
          nearbyFound = true;
          break;
        }

        if (nearbyFound == false) {
          continue;
        }
      }

      ES::V3d Z = X.cross(hintY);
      Z.normalize();

      ES::V3d Y = Z.cross(X);
      Y.normalize();

      ES::M3d R;
      R.row(0) = X;
      R.row(1) = Y;
      R.row(2) = Z;

      data->fiberAxesRest.block<3, 3>(0, ele * 3) = R;
      axisComputed[ele] = 1;
    }

    int numAxes = std::count_if(axisComputed.begin(), axisComputed.end(), [](int val) { return val == 1; });
    SPDLOG_LOGGER_INFO(pgo::Logging::lgr(), "# axes computed:{}", numAxes);
  }

  data->vertexFiberAxesRest.resize(3, data->nvtx * 3);
  for (int vi = 0; vi < data->nvtx; vi++) {
    ES::V3d hintY = data->fiberAxesRest.col(vertexNearbyElements[vi][0] * 3 + 1);
    ES::V3d X = data->vertexFiberDirections.segment<3>(vi * 3);
    X.normalize();

    ES::V3d Z = X.cross(hintY);
    Z.normalize();

    ES::V3d Y = Z.cross(X);
    Y.normalize();

    ES::M3d R;
    R.row(0) = X;
    R.row(1) = Y;
    R.row(2) = Z;

    data->vertexFiberAxesRest.block<3, 3>(0, vi * 3) = R;
  }

  data->fiberAxes = data->fiberAxesRest;
  data->vertexFiberAxes = data->vertexFiberAxesRest;

  SPDLOG_LOGGER_INFO(pgo::Logging::lgr(), "Initializing element material model...");

  if (elasticMaterialType == DeformationModelElasticMaterial::HILL_STABLE_NEO ||
    elasticMaterialType == DeformationModelElasticMaterial::HILL_STVK ||
    elasticMaterialType == DeformationModelElasticMaterial::HILL_STVK_VOL) {
    data->hillTypeMaterials.assign(data->nele, nullptr);
  }

  if (elasticMaterialType == DeformationModelElasticMaterial::HILL_STABLE_NEO ||
    elasticMaterialType == DeformationModelElasticMaterial::STABLE_NEO) {
    data->stableNeoHookeanMaterials.assign(data->nele, nullptr);
  }

  if (elasticMaterialType == DeformationModelElasticMaterial::HILL_STVK ||
    elasticMaterialType == DeformationModelElasticMaterial::HILL_STVK_VOL ||
    elasticMaterialType == DeformationModelElasticMaterial::STVK ||
    elasticMaterialType == DeformationModelElasticMaterial::STVK_VOL) {
    data->invariantModels.assign(data->nele, nullptr);
    data->invariantBasedMaterials.assign(data->nele, nullptr);
  }

  if (elasticMaterialType == DeformationModelElasticMaterial::HILL_STVK_VOL ||
    elasticMaterialType == DeformationModelElasticMaterial::VOLUME ||
    elasticMaterialType == DeformationModelElasticMaterial::STVK_VOL) {
    data->volumeMaterials.assign(data->nele, nullptr);
  }

  if (elasticMaterialType == DeformationModelElasticMaterial::HILL_STABLE_NEO ||
    elasticMaterialType == DeformationModelElasticMaterial::HILL_STVK ||
    elasticMaterialType == DeformationModelElasticMaterial::STVK_VOL) {
    data->combined2Materials.assign(data->nele, nullptr);
  }

  if (elasticMaterialType == DeformationModelElasticMaterial::HILL_STVK_VOL) {
    data->combined3Materials.assign(data->nele, nullptr);
  }

  // if (elasticMaterialType == DeformationModelElasticMaterial::CLOTH_STVK) {
  //   data->clothStVKMaterials.assign(data->nele, nullptr);
  // }

  // if (elasticMaterialType == DeformationModelElasticMaterial::KOITER_SHELL) {
  //   data->koiterShellMaterials.assign(data->nele, nullptr);
  // }

  // if (elasticMaterialType == DeformationModelElasticMaterial::CLOTH_BENDING_QUADRATIC) {
  //   data->quadratic1D.assign(data->nele, nullptr);
  // }

  // if (elasticMaterialType == DeformationModelElasticMaterial::CLOTH_BENDING_CUBIC_SPLINE) {
  //   data->spline1D.assign(data->nele, nullptr);
  // }

  // if (elasticMaterialType == DeformationModelElasticMaterial::LOWER_INVARIENT_CUBIC_SPLINE) {
  //   data->lowerInvariantMaterials.assign(data->nele, nullptr);
  //   data->lowerInvariantSplineModels.assign(data->nele, nullptr);
  // }

  if (elasticMaterialType == DeformationModelElasticMaterial::LINEAR) {
    data->linearMaterials.assign(data->nele, nullptr);
  }

  data->elementFEMs.assign(data->nele, nullptr);
  data->elementMaterials.assign(data->nele, nullptr);

  if (plasticModelType == DeformationModelPlasticMaterial::VOLUMETRIC_DOF0) {
    data->plasticVolConstant.assign(data->nele, nullptr);
  }
  else if (plasticModelType == DeformationModelPlasticMaterial::VOLUMETRIC_DOF3) {
    data->plasticVol3DOF.assign(data->nele, nullptr);
  }
  else if (plasticModelType == DeformationModelPlasticMaterial::VOLUMETRIC_DOF6) {
    data->plasticVol6DOF.assign(data->nele, nullptr);
  }
  // else if (plasticModelType == DeformationModelPlasticMaterial::SURFACE_DOF3) {
  //   data->plasticSurface3DOF.assign(data->nele, nullptr);
  // }
  // else if (plasticModelType == DeformationModelPlasticMaterial::SURFACE_DOF0) {
  //   data->plasticSurfaceConstant.assign(data->nele, nullptr);
  // }
  // else if (plasticModelType == DeformationModelPlasticMaterial::BEND_DOF0) {
  //   data->plasticBendConstant.assign(data->nele, nullptr);
  // }
  // else if (plasticModelType == DeformationModelPlasticMaterial::BEND_DOF1) {
  //   data->plasticBend1DOF.assign(data->nele, nullptr);
  // }
  // else if (plasticModelType == DeformationModelPlasticMaterial::SHELL_DOF6) {
  //   data->plasticKoiterShell6DOF.assign(data->nele, nullptr);
  // }

  if (muscleParameterConfigFilename && strlen(muscleParameterConfigFilename)) {
    ConfigFileJSON jconfig;
    if (jconfig.open(muscleParameterConfigFilename)) {
      if (jconfig.exist("E-active"))
        data->Eact = jconfig.getDouble("E-active");

      if (jconfig.exist("E-passive"))
        data->Epass = jconfig.getDouble("E-passive");

      if (jconfig.exist("E-tendon"))
        data->Etendon = jconfig.getDouble("E-tendon");

      if (jconfig.exist("nu"))
        data->nu = jconfig.getDouble("nu");

      if (jconfig.exist("gamma"))
        data->gamma = jconfig.getDouble("gamma");

      if (jconfig.exist("lo"))
        data->lo = jconfig.getDouble("lo");

      if (jconfig.exist("compression"))
        data->compressionRatio = jconfig.getDouble("compression");

      if (jconfig.exist("h"))
        data->h = jconfig.getDouble("h");

      if (jconfig.exist("spd"))
        data->enforceSPD = jconfig.getInt("spd");

      if (jconfig.exist("use")) {
        data->hasUserDefinedMaterial = jconfig.getInt("use");
      }
    }
  }

  std::stringstream ss;
  ss << "\n"
     << "E active : " << data->Eact << '\n'
     << "E passive: " << data->Epass << '\n'
     << "E tendon : " << data->Etendon << '\n'
     << "nu       : " << data->nu << '\n'
     << "gamma    : " << data->gamma << '\n'
     << "lo       : " << data->lo << '\n'
     << "h        : " << data->h << '\n'
     << "SPD      : " << data->enforceSPD << '\n'
     << "compress : " << data->compressionRatio << '\n'
     << "Use json : " << data->hasUserDefinedMaterial;

  SPDLOG_LOGGER_INFO(pgo::Logging::lgr(), ss.str());

  // for (int ele = 0; ele < data->nele; ele++) {
  tbb::parallel_for(
    0, data->nele, [&](int ele) {
      double mu = -1, lambda = -1, E = -1, nu = -1, h = -1;
      if (data->hasUserDefinedMaterial == 0) {
        const SimulationMeshENuMaterial *matBase;
        if (const SimulationMeshENuhMaterial *mat = dynamic_cast<const SimulationMeshENuhMaterial *>(data->simulationMesh->getElementMaterial(ele))) {
          matBase = mat;
          h = mat->geth();
        }
        else if (const SimulationMeshENuMaterial *mat = dynamic_cast<const SimulationMeshENuMaterial *>(data->simulationMesh->getElementMaterial(ele))) {
          matBase = mat;
        }
        else {
          SPDLOG_LOGGER_ERROR(pgo::Logging::lgr(), "Unsupported material.");
          throw std::logic_error("Unsupported material");
        }

        mu = matBase->getMuLame();
        lambda = matBase->getLambdaLame();
        E = matBase->getE();
        nu = matBase->getNu();
      }
      else {
        SimulationMeshENuhMaterial mat(data->Epass, data->nu, data->h);
        mu = mat.getMuLame();
        lambda = mat.getLambdaLame();
        E = data->Epass;
        nu = data->nu;
        h = data->h;
      }

      double coeffJ = 1;
      if (data->compressionRatio > 0)
        coeffJ = data->compressionRatio * E / (1.0 - 2.0 * nu);
      else
        coeffJ = E / (1.0 - 2.0 * nu);

      ES::V3d dir = data->fiberAxesRest.block<3, 3>(0, ele * 3).row(0);
      if (elasticMaterialType == DeformationModelElasticMaterial::HILL_STABLE_NEO) {
        data->stableNeoHookeanMaterials[ele] = new ElasticModelStableNeoHookeanMaterial(mu, lambda);
        data->stableNeoHookeanMaterials[ele]->enforceSPD(data->enforceSPD);

        data->hillTypeMaterials[ele] = new ElasticModelHillTypeMaterial(data->gamma, data->Eact, data->lo, dir.data());
        data->combined2Materials[ele] = new ElasticModelCombinedMaterial<2>(data->stableNeoHookeanMaterials[ele], data->hillTypeMaterials[ele]);

        data->elementMaterials[ele] = data->combined2Materials[ele];
      }
      else if (elasticMaterialType == DeformationModelElasticMaterial::LINEAR) {
        data->linearMaterials[ele] = new ElasticModelLinearMaterial(mu, lambda);
        data->elementMaterials[ele] = data->linearMaterials[ele];
      }
      else if (elasticMaterialType == DeformationModelElasticMaterial::HILL_STVK) {
        data->invariantModels[ele] = new InvariantBasedMaterialStVK(E, nu, data->compressionRatio);

        data->invariantBasedMaterials[ele] = new ElasticModelInvariantBasedMaterial(data->invariantModels[ele]);
        data->invariantBasedMaterials[ele]->enforceSPD(data->enforceSPD);

        data->hillTypeMaterials[ele] = new ElasticModelHillTypeMaterial(data->gamma, data->Eact, data->lo, dir.data());
        data->combined2Materials[ele] = new ElasticModelCombinedMaterial<2>(data->invariantBasedMaterials[ele], data->hillTypeMaterials[ele]);

        data->elementMaterials[ele] = data->combined2Materials[ele];
      }
      else if (elasticMaterialType == DeformationModelElasticMaterial::HILL_STVK_VOL) {
        data->invariantModels[ele] = new InvariantBasedMaterialStVK(E, nu, data->compressionRatio);

        data->invariantBasedMaterials[ele] = new ElasticModelInvariantBasedMaterial(data->invariantModels[ele]);
        data->invariantBasedMaterials[ele]->enforceSPD(data->enforceSPD);

        data->hillTypeMaterials[ele] = new ElasticModelHillTypeMaterial(data->gamma, data->Eact, data->lo, dir.data());

        data->volumeMaterials[ele] = new ElasticModelVolumeMaterial(coeffJ);

        data->combined3Materials[ele] = new ElasticModelCombinedMaterial<3>(data->invariantBasedMaterials[ele], data->hillTypeMaterials[ele], data->volumeMaterials[ele]);

        data->elementMaterials[ele] = data->combined3Materials[ele];
      }
      else if (elasticMaterialType == DeformationModelElasticMaterial::STABLE_NEO) {
        data->stableNeoHookeanMaterials[ele] = new ElasticModelStableNeoHookeanMaterial(mu, lambda);
        data->stableNeoHookeanMaterials[ele]->enforceSPD(data->enforceSPD);
        // data->stableNeoHookeanMaterials[ele]->enforceSPD(true);

        data->elementMaterials[ele] = data->stableNeoHookeanMaterials[ele];
      }
      else if (elasticMaterialType == DeformationModelElasticMaterial::STVK) {
        data->invariantModels[ele] = new InvariantBasedMaterialStVK(E, nu, data->compressionRatio);

        data->invariantBasedMaterials[ele] = new ElasticModelInvariantBasedMaterial(data->invariantModels[ele]);
        data->invariantBasedMaterials[ele]->enforceSPD(data->enforceSPD);

        data->elementMaterials[ele] = data->invariantBasedMaterials[ele];
      }
      else if (elasticMaterialType == DeformationModelElasticMaterial::STVK_VOL) {
        data->invariantModels[ele] = new InvariantBasedMaterialStVK(E, nu, data->compressionRatio);

        data->invariantBasedMaterials[ele] = new ElasticModelInvariantBasedMaterial(data->invariantModels[ele]);
        data->invariantBasedMaterials[ele]->enforceSPD(data->enforceSPD);

        data->volumeMaterials[ele] = new ElasticModelVolumeMaterial(coeffJ);

        data->combined2Materials[ele] = new ElasticModelCombinedMaterial<2>(data->invariantBasedMaterials[ele], data->volumeMaterials[ele]);

        data->elementMaterials[ele] = data->combined2Materials[ele];
      }
      else if (elasticMaterialType == DeformationModelElasticMaterial::VOLUME) {
        data->volumeMaterials[ele] = new ElasticModelVolumeMaterial(coeffJ);

        data->elementMaterials[ele] = data->volumeMaterials[ele];
      }
      // else if (elasticMaterialType == DeformationModelElasticMaterial::CLOTH_STVK) {
      //   data->clothStVKMaterials[ele] = new ElasticModelTensileMaterial(E, nu);
      //   data->elementMaterials[ele] = data->clothStVKMaterials[ele];
      // }
      // else if (elasticMaterialType == DeformationModelElasticMaterial::KOITER_SHELL) {
      //   data->koiterShellMaterials[ele] = new ElasticModelKoiterShellMaterial(E, nu);
      //   data->elementMaterials[ele] = data->koiterShellMaterials[ele];
      // }
      // else if (elasticMaterialType == DeformationModelElasticMaterial::CLOTH_BENDING_QUADRATIC) {
      //   PGO_ALOG(h >= 0);
      //   double coeff = E * h * h * h / 24.0 / (1.0 - nu * nu);
      //   data->quadratic1D[ele] = new ElasticModel1DQuadratic(coeff);
      //   data->elementMaterials[ele] = data->quadratic1D[ele];
      // }
      // else if (elasticMaterialType == DeformationModelElasticMaterial::CLOTH_BENDING_CUBIC_SPLINE) {
      //   PGO_ALOG(h >= 0);
      //   double coeff = E * h * h * h / 24.0 / (1.0 - nu * nu);
      //   data->spline1D[ele] = new ElasticModel1DCubicSpline(coeff, 5, -M_PI, M_PI);
      //   data->elementMaterials[ele] = data->spline1D[ele];
      // }
      // else if (elasticMaterialType == DeformationModelElasticMaterial::LOWER_INVARIENT_CUBIC_SPLINE) {
      //   data->lowerInvariantSplineModels[ele] = new LowerInvariantBasedMaterialSpline(data->nSplinePt, data->x, data->y1, data->y2, data->y3);
      //   data->lowerInvariantMaterials[ele] = new ElasticModelLowerInvariantBasedMaterial(data->lowerInvariantSplineModels[ele]);
      //   data->elementMaterials[ele] = data->lowerInvariantMaterials[ele];
      // }
      else {
        throw std::runtime_error("unknown elastic model");
      }

      // Plastic Model Type
      PlasticModel *pm;
      if (plasticModelType == DeformationModelPlasticMaterial::VOLUMETRIC_DOF0) {
        ES::M3d I = ES::M3d::Identity();
        data->plasticVolConstant[ele] = new PlasticModel3DConstant(I.data());
        pm = data->plasticVolConstant[ele];
      }
      else if (plasticModelType == DeformationModelPlasticMaterial::VOLUMETRIC_DOF3) {
        data->plasticVol3DOF[ele] = new PlasticModel3D3DOF(data->fiberAxesRest.data() + ele * 9);
        pm = data->plasticVol3DOF[ele];
      }
      else if (plasticModelType == DeformationModelPlasticMaterial::VOLUMETRIC_DOF6) {
        data->plasticVol6DOF[ele] = new PlasticModel3D6DOF();
        pm = data->plasticVol6DOF[ele];
      }
      // else if (plasticModelType == DeformationModelPlasticMaterial::SURFACE_DOF0) {
      //   ES::M2d I = ES::M2d::Identity();
      //   data->plasticSurfaceConstant[ele] = new PlasticModel2DConstant(I.data());
      //   pm = data->plasticSurfaceConstant[ele];
      // }
      // else if (plasticModelType == DeformationModelPlasticMaterial::SURFACE_DOF3) {
      //   ES::M2d I = ES::M2d::Identity();
      //   data->plasticSurface3DOF[ele] = new PlasticModel2D3DOF(I.data());
      //   pm = data->plasticSurface3DOF[ele];
      // }
      // else if (plasticModelType == DeformationModelPlasticMaterial::SHELL_DOF6) {
      //   data->plasticKoiterShell6DOF[ele] = new PlasticModelKoiterShell6DOF();
      //   pm = data->plasticKoiterShell6DOF[ele];
      // }
      // else if (plasticModelType == DeformationModelPlasticMaterial::BEND_DOF0) {
      //   data->plasticBendConstant[ele] = new PlasticModelBendConstant(0);
      //   pm = data->plasticBendConstant[ele];
      // }
      // else if (plasticModelType == DeformationModelPlasticMaterial::BEND_DOF1) {
      //   data->plasticBend1DOF[ele] = new PlasticModelBend1DOF(0);
      //   pm = data->plasticBend1DOF[ele];
      // }
      else {
        throw std::runtime_error("unknown plastic model");
      }

      // Mesh Element Type
      if (data->simulationMesh->getElementType() == SimulationMeshType::TET) {
        ES::V12d restPosition;
        for (int j = 0; j < 4; j++) {
          ES::V3d p;
          data->simulationMesh->getVertex(ele, j, p.data());
          restPosition.segment<3>(j * 3) = p;
        }

        data->elementFEMs[ele] = new TetMeshDeformationModel(
          restPosition.data(), restPosition.data() + 3, restPosition.data() + 6, restPosition.data() + 9,
          data->elementMaterials[ele], pm);
      }
      // else if (data->simulationMesh->getElementType() == SimulationMeshType::TRIANGLE) {
      //   ES::V9d restPosition;
      //   for (int j = 0; j < 3; j++) {
      //     ES::V3d p;
      //     data->simulationMesh->getVertex(ele, j, p.data());
      //     restPosition.segment<3>(j * 3) = p;
      //   }

      //   if (data->simulationMesh->hasElementUV()) {
      //     ES::V2d restUVs[3];
      //     for (int j = 0; j < 3; j++) {
      //       data->simulationMesh->getElementUV(ele, j, restUVs[j].data());
      //     }

      //     data->elementFEMs[ele] = new TriangleMeshDeformationModel(
      //       restPosition.data(), restPosition.data() + 3, restPosition.data() + 6,
      //       data->elementMaterials[ele], pm,
      //       restUVs[0].data(), restUVs[1].data(), restUVs[2].data());
      //   }
      //   else {
      //     data->elementFEMs[ele] = new TriangleMeshDeformationModel(
      //       restPosition.data(), restPosition.data() + 3, restPosition.data() + 6,
      //       data->elementMaterials[ele], pm);
      //   }
      // }
      // else if (data->simulationMesh->getElementType() == SimulationMeshType::SHELL) {
      //   ES::V18d restPosition;
      //   for (int j = 0; j < 6; j++) {
      //     ES::V3d p;
      //     data->simulationMesh->getVertex(ele, j, p.data());
      //     restPosition.segment<3>(3 * j) = p;
      //   }

      //   data->elementFEMs[ele] = new KoiterShellDeformationModel(restPosition.data(), restPosition.data() + 3, restPosition.data() + 6, restPosition.data() + 9, restPosition.data() + 12,
      //     restPosition.data() + 15, data->elementMaterials[ele], pm, h);
      // }
      // else if (data->simulationMesh->getElementType() == SimulationMeshType::EDGE_QUAD) {
      //   ES::V12d restPosition;
      //   for (int j = 0; j < 4; j++) {
      //     ES::V3d p;
      //     data->simulationMesh->getVertex(ele, j, p.data());
      //     restPosition.segment<3>(j * 3) = p;
      //   }

      //   int useTangent = 0;
      //   if (elasticMaterialType == DeformationModelElasticMaterial::CLOTH_BENDING_QUADRATIC)
      //     useTangent = 1;

      //   data->elementFEMs[ele] = new EdgeDeformationModel(
      //     restPosition.data(), restPosition.data() + 3, restPosition.data() + 6, restPosition.data() + 9,
      //     data->elementMaterials[ele], pm, useTangent);
      // }
      else {
        throw std::runtime_error("only tet mesh is supported");
      }
    },
    tbb::static_partitioner());
}

const DeformationModel *DeformationModelManager::getDeformationModel(int eleID) const
{
  return data->elementFEMs[eleID];
}

const SimulationMesh *DeformationModelManager::getMesh() const
{
  return data->simulationMesh;
}

void DeformationModelManager::updateMeshRigidTransformation(const double R[9])
{
  data->globalRotation = Eigen::Map<const ES::M3d>(R);
  tbb::parallel_for(
    0, (int)data->fiberAxes.cols() / 3, [this](int i) {
      data->fiberAxes.block<3, 3>(0, i * 3) = data->fiberAxesRest.block<3, 3>(0, i * 3) * data->globalRotation.transpose();
    },
    tbb::static_partitioner());

  tbb::parallel_for(
    0, (int)data->vertexFiberAxes.cols() / 3, [this](int i) {
      data->vertexFiberAxes.block<3, 3>(0, i * 3) = data->vertexFiberAxesRest.block<3, 3>(0, i * 3) * data->globalRotation.transpose();
    },
    tbb::static_partitioner());

  tbb::parallel_for(
    0, data->nele, [this](int ele) {
      if (data->plasticVol3DOF[ele])
        data->plasticVol3DOF[ele]->setR(data->fiberAxes.data() + ele * 9);
    },
    tbb::static_partitioner());
}

void DeformationModelManager::getVertexAlignedMatrix(int id, double R[9]) const
{
  (Eigen::Map<ES::M3d>(R)) = data->vertexFiberAxes.block<3, 3>(0, id * 3);
}

void DeformationModelManager::getElementAlignedMatrix(int id, double R[9]) const
{
  if (data->plasticVol6DOF[id] || data->plasticVolConstant[id]) {
    (Eigen::Map<ES::M3d>(R)) = ES::M3d::Identity();
  }
  else {
    (Eigen::Map<ES::M3d>(R)) = data->fiberAxes.block<3, 3>(0, id * 3);
  }
}

void DeformationModelManager::setElementAlignedMatrix(int id, double R[9])
{
  data->fiberAxesRest.block<3, 3>(0, id * 3) = Eigen::Map<ES::M3d>(R);
  data->fiberAxes.block<3, 3>(0, id * 3) = data->fiberAxesRest.block<3, 3>(0, id * 3) * data->globalRotation.transpose();

  if (data->plasticVol3DOF[id] == nullptr)
    return;

  data->plasticVol3DOF[id]->setR(data->fiberAxes.data() + id * 9);
}

void DeformationModelManager::getMaterialParameters(double &EActive, double &EPassive, double &nu, double &gamma, double &lo) const
{
  EActive = data->Eact;
  EPassive = data->Epass;
  nu = data->nu;
  gamma = data->gamma;
  lo = data->lo;
}

int DeformationModelManager::saveMaterialParameters(const char *filename) const
{
  nlohmann::json jconfig;

  jconfig["E-active"] = data->Eact;
  jconfig["E-passive"] = data->Epass;
  jconfig["E-tendon"] = data->Etendon;
  jconfig["nu"] = data->nu;
  jconfig["gamma"] = data->gamma;
  jconfig["lo"] = data->lo;
  jconfig["spd"] = data->enforceSPD;
  jconfig["compression"] = data->compressionRatio;
  jconfig["use"] = data->hasUserDefinedMaterial;

  if (std::ofstream(filename) << jconfig.dump(2))
    return 0;
  else
    return 1;
}

int DeformationModelManager::saveDefaultMaterialParameters(const char *filename)
{
  nlohmann::json jconfig;

  jconfig["E-active"] = 0.1e6;
  jconfig["E-passive"] = 6e3;
  jconfig["E-tendon"] = 40e3;
  jconfig["nu"] = 0.49;
  jconfig["gamma"] = 1.0;
  jconfig["lo"] = 0.6;
  jconfig["spd"] = 0;
  jconfig["compression"] = 0;
  jconfig["use"] = 0;

  if (std::ofstream(filename) << jconfig.dump(2))
    return 0;
  else
    return 1;
}

int DeformationModelManager::getNumPlasticParameters() const
{
  return data->numPlasticParams;
}

int DeformationModelManager::getNumElasticParameters() const
{
  return data->elementMaterials[0]->getNumParameters();
}
