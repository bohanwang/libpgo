/*
author: Bohan Wang
copyright to USC, MIT, NUS
*/

#include "deformationModelManager.h"

#include "deformationModel.h"
#include "tetMeshDeformationModel.h"
#include "koiterDeformationModel.h"

#include "simulationMesh.h"

#include "elasticModel.h"
#include "elasticModel3DDeformationGradient.h"
#include "elasticModelCombinedMaterial.h"
#include "elasticModelHillTypeMaterial.h"
#include "elasticModelInvariantBasedMaterial.h"
#include "elasticModelStableNeoHookeanMaterial.h"
#include "elasticModelVolumeMaterial.h"
#include "invariantBasedMaterialStVK.h"
#include "elasticModelLinearMaterial.h"
#include "elasticModel3DSTVKMaterial.h"
#include "elasticModel3DMooneyRivlin.h"

#include "elasticModel2DFundamentalForms.h"
#include "elasticModel2DFundamentalFormsFabric.h"
#include "elasticModel2DFundamentalFormsSTVK.h"

#include "plasticModel.h"
#include "plasticModel3DDeformationGradient.h"
#include "plasticModel3D3DOF.h"
#include "plasticModel3D6DOF.h"
#include "plasticModel3DConstant.h"

#include "plasticModel2DFundamentalForms.h"
#include "plasticModel2DFundamentalFormsUniformStretch.h"

#include "pgoLogging.h"
#include "EigenSupport.h"

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

  // element elastic material
  // volumetric elastic material
  std::vector<ElasticModelStableNeoHookeanMaterial *> stableNeoHookeanMaterials;
  std::vector<ElasticModelLinearMaterial *> linearMaterials;
  std::vector<ElasticModelHillTypeMaterial *> hillTypeMaterials;
  std::vector<ElasticModelInvariantBasedMaterial *> invariantBasedMaterials;
  std::vector<ElasticModelVolumeMaterial *> volumeMaterials;
  std::vector<ElasticModel3DSTVKMaterial *> stvkMaterials;
  std::vector<ElasticModel3DMooneyRivlin *> mooneyRivlinMaterials;
  std::vector<ElasticModelCombinedMaterial<2> *> combined2Materials;
  std::vector<ElasticModelCombinedMaterial<3> *> combined3Materials;

  // shell elastic material
  std::vector<ElasticModel2DFundamentalFormsFabric *> shellFabricMaterials;
  std::vector<ElasticModel2DFundamentalFormsSTVK *> shellSTVKMaterials;

  std::vector<ElasticModel *> elementMaterials;
  std::vector<InvariantBasedMaterial *> invariantModels;

  // plastic models
  std::vector<PlasticModel3DConstant *> plasticVolConstant;
  std::vector<PlasticModel3D3DOF *> plasticVol3DOF;
  std::vector<PlasticModel3D6DOF *> plasticVol6DOF;

  // plastic model shell
  std::vector<PlasticModel2DFundamentalForms *> plasticShellConstant;
  std::vector<PlasticModel2DFundamentalFormsUniformStretch *> plasticShellUniformStretch;

  ES::VXd fiberDirections;
  ES::VXd vertexFiberDirections;
  ES::M3Xd fiberAxesRest, vertexFiberAxesRest;
  ES::M3Xd fiberAxes, vertexFiberAxes;
  ES::M3d globalRotation;

  int numPlasticParams;
  int nele;
  int nvtx;
  int enforceSPD = 0;

  void computeFiberAxes();
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

  for (auto ptr : volumeMaterials)
    if (ptr)
      delete ptr;

  for (auto ptr : combined2Materials)
    if (ptr)
      delete ptr;

  for (auto ptr : combined3Materials)
    if (ptr)
      delete ptr;

  for (auto ptr : stvkMaterials)
    if (ptr)
      delete ptr;

  for (auto ptr : mooneyRivlinMaterials)
    if (ptr)
      delete ptr;

  for (auto ptr : invariantModels)
    if (ptr)
      delete ptr;

  for (auto ptr : shellFabricMaterials)
    if (ptr)
      delete ptr;

  for (auto ptr : shellSTVKMaterials)
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

  for (auto ptr : plasticShellConstant)
    if (ptr)
      delete ptr;

  for (auto ptr : plasticShellUniformStretch)
    if (ptr)
      delete ptr;
}

void DeformationModelManagerImpl::computeFiberAxes()
{
  SPDLOG_LOGGER_INFO(pgo::Logging::lgr(), "Computing the fiber alignment transformation...");

  std::vector<std::vector<int>> vertexNearbyElements(nvtx);
  std::vector<std::vector<int>> elementNearbyElements(nele);

  for (int i = 0; i < nele; i++) {
    for (int j = 0; j < simulationMesh->getNumElementVertices(); j++) {
      vertexNearbyElements[simulationMesh->getVertexIndex(i, j)].push_back(i);
    }
  }

  for (int i = 0; i < nele; i++) {
    for (int j = 0; j < simulationMesh->getNumElementVertices(); j++) {
      const auto &eles = vertexNearbyElements[simulationMesh->getVertexIndex(i, j)];
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

  fiberAxesRest.resize(3, nele * 3);

  ES::V3d guideVector(0, 0, 1);
  std::vector<int> axisComputed(nele, 0);
  while (1) {
    int allComputed = 1;
    for (int f : axisComputed)
      allComputed &= f;

    if (allComputed)
      break;

    for (int ele = 0; ele < nele; ele++) {
      if (axisComputed[ele] == 1)
        continue;

      ES::V3d hintY = guideVector;
      ES::V3d X = fiberDirections.segment<3>(ele * 3);
      X.normalize();

      // if it is the degenerated case
      if (std::abs(X.dot(hintY)) > 1 - 1e-6) {
        // we first search nearby elements
        bool nearbyFound = false;
        for (int elej : elementNearbyElements[ele]) {
          if (axisComputed[elej] == 0) {
            continue;
          }

          // fetch the computed direction
          ES::M3d R = fiberAxesRest.block<3, 3>(0, elej * 3);
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

      fiberAxesRest.block<3, 3>(0, ele * 3) = R;
      axisComputed[ele] = 1;
    }

    int numAxes = std::count_if(axisComputed.begin(), axisComputed.end(), [](int val) { return val == 1; });
    SPDLOG_LOGGER_INFO(pgo::Logging::lgr(), "# axes computed:{}", numAxes);
  }

  vertexFiberAxesRest.resize(3, nvtx * 3);
  for (int vi = 0; vi < nvtx; vi++) {
    ES::V3d hintY = fiberAxesRest.col(vertexNearbyElements[vi][0] * 3 + 1);
    ES::V3d X = vertexFiberDirections.segment<3>(vi * 3);
    X.normalize();

    ES::V3d Z = X.cross(hintY);
    Z.normalize();

    ES::V3d Y = Z.cross(X);
    Y.normalize();

    ES::M3d R;
    R.row(0) = X;
    R.row(1) = Y;
    R.row(2) = Z;

    vertexFiberAxesRest.block<3, 3>(0, vi * 3) = R;
  }

  fiberAxes = fiberAxesRest;
  vertexFiberAxes = vertexFiberAxesRest;
}

std::map<DeformationModelPlasticMaterial, int> numPlasticDOFs{
  { DeformationModelPlasticMaterial::VOLUMETRIC_DOF0, 0 },
  { DeformationModelPlasticMaterial::VOLUMETRIC_DOF3, 3 },
  { DeformationModelPlasticMaterial::VOLUMETRIC_DOF6, 6 },

  { DeformationModelPlasticMaterial::SHELL_FF_DOF0, 0 },
  { DeformationModelPlasticMaterial::SHELL_FF_DOF1, 1 },

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
    data->fiberDirections.setZero(0);
  }

  if (vertexFiberDirections) {
    data->vertexFiberDirections = Eigen::Map<const ES::VXd>(vertexFiberDirections, data->nvtx * 3);
  }
  else {
    data->vertexFiberDirections.setZero(0);
  }
}

void DeformationModelManager::init(DeformationModelPlasticMaterial plasticModelType, DeformationModelElasticMaterial elasticMaterialType, int enforceSPD)
{
  data->enforceSPD = enforceSPD;

  SPDLOG_LOGGER_INFO(pgo::Logging::lgr(), "Computing the type of each element...");

  if (data->fiberDirections.size() || data->vertexFiberDirections.size()) {
    data->computeFiberAxes();
  }

  auto it = numPlasticDOFs.find(plasticModelType);
  PGO_ALOG(it != numPlasticDOFs.end());
  data->numPlasticParams = it->second;
  data->globalRotation = ES::M3d::Identity();

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
    elasticMaterialType == DeformationModelElasticMaterial::INV_STVK ||
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

  if (elasticMaterialType == DeformationModelElasticMaterial::LINEAR) {
    data->linearMaterials.assign(data->nele, nullptr);
  }

  if (elasticMaterialType == DeformationModelElasticMaterial::STVK) {
    data->stvkMaterials.assign(data->nele, nullptr);
  }

  if (elasticMaterialType == DeformationModelElasticMaterial::MOONEY_RIVLIN) {
    data->mooneyRivlinMaterials.assign(data->nele, nullptr);
  }

  if (elasticMaterialType == DeformationModelElasticMaterial::KOITER_FABRIC) {
    data->shellFabricMaterials.assign(data->nele, nullptr);
  }

    if (elasticMaterialType == DeformationModelElasticMaterial::KOITER_STVK) {
    data->shellSTVKMaterials.assign(data->nele, nullptr);
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
  else if (plasticModelType == DeformationModelPlasticMaterial::SHELL_FF_DOF0) {
    data->plasticShellConstant.assign(data->nele, nullptr);
  }
  else if (plasticModelType == DeformationModelPlasticMaterial::SHELL_FF_DOF1) {
    data->plasticShellUniformStretch.assign(data->nele, nullptr);
  }

  // for (int ele = 0; ele < data->nele; ele++) {
  tbb::parallel_for(
    0, data->nele, [&](int ele) {
      ES::V3d dir(1, 0, 0);
      if (data->fiberAxesRest.size() > 0) {
        dir = data->fiberAxesRest.block<3, 3>(0, ele * 3).row(0);
      }

      if (elasticMaterialType == DeformationModelElasticMaterial::HILL_STABLE_NEO) {
        const SimulationMeshENuMaterial *mat = dynamic_cast<const SimulationMeshENuMaterial *>(data->simulationMesh->getElementMaterial(ele, 0));
        PGO_ALOG(mat != nullptr);

        double mu = mat->getMuLame();
        double lambda = mat->getLambdaLame();

        data->stableNeoHookeanMaterials[ele] = new ElasticModelStableNeoHookeanMaterial(mu, lambda);
        data->stableNeoHookeanMaterials[ele]->enforceSPD(data->enforceSPD);

        const SimulationMeshHillMaterial *hillMat = dynamic_cast<const SimulationMeshHillMaterial *>(data->simulationMesh->getElementMaterial(ele, 1));
        PGO_ALOG(hillMat != nullptr);
        double Eact = hillMat->getEact();
        double gamma = hillMat->getGamma();
        double lo = hillMat->getLo();

        data->hillTypeMaterials[ele] = new ElasticModelHillTypeMaterial(gamma, Eact, lo, dir.data());
        data->combined2Materials[ele] = new ElasticModelCombinedMaterial<2>(data->stableNeoHookeanMaterials[ele], data->hillTypeMaterials[ele]);

        data->elementMaterials[ele] = data->combined2Materials[ele];
      }
      else if (elasticMaterialType == DeformationModelElasticMaterial::LINEAR) {
        const SimulationMeshENuMaterial *mat = dynamic_cast<const SimulationMeshENuMaterial *>(data->simulationMesh->getElementMaterial(ele, 0));
        PGO_ALOG(mat != nullptr);

        double mu = mat->getMuLame();
        double lambda = mat->getLambdaLame();

        data->linearMaterials[ele] = new ElasticModelLinearMaterial(mu, lambda);
        data->elementMaterials[ele] = data->linearMaterials[ele];
      }
      else if (elasticMaterialType == DeformationModelElasticMaterial::HILL_STVK) {
        const SimulationMeshENuMaterial *mat = dynamic_cast<const SimulationMeshENuMaterial *>(data->simulationMesh->getElementMaterial(ele, 0));
        PGO_ALOG(mat != nullptr);

        double E = mat->getE();
        double nu = mat->getNu();
        double compressionRatio = mat->getCompressionRatio();

        data->invariantModels[ele] = new InvariantBasedMaterialStVK(E, nu, compressionRatio);

        data->invariantBasedMaterials[ele] = new ElasticModelInvariantBasedMaterial(data->invariantModels[ele]);
        data->invariantBasedMaterials[ele]->enforceSPD(data->enforceSPD);

        const SimulationMeshHillMaterial *hillMat = dynamic_cast<const SimulationMeshHillMaterial *>(data->simulationMesh->getElementMaterial(ele, 1));
        PGO_ALOG(hillMat != nullptr);
        double Eact = hillMat->getEact();
        double gamma = hillMat->getGamma();
        double lo = hillMat->getLo();

        data->hillTypeMaterials[ele] = new ElasticModelHillTypeMaterial(gamma, Eact, lo, dir.data());
        data->combined2Materials[ele] = new ElasticModelCombinedMaterial<2>(data->invariantBasedMaterials[ele], data->hillTypeMaterials[ele]);

        data->elementMaterials[ele] = data->combined2Materials[ele];
      }
      else if (elasticMaterialType == DeformationModelElasticMaterial::HILL_STVK_VOL) {
        const SimulationMeshENuMaterial *mat = dynamic_cast<const SimulationMeshENuMaterial *>(data->simulationMesh->getElementMaterial(ele, 0));
        PGO_ALOG(mat != nullptr);

        double E = mat->getE();
        double nu = mat->getNu();
        double compressionRatio = mat->getCompressionRatio();

        data->invariantModels[ele] = new InvariantBasedMaterialStVK(E, nu, compressionRatio);

        data->invariantBasedMaterials[ele] = new ElasticModelInvariantBasedMaterial(data->invariantModels[ele]);
        data->invariantBasedMaterials[ele]->enforceSPD(data->enforceSPD);

        const SimulationMeshHillMaterial *hillMat = dynamic_cast<const SimulationMeshHillMaterial *>(data->simulationMesh->getElementMaterial(ele, 1));
        PGO_ALOG(hillMat != nullptr);
        double Eact = hillMat->getEact();
        double gamma = hillMat->getGamma();
        double lo = hillMat->getLo();

        data->hillTypeMaterials[ele] = new ElasticModelHillTypeMaterial(gamma, Eact, lo, dir.data());

        data->volumeMaterials[ele] = new ElasticModelVolumeMaterial(compressionRatio);

        data->combined3Materials[ele] = new ElasticModelCombinedMaterial<3>(data->invariantBasedMaterials[ele], data->hillTypeMaterials[ele], data->volumeMaterials[ele]);

        data->elementMaterials[ele] = data->combined3Materials[ele];
      }
      else if (elasticMaterialType == DeformationModelElasticMaterial::STABLE_NEO) {
        const SimulationMeshENuMaterial *mat = dynamic_cast<const SimulationMeshENuMaterial *>(data->simulationMesh->getElementMaterial(ele, 0));
        PGO_ALOG(mat != nullptr);

        double mu = mat->getMuLame();
        double lambda = mat->getLambdaLame();

        data->stableNeoHookeanMaterials[ele] = new ElasticModelStableNeoHookeanMaterial(mu, lambda);
        data->stableNeoHookeanMaterials[ele]->enforceSPD(data->enforceSPD);
        // data->stableNeoHookeanMaterials[ele]->enforceSPD(true);

        data->elementMaterials[ele] = data->stableNeoHookeanMaterials[ele];
      }
      else if (elasticMaterialType == DeformationModelElasticMaterial::INV_STVK) {
        const SimulationMeshENuMaterial *mat = dynamic_cast<const SimulationMeshENuMaterial *>(data->simulationMesh->getElementMaterial(ele, 0));
        PGO_ALOG(mat != nullptr);

        double E = mat->getE();
        double nu = mat->getNu();
        double compressionRatio = mat->getCompressionRatio();

        data->invariantModels[ele] = new InvariantBasedMaterialStVK(E, nu, compressionRatio);

        data->invariantBasedMaterials[ele] = new ElasticModelInvariantBasedMaterial(data->invariantModels[ele]);
        data->invariantBasedMaterials[ele]->enforceSPD(data->enforceSPD);

        data->elementMaterials[ele] = data->invariantBasedMaterials[ele];
      }
      else if (elasticMaterialType == DeformationModelElasticMaterial::STVK_VOL) {
        const SimulationMeshENuMaterial *mat = dynamic_cast<const SimulationMeshENuMaterial *>(data->simulationMesh->getElementMaterial(ele, 0));
        PGO_ALOG(mat != nullptr);

        double E = mat->getE();
        double nu = mat->getNu();
        double compressionRatio = mat->getCompressionRatio();

        data->invariantModels[ele] = new InvariantBasedMaterialStVK(E, nu, compressionRatio);

        data->invariantBasedMaterials[ele] = new ElasticModelInvariantBasedMaterial(data->invariantModels[ele]);
        data->invariantBasedMaterials[ele]->enforceSPD(data->enforceSPD);

        data->volumeMaterials[ele] = new ElasticModelVolumeMaterial(compressionRatio);

        data->combined2Materials[ele] = new ElasticModelCombinedMaterial<2>(data->invariantBasedMaterials[ele], data->volumeMaterials[ele]);

        data->elementMaterials[ele] = data->combined2Materials[ele];
      }
      else if (elasticMaterialType == DeformationModelElasticMaterial::VOLUME) {
        const SimulationMeshENuMaterial *mat = dynamic_cast<const SimulationMeshENuMaterial *>(data->simulationMesh->getElementMaterial(ele, 0));
        PGO_ALOG(mat != nullptr);

        double compressionRatio = mat->getCompressionRatio();

        data->volumeMaterials[ele] = new ElasticModelVolumeMaterial(compressionRatio);

        data->elementMaterials[ele] = data->volumeMaterials[ele];
      }
      else if (elasticMaterialType == DeformationModelElasticMaterial::STVK) {
        const SimulationMeshENuMaterial *mat = dynamic_cast<const SimulationMeshENuMaterial *>(data->simulationMesh->getElementMaterial(ele, 0));
        PGO_ALOG(mat != nullptr);

        double mu = mat->getMuLame();
        double lambda = mat->getLambdaLame();

        data->stvkMaterials[ele] = new ElasticModel3DSTVKMaterial(mu, lambda);
        data->elementMaterials[ele] = data->stvkMaterials[ele];
      }
      else if (elasticMaterialType == DeformationModelElasticMaterial::MOONEY_RIVLIN) {
        const SimulationMeshMooneyRivlinMaterial *mat = dynamic_cast<const SimulationMeshMooneyRivlinMaterial *>(data->simulationMesh->getElementMaterial(ele, 0));
        PGO_ALOG(mat != nullptr);

        int N = mat->getN();
        int M = mat->getM();

        data->mooneyRivlinMaterials[ele] = new ElasticModel3DMooneyRivlin(N, mat->getC(), M, mat->getD());
        data->elementMaterials[ele] = data->mooneyRivlinMaterials[ele];
      }
      else if (elasticMaterialType == DeformationModelElasticMaterial::KOITER_FABRIC) {
        // const SimulationMeshKoiterFabricMaterial *mat = dynamic_cast<const SimulationMeshKoiterFabricMaterial *>(data->simulationMesh->getElementMaterial(ele, 0));
        // PGO_ALOG(mat != nullptr);
        // double E1 = mat->getE1();
        // double E2 = mat->getE2();
        // double G12 = mat->getG12();
        // double nu12 = mat->getNu12();
        // double nu21 = mat->getNu21();
        // double bendingE1 = mat->getBendingE1();
        // double bendingE2 = mat->getBendingE2();
        // double bendingG12 = mat->getBendingG12();
        ES::V2d dir0(1, 0), dir1(0, 1);
        data->shellFabricMaterials[ele] = new ElasticModel2DFundamentalFormsFabric(dir0, dir1);
        data->elementMaterials[ele] = data->shellFabricMaterials[ele];
      }
      else if (elasticMaterialType == DeformationModelElasticMaterial::KOITER_STVK) {        
        data->shellSTVKMaterials[ele] = new ElasticModel2DFundamentalFormsSTVK;
        data->elementMaterials[ele] = data->shellSTVKMaterials[ele];
      }
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
      else if (plasticModelType == DeformationModelPlasticMaterial::SHELL_FF_DOF0) {
        data->plasticShellConstant[ele] = new PlasticModel2DFundamentalForms();
        pm = data->plasticShellConstant[ele];
      }
      else if (plasticModelType == DeformationModelPlasticMaterial::SHELL_FF_DOF1) {
        data->plasticShellUniformStretch[ele] = new PlasticModel2DFundamentalFormsUniformStretch();
        pm = data->plasticShellUniformStretch[ele];
      }
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
      else if (data->simulationMesh->getElementType() == SimulationMeshType::SHELL) {
        ES::V18d restPosition;
        for (int j = 0; j < 6; j++) {
          ES::V3d p;
          if (data->simulationMesh->getVertexIndex(ele, j) < 0) {
            p = ES::V3d(-10496, -10496, -10496);
          }
          else {
            data->simulationMesh->getVertex(ele, j, p.data());
          }
          restPosition.segment<3>(3 * j) = p;
        }

        double h = 0;
        const SimulationMeshENuhMaterial *mat = dynamic_cast<const SimulationMeshENuhMaterial *>(data->simulationMesh->getElementMaterial(ele, 0));
        if (mat == nullptr) {
          const SimulationMeshMooneyRivlinhMaterial *mat = dynamic_cast<const SimulationMeshMooneyRivlinhMaterial *>(data->simulationMesh->getElementMaterial(ele, 0));
          PGO_ALOG(mat != nullptr);

          h = mat->geth();
        }
        else {
          h = mat->geth();
        }

        if (elasticMaterialType == DeformationModelElasticMaterial::STVK ||
          elasticMaterialType == DeformationModelElasticMaterial::LINEAR ||
          elasticMaterialType == DeformationModelElasticMaterial::MOONEY_RIVLIN) {
            throw std::logic_error("unsupported elastic material for shell element");
         
        }
        else if (elasticMaterialType == DeformationModelElasticMaterial::KOITER_FABRIC ||
          elasticMaterialType == DeformationModelElasticMaterial::KOITER_STVK) {
          data->elementFEMs[ele] = new KoiterDeformationModel(restPosition.data(), restPosition.data() + 3, restPosition.data() + 6, restPosition.data() + 9, restPosition.data() + 12,
            restPosition.data() + 15, data->elementMaterials[ele], pm, data->enforceSPD);
        }
        else {
          throw std::logic_error("unsupported elastic material for shell element");
          {};
        }
      }
      else {
        throw std::logic_error("unknown mesh element type");
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

int DeformationModelManager::getNumPlasticParameters() const
{
  return data->numPlasticParams;
}

int DeformationModelManager::getNumElasticParameters() const
{
  return data->elementMaterials[0]->getNumParameters();
}
