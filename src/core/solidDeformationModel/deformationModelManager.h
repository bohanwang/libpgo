/*
author: Bohan Wang
copyright to USC, MIT, NUS
*/

#pragma once

namespace pgo
{
namespace SolidDeformationModel
{
class SimulationMesh;
class DeformationModel;
class DeformationModelManagerImpl;

enum class DeformationModelElasticMaterial
{
  STABLE_NEO,
  STVK_VOL,
  INV_STVK,
  LINEAR,
  VOLUME,

  HILL_STABLE_NEO,
  HILL_STVK_VOL,
  HILL_STVK,

  STVK,
  MOONEY_RIVLIN,

  KOITER_FABRIC,
  KOITER_STVK,
};

enum class DeformationModelPlasticMaterial
{
  VOLUMETRIC_DOF0 = 0,
  VOLUMETRIC_DOF3 = 1,
  VOLUMETRIC_DOF6 = 2,

  SHELL_FF_DOF0 = 3,
  SHELL_FF_DOF1 = 4,
};

class DeformationModelManager
{
public:
  DeformationModelManager();
  ~DeformationModelManager();

  void setMesh(const SimulationMesh *simulationMesh, const double *elementFiberDirections = nullptr, const double *vertexFiberDirections = nullptr);
  void init(DeformationModelPlasticMaterial plasticModelType, DeformationModelElasticMaterial elasticMaterialType, int enforceSPD);
  void updateMeshRigidTransformation(const double R[9]);

  int getNumPlasticParameters() const;
  int getNumElasticParameters() const;
  const SimulationMesh *getMesh() const;

  void setElementAlignedMatrix(int id, double R[9]);
  void getElementAlignedMatrix(int id, double R[9]) const;
  void getVertexAlignedMatrix(int id, double R[9]) const;

  const DeformationModel *getDeformationModel(int eleID) const;

protected:
  DeformationModelManagerImpl *data;
};

}  // namespace SolidDeformationModel
}  // namespace pgo
