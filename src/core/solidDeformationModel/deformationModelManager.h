/*
author: Bohan Wang
copyright to USC,MIT
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
  STVK,
  LINEAR,
  VOLUME,

  HILL_STABLE_NEO,
  HILL_STVK_VOL,
  HILL_STVK,

  // LOWER_INVARIENT_CUBIC_SPLINE,

  // CLOTH_STVK,
  // CLOTH_BENDING_QUADRATIC,
  // CLOTH_BENDING_CUBIC_SPLINE,
  // KOITER_SHELL,
};

enum class DeformationModelPlasticMaterial
{
  VOLUMETRIC_DOF0 = 0,
  VOLUMETRIC_DOF3 = 1,
  VOLUMETRIC_DOF6 = 2,

  // SURFACE_DOF0 = 3,
  // SURFACE_DOF3 = 4,
  // BEND_DOF0 = 5,
  // BEND_DOF1 = 6,

  // SHELL_DOF6 = 7,
};

class DeformationModelManager
{
public:
  DeformationModelManager();
  ~DeformationModelManager();

  void setMesh(const SimulationMesh *simulationMesh, const double *elementFiberDirections = nullptr, const double *vertexFiberDirections = nullptr);
  void setSplineFunction(int numPts, double *x, double *y1, double *y2, double *y3);

  void init(DeformationModelPlasticMaterial plasticModelType, DeformationModelElasticMaterial elasticMaterialType, const char *muscleParameterConfigFilename = nullptr);
  void updateMeshRigidTransformation(const double R[9]);

  int getNumPlasticParameters() const;
  int getNumElasticParameters() const;
  const SimulationMesh *getMesh() const;

  void setElementAlignedMatrix(int id, double R[9]);
  void getElementAlignedMatrix(int id, double R[9]) const;
  void getVertexAlignedMatrix(int id, double R[9]) const;
  const DeformationModel *getDeformationModel(int eleID) const;

  void getMaterialParameters(double &EActive, double &EPassive, double &nu, double &gamma, double &lo) const;
  int saveMaterialParameters(const char *filename) const;
  static int saveDefaultMaterialParameters(const char *filename);

protected:
  DeformationModelManagerImpl *data;
};

}  // namespace SolidDeformationModel
}  // namespace pgo
