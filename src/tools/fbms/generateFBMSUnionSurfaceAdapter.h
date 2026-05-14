#pragma once

#include "generateFBMSUnionSurfaceOptions.h"
#include "generateFBMSUnionSurfaceVolume.h"
#include "EigenSupport.h"
#include "triMeshGeo.h"

namespace pgo::Tools::FBMSUnionSurface
{

// ---- Sphere parameters (was generateFBMSUnionSurfaceSphere.h) ----

struct SphereParameters
{
  EigenSupport::V3d center = EigenSupport::V3d::Zero();
  double radius = 0.0;
  bool fromHeader = false;
};

SphereParameters loadSphereParameters(const std::string &spherePath, const Mesh::TriMeshGeo &sphereMesh);
int projectBoundaryVerticesToSphere(Mesh::TriMeshGeo &mesh, const SphereParameters &sphere);

// ---- Field operations (was generateFBMSUnionSurfaceField.h) ----

void computeUnionBBox(const Mesh::TriMeshGeo &fbmsMesh, const Mesh::TriMeshGeo &sphereMesh,
  double fbmsThickness, double sphereThickness, double paddingRatio,
  EigenSupport::V3d &bmin, EigenSupport::V3d &bmax);

void computeFBMSDistance(const Mesh::TriMeshGeo &fbmsMesh,
  const EigenSupport::V3d &bmin, const EigenSupport::V3d &bmax, int resolution,
  EigenSupport::VXd &fbmsDistance);

void assembleUnionField(const EigenSupport::VXd &fbmsDistance, const SphereParameters &sphere,
  const EigenSupport::V3d &bmin, const EigenSupport::V3d &bmax, int resolution,
  double fbmsThickness, double sphereThickness, bool enableTruncating, const std::string &debugFieldMode,
  EigenSupport::VXd &unionField, double &fieldMin, double &fieldMax);

double selectIsoOffset(const std::string &isoOffsetMode, double fixedIsoOffset,
  const EigenSupport::VXd &field, double fbmsThickness, double sphereThickness,
  const EigenSupport::V3d &bmin, const EigenSupport::V3d &bmax, int resolution);

// ---- Extraction operations (was generateFBMSUnionSurfaceExtraction.h) ----

void extractSurfaceMarchingCubes(const EigenSupport::V3d &bmin, const EigenSupport::V3d &bmax,
  int resolution, const EigenSupport::VXd &field, double isoOffset, Mesh::TriMeshGeo &outMesh);

void extractSurfaceOpenVDB(const Mesh::TriMeshGeo &fbmsMesh, const SphereParameters &sphere,
  double fbmsThickness, double sphereThickness, bool enableTruncating, const std::string &debugFieldMode,
  double voxelSize, double halfWidth, double adaptivity, int smoothSteps,
  Mesh::TriMeshGeo &outMesh);

ThicknessSearchResult findThicknessForVolumeBudgetOpenVDB(
  const Mesh::TriMeshGeo &fbmsMesh, const SphereParameters &sphere,
  double sphereThickness, bool enableTruncating,
  double voxelSize, double halfWidth, double adaptivity, int smoothSteps,
  double targetVolume, double tLo, double tHi, int maxIterations, double relTol);

// ---- Volume operations (was generateFBMSUnionSurfaceVolume.h) ----

double computeMeshVolume(const Mesh::TriMeshGeo &mesh);

ThicknessSearchResult findThicknessForVolumeBudget(
  const EigenSupport::VXd &fbmsDistance, const SphereParameters &sphere,
  const EigenSupport::V3d &bmin, const EigenSupport::V3d &bmax, int resolution,
  double sphereThickness, bool enableTruncating, const std::string &isoOffsetMode, double fixedIsoOffset,
  double targetVolume,
  double tLo, double tHi, int maxIterations, double relTol);

}  // namespace pgo::Tools::FBMSUnionSurface
