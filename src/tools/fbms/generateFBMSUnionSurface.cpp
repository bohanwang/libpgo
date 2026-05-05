#include "generateFBMSUnionSurfaceExtraction.h"
#include "generateFBMSUnionSurfaceField.h"
#include "generateFBMSUnionSurfaceIO.h"
#include "generateFBMSUnionSurfaceOptions.h"
#include "generateFBMSUnionSurfaceSphere.h"
#include "generateFBMSUnionSurfaceVolume.h"

#include "boundingBox.h"
#include "pgoLogging.h"

#include <algorithm>
#include <exception>
#include <iostream>

namespace fbms_union = pgo::Tools::FBMSUnionSurface;

int main(int argc, char *argv[])
{
  try {
    const fbms_union::Options options = fbms_union::parseOptions(argc, argv);
    fbms_union::validateOptions(options);

    pgo::Logging::init();

    pgo::Mesh::TriMeshGeo fbmsMesh;
    pgo::Mesh::TriMeshGeo sphereMesh;
    fbms_union::loadMeshOrThrow(options.fbmsPath, fbmsMesh, "FBMS");
    fbms_union::loadMeshOrThrow(options.spherePath, sphereMesh, "sphere");

    if (fbmsMesh.numTriangles() == 0)
      throw std::runtime_error("FBMS mesh has no triangles: " + options.fbmsPath);

    const fbms_union::SphereParameters sphere = fbms_union::loadSphereParameters(options.spherePath, sphereMesh);
    std::cout << "Sphere parameter source: " << (sphere.fromHeader ? "header" : "bbox fallback") << std::endl;
    fbms_union::printVector("Sphere center", sphere.center);
    std::cout << "Sphere radius = " << sphere.radius << std::endl;

    if (options.projectFBMSBoundaryToSphere) {
      const int projectedVertices = fbms_union::projectBoundaryVerticesToSphere(fbmsMesh, sphere);
      std::cout << "Projected FBMS boundary vertices to sphere = " << projectedVertices << std::endl;
    }

    if (options.enableTruncating)
      std::cout << "Truncation enabled: FBMS clipped to ball of radius "
                << sphere.radius << std::endl;

    pgo::Mesh::TriMeshGeo rawSurface;
    pgo::EigenSupport::V3d bmin, bmax;

    if (options.fbmsThickness > 0.0) {
      fbms_union::computeUnionBBox(fbmsMesh, sphereMesh, options.fbmsThickness,
        options.sphereThickness, options.paddingRatio, bmin, bmax);

      if (options.extractionBackend == "openvdb") {
        const double derivedVoxelSize = (bmax - bmin).maxCoeff() / static_cast<double>(options.resolution - 1);
        const double voxelSize = options.vdbVoxelSize > 0.0 ? options.vdbVoxelSize : derivedVoxelSize;
        fbms_union::extractSurfaceOpenVDB(fbmsMesh, sphere,
          options.fbmsThickness, options.sphereThickness, options.enableTruncating, options.debugFieldMode,
          voxelSize, options.vdbHalfWidth, options.vdbAdaptivity, options.vdbSmoothSteps,
          rawSurface);

        std::cout << "Extraction backend = openvdb" << std::endl;
        std::cout << "Reference grid resolution = " << options.resolution << std::endl;
        fbms_union::printVector("Reference bbox min", bmin);
        fbms_union::printVector("Reference bbox max", bmax);
        std::cout << "OpenVDB voxel size = " << voxelSize << std::endl;
        std::cout << "OpenVDB half width = " << options.vdbHalfWidth << std::endl;
        std::cout << "OpenVDB adaptivity = " << options.vdbAdaptivity << std::endl;
        std::cout << "OpenVDB smooth steps = " << options.vdbSmoothSteps << std::endl;
      }
      else {
        pgo::EigenSupport::VXd fbmsDistance;
        fbms_union::computeFBMSDistance(fbmsMesh, bmin, bmax, options.resolution, fbmsDistance);

        pgo::EigenSupport::VXd unionField;
        double fieldMin = 0.0;
        double fieldMax = 0.0;
        fbms_union::assembleUnionField(fbmsDistance, sphere, bmin, bmax, options.resolution,
          options.fbmsThickness, options.sphereThickness, options.enableTruncating, options.debugFieldMode,
          unionField, fieldMin, fieldMax);

        const double isoOffset = fbms_union::selectIsoOffset(options.isoOffsetMode, options.isoOffset,
          unionField, options.fbmsThickness, options.sphereThickness, bmin, bmax, options.resolution);

        if (!(fieldMin <= isoOffset && fieldMax >= isoOffset))
          throw std::runtime_error("Union field does not cross the selected isosurface");

        fbms_union::extractSurfaceMarchingCubes(bmin, bmax, options.resolution, unionField, isoOffset, rawSurface);

        std::cout << "Extraction backend = marching-cubes" << std::endl;
        std::cout << "Grid resolution = " << options.resolution << std::endl;
        fbms_union::printVector("BBox min", bmin);
        fbms_union::printVector("BBox max", bmax);
        std::cout << "Field min = " << fieldMin << std::endl;
        std::cout << "Field max = " << fieldMax << std::endl;
        std::cout << "Iso offset mode = " << options.isoOffsetMode << std::endl;
        std::cout << "Selected iso offset = " << isoOffset << std::endl;
      }
    }
    else {
      // Volume-budget mode: |options.fbmsThickness| = target volume.
      const double targetVolume = -options.fbmsThickness;

      const pgo::Mesh::BoundingBox fbmsBBox(fbmsMesh.positions());
      const pgo::Mesh::BoundingBox sphereBBox(sphereMesh.positions());
      const pgo::EigenSupport::V3d baseSides =
        fbmsBBox.bmax().cwiseMax(sphereBBox.bmax())
        - fbmsBBox.bmin().cwiseMin(sphereBBox.bmin());
      const double baseDiag = baseSides.norm();

      // Search ceiling. With truncating, volume saturates once t exceeds the
      // ball diameter; without truncating, t can grow until the union fills
      // a sizable fraction of the union bbox.
      const double tHiCap = options.enableTruncating
        ? std::max(2.0 * sphere.radius, 4.0 * options.sphereThickness)
        : std::max(0.5 * baseDiag, 4.0 * options.sphereThickness);

      // Size the grid so the precomputed FBMS UDF is valid across the whole
      // search range. With truncating, the meshable region is clipped to the
      // ball, so we only need padding for the sphere shell. Otherwise we must
      // accommodate the search ceiling.
      const double bboxThickness = options.enableTruncating ? options.sphereThickness : tHiCap;
      fbms_union::computeUnionBBox(fbmsMesh, sphereMesh, bboxThickness,
        options.sphereThickness, options.paddingRatio, bmin, bmax);

      const pgo::EigenSupport::V3d voxel = (bmax - bmin) / static_cast<double>(options.resolution - 1);
      const double derivedVoxelSize = voxel.maxCoeff();
      const double openVDBVoxelSize = options.vdbVoxelSize > 0.0 ? options.vdbVoxelSize : derivedVoxelSize;
      const double tLo = 2.0 * (options.extractionBackend == "openvdb" ? openVDBVoxelSize : derivedVoxelSize);
      double searchHi = tHiCap;
      if (options.extractionBackend == "openvdb")
        searchHi = std::min(tHiCap, std::max(0.06, 2.0 * options.sphereThickness));

      std::cout << "[volume-search] target volume = " << targetVolume
                << ", search range t in [" << tLo << ", " << searchHi << "]" << std::endl;

      fbms_union::ThicknessSearchResult result;
      if (options.extractionBackend == "openvdb") {
        while (true) {
          result = fbms_union::findThicknessForVolumeBudgetOpenVDB(
            fbmsMesh, sphere, options.sphereThickness, options.enableTruncating,
            openVDBVoxelSize, options.vdbHalfWidth, options.vdbAdaptivity, options.vdbSmoothSteps,
            targetVolume, tLo, searchHi, /*maxIterations=*/20, /*relTol=*/1e-3);
          if (!result.budgetUnreachedAtMax || searchHi >= tHiCap)
            break;
          const double nextSearchHi = std::min(tHiCap, 2.0 * searchHi);
          std::cout << "[volume-search] OpenVDB search range did not exhaust budget; expanding tHi from "
                    << searchHi << " to " << nextSearchHi << std::endl;
          searchHi = nextSearchHi;
        }
      }
      else {
        pgo::EigenSupport::VXd fbmsDistance;
        fbms_union::computeFBMSDistance(fbmsMesh, bmin, bmax, options.resolution, fbmsDistance);

        result = fbms_union::findThicknessForVolumeBudget(
          fbmsDistance, sphere, bmin, bmax, options.resolution,
          options.sphereThickness, options.enableTruncating, options.isoOffsetMode, options.isoOffset, targetVolume,
          tLo, searchHi, /*maxIterations=*/20, /*relTol=*/1e-3);
      }

      if (result.budgetExceededAtMin)
        std::cout << "[volume-search] WARNING: even t=" << result.thickness
                  << " produces volume " << result.volume
                  << " > target " << targetVolume
                  << "; saving the minimum-thickness mesh anyway." << std::endl;
      else if (result.budgetUnreachedAtMax)
        std::cout << "[volume-search] WARNING: top of search range t=" << result.thickness
                  << " still has volume " << result.volume
                  << " <= target " << targetVolume
                  << "; budget is not exhausted." << std::endl;

      std::cout << "[volume-search] selected thickness = " << result.thickness
                << " volume = " << result.volume
                << " iterations = " << result.iterations << std::endl;

      rawSurface = result.mesh;

      if (options.extractionBackend == "openvdb") {
        std::cout << "Extraction backend = openvdb" << std::endl;
        std::cout << "Reference grid resolution = " << options.resolution << std::endl;
        fbms_union::printVector("Reference bbox min", bmin);
        fbms_union::printVector("Reference bbox max", bmax);
        std::cout << "OpenVDB voxel size = " << openVDBVoxelSize << std::endl;
        std::cout << "OpenVDB half width = " << options.vdbHalfWidth << std::endl;
        std::cout << "OpenVDB adaptivity = " << options.vdbAdaptivity << std::endl;
        std::cout << "OpenVDB smooth steps = " << options.vdbSmoothSteps << std::endl;
      }
      else {
        std::cout << "Extraction backend = marching-cubes" << std::endl;
        std::cout << "Grid resolution = " << options.resolution << std::endl;
        fbms_union::printVector("BBox min", bmin);
        fbms_union::printVector("BBox max", bmax);
      }
    }

    std::cout << "Raw vertex count = " << rawSurface.numVertices() << std::endl;
    std::cout << "Raw face count = " << rawSurface.numTriangles() << std::endl;

    if (rawSurface.numVertices() == 0 || rawSurface.numTriangles() == 0)
      throw std::runtime_error("Surface extraction produced an empty raw surface");

    fbms_union::applySmallComponentFilterIfRequested(rawSurface, options);
    fbms_union::saveSurfaceOrThrow(rawSurface, options.outputSurfacePath);
    return 0;
  }
  catch (const std::exception &err) {
    std::cerr << "Error: " << err.what() << std::endl;
    return 1;
  }
}
