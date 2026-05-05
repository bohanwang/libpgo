#include "generateFBMSUnionSurfaceVolume.h"

#include "generateFBMSUnionSurfaceExtraction.h"
#include "generateFBMSUnionSurfaceField.h"

#include <cmath>
#include <iostream>

namespace pgo::Tools::FBMSUnionSurface
{

double computeMeshVolume(const Mesh::TriMeshGeo &mesh)
{
  // Signed volume of a closed polyhedron via the divergence theorem on triangles:
  // V = (1/6) * sum_t (a . (b x c)). Sign depends on triangle winding; absolute
  // value gives the enclosed volume for a watertight mesh.
  double v6 = 0.0;
  const int numTri = mesh.numTriangles();
  for (int t = 0; t < numTri; ++t) {
    const Vec3d &a = mesh.pos(t, 0);
    const Vec3d &b = mesh.pos(t, 1);
    const Vec3d &c = mesh.pos(t, 2);
    v6 += a.dot(b.cross(c));
  }
  return std::abs(v6) / 6.0;
}

ThicknessSearchResult findThicknessForVolumeBudget(
  const EigenSupport::VXd &fbmsDistance, const SphereParameters &sphere,
  const EigenSupport::V3d &bmin, const EigenSupport::V3d &bmax, int resolution,
  double sphereThickness, bool enableTruncating, const std::string &isoOffsetMode, double fixedIsoOffset,
  double targetVolume,
  double tLo, double tHi, int maxIterations, double relTol)
{
  EigenSupport::VXd unionField;
  double fmin = 0.0;
  double fmax = 0.0;

  auto evalAt = [&](double t, Mesh::TriMeshGeo &outMesh) -> double {
    assembleUnionField(fbmsDistance, sphere, bmin, bmax, resolution,
      t, sphereThickness, enableTruncating, "union", unionField, fmin, fmax);
    const double isoOffset = selectIsoOffset(isoOffsetMode, fixedIsoOffset,
      unionField, t, sphereThickness, bmin, bmax, resolution);
    if (!(fmin <= isoOffset && fmax >= isoOffset)) {
      // No selected isosurface crossing: nothing to mesh; treat as volume 0.
      outMesh = Mesh::TriMeshGeo();
      return 0.0;
    }
    extractSurfaceMarchingCubes(bmin, bmax, resolution, unionField, isoOffset, outMesh);
    if (outMesh.numTriangles() == 0)
      return 0.0;
    return computeMeshVolume(outMesh);
  };

  ThicknessSearchResult result;
  Mesh::TriMeshGeo meshLo, meshHi, meshMid;

  const double vLo = evalAt(tLo, meshLo);
  std::cout << "[volume-search] tLo=" << tLo << " vol=" << vLo << " (target=" << targetVolume << ")" << std::endl;

  if (vLo > targetVolume) {
    result.budgetExceededAtMin = true;
    result.thickness = tLo;
    result.volume = vLo;
    result.mesh = std::move(meshLo);
    return result;
  }

  const double vHi = evalAt(tHi, meshHi);
  std::cout << "[volume-search] tHi=" << tHi << " vol=" << vHi << std::endl;

  if (vHi <= targetVolume) {
    result.budgetUnreachedAtMax = true;
    result.thickness = tHi;
    result.volume = vHi;
    result.mesh = std::move(meshHi);
    return result;
  }

  result.thickness = tLo;
  result.volume = vLo;
  result.mesh = meshLo;

  for (int it = 0; it < maxIterations; ++it) {
    const double mid = 0.5 * (tLo + tHi);
    const double vMid = evalAt(mid, meshMid);
    std::cout << "[volume-search] iter=" << it << " t=" << mid << " vol=" << vMid << std::endl;

    if (vMid <= targetVolume) {
      tLo = mid;
      result.thickness = mid;
      result.volume = vMid;
      result.mesh = meshMid;
    }
    else {
      tHi = mid;
    }
    result.iterations = it + 1;
    if (tHi <= 0.0 || (tHi - tLo) / tHi < relTol)
      break;
  }
  return result;
}

}  // namespace pgo::Tools::FBMSUnionSurface
