/*
author: Bohan Wang
copyright to USC,MIT,NUS
*/

#include "tetMeshOccupation.h"

#include "windingNumberTree.h"
#include "pgoLogging.h"
#include "triMeshPseudoNormal.h"
#include "boundingVolumeTree.h"

#include <tbb/parallel_for.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/cache_aligned_allocator.h>

#include <random>
#include <vector>
#include <atomic>

void pgo::SolidDeformationModel::computeTetMeshOccupation(int numTetVertices, const Vec3d *tetVertices, int numTets, const Vec4i *tets,
  int numSurfaceVertices, const Vec3d *surfaceVertices, int numTriangles, const Vec3i *triangles,
  double minThreshold, int sampleCount, double *weights)
{
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<double> randomNumber(0.0, 1.0);

  Mesh::TriMeshRef surfaceMeshRef(numSurfaceVertices, surfaceVertices, numTriangles, triangles);
  Mesh::TriMeshBVTree surfaceMeshBVTree;
  surfaceMeshBVTree.buildByInertiaPartition(surfaceMeshRef);

  Mesh::TriMeshPseudoNormal surfaceMeshNormal;
  surfaceMeshNormal.buildPseudoNormals(surfaceMeshRef);
  // WindingNumberTree windingNumberTree;
  // windingNumberTree.build(surfaceMeshRef);

  std::atomic<int> numFinished(0);

  SPDLOG_LOGGER_INFO(pgo::Logging::lgr(), "#tets: {}", numTets);

  tbb::parallel_for(0, numTets, [&](int ei) {
    int insideCounter = 0;
    for (int si = 0; si < sampleCount; si++) {
      double w[4];
      /*
      w[0] = randomNumber(gen);
      w[1] = (1 - w[0]) * randomNumber(gen);
      w[2] = (1 - w[0] - w[1]) * randomNumber(gen);
      w[3] = 1 - w[0] - w[1] - w[2];
      */

      double s = randomNumber(gen);
      double t = randomNumber(gen);
      double u = randomNumber(gen);
      if (s + t > 1.0) {  // cut'n fold the cube into a prism

        s = 1.0 - s;
        t = 1.0 - t;
      }

      if (t + u > 1.0) {  // cut'n fold the prism into a tetrahedron

        double tmp = u;
        u = 1.0 - s - t;
        t = 1.0 - tmp;
      }
      else if (s + t + u > 1.0) {
        double tmp = u;
        u = s + t + u - 1.0;
        s = 1 - t - tmp;
      }

      w[3] = 1 - s - t - u;  // a,s,t,u are the barycentric coordinates of the random point.
      w[2] = u;
      w[1] = t;
      w[0] = s;

      Vec3d p;
      p.setZero();
      for (int vi = 0; vi < 4; vi++) {
        p += w[vi] * tetVertices[tets[ei][vi]];
      }

      auto ret = surfaceMeshBVTree.closestTriangleQuery(surfaceMeshRef, p);
      Vec3d diff = p - ret.closestPosition;
      Vec3d n = surfaceMeshNormal.getPseudoNormal(surfaceMeshRef.triangles(), ret.triID, ret.feature);

      if (diff.dot(n) < 0) {
        insideCounter++;
      }
      // if (windingNumberTree.windingNumber(surfaceMeshRef, p) > 0.5) {
      //   insideCounter++;
      // }
    }

    weights[ei] = std::max((double)insideCounter / sampleCount, minThreshold);

    int cur = numFinished.fetch_add(1, std::memory_order_relaxed);

    if (cur % 1000 == 0) {
      std::cout << cur << ' ' << std::flush;
    }
  });

  std::cout << std::endl;
}

void pgo::SolidDeformationModel::computeTetMeshOccupation(int numTetVertices, const Vec3d *tetVertices, int numTets, const Vec4i *tets,
  const Vec3d bbIn[2], double minThreshold, int sampleCount, double *weights)
{
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<double> randomNumber(0.0, 1.0);
  std::atomic<int> numFinished(0);

  Mesh::BoundingBox bb(bbIn[0], bbIn[1]);
  SPDLOG_LOGGER_INFO(pgo::Logging::lgr(), "#tets: {}", numTets);

  tbb::parallel_for(0, numTets, [&](int ei) {
    int insideCounter = 0;
    for (int si = 0; si < sampleCount; si++) {
      double w[4];
      /*
      w[0] = randomNumber(gen);
      w[1] = (1 - w[0]) * randomNumber(gen);
      w[2] = (1 - w[0] - w[1]) * randomNumber(gen);
      w[3] = 1 - w[0] - w[1] - w[2];
      */

      double s = randomNumber(gen);
      double t = randomNumber(gen);
      double u = randomNumber(gen);
      if (s + t > 1.0) {  // cut'n fold the cube into a prism

        s = 1.0 - s;
        t = 1.0 - t;
      }

      if (t + u > 1.0) {  // cut'n fold the prism into a tetrahedron

        double tmp = u;
        u = 1.0 - s - t;
        t = 1.0 - tmp;
      }
      else if (s + t + u > 1.0) {
        double tmp = u;
        u = s + t + u - 1.0;
        s = 1 - t - tmp;
      }

      w[3] = 1 - s - t - u;  // a,s,t,u are the barycentric coordinates of the random point.
      w[2] = u;
      w[1] = t;
      w[0] = s;

      Vec3d p;
      p.setZero();
      for (int vi = 0; vi < 4; vi++) {
        p += w[vi] * tetVertices[tets[ei][vi]];
      }

      if (bb.checkInside(p)) {
        insideCounter++;
      }
    }

    weights[ei] = std::max((double)insideCounter / sampleCount, minThreshold);

    int cur = numFinished.fetch_add(1, std::memory_order_relaxed);

    if (cur % 1000 == 0) {
      std::cout << cur << ' ' << std::flush;
    }
  });

  std::cout << std::endl;
}