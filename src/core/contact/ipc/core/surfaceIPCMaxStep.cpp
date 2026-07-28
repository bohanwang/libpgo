/*
copyright to Bohan Wang
*/

#include "surfaceIPCMaxStep.h"

#include "../broadPhase/spatialHashGrid.h"
#include "../geometry/ipcCCD.h"

#include <tbb/blocked_range.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>

#include <algorithm>
#include <functional>
#include <vector>

namespace pgo {
namespace Contact {
namespace CIPC {
using namespace pgo::EigenSupport;

double SurfaceIPCMaxStep::compute(
  const SurfaceIPCTopology &topology,
  EigenSupport::ConstRefVecXd x,
  EigenSupport::ConstRefVecXd dx,
  double dhat,
  double slackness) const
{
  const VXd pos = VXd(x);

  auto getV = [&](int i) -> V3d {
    return pos.segment<3>(3 * i);
  };
  auto getdV = [&](int i) -> V3d {
    return dx.segment<3>(3 * i);
  };

  int nTri = (int)topology.triangles.size();
  int nEdge = (int)topology.edges.size();

  std::vector<SpatialHashGrid::AABB> vertBox(topology.numVerts);
  std::vector<SpatialHashGrid::AABB> triBox(nTri);
  std::vector<SpatialHashGrid::AABB> edgeBox(nEdge);
  double cellSize = 0.0;

  {
    // Build swept-volume AABBs (parallel)
    tbb::parallel_for(tbb::blocked_range<int>(0, topology.numVerts),
      [&](const tbb::blocked_range<int> &r) {
        for (int vi = r.begin(); vi < r.end(); ++vi) {
          V3d p0 = getV(vi), p1 = p0 + getdV(vi);
          vertBox[vi].init(p0, 0.0);
          vertBox[vi].expand(p1);
        }
      });

    tbb::parallel_for(tbb::blocked_range<int>(0, nTri),
      [&](const tbb::blocked_range<int> &r) {
        for (int fi = r.begin(); fi < r.end(); ++fi) {
          auto &tri = topology.triangles[fi];
          V3d v0 = getV(tri[0]), v1 = getV(tri[1]), v2 = getV(tri[2]);
          V3d d0 = getdV(tri[0]), d1 = getdV(tri[1]), d2 = getdV(tri[2]);
          triBox[fi].init(v0, 0.0);
          triBox[fi].expand(v1);
          triBox[fi].expand(v2);
          triBox[fi].expand(v0 + d0);
          triBox[fi].expand(v1 + d1);
          triBox[fi].expand(v2 + d2);
        }
      });

    tbb::parallel_for(tbb::blocked_range<int>(0, nEdge),
      [&](const tbb::blocked_range<int> &r) {
        for (int ei = r.begin(); ei < r.end(); ++ei) {
          V3d a0 = getV(topology.edges[ei][0]), a1 = getV(topology.edges[ei][1]);
          V3d da0 = getdV(topology.edges[ei][0]), da1 = getdV(topology.edges[ei][1]);
          edgeBox[ei].init(a0, 0.0);
          edgeBox[ei].expand(a1);
          edgeBox[ei].expand(a0 + da0);
          edgeBox[ei].expand(a1 + da1);
        }
      });

    // Cell size: average swept-AABB diagonal (parallel reduction)
    double avgBoxDiag = tbb::parallel_reduce(
      tbb::blocked_range<int>(0, nTri), 0.0,
      [&](const tbb::blocked_range<int> &r, double sum) {
        for (int fi = r.begin(); fi < r.end(); ++fi)
          sum += (triBox[fi].hi - triBox[fi].lo).norm();
        return sum;
      },
      std::plus<double>());
    cellSize = nTri > 0 ? std::max(avgBoxDiag / nTri, 1e-6) : std::max(1e-6, dhat);
  }

  double alpha = 1.0;

  // --- PT CCD: insert triangles (serial), query with vertices (parallel reduce) ---
  {
    SpatialHashGrid triHash(nTri);
    triHash.setCellSize(cellSize);
    for (int fi = 0; fi < nTri; ++fi)
      triHash.insert(triBox[fi], fi);

    tbb::enumerable_thread_specific<std::vector<int>> tls_visited(
      [nTri]() { return std::vector<int>(nTri, 0); });
    tbb::enumerable_thread_specific<std::vector<int>> tls_candidates;

    // {
    //   double alpha = 1.0;
    //   auto &visited = tls_visited.local();
    //   auto &candidates = tls_candidates.local();

    //   for (int vi = 0; vi < topology.numVerts; ++vi) {
    //     candidates.clear();
    //     triHash.query(vertBox[vi], -1, visited, vi + 1, candidates);

    //     for (int fi : candidates) {
    //       auto &tri = topology.triangles[fi];
    //       if (vi == tri[0] || vi == tri[1] || vi == tri[2])
    //         continue;
    //       if (!vertBox[vi].overlaps(triBox[fi]))
    //         continue;

    //       V3d p = getV(vi), dp = getdV(vi);
    //       V3d t0 = getV(tri[0]), dt0 = getdV(tri[0]);
    //       V3d t1 = getV(tri[1]), dt1 = getdV(tri[1]);
    //       V3d t2 = getV(tri[2]), dt2 = getdV(tri[2]);

    //       double toi = ccd::pointTriangleCCD(p, t0, t1, t2,
    //         dp, dt0, dt1, dt2,
    //         0.0, alpha);

    //       if (toi < 1e-3) {
    //         toi = ccd::pointTriangleCCD(p, t0, t1, t2,
    //           dp, dt0, dt1, dt2,
    //           0.0, 1.0);

    //         V3d q = p + dp;
    //         V3d s0 = t0 + dt0;
    //         V3d s1 = t1 + dt1;
    //         V3d s2 = t2 + dt2;

    //         std::cout << p[0] << ',' << p[1] << ',' << p[2] << std::endl;
    //         std::cout << q[0] << ',' << q[1] << ',' << q[2] << std::endl;

    //         std::cout << t0[0] << ',' << t0[1] << ',' << t0[2] << std::endl;
    //         std::cout << t1[0] << ',' << t1[1] << ',' << t1[2] << std::endl;
    //         std::cout << t2[0] << ',' << t2[1] << ',' << t2[2] << std::endl;

    //         std::cout << s0[0] << ',' << s0[1] << ',' << s0[2] << std::endl;
    //         std::cout << s1[0] << ',' << s1[1] << ',' << s1[2] << std::endl;
    //         std::cout << s2[0] << ',' << s2[1] << ',' << s2[2] << std::endl;

    //         std::ofstream outfile("/home/user/code/libpgo-private/data/deb.obj");
    //         V3d e(1e-5, 1e-5, 1e-5);

    //         // Triangle: p, p+dp, p+e
    //         auto wv = [&](const V3d &v) {
    //           outfile << "v " << v[0] << " " << v[1] << " " << v[2] << "\n";
    //         };
    //         wv(p);       // 1
    //         wv(p + dp);  // 2
    //         wv(p + e);   // 3
    //         outfile << "f 1 2 3\n";

    //         // Prism: base t0,t1,t2 / top t0+dt0,t1+dt1,t2+dt2
    //         wv(t0);                    // 4
    //         wv(t1);                    // 5
    //         wv(t2);                    // 6
    //         wv(t0 + dt0);              // 7
    //         wv(t1 + dt1);              // 8
    //         wv(t2 + dt2);              // 9
    //         outfile << "f 4 5 6\n";    // bottom
    //         outfile << "f 7 9 8\n";    // top (flipped for outward normal)
    //         outfile << "f 4 5 8 7\n";  // side 01
    //         outfile << "f 5 6 9 8\n";  // side 12
    //         outfile << "f 6 4 7 9\n";  // side 20
    //       }

    //       if (toi < alpha)
    //         alpha = toi * slackness;
    //     }
    //   }
    //   std::cout << "!" << alpha << std::endl;
    // }

    alpha = tbb::parallel_reduce(
      tbb::blocked_range<int>(0, topology.numVerts), 1.0,
      [&](const tbb::blocked_range<int> &range, double localAlpha) {
        auto &visited = tls_visited.local();
        auto &candidates = tls_candidates.local();

        for (int vi = range.begin(); vi < range.end(); ++vi) {
          candidates.clear();
          triHash.query(vertBox[vi], -1, visited, vi + 1, candidates);

          for (int fi : candidates) {
            auto &tri = topology.triangles[fi];
            if (vi == tri[0] || vi == tri[1] || vi == tri[2])
              continue;
            if (!vertBox[vi].overlaps(triBox[fi]))
              continue;

            V3d p = getV(vi), dp = getdV(vi);
            V3d t0 = getV(tri[0]), dt0 = getdV(tri[0]);
            V3d t1 = getV(tri[1]), dt1 = getdV(tri[1]);
            V3d t2 = getV(tri[2]), dt2 = getdV(tri[2]);

            double toi = ccd::pointTriangleCCD(p, t0, t1, t2,
              dp, dt0, dt1, dt2,
              0.0, localAlpha);
            if (toi < localAlpha)
              localAlpha = toi * slackness;
          }
        }
        return localAlpha;
      },
      [](double a, double b) { return std::min(a, b); });
  }

  // --- EE CCD: insert edges (serial), query with edges (parallel reduce) ---
  {
    SpatialHashGrid edgeHash(nEdge);
    edgeHash.setCellSize(cellSize);
    for (int ei = 0; ei < nEdge; ++ei)
      edgeHash.insert(edgeBox[ei], ei);

    tbb::enumerable_thread_specific<std::vector<int>> tls_visited(
      [nEdge]() { return std::vector<int>(nEdge, 0); });
    tbb::enumerable_thread_specific<std::vector<int>> tls_candidates;

    alpha = tbb::parallel_reduce(
      tbb::blocked_range<int>(0, nEdge), alpha,
      [&](const tbb::blocked_range<int> &range, double localAlpha) {
        auto &visited = tls_visited.local();
        auto &candidates = tls_candidates.local();

        for (int ei = range.begin(); ei < range.end(); ++ei) {
          candidates.clear();
          edgeHash.query(edgeBox[ei], ei, visited, ei + 1, candidates);

          int a0 = topology.edges[ei][0], a1 = topology.edges[ei][1];
          for (int ej : candidates) {
            if (ej <= ei)
              continue;

            int b0 = topology.edges[ej][0], b1 = topology.edges[ej][1];
            if (a0 == b0 || a0 == b1 || a1 == b0 || a1 == b1)
              continue;
            if (!edgeBox[ei].overlaps(edgeBox[ej]))
              continue;

            V3d va0 = getV(a0), va1 = getV(a1);
            V3d vb0 = getV(b0), vb1 = getV(b1);
            V3d da0 = getdV(a0), da1 = getdV(a1);
            V3d db0 = getdV(b0), db1 = getdV(b1);

            double toi = ccd::edgeEdgeCCD(va0, va1, vb0, vb1,
              da0, da1, db0, db1,
              0.0, localAlpha);
            if (toi < localAlpha)
              localAlpha = toi * slackness;
          }
        }
        return localAlpha;
      },
      [](double a, double b) { return std::min(a, b); });
  }

  return alpha;
}

}  // namespace CIPC
}  // namespace Contact
}  // namespace pgo
