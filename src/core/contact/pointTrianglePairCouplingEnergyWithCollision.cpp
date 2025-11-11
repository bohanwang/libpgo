#include "pointTrianglePairCouplingEnergyWithCollision.h"

#include "triMeshGeo.h"
#include "triMeshNeighbor.h"
#include "pgoLogging.h"
#include "basicAlgorithms.h"
#include "triangle.h"
#include "geometryQuery.h"
#include "triMeshPseudoNormal.h"
#include "EigenSupport.h"
#include "automaticDifferentiation_autodiff.h"

#include <tbb/parallel_for.h>
#include <tbb/cache_aligned_allocator.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/spin_mutex.h>
#include <tbb/concurrent_vector.h>
#include <tbb/concurrent_set.h>

#include <numeric>
#include <iomanip>

using namespace pgo;
using namespace pgo::Contact;
using namespace pgo::Mesh;
using namespace pgo::BasicAlgorithms;

namespace ES = pgo::EigenSupport;

namespace pgo::Contact
{
class PointTrianglePairCouplingEnergyWithCollisionBuffer
{
public:
  PointTrianglePairCouplingEnergyWithCollisionBuffer():
    energyBuffer(0.0), energyBufferFC(0.0), energyBufferFF(0.0) {}

  tbb::enumerable_thread_specific<double> energyBuffer;
  tbb::enumerable_thread_specific<double> energyBufferFF;
  tbb::enumerable_thread_specific<double> energyBufferFC;
  std::vector<tbb::spin_mutex> entryLock;

  std::vector<TriMeshGeo> surfaceMeshesRuntime;
  std::vector<TriMeshPseudoNormal> surfaceNormals;
};
}  // namespace pgo::Contact

using dual2nd = NonlinearOptimization::AutomaticDifferentiation_autodiff::dual2nd;
using dualV3d = Eigen::Matrix<dual2nd, 3, 1>;
using hclock = std::chrono::high_resolution_clock;

static dual2nd relativeDistanceOnFixedDirections(const dual2nd x[12], const ES::V3d &n, const ES::V3d &w)
{
  dualV3d v0(x[0], x[1], x[2]);
  dualV3d v1(x[3], x[4], x[5]);
  dualV3d v2(x[6], x[7], x[8]);
  dualV3d v3(x[9], x[10], x[11]);

  dualV3d diff = v0 - v1 * w[0] - v2 * w[1] - v3 * w[2];
  dual2nd diffScale = diff.dot(n) - 1e-6;

  return (diffScale * diffScale) * 0.5;
}

PointTrianglePairCouplingEnergyWithCollision::PointTrianglePairCouplingEnergyWithCollision(
  int numPairs_, int numObjects_,
  const std::array<int, 4> *const objectIDs_,
  const std::array<int, 4> *const pointTrianglePairs_,
  const int *const triangleIDs_,
  const int *const objectDOFOffsets_,
  const TriMeshRef *const surfaceMeshes_,
  const Mesh::TriMeshNeighbor *const surfaceMeshNeighbors_,
  const std::vector<int> *const vertexIDToSampleIDs_,
  const ES::SpMatD *const sampleEmbeddingWeights_,
  const std::vector<double> *const sampleWeights_):
  numPairs(numPairs_),
  numObjects(numObjects_),
  objectIDs(objectIDs_), pointTrianglePairs(pointTrianglePairs_), triangleIDs(triangleIDs_),
  objectDOFOffsets(objectDOFOffsets_), surfaceMeshesRest(surfaceMeshes_),
  vertexIDToSampleIDs(vertexIDToSampleIDs_),
  sampleEmbeddingWeights(sampleEmbeddingWeights_), sampleWeights(sampleWeights_)
{
  // tbb::concurrent_vector<ES::TripletD> entriesCC;
  // std::cout << 1 << std::endl;
  // entries.reserve(numPairs * 1000ull);
  tbb::concurrent_set<uint64_t> entriesCC;

  // tbb::enumerable_thread_specific<std::vector<ES::TripletD, tbb::cache_aligned_allocator<ES::TripletD>>> entriesTLS;

  neighboringTriangles.assign(numPairs, std::vector<int>());

  // for (int pi = 0; pi < numPairs; pi++) {
  tbb::parallel_for(0, numPairs, [&](int pi) {
#if 0
    auto &localBuf = entriesTLS.local();
    if (localBuf.size() == 0ull)
      localBuf.reserve(1 << 20);
#else
    auto &localBuf = entriesCC;
#endif
    for (int vi = 0; vi < 4; vi++) {
      int vidi = pointTrianglePairs[pi][vi];
      int mi = objectIDs[pi][vi];

      for (int vj = 0; vj < 4; vj++) {
        int vidj = pointTrianglePairs[pi][vj];
        int mj = objectIDs[pi][vj];

        for (ES::SpMatD::InnerIterator it_i(sampleEmbeddingWeights[mi], vidi * 3); it_i; ++it_i) {
          int actual_vidi = (int)it_i.col() / 3;

          for (ES::SpMatD::InnerIterator it_j(sampleEmbeddingWeights[mj], vidj * 3); it_j; ++it_j) {
            int actual_vidj = (int)it_j.col() / 3;

            for (int i1 = 0; i1 < 3; i1++) {
              for (int i2 = 0; i2 < 3; i2++) {
                int row = objectDOFOffsets[mi] + actual_vidi * 3 + i1;
                int col = objectDOFOffsets[mj] + actual_vidj * 3 + i2;
                // entries.emplace_back(actual_vidi * 3 + i1, actual_vidj * 3 + i2, 0);
                // localBuf.emplace_back(
                //  objectDOFOffsets[mi] + actual_vidi * 3 + i1,
                //  objectDOFOffsets[mj] + actual_vidj * 3 + i2,
                //  1);

                uint64_t rc = (uint64_t)row << 32ull | (uint64_t)col;
                localBuf.emplace(rc);
              }
            }
          }
        }
      }  // vj
    }  // vi

    // for the triangle, we get the neighbors
    int triObjID = getTriangleObjectID(pi);
    int triID = triangleIDs[pi];
    ES::V3i triVtxID = surfaceMeshesRest[triObjID].tri(triID);

    neighboringTriangles[pi].reserve(50);

    for (int vi = 0; vi < 3; vi++) {
      const auto &neighbors = surfaceMeshNeighbors_[triObjID].getVtxNearbyTriangles(triVtxID[vi]);
      for (int ni : neighbors)
        if (ni >= 0)
          neighboringTriangles[pi].emplace_back(ni);
    }

    int count = (int)neighboringTriangles[pi].size();
    for (int ti = 0; ti < count; ti++) {
      ES::V3i n = surfaceMeshNeighbors_[triObjID].getTriangleNeighbors(neighboringTriangles[pi][ti]);
      for (int j = 0; j < 3; j++) {
        int ni = n[j];
        if (ni >= 0)
          neighboringTriangles[pi].emplace_back(ni);
      }
    }

    sortAndDeduplicateWithErase(neighboringTriangles[pi]);

    // put the direct triangle at the first of the array
    for (size_t i = 0; i < neighboringTriangles[pi].size(); i++) {
      if (neighboringTriangles[pi][i] == triID) {
        if (i != 0)
          std::swap(neighboringTriangles[pi][0], neighboringTriangles[pi][i]);

        break;
      }
    }
  });  // contacted pair

  std::vector<ES::TripletD> entries;
#if 0
  size_t totalEntries = 0;
  for (auto it = entriesTLS.begin(); it != entriesTLS.end(); ++it) {
    totalEntries += it->size();
  }
  entries.reserve(totalEntries);
  std::cout << 3 << std::endl;

  for (auto it = entriesTLS.begin(); it != entriesTLS.end(); ++it) {
    entries.insert(entries.end(), it->begin(), it->end());
  }

  std::cout << "entries:" << entries.size() << std::endl;
#else
  // std::cout << "entries before:" << entriesCC.size() << std::endl;

  // ankerl::unordered_dense::set<uint64_t> entriesSet;
  // for (auto it = entriesCC.begin(); it != entriesCC.end(); ++it) {
  //   uint64_t rc = (uint64_t)it->row() << 32ull | (uint64_t)it->col();
  //   entriesSet.emplace(rc);
  // }
  std::cout << "entries after:" << entriesCC.size() << std::endl;

  entries.reserve(entriesCC.size());

  for (const auto &rc : entriesCC) {
    int row = (int)((rc >> 32ull) & 0xffffffffull);
    int col = (int)(rc & 0xffffffffull);
    entries.emplace_back(row, col, 1.0);
  }
#endif

  hessTemplate.resize(objectDOFOffsets[numObjects], objectDOFOffsets[numObjects]);
  hessTemplate.setFromTriplets(entries.begin(), entries.end());

  std::cout << "#nonzeros: " << hessTemplate.nonZeros() << std::endl;

  allDOFs.resize(objectDOFOffsets[numObjects]);
  // allDOFs.resize(sampleEmbeddingWeights[0].rows());
  std::iota(allDOFs.begin(), allDOFs.end(), 0);

  ES::VXd x(objectDOFOffsets[numObjects]);
  x.setZero();

  // ES::SpMatD h0 = hessTemplate;
  ES::buildEntryMap(hessTemplate, entryMap);
  // std::cout << (entryMap.find(std::make_pair(156549, 1995)) != entryMap.end()) << std::endl;

  std::cout << "done" << std::endl;
}

std::tuple<bool, int> PointTrianglePairCouplingEnergyWithCollision::checkContact(ES::ConstRefVecXd x, int pi, ES::V12d &xlocal,
  int checkNeighboringTriangles) const
{
  xlocal.setZero();
  int pointID = pointTrianglePairs[pi][0];
  int pointObjID = getPointObjectID(pi);
  xlocal.segment<3>(0) = computePosition(x, pointObjID, pointID);

  if (checkNeighboringTriangles) {
    double minDistToTriangle = 1e100;
    int closestTriID = neighboringTriangles[pi][0];
    int closestFeature = -1;
    int triObjID = getTriangleObjectID(pi);
    ES::V3d closestPos;
    std::array<ES::V3d, 3> closestTrip;

    for (int trii = 0; trii < (int)neighboringTriangles[pi].size(); trii++) {
      int triID = neighboringTriangles[pi][trii];
      ES::V3i triVtxID = buf->surfaceMeshesRuntime[triObjID].tri(triID);
      std::array<ES::V3d, 3> trip{ ES::V3d::Zero(), ES::V3d::Zero(), ES::V3d::Zero() };

      // compute position
      for (int vi = 0; vi < 3; vi++) {
        int vid = triVtxID[vi];
        int sampleID = vertexIDToSampleIDs[triObjID][vid];
        trip[vi] = computePosition(x, triObjID, sampleID);
      }

      int feature;
      ES::V3d pt, w;
      double dist = getSquaredDistanceToTriangle(asVec3d(xlocal.data()),
        trip[0], trip[1], trip[2], feature, pt, w);

      if (dist < minDistToTriangle) {
        closestTrip = trip;

        minDistToTriangle = dist;
        closestTriID = neighboringTriangles[pi][trii];
        closestFeature = feature;
        closestPos = pt;
      }

      if (trii == 0) {
        xlocal.segment<3>(3) = trip[0];
        xlocal.segment<3>(6) = trip[1];
        xlocal.segment<3>(9) = trip[2];
      }
    }

    ES::V3d n = buf->surfaceNormals[triObjID].getPseudoNormal(buf->surfaceMeshesRuntime[triObjID].triangles().data(), closestTriID, closestFeature);
    ES::V3d dir = xlocal.segment<3>(0) - closestPos;

    return std::make_tuple(dir.dot(n) < 0, neighboringTriangles[pi][0]);
  }
  else {
    int triObjID = getTriangleObjectID(pi);
    int triID = neighboringTriangles[pi][0];
    ES::V3i triVtxID = surfaceMeshesRest[triObjID].tri(triID);
    // compute position
    for (int vi = 0; vi < 3; vi++) {
      int vid = triVtxID[vi];
      int sampleID = vertexIDToSampleIDs[triObjID][vid];
      xlocal.segment<3>(vi * 3 + 3) = computePosition(x, triObjID, sampleID);
      PGO_ALOG(pointTrianglePairs[pi][vi + 1] == sampleID);
    }

    ES::V3d n = (xlocal.segment<3>(6) - xlocal.segment<3>(3)).cross(xlocal.segment<3>(9) - xlocal.segment<3>(3));
    n.normalize();

    // ES::V3d p = xlocal.segment<3>(3) * barycentricWeights[pi][0] +
    //   xlocal.segment<3>(6) * barycentricWeights[pi][1] +
    //   xlocal.segment<3>(9) * barycentricWeights[pi][2];

    // return std::make_tuple(true, neighboringTriangles[pi][0]);
    return std::make_tuple((xlocal.segment<3>(0) - xlocal.segment<3>(3)).dot(n) <= 0, neighboringTriangles[pi][0]);
  }
}

double PointTrianglePairCouplingEnergyWithCollision::frictionPotential(const ES::V12d &x, const ES::V12d &xlast,
  const ES::V3d &n, const ES::V3d &w, double eps, double timestep) const
{
  // compute velocity v = (x_t+1 - x_t) / h
  // (v0 - vbc) - (v0 - vbc).n * n
  // w xdiff - n^T(w xdiff)n
  ES::V12d xdiff = x - xlast;

  Eigen::Matrix<double, 3, 12> W;
  ES::M3d InnT = ES::M3d::Identity() - ES::tensorProduct(n, n);
  W.block<3, 3>(0, 0) = InnT;
  W.block<3, 3>(0, 3) = InnT * -w[0];
  W.block<3, 3>(0, 6) = InnT * -w[1];
  W.block<3, 3>(0, 9) = InnT * -w[2];

  ES::V3d r = W * xdiff;

  // compute relative velocity
  /*ES::V3d vel0 = xdiff.segment<3>(0);
  ES::V3d vel1 = xdiff.segment<3>(3) * w[0] + xdiff.segment<3>(6) * w[1] + xdiff.segment<3>(9) * w[2];
  ES::V3d veldiff = vel0 - vel1;
  ES::V3d veldiff_plane = veldiff - n * n.dot(veldiff);

  ALOG((veldiff_plane2 - veldiff_plane).norm() < 1e-8);*/

  // compute energy
  double d = r.norm();

  if (d < timestep * eps)
    return (d * d * d) * (-1.0 / (3.0 * eps * timestep * eps * timestep)) + (d * d) / (eps * timestep) + eps * timestep / 3;
  // when velmag == timestep * eps
  // the function value is timestep * eps
  // the derivative is 1
  else
    return d;
}

void PointTrianglePairCouplingEnergyWithCollision::frictionPotentialGradient(const ES::V12d &x, const ES::V12d &xlast,
  const ES::V3d &n, const ES::V3d &w, double eps, double timestep, ES::V12d &grad) const
{
  // compute volecity v = (x_t+1 - x_t) / h
  ES::V12d xdiff = x - xlast;

  Eigen::Matrix<double, 3, 12> W;
  ES::M3d InnT = ES::M3d::Identity() - ES::tensorProduct(n, n);
  W.block<3, 3>(0, 0) = InnT;
  W.block<3, 3>(0, 3) = InnT * -w[0];
  W.block<3, 3>(0, 6) = InnT * -w[1];
  W.block<3, 3>(0, 9) = InnT * -w[2];

  ES::V3d r = W * xdiff;

  // compute energy
  double d = r.norm();
  double c1 = 0;

  if (d < timestep * eps)
    c1 = (-d / (eps * timestep * eps * timestep) + 2 / (eps * timestep));
  else
    c1 = 1.0 / d;

  ES::V3d grad_r = c1 * r;
  if (grad_r.norm() > 1.1) {
    std::cerr << "grad_r > 1" << std::endl;
    std::cerr << "grad_r: " << std::setprecision(16) << grad_r[0] << ',' << grad_r[1] << ',' << grad_r[2] << '\n'
              << "|grad_r|: " << grad_r.norm() << '\n'
              << eps << ',' << timestep << '\n'
              << c1 << '\n'
              << "r: " << r[0] << ',' << r[1] << ',' << r[2] << '\n'
              << "|r|: " << r.norm() << std::endl;

    exit(1);
  }

  grad = W.transpose() * grad_r;
}

void PointTrianglePairCouplingEnergyWithCollision::frictionPotentialHessian(const ES::V12d &x, const ES::V12d &xlast,
  const ES::V3d &n, const ES::V3d &w, double eps, double timestep, ES::M12d &hess) const
{
  Eigen::Matrix<double, 3, 12> W;
  ES::M3d InnT = ES::M3d::Identity() - ES::tensorProduct(n, n);
  W.block<3, 3>(0, 0) = InnT;
  W.block<3, 3>(0, 3) = InnT * -w[0];
  W.block<3, 3>(0, 6) = InnT * -w[1];
  W.block<3, 3>(0, 9) = InnT * -w[2];

  ES::V12d xdiff = x - xlast;
  ES::V3d r = W * xdiff;

  double d = r.norm();

  ES::M3d hess_r;
  if (d < timestep * eps) {
    if (d < 1e-16) {
      hess_r.setZero();
    }
    else {
      hess_r = -ES::tensorProduct(r, r) / d / (eps * eps * timestep * timestep);
    }

    hess_r += ES::M3d::Identity() * (2 / (timestep * eps) - d / (eps * eps * timestep * timestep));
  }
  else {
    // grad = r (rTr)^-0.5
    // hess = I (rTr)^-0.5 + r * -(rTr)^-1.5 rT
    hess_r = ES::M3d::Identity() / d - ES::tensorProduct(r, r) / (d * d * d);
  }

  if (1) {
    Eigen::SelfAdjointEigenSolver<ES::M3d> eig(hess_r);
    ES::V3d eigv = eig.eigenvalues().cwiseMax(0.0);
    ES::M3d eigb = eig.eigenvectors();
    hess_r = eigb * eigv.asDiagonal() * eigb.transpose();
  }

  hess = W.transpose() * hess_r * W;
}

double PointTrianglePairCouplingEnergyWithCollision::func(ES::ConstRefVecXd x) const
{
  for (auto it = buf->energyBuffer.begin(); it != buf->energyBuffer.end(); ++it)
    *it = 0;

  for (auto it = buf->energyBufferFF.begin(); it != buf->energyBufferFF.end(); ++it)
    *it = 0;

  for (auto it = buf->energyBufferFC.begin(); it != buf->energyBufferFC.end(); ++it)
    *it = 0;

  tbb::parallel_for(0, numObjects, [&](int mi) {
    tbb::parallel_for(0, buf->surfaceMeshesRuntime[mi].numVertices(), [&](int vi) {
      int sampleID = vertexIDToSampleIDs[mi][vi];
      ES::V3d p = computePosition(x, mi, sampleID);
      buf->surfaceMeshesRuntime[mi].pos(vi) = p;
    });
    buf->surfaceNormals[mi].updateVertexPositions(buf->surfaceMeshesRuntime[mi]);
  });

  // ES::VXd p(x.size());
  // for (int i = 0; i < p.size() / 3; i++) {
  //   ES::V3d p3;
  //   toPosFunc(x.segment<3>(i * 3), p3, i * 3);
  //   p.segment<3>(i * 3) = p3;
  // }

  // ES::VXd s(sampleEmbeddingWeights[0].rows());
  // ES::Mv(sampleEmbeddingWeights[0], p, s);

  auto localEnergy = [&](int pi) {
    if (contactStatus[pi] == 0)
      return;

    ES::V12d xlocal;
    bool inContact = false;

    inContact = std::get<0>(checkContact(x, pi, xlocal, 1));

    if (inContact == false)
      return;

    // for (int vi = 0; vi < 4; vi++) {
    //   xlocal.segment<3>(vi * 3) = s.segment<3>(pointTrianglePairs[pi][vi] * 3);
    // }

    double eng = 0;
    double eng_fc = 0, eng_ff = 0;
    eng = (xlocal.dot(hessianBlocks[pi] * xlocal) * 0.5 + gradientBlocks[pi].dot(xlocal)) * coeff;

    eng_fc = eng;

    // friction forces
    // const dual2nd x[12], const ES::V12d &xlast, const ES::V3d &n, const ES::V3d &w, double eps, double timestep
    if (frictionCoeff > 0 && frictionContactForceMag > 0) {
      ES::V12d fnow = hessianBlocks[pi] * xlocal + gradientBlocks[pi];
      double fc = fnow.head<3>().norm() * coeff;
      // double fc = contactForceMags[pi] * coeff;

      ES::V12d xlast;
      for (int j = 0; j < 4; j++) {
        xlast.segment<3>(j * 3) = computeLastPosition(x, objectIDs[pi][j], pointTrianglePairs[pi][j]);
      }

      double frictionEnergy = frictionPotential(xlocal, xlast, normals[pi], barycentricWeights[pi], eps, timestep);
      frictionEnergy *= frictionContactForceMag * frictionCoeff * fc;
      eng += frictionEnergy;
      eng_ff = frictionEnergy;

      // ES::V12d xdiff = xlocal - xlast;
      // if (xdiff.norm() > 1e-4) {
      //   LGI << "pi: " << pi << '\n'
      //       << "xlast:\n"
      //       << xlast << '\n';

      //  for (double alpha = 1.0; alpha < 3.0; alpha += 0.1) {
      //    ES::V12d xcur = xdiff * alpha + xlast;
      //    double ee = frictionPotential(xcur, xlast, normals[pi], barycentricWeights[pi], eps, timestep);
      //    LGI
      //        // << "xlast:\n"
      //        // << xlast << '\n'
      //        // << "xcur:\n"
      //        // << xlocal << '\n'
      //        << "alpha: " << alpha << '\n'
      //        << "Eng:\n"
      //        << ee;
      //  }
      //}
    }

    if (sampleWeights) {
      eng *= sampleWeights[objectIDs[pi][0]][pointTrianglePairs[pi][0]];
      eng_ff *= sampleWeights[objectIDs[pi][0]][pointTrianglePairs[pi][0]];
      eng_fc *= sampleWeights[objectIDs[pi][0]][pointTrianglePairs[pi][0]];
    }

    buf->energyBuffer.local() += eng;
    buf->energyBufferFF.local() += eng_ff;
    buf->energyBufferFC.local() += eng_fc;
  };

  // for (int pi = 0; pi < numPairs; pi++) {
  tbb::parallel_for(0, numPairs, [&](int pi) {
    localEnergy(pi);
  });
  // }

  double energyAll = std::accumulate(buf->energyBuffer.begin(), buf->energyBuffer.end(), 0.0);

  eng_fn = std::accumulate(buf->energyBufferFC.begin(), buf->energyBufferFC.end(), 0.0);
  eng_ff = std::accumulate(buf->energyBufferFF.begin(), buf->energyBufferFF.end(), 0.0);

  return energyAll;
}

// p = W x
// E = 1/2 pT H p + pT b
// dg/dp = H p + b
// dg/dx = WT (Hp + b)

void PointTrianglePairCouplingEnergyWithCollision::gradient(ES::ConstRefVecXd x, ES::RefVecXd grad) const
{
  grad.setZero();
  fn.setZero(grad.size());
  ff.setZero(grad.size());

  tbb::parallel_for(0, numObjects, [&](int mi) {
    tbb::parallel_for(0, buf->surfaceMeshesRuntime[mi].numVertices(), [&](int vi) {
      int sampleID = vertexIDToSampleIDs[mi][vi];
      ES::V3d p = computePosition(x, mi, sampleID);
      buf->surfaceMeshesRuntime[mi].pos(vi) = p;
    });
    buf->surfaceNormals[mi].updateVertexPositions(buf->surfaceMeshesRuntime[mi]);
  });

  auto localGradFunc = [&](int pi) {
    if (contactStatus[pi] == 0)
      return;

    ES::V12d xlocal;
    bool inContact = false;
    std::array<int, 4> sids = pointTrianglePairs[pi];

    inContact = std::get<0>(checkContact(x, pi, xlocal, 1));
    if (inContact == false)
      return;

    ES::V12d dE_dsamplex, fnlocal, fflocal;
    dE_dsamplex = hessianBlocks[pi] * xlocal + gradientBlocks[pi];

    dE_dsamplex *= coeff;
    fnlocal = dE_dsamplex;
    // std::cout << pi << "f:" << dE_dsamplex.head<3>().norm();

    fflocal.setZero();

    // friction
    if (frictionContactForceMag > 0 && frictionCoeff > 0) {
      ES::V12d fnow = hessianBlocks[pi] * xlocal + gradientBlocks[pi];
      double fc = fnow.head<3>().norm() * coeff;
      // double fc = contactForceMags[pi] * coeff;

      ES::V12d xlast;
      for (int j = 0; j < 4; j++) {
        xlast.segment<3>(j * 3) = computeLastPosition(x, objectIDs[pi][j], pointTrianglePairs[pi][j]);
      }

      ES::V12d gradLocal;
      frictionPotentialGradient(xlocal, xlast, normals[pi], barycentricWeights[pi], eps, timestep, gradLocal);
      gradLocal *= frictionContactForceMag * frictionCoeff * fc;

      fflocal = gradLocal;
      dE_dsamplex += gradLocal;
      // std::cout << " ff:" << gradLocal.head<3>().norm() << '\n';
    }

    if (sampleWeights) {
      dE_dsamplex *= sampleWeights[objectIDs[pi][0]][pointTrianglePairs[pi][0]];
      fflocal *= sampleWeights[objectIDs[pi][0]][pointTrianglePairs[pi][0]];
      fnlocal *= sampleWeights[objectIDs[pi][0]][pointTrianglePairs[pi][0]];
    }

    for (int vi = 0; vi < 4; vi++) {
      int local_vid = sids[vi];
      int mi = objectIDs[pi][vi];
      ES::V3d grad_trix_i = dE_dsamplex.segment<3>(vi * 3);

      ES::V3d fn_trix_i = fnlocal.segment<3>(vi * 3);
      ES::V3d ff_trix_i = fflocal.segment<3>(vi * 3);

      for (ES::SpMatD::InnerIterator it(sampleEmbeddingWeights[mi], local_vid * 3); it; ++it) {
        int globalStart = objectDOFOffsets[mi] + (int)it.col();
        ES::V3d gradx_j = grad_trix_i * it.value();

        ES::V3d fn_j = fn_trix_i * it.value();
        ES::V3d ff_j = ff_trix_i * it.value();
        // LGI << "\n"
        //     << sampleEmbeddingWeights[mi].block(it.row(), it.col(), 3, 3) << '\n';

        buf->entryLock[globalStart / 3].lock();

        grad.segment<3>(globalStart) += gradx_j;
        fn.segment<3>(globalStart) += fn_j;
        ff.segment<3>(globalStart) += ff_j;

        buf->entryLock[globalStart / 3].unlock();
      }
      // LGI << "zzzz";
    }
  };

  tbb::parallel_for(0, numPairs, [&](int pi) {
    // for (int pi = 0; pi < numPairs; pi++) {
    localGradFunc(pi);
  });
}

void PointTrianglePairCouplingEnergyWithCollision::hessian(ES::ConstRefVecXd x, ES::SpMatD &hess) const
{
  // LGI << "#nonzeros: " << hess.nonZeros();

  memset(hess.valuePtr(), 0, sizeof(double) * hess.nonZeros());

  auto localHessian = [&](int pi) {
    if (contactStatus[pi] == 0)
      return;

    ES::V12d xlocal;
    bool inContact = false;
    std::array<int, 4> sids = pointTrianglePairs[pi];

    inContact = std::get<0>(checkContact(x, pi, xlocal, 1));

    if (inContact == false)
      return;

    ES::M12d K;

    K = hessianBlocks[pi];
    K *= coeff;

    // friction
    if (frictionContactForceMag > 0 && frictionCoeff > 0) {
      ES::V12d fnow = hessianBlocks[pi] * xlocal + gradientBlocks[pi];
      double fc = fnow.head<3>().norm() * coeff;
      // double fc = contactForceMags[pi] * coeff;

      ES::V12d xlast;
      for (int j = 0; j < 4; j++) {
        xlast.segment<3>(j * 3) = computeLastPosition(x, objectIDs[pi][j], pointTrianglePairs[pi][j]);
      }

      ES::M12d K1;
      frictionPotentialHessian(xlocal, xlast, normals[pi], barycentricWeights[pi], eps, timestep, K1);
      K1 *= frictionContactForceMag * frictionCoeff * fc;
      K += K1;
    }

    if (sampleWeights)
      K *= sampleWeights[objectIDs[pi][0]][pointTrianglePairs[pi][0]];

    for (int vi = 0; vi < 4; vi++) {
      int vidi = sids[vi];
      int mi = objectIDs[pi][vi];

      for (int vj = 0; vj < 4; vj++) {
        int vidj = sids[vj];
        int mj = objectIDs[pi][vj];

        ES::M3d Kij = K.block<3, 3>(vi * 3, vj * 3);

        for (ES::SpMatD::InnerIterator it_i(sampleEmbeddingWeights[mi], vidi * 3); it_i; ++it_i) {
          int actual_vidi = (int)it_i.col() / 3;

          for (ES::SpMatD::InnerIterator it_j(sampleEmbeddingWeights[mj], vidj * 3); it_j; ++it_j) {
            int actual_vidj = (int)it_j.col() / 3;

            ES::M3d WiTKijWj = Kij * it_i.value() * it_j.value();

            for (int ii = 0; ii < 3; ii++) {
              for (int jj = 0; jj < 3; jj++) {
                // std::ptrdiff_t offset = ES::FindEntryOffset(hess, objectDOFOffsets[mi] + actual_vidi * 3 + ii, objectDOFOffsets[mj] + actual_vidj * 3 + jj);
                // std::ptrdiff_t offset1 = ES::FindEntryOffset(hessTemplate, objectDOFOffsets[mi] + actual_vidi * 3 + ii, objectDOFOffsets[mj] + actual_vidj * 3 + jj);

                auto itt = entryMap.find(std::make_pair(objectDOFOffsets[mi] + actual_vidi * 3 + ii, objectDOFOffsets[mj] + actual_vidj * 3 + jj));
                std::ptrdiff_t offset = -1;
                if (itt != entryMap.end()) {
                  offset = itt->second;
                }

                if (offset < 0) {
                  // std::cerr << "-";
                  std::cerr << "Critical error: entry not found." << std::endl;
                  std::cout << pi << ',' << (objectDOFOffsets[mi] + actual_vidi * 3 + ii) << ',' << (objectDOFOffsets[mj] + actual_vidj * 3 + jj) << std::endl;
                  abort();
                }
                else {
                  buf->entryLock[objectDOFOffsets[mj] + actual_vidj * 3 + jj].lock();

                  hess.valuePtr()[offset] += WiTKijWj(ii, jj);

                  buf->entryLock[objectDOFOffsets[mj] + actual_vidj * 3 + jj].unlock();
                }
              }
            }
          }
        }
      }  // vj
    }
  };

  // for (int pi = 0; pi < numPairs; pi++) {
  tbb::parallel_for(0, numPairs, [&](int pi) {
    localHessian(pi);
  });

  // std::cerr << std::endl;
}

PointTrianglePairCouplingEnergyWithCollisionBuffer *PointTrianglePairCouplingEnergyWithCollision::allocateBuffer() const
{
  PointTrianglePairCouplingEnergyWithCollisionBuffer *b = new PointTrianglePairCouplingEnergyWithCollisionBuffer;

  // b->entryLock = std::vector<tbb::spin_mutex>(std::max(hessTemplate.nonZeros(), (ES::IDX)objectDOFOffsets[numObjects]));
  b->entryLock = std::vector<tbb::spin_mutex>((ES::IDX)objectDOFOffsets[numObjects]);
  b->surfaceNormals.resize(numObjects);
  b->surfaceMeshesRuntime.resize(numObjects);

  for (size_t i = 0; i < b->surfaceNormals.size(); i++) {
    b->surfaceMeshesRuntime[i] = TriMeshGeo(surfaceMeshesRest[i]);
    b->surfaceNormals[i].buildPseudoNormals(surfaceMeshesRest[i]);
  }

  return b;
}

void PointTrianglePairCouplingEnergyWithCollision::freeBuffer(PointTrianglePairCouplingEnergyWithCollisionBuffer *buf) const
{
  delete buf;
}

void PointTrianglePairCouplingEnergyWithCollision::computeClosestPosition(const double *x_)
{
  if (hessianBlocks.size() == 0) {
    hessianBlocks.assign(numPairs, ES::M12d::Zero());
  }

  if (gradientBlocks.size() == 0) {
    gradientBlocks.assign(numPairs, ES::V12d::Zero());
  }

  if (barycentricWeights.size() == 0) {
    barycentricWeights.assign(numPairs, ES::V3d::Zero());
  }

  if (normals.size() == 0) {
    normals.assign(numPairs, ES::V3d::Zero());
  }

  if (closestPositions.size() == 0) {
    closestPositions.assign(numPairs, ES::V3d::Zero());
  }

  if (contactStatus.size() == 0) {
    contactStatus.assign(numPairs, 1);
  }

  if (contactDistances.size() == 0) {
    contactDistances.assign(numPairs, 0);
  }

  if (contactForceMags.size() == 0) {
    contactForceMags.assign(numPairs, 0);
  }

  ES::Mp<const ES::VXd> x(x_, objectDOFOffsets[numObjects]);

  // std::vector<ES::V3d> ppp;
  // std::vector<ES::V3i> ttt;
  ////std::vector<int> depidx;
  // for (int pi = 0; pi < numPairs; pi++) {
  //   //  LG_ << "pi: " << pi << '\n';
  //   //  for (int j = 0; j < 4; j++) {
  //   //    int sid = pointTrianglePairs[pi][j];
  //   //    for (ES::SpMatD::InnerIterator it(sampleEmbeddingWeights[0], sid * 3); it; ++it) {
  //   //      depidx.push_back((int)it.col() / 3);
  //   //      LG_ << it.col() / 3 << ',';
  //   //    }
  //   //    LG_ << '\n';
  //   //  }
  //   //}
  //   //sortAndDeduplicate(depidx);
  //   ES::V12d xlocal;
  //   checkContact(x, pi, xlocal, 0, 0);

  //  int base = (int)ppp.size();

  //  ppp.emplace_back(xlocal.data());
  //  ppp.emplace_back(xlocal.data() + 3);
  //  ppp.emplace_back(xlocal.data() + 6);
  //  ppp.emplace_back(xlocal.data() + 9);

  //  ttt.emplace_back(base, base + 1, base + 2);
  //  ttt.emplace_back(base, base + 2, base + 3);
  //  ttt.emplace_back(base, base + 3, base + 1);
  //  ttt.emplace_back(base + 1, base + 2, base + 3);
  //}

  // ObjMesh(ppp, ttt).save("zzz.obj");

  // for (int pi = 0; pi < numPairs; pi++) {
  tbb::parallel_for(0, numPairs, [&](int pi) {
    ES::V12d xlocal;
    bool inContact = true;
    inContact = std::get<0>(checkContact(x, pi, xlocal, 1));
    contactStatus[pi] = 1;  // inContact ? 1 : 0;

    ES::V3d va = xlocal.segment<3>(3);
    ES::V3d vb = xlocal.segment<3>(6);
    ES::V3d vc = xlocal.segment<3>(9);
    ES::V3d v0 = xlocal.segment<3>(0);
    int feature;
    ES::V3d w, pt;
    double d2 = getSquaredDistanceToTriangle(v0, va, vb, vc, feature, pt, w);

    barycentricWeights[pi] = ES::V3d(w[0], w[1], w[2]);
    closestPositions[pi] = ES::V3d(pt[0], pt[1], pt[2]);

    normals[pi] = (xlocal.segment<3>(6) - xlocal.segment<3>(3)).cross(xlocal.segment<3>(9) - xlocal.segment<3>(3));
    normals[pi].normalize();

    ES::M12d K;
    ES::V12d f;
    NonlinearOptimization::AutomaticDifferentiation_autodiff::computeHessian<12>(
      relativeDistanceOnFixedDirections, xlocal, K, normals[pi], barycentricWeights[pi]);

    ES::V12d zero12 = ES::V12d::Zero();
    NonlinearOptimization::AutomaticDifferentiation_autodiff::computeGradient<12>(
      relativeDistanceOnFixedDirections, zero12, f, normals[pi], barycentricWeights[pi]);

    double eng = NonlinearOptimization::AutomaticDifferentiation_autodiff::computeEnergy<12>(
      relativeDistanceOnFixedDirections, xlocal, normals[pi], barycentricWeights[pi]);

    hessianBlocks[pi] = K;
    gradientBlocks[pi] = f;

    ES::V12d fc = K * xlocal + f;
    contactForceMags[pi] = f.segment<3>(0).norm();

    // double f1 = f.segment<3>(0).dot(normals[pi]);
    // LGI << contactForceMags[pi] << ',' << f1;

    // double e1 = (K * xlocal).dot(xlocal) * 0.5 + f.dot(xlocal);
    // double e2 = eng;
    // LGI << e1 << ',' << e2;
  });
}

ES::V3d PointTrianglePairCouplingEnergyWithCollision::computePosition(ES::ConstRefVecXd x, int objID, int vid) const
{
  ES::V3d p;
  p.setZero();

  for (ES::SpMatD::InnerIterator it(sampleEmbeddingWeights[objID], vid * 3); it; ++it) {
    int offset = objectDOFOffsets[objID] + (int)it.col();
    ES::V3d x3 = x.segment<3>(offset);
    ES::V3d p3;
    toPosFunc(x3, p3, offset);

    p += p3 * it.value();
  }

  return p;
}

ES::V3d PointTrianglePairCouplingEnergyWithCollision::computeLastPosition(ES::ConstRefVecXd x, int objID, int vid) const
{
  ES::V3d p;
  p.setZero();

  for (ES::SpMatD::InnerIterator it(sampleEmbeddingWeights[objID], vid * 3); it; ++it) {
    int offset = objectDOFOffsets[objID] + (int)it.col();
    ES::V3d x3 = x.segment<3>(offset);
    ES::V3d p3;
    toLastPosFunc(x3, p3, offset);

    p += p3 * it.value();
  }

  return p;
}
