#include "surfaceSmoothnessAbsoluteMeanCurvature.h"
#include "pgoLogging.h"
#include "triMeshNeighbor.h"
#include "basicAlgorithms.h"
#include "EigenSupport.h"

#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>
#include <tbb/spin_mutex.h>

#include <numeric>

using namespace pgo::PredefinedPotentialEnergies;
namespace ES = pgo::EigenSupport;

namespace pgo
{
namespace PredefinedPotentialEnergies
{
struct SurfaceSmoothnessAbsoluteMeanCurvatureBuf
{
  SurfaceSmoothnessAbsoluteMeanCurvatureBuf(int numV):
    locks(numV) {}
  std::vector<tbb::spin_mutex> locks;
};
}  // namespace PredefinedPotentialEnergies
}  // namespace pgo

SurfaceSmoothnessAbsoluteMeanCurvature::SurfaceSmoothnessAbsoluteMeanCurvature(const EigenSupport::VXd &restp, const Mesh::TriMeshGeo &meshIn, double eps):
  restPositions(restp), mesh(meshIn)
{
  std::vector<std::vector<int>> vertexNeigboringVertices(mesh.numVertices());
  std::vector<std::vector<int>> vertexNeigboringTriangles(mesh.numVertices());
  try {
    Mesh::TriMeshNeighbor neighbor(mesh);

    // for (int vi = 0; vi < surfaceIn.numVertices(); vi++) {
    tbb::parallel_for(0, mesh.numVertices(), [&](int vi) {
      vertexNeigboringTriangles[vi] = neighbor.getVtxNearbyTriangles(vi);
      vertexNeigboringVertices[vi] = neighbor.getVtxNearbyVertices(vi, mesh);
      BasicAlgorithms::sortAndDeduplicate(vertexNeigboringVertices[vi]);
      auto it = std::lower_bound(vertexNeigboringVertices[vi].begin(), vertexNeigboringVertices[vi].end(), vi);
      if (it != vertexNeigboringVertices[vi].end() && *it == vi) {
        vertexNeigboringVertices[vi].erase(it);
      }
    });
  }
  catch (std::exception &) {
    SPDLOG_LOGGER_ERROR(pgo::Logging::lgr(), "Input mesh is not manifold.");
    throw std::invalid_argument("Input mesh is not manifold.");
  }

  using EdgeIndex = std::pair<int, int>;
  std::map<EdgeIndex, std::array<int, 3>> edgeTriangles;

  for (int tri = 0; tri < mesh.numTriangles(); tri++) {
    for (int ei = 0; ei < 3; ei++) {
      int v0 = mesh.tri(tri)[ei];
      int v1 = mesh.tri(tri)[(ei + 1) % 3];

      if (v0 > v1)
        std::swap(v0, v1);

      EdgeIndex eidx(v0, v1);
      auto iter = edgeTriangles.find(eidx);

      // if it exist
      if (iter != edgeTriangles.end()) {
        if (iter->second[2] >= 2) {
          SPDLOG_LOGGER_ERROR(pgo::Logging::lgr(), "Input mesh is not manifold.");
          throw std::invalid_argument("Input mesh is not manifold.");
        }

        iter->second[1] = tri;
        iter->second[2]++;
      }
      else {
        std::array<int, 3> idx = { tri, -1, 1 };
        edgeTriangles.emplace(eidx, idx);
      }
    }
  }

  for (const auto &pr : edgeTriangles) {
    int v0 = pr.first.first;
    int v1 = pr.first.second;

    if (pr.second[2] < 2)
      continue;

    int vidx = 0;
    for (; vidx < 3; vidx++) {
      if (mesh.tri(pr.second[0])[vidx] == v0)
        break;
    }

    if (vidx >= 3) {
      SPDLOG_LOGGER_ERROR(pgo::Logging::lgr(), "Input mesh is ill.");
      throw std::invalid_argument("Input mesh is ill.");
    }

    int vidxNext = (vidx + 1) % 3;
    int v2 = -1, v3 = -1;

    // v2 belongs to tri0
    if (mesh.tri(pr.second[0])[vidxNext] == v1) {
      v2 = mesh.tri(pr.second[0])[(vidxNext + 1) % 3];
      // find the index for v3
      uint64_t vAll = (uint64_t)mesh.tri(pr.second[1])[0] ^ (uint64_t)mesh.tri(pr.second[1])[1] ^ (uint64_t)mesh.tri(pr.second[1])[2];
      v3 = (int)(vAll ^ (uint64_t)v0 ^ (uint64_t)v1);
    }
    else {
      vidxNext = 0;
      for (; vidxNext < 3; vidxNext++) {
        if (mesh.tri(pr.second[1])[vidxNext] == v1)
          break;
      }

      if (!(vidxNext < 3)) {
        SPDLOG_LOGGER_ERROR(pgo::Logging::lgr(), "Input mesh is ill.");
        throw std::invalid_argument("Input mesh is ill.");
      }

      v2 = mesh.tri(pr.second[1])[(vidxNext + 1) % 3];
      // find the index for v3
      uint64_t vAll = (uint64_t)mesh.tri(pr.second[0])[0] ^ (uint64_t)mesh.tri(pr.second[0])[1] ^ (uint64_t)mesh.tri(pr.second[0])[2];
      v3 = (int)(vAll ^ (uint64_t)v0 ^ (uint64_t)v1);
    }
    if (!(v2 >= 0 && v3 >= 0)) {
      SPDLOG_LOGGER_ERROR(pgo::Logging::lgr(), "Input mesh is ill.");
      throw std::invalid_argument("Input mesh is ill.");
    }

    std::array<int, 4> quad = { v0, v1, v2, v3 };
    surfaceQuads.emplace_back(quad);
  }

  elementWeights.resize(surfaceQuads.size());
  elementVertexWeights.resize(surfaceQuads.size());
  restValues.resize(surfaceQuads.size());
  updateRestInfo();

  std::vector<ES::TripletD> entries;
  for (int ei = 0; ei < (int)surfaceQuads.size(); ei++) {
    for (int vi = 0; vi < 4; vi++) {
      int vidx_i = surfaceQuads[ei][vi];
      for (int vj = 0; vj < 4; vj++) {
        int vidx_j = surfaceQuads[ei][vj];

        for (int r = 0; r < 3; r++) {
          for (int c = 0; c < 3; c++) {
            entries.emplace_back(vidx_i * 3 + r, vidx_j * 3 + c, 1.0);
          }
        }
      }
    }
  }

  hessTemplate.resize(restPositions.size(), restPositions.size());
  hessTemplate.setFromTriplets(entries.begin(), entries.end());

  allDOFs.resize(restPositions.size());
  std::iota(allDOFs.begin(), allDOFs.end(), 0);

  elementOffsets.resize(surfaceQuads.size());
  for (int ei = 0; ei < (int)surfaceQuads.size(); ei++) {
    for (int vi = 0; vi < 4; vi++) {
      int vidx_i = surfaceQuads[ei][vi];
      for (int vj = 0; vj < 4; vj++) {
        int vidx_j = surfaceQuads[ei][vj];

        for (int r = 0; r < 3; r++) {
          for (int c = 0; c < 3; c++) {
            elementOffsets[ei](vi * 3 + r, vj * 3 + c) = ES::findEntryOffset(hessTemplate, vidx_i * 3 + r, vidx_j * 3 + c);
            PGO_ALOG(elementOffsets[ei](vi * 3 + r, vj * 3 + c) >= 0);
          }
        }
      }
    }
  }

  buf = std::make_shared<SurfaceSmoothnessAbsoluteMeanCurvatureBuf>(mesh.numVertices());
}

SurfaceSmoothnessAbsoluteMeanCurvature::~SurfaceSmoothnessAbsoluteMeanCurvature()
{
}

void SurfaceSmoothnessAbsoluteMeanCurvature::setDOFs(const std::vector<int> &dofs)
{
  PGO_ALOG((int)dofs.size() == (int)allDOFs.size());
  allDOFs = dofs;
}

void SurfaceSmoothnessAbsoluteMeanCurvature::updateRestInfo()
{
  tbb::parallel_for(0, (int)surfaceQuads.size(), [&](int ei) {
    Eigen::Matrix<double, 3, 4> p;
    for (int i = 0; i < 4; ++i) {
      p.col(i) = restPositions.segment<3>(surfaceQuads[ei][i] * 3);
    }

    double l01 = (p.col(0) - p.col(1)).norm();
    double l02 = (p.col(0) - p.col(2)).norm();
    double l12 = (p.col(1) - p.col(2)).norm();
    double r0 = 0.5 * (l01 + l02 + l12);
    double A0 = std::sqrt(r0 * (r0 - l01) * (r0 - l02) * (r0 - l12));
    double l03 = (p.col(0) - p.col(3)).norm();
    double l13 = (p.col(1) - p.col(3)).norm();
    double r1 = 0.5 * (l01 + l03 + l13);
    double A1 = std::sqrt(r1 * (r1 - l01) * (r1 - l03) * (r1 - l13));

    elementWeights[ei] = std::sqrt(3.0 / (A0 + A1));

    double cot02 = ((l01 * l01) - (l02 * l02) + (l12 * l12)) / (4.0 * A0);
    double cot12 = ((l01 * l01) + (l02 * l02) - (l12 * l12)) / (4.0 * A0);
    double cot03 = ((l01 * l01) - (l03 * l03) + (l13 * l13)) / (4.0 * A1);
    double cot13 = ((l01 * l01) + (l03 * l03) - (l13 * l13)) / (4.0 * A1);

    elementVertexWeights[ei].block<3, 3>(0, 0) = ES::M3d::Identity() * (cot02 + cot03);
    elementVertexWeights[ei].block<3, 3>(0, 3) = ES::M3d::Identity() * (cot12 + cot13);
    elementVertexWeights[ei].block<3, 3>(0, 6) = ES::M3d::Identity() * -(cot02 + cot12);
    elementVertexWeights[ei].block<3, 3>(0, 9) = ES::M3d::Identity() * -(cot03 + cot13);

    restValues[ei] = (elementVertexWeights[ei] * ES::Mp<ES::V12d>(p.data())).norm();
  });
}

// ||norm(x) - norm_rest ||^2
double SurfaceSmoothnessAbsoluteMeanCurvature::func(EigenSupport::ConstRefVecXd x) const
{
  double energyAll = tbb::parallel_reduce(
    tbb::blocked_range<int>(0, (int)surfaceQuads.size()), 0.0,
    [&](const tbb::blocked_range<int> &r, double init) -> double {
      for (int ei = r.begin(); ei != r.end(); ++ei) {
        ES::V12d p;
        for (int i = 0; i < 4; ++i) {
          p.segment<3>(i * 3) = x.segment<3>(surfaceQuads[ei][i] * 3);
        }

        ES::V3d Hi = elementVertexWeights[ei] * p;
        double v_norm = Hi.norm();
        init += (v_norm - restValues[ei]) * (v_norm - restValues[ei]) * elementWeights[ei];
      }

      return init;
    },
    [](double x, double y) -> double {
      return x + y;
    });

  return energyAll * 0.5;
}

// E = 1/2 \sum mi || |Hi(x)| - |Hi_0| ||^2
// dE/dHi = mi (|Hi(x)| - |Hi_0|) Hi/|Hi|
// dH/dx = W
void SurfaceSmoothnessAbsoluteMeanCurvature::gradient(EigenSupport::ConstRefVecXd x, EigenSupport::RefVecXd grad) const
{
  grad.setZero();

  // for (int ei = 0; ei < (int)surfaceQuads.size(); ei++) {
  tbb::parallel_for(0, (int)surfaceQuads.size(), [&](int ei) {
    ES::V12d p;
    for (int i = 0; i < 4; ++i) {
      p.segment<3>(i * 3) = x.segment<3>(surfaceQuads[ei][i] * 3);
    }

    ES::V3d Hi = elementVertexWeights[ei] * p;
    double curHNorm = std::max(Hi.norm(), eps);
    ES::V3d dEdHi = elementWeights[ei] * (curHNorm - restValues[ei]) / curHNorm * Hi;
    ES::V12d dEdx = elementVertexWeights[ei].transpose() * dEdHi;
    for (int i = 0; i < 4; i++) {
      buf->locks[surfaceQuads[ei][i]].lock();

      grad.segment<3>(surfaceQuads[ei][i] * 3) += dEdx.segment<3>(i * 3);

      buf->locks[surfaceQuads[ei][i]].unlock();
    }
  });
}

// E = 1/2 \sum mi || |Hi(x)| - |Hi_0| ||^2
// dE/dHi = mi (|Hi(x)| - |Hi_0|) Hi/|Hi|
// dH/dx = M^-1L

// d2E/dx2 = d(dE/dH dH/dx)/dx = dH/dx d2E/dH2 dH/dx
void SurfaceSmoothnessAbsoluteMeanCurvature::hessian(EigenSupport::ConstRefVecXd x, EigenSupport::SpMatD &hess) const
{
  memset(hess.valuePtr(), 0, sizeof(double) * hess.nonZeros());

  // for (int ei = 0; ei < (int)surfaceQuads.size(); ei++) {
  tbb::parallel_for(0, (int)surfaceQuads.size(), [&](int ei) {
    ES::V12d p;
    for (int i = 0; i < 4; ++i) {
      p.segment<3>(i * 3) = x.segment<3>(surfaceQuads[ei][i] * 3);
    }

    ES::V3d Hi = elementVertexWeights[ei] * p;
    double curHNorm = std::max(Hi.norm(), eps);

    ES::V3d dzdH;
    ES::M3d d2zdH2;

    // dz/dH = 0.5 * (H^T H)^-0.5 * 2 H  = H / (HT H)^0.5
    // d2z/dH2 = 0.5 * (-0.5) (H^T H)^(-1.5) 2H * 2H + 0.5 * (H^T H)^(-0.5) * 2.0
    //         = - HHT / (HTH)^1.5 + 1 / (HTH)^0.5
    dzdH = Hi / curHNorm;
    d2zdH2 = ES::tensorProduct(Hi, Hi) / (-curHNorm * curHNorm * curHNorm) + ES::M3d::Identity() / curHNorm;

    ES::M3d d2EdH2;
    // dE/dz = mi (z - c0) * dz/dH
    // d2E/dz2 = mi (dz/dH * dz/dH + (z - c0) * d2z/dH2)
    d2EdH2 = (ES::tensorProduct(dzdH, dzdH) + d2zdH2 * (curHNorm - restValues[ei])) * elementWeights[ei];

    ES::M12d hessLocal;
    hessLocal = elementVertexWeights[ei].transpose() * d2EdH2 * elementVertexWeights[ei];

    for (int vi = 0; vi < 4; vi++) {
      int vidx_i = surfaceQuads[ei][vi];

      for (int vj = 0; vj < 4; vj++) {
        int vidx_j = surfaceQuads[ei][vj];

        buf->locks[vidx_i].lock();

        for (int r = 0; r < 3; r++) {
          for (int c = 0; c < 3; c++) {
            std::ptrdiff_t offset = elementOffsets[ei](vi * 3 + r, vj * 3 + c);
            hess.valuePtr()[offset] += hessLocal(vi * 3 + r, vj * 3 + c);
          }
        }

        buf->locks[vidx_i].unlock();
      }
    }
  });
}