#include "triangleAffineToPositionEnergy.h"

#include "triMeshNeighbor.h"
#include "geometryQuery.h"
#include "EigenSupport.h"

#include <numeric>

using namespace pgo::PredefinedPotentialEnergies;
namespace ES = pgo::EigenSupport;

namespace pgo::PredefinedPotentialEnergies
{
struct TriangleAffineToPositionEnergyBuf
{
  ES::SpMatD hess_orig, hess_orig11, hess_orig12, hess_orig22;
  ES::SpMatD H11C, CTH11C, CTH12, H21C;

  ES::SpMatI hess_orig11Map, hess_orig12Map, hess_orig22Map, H21CTransposeMap;
  ES::SpMatI CTH11CMap, CTH12Map, H21CMap, H22Map;

  ES::SymbolicMmData *mul_hc = nullptr, *mul_cthc = nullptr, *mul_cth = nullptr;

  ES::VXd grad_orig;
  ES::VXd A_orig;
};
}  // namespace pgo::PredefinedPotentialEnergies

TriangleAffineToPositionEnergy::TriangleAffineToPositionEnergy(int totalNumDOFs,
  const Mesh::TriMeshGeo &mesh, std::shared_ptr<const PotentialEnergy> origE):
  nAll(totalNumDOFs),
  inputMesh(mesh), originalEnergy(origE)
{
  buf = std::make_shared<TriangleAffineToPositionEnergyBuf>();
  originalEnergy->createHessian(buf->hess_orig);
  buf->grad_orig.setZero(originalEnergy->getNumDOFs());
  buf->A_orig.setZero(originalEnergy->getNumDOFs());

  computeWA();

  if (totalNumDOFs == (int)WA.cols()) {
    nRestDOFs = 0;

    // Horig \in R^{3n x 3n}
    // AToX \in R^{3n x 12n}
    ES::symbolicMm(buf->hess_orig, WA, buf->H11C, &buf->mul_hc);
    ES::symbolicMm(WA, buf->H11C, hessTemplate, &buf->mul_cthc, 1);

    allDOFs.resize(WA.cols());
    std::iota(allDOFs.begin(), allDOFs.end(), 0);
  }
  else {
    // have not debugged this part
    abort();
  }
}

TriangleAffineToPositionEnergy::~TriangleAffineToPositionEnergy()
{
  ES::destroySymbolicMmData(buf->mul_hc);
  ES::destroySymbolicMmData(buf->mul_cthc);
  ES::destroySymbolicMmData(buf->mul_cth);
}

void TriangleAffineToPositionEnergy::computeWA()
{
  WA.setZero();

  // x = \sum wi Ai x
  Mesh::TriMeshNeighbor meshNeighbor(inputMesh);
  vertexNeighboringTriangleWeights.resize(inputMesh.numVertices());
  for (int vi = 0; vi < inputMesh.numVertices(); vi++) {
    const auto &nTris = meshNeighbor.getVtxNearbyTriangles(vi);
    vertexNeighboringTriangleWeights[vi].resize(nTris.size());

    for (int j = 0; j < (int)nTris.size(); j++) {
      std::get<0>(vertexNeighboringTriangleWeights[vi][j]) = nTris[j];
      for (int k = 0; k < 3; k++) {
        if (inputMesh.triVtxID(nTris[j], k) == vi) {
          std::get<1>(vertexNeighboringTriangleWeights[vi][j]) = Mesh::getTriangleAngle(inputMesh.pos(nTris[j], k), inputMesh.pos(nTris[j], (k + 1) % 3), inputMesh.pos(nTris[j], (k + 2) % 3));
          break;
        }
      }
    }

    double sum = 0.0;
    for (int j = 0; j < (int)nTris.size(); j++) {
      sum += std::get<1>(vertexNeighboringTriangleWeights[vi][j]);
    }

    for (int j = 0; j < (int)nTris.size(); j++) {
      std::get<1>(vertexNeighboringTriangleWeights[vi][j]) /= sum;
    }
  }

  entries.reserve(inputMesh.numVertices() * 12 * 10);
  for (int vi = 0; vi < inputMesh.numVertices(); vi++) {
    for (int ni = 0; ni < (int)vertexNeighboringTriangleWeights[vi].size(); ni++) {
      for (int k = 0; k < 12; k++) {
        entries.emplace_back(vi * 12 + k, std::get<0>(vertexNeighboringTriangleWeights[vi][ni]) * 12 + k, std::get<1>(vertexNeighboringTriangleWeights[vi][ni]));
      }
    }
  }

  WA.resize(inputMesh.numVertices() * 12, inputMesh.numTriangles() * 12);
  WA.setFromTriplets(entries.begin(), entries.end());
}

void TriangleAffineToPositionEnergy::compute_A_orig(pgo::EigenSupport::ConstRefVecXd A_tri, pgo::EigenSupport::RefVecXd A_orig) const
{
  // x = C A
  ES::mv(WA, A_tri.head(WA.cols()), A_orig.head(WA.rows()), 0);
  if (nRestDOFs) {
    A_orig.tail(nRestDOFs) = A_tri.tail(nRestDOFs);
  }
}

double TriangleAffineToPositionEnergy::func(pgo::EigenSupport::ConstRefVecXd A) const
{
  compute_A_orig(A, buf->A_orig);

  return originalEnergy->func(buf->A_orig);
}

void TriangleAffineToPositionEnergy::gradient(pgo::EigenSupport::ConstRefVecXd A, pgo::EigenSupport::RefVecXd grad) const
{
  compute_A_orig(A, buf->A_orig);

  buf->grad_orig.setZero();
  originalEnergy->gradient(buf->A_orig, buf->grad_orig);

  ES::mv(WA, buf->grad_orig.head(WA.rows()), grad.head(WA.cols()), 1);

  if (nRestDOFs) {
    grad.tail(nRestDOFs) = buf->grad_orig.tail(nRestDOFs);
  }
}

void TriangleAffineToPositionEnergy::hessian(pgo::EigenSupport::ConstRefVecXd A, pgo::EigenSupport::SpMatD &hess) const
{
  compute_A_orig(A, buf->A_orig);

  std::memset(buf->hess_orig.valuePtr(), 0, sizeof(double) * buf->hess_orig.nonZeros());
  originalEnergy->hessian(buf->A_orig, buf->hess_orig);

  if (nRestDOFs == 0) {
    ES::mm(buf->hess_orig, WA, buf->mul_hc, buf->H11C, 0);
    ES::mm(WA, buf->H11C, buf->mul_cthc, hess, 1);
  }
  else {
    throw std::runtime_error("not done");
  }
}
