/*
author: Bohan Wang
copyright to USC,MIT
*/

#include "generateTetMeshMatrix.h"

#include "tetMeshGeo.h"
#include "pgoLogging.h"
#include "geometryQuery.h"

#include <tbb/parallel_for.h>
#include <tbb/concurrent_vector.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/cache_aligned_allocator.h>

using namespace pgo;
using namespace pgo::SolidDeformationModel;

namespace pgo
{
namespace SolidDeformationModel
{
void generateElementMatrixEntries(const Mesh::TetMeshRef &tetMesh, int tetID, double m[12])
{
  // grad is constant inside a tet
  Vec3d vtx[4];
  for (int i = 0; i < 4; i++)
    vtx[i] = tetMesh.pos(tetID, i);

  // form M =
  // [b - a]
  // [c - a]
  // [d - a]

  Mat3d M = asMat3d(vtx[1] - vtx[0], vtx[2] - vtx[0], vtx[3] - vtx[0]);
  Mat3d MInvT = M.inverse().transpose();

  // the 12x1 column vector m is seen as [m0, m1, m2, m3], where mi is a 3x1 column vector
  Vec3d m0;
  m0.setZero();

  for (int i = 0; i < 3; i++) {
    Vec3d r = MInvT.row(i);
    (Eigen::Map<Vec3d>(m + 3 * (i + 1))) = r;

    // MInvT[i].convertToArray(m + 3 * (i + 1));  // assign m{i+1} to m
    m0 -= r;
  }
  (Eigen::Map<Vec3d>(m)) = m0;
}
}  // namespace SolidDeformationModel
}  // namespace pgo

void TetMeshMatrix::generateElementGradientMatrix(const Mesh::TetMeshRef &tetMesh, int tetID, ES::M9x12d &G)
{
  double m[12];
  generateElementMatrixEntries(tetMesh, tetID, m);

  // G is 9 x 12:
  //        [ m0       m1       m2       m3       ]
  //  G =   [    m0       m1       m2       m3    ]
  //        [       m0       m1       m2       m3 ]
  //
  // where mi are 3-vectors

  G.setZero();

  for (int vtx = 0; vtx < 4; vtx++)
    for (int dof = 0; dof < 3; dof++)
      for (int j = 0; j < 3; j++)
        G(3 * j + dof, 3 * vtx + dof) = m[3 * vtx + j];
}

void TetMeshMatrix::generateGradientMatrix(const Mesh::TetMeshRef &tetMesh, ES::SpMatD &G)
{
  // G is 9 #tets x 3 #vertices
  // std::vector<ES::TripletD> entries;
  // tbb::enumerable_thread_specific<std::vector<ES::TripletD, tbb::cache_aligned_allocator<ES::TripletD>>> entriesTLS;
  tbb::concurrent_vector<ES::TripletD> entries;

  // for (int tetID = 0; tetID < tetMesh.getNumElements(); tetID++) {
  tbb::parallel_for(
    0, tetMesh.numTets(), [&](int tetID) {
      // auto &entriesBuf = entriesTLS.local();

      double m[12];
      generateElementMatrixEntries(tetMesh, tetID, m);

      //// write dFduPacked in place
      // write m in place
      for (int vtx = 0; vtx < 4; vtx++)
        for (int dof = 0; dof < 3; dof++) {
          int column = 3 * tetMesh.tetVtxID(tetID, vtx) + dof;
          for (int j = 0; j < 3; j++) {
            int row = 9 * tetID + 3 * j + dof;
            // int row = 9 * tetID + 3 * dof + j;
            // double entry = dFduPacked[ELT(9, 3 * dof + j, vtx)];
            double entry = m[3 * vtx + j];
            entries.emplace_back(row, column, entry);
          }
        }
    },
    tbb::static_partitioner());

  // for (auto it = entriesTLS.begin(); it != entriesTLS.end(); ++it) {
  //   entries.insert(entries.end(), it->begin(), it->end());
  // }

  G.resize(tetMesh.numTets() * 9, tetMesh.numVertices() * 3);
  G.setFromTriplets(entries.begin(), entries.end());
}

void TetMeshMatrix::generateBasicElementLaplacianMatrix(const Mesh::TetMeshRef &tetMesh, ES::SpMatD &L, int faceNeighbor, int scale)
{
  int numElements = tetMesh.numTets();
  int numElementVertices = 4;
  int numVertices = tetMesh.numVertices();

  // build elements that neighbor each vertex
  std::vector<std::vector<int>> vertexNeighbors(numVertices);
  for (int el = 0; el < numElements; el++) {
    for (int vtxIdx = 0; vtxIdx < numElementVertices; vtxIdx++) {
      int vertexIndex = tetMesh.tetVtxID(el, vtxIdx);
      vertexNeighbors[vertexIndex].push_back(el);
    }
  }

  tbb::parallel_for((size_t)0, vertexNeighbors.size(), [&](size_t vi) {
    std::sort(vertexNeighbors[vi].begin(), vertexNeighbors[vi].end());
    auto itt = std::unique(vertexNeighbors[vi].begin(), vertexNeighbors[vi].end());
    vertexNeighbors[vi].erase(itt, vertexNeighbors[vi].end());
  },
    tbb::static_partitioner());

  // build elements that neighbor each element, and assemble L
  tbb::enumerable_thread_specific<std::vector<int, tbb::cache_aligned_allocator<int>>> elementNeighborsTLS;
  tbb::concurrent_vector<ES::TripletD> entries;

  // for (int el = 0; el < numElements; el++)
  tbb::parallel_for(
    0, numElements, [&](int el) {
      auto &elementNeighbors = elementNeighborsTLS.local();

      elementNeighbors.clear();
      for (int vtxIdx = 0; vtxIdx < numElementVertices; vtxIdx++) {
        int vertexIndex = tetMesh.tetVtxID(el, vtxIdx);
        elementNeighbors.insert(elementNeighbors.end(), vertexNeighbors[vertexIndex].begin(), vertexNeighbors[vertexIndex].end());
      }

      std::sort(elementNeighbors.begin(), elementNeighbors.end());
      auto itt = std::unique(elementNeighbors.begin(), elementNeighbors.end());
      elementNeighbors.erase(itt, elementNeighbors.end());

      itt = std::find(elementNeighbors.begin(), elementNeighbors.end(), el);
      if (itt != elementNeighbors.end()) {
        elementNeighbors.erase(itt);
      }

      if (faceNeighbor) {
        int eleVtxIDs[4];
        for (int vtxIdx = 0; vtxIdx < numElementVertices; vtxIdx++) {
          eleVtxIDs[vtxIdx] = tetMesh.tetVtxID(el, vtxIdx);
        }

        int inc = 0;
        std::array<int, 8> eleNeighbor;
        for (int i = 0; i < (int)elementNeighbors.size(); i++) {
          if (elementNeighbors[i] == el)
            continue;

          int eleVtxIDs2[4];
          for (int vtxIdx = 0; vtxIdx < numElementVertices; vtxIdx++) {
            eleVtxIDs2[vtxIdx] = tetMesh.tetVtxID(elementNeighbors[i], vtxIdx);
          }

          int count = 0;
          for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
              if (eleVtxIDs[j] == eleVtxIDs2[k]) {
                count++;
                break;
              }
            }
          }

          if (count == 3) {
            eleNeighbor[inc++] = elementNeighbors[i];
            if (inc >= 8)
              break;
          }
        }
        PGO_ALOG(inc <= 4);

        elementNeighbors.clear();
        elementNeighbors.assign(eleNeighbor.data(), eleNeighbor.data() + inc);
      }

      int c = 0;
      for (int otherEl : elementNeighbors) {
        if (otherEl != el) {
          entries.emplace_back(el, otherEl, -1.0);
          c++;
        }
      }

      entries.emplace_back(el, el, (double)c);
    });

  L.resize(numElements, numElements);
  L.setFromTriplets(entries.begin(), entries.end());

  if (scale) {
    for (ES::IDX rowi = 0; rowi < L.rows(); rowi++) {
      double s = L.coeff(rowi, rowi);
      for (ES::SpMatD::InnerIterator it(L, rowi); it; ++it) {
        it.valueRef() /= s;
      }
    }
  }
}

void TetMeshMatrix::generateVertexLapacianMatrix(const Mesh::TetMeshRef &tetMesh, ES::SpMatD &L)
{
  // can be done using libigl
  SPDLOG_LOGGER_CRITICAL(pgo::Logging::lgr(), "Not implemented yet.");
}