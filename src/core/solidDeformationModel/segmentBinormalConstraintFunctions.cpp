/*
author: Bohan Wang
copyright to USC,MIT
*/

#include "segmentBinormalConstraintFunctions.h"

#include "pgoLogging.h"
#include "EigenSupport.h"

#include <tbb/parallel_for.h>
#include <tbb/partitioner.h>

using namespace pgo;
using namespace pgo::SolidDeformationModel;
namespace ES = pgo::EigenSupport;

SegmentBinormalConstraintFunctions::SegmentBinormalConstraintFunctions(int nAll, int pStart, int sStart, int numSegments, const double *positions):
  NonlinearOptimization::ConstraintFunctions(nAll), nSeg(numSegments), positionDOFStart(pStart), segmentDOFStart(sStart)
{
  if (positions) {
    restPositions = ES::Mp<const ES::VXd>(positions, numSegments * 3 + 3);
  }

  std::vector<ES::TripletD> entries;
  // each segment has two constraints
  // perpendicular constraint
  // unit vector constraint

  for (int ei = 0; ei < nSeg; ei++) {
    // unit vector constraint
    entries.emplace_back(ei * 2, segmentDOFStart + ei * 3, 1);
    entries.emplace_back(ei * 2, segmentDOFStart + ei * 3 + 1, 1);
    entries.emplace_back(ei * 2, segmentDOFStart + ei * 3 + 2, 1);

    // perpendicular constraint
    if (restPositions.size()) {
      entries.emplace_back(ei * 2 + 1, segmentDOFStart + ei * 3, 1);
      entries.emplace_back(ei * 2 + 1, segmentDOFStart + ei * 3 + 1, 1);
      entries.emplace_back(ei * 2 + 1, segmentDOFStart + ei * 3 + 2, 1);
    }
    else {
      entries.emplace_back(ei * 2 + 1, positionDOFStart + ei * 3, 1);
      entries.emplace_back(ei * 2 + 1, positionDOFStart + ei * 3 + 1, 1);
      entries.emplace_back(ei * 2 + 1, positionDOFStart + ei * 3 + 2, 1);

      entries.emplace_back(ei * 2 + 1, positionDOFStart + ei * 3 + 3, 1);
      entries.emplace_back(ei * 2 + 1, positionDOFStart + ei * 3 + 4, 1);
      entries.emplace_back(ei * 2 + 1, positionDOFStart + ei * 3 + 5, 1);

      entries.emplace_back(ei * 2 + 1, segmentDOFStart + ei * 3, 1);
      entries.emplace_back(ei * 2 + 1, segmentDOFStart + ei * 3 + 1, 1);
      entries.emplace_back(ei * 2 + 1, segmentDOFStart + ei * 3 + 2, 1);
    }
  }
  jacobianTemplate.resize(nSeg * 2, nAll);
  jacobianTemplate.setFromTriplets(entries.begin(), entries.end());

  jacobianIndices.assign(nSeg, JacIndex::Constant(-1));
  for (int ei = 0; ei < nSeg; ei++) {
    JacIndex jacIdx;

    jacIdx[0] = ES::findEntryOffset(jacobianTemplate, ei * 2, segmentDOFStart + ei * 3);
    PGO_ALOG(jacIdx[0] >= 0);

    if (restPositions.size() == 0) {
      jacIdx[1] = ES::findEntryOffset(jacobianTemplate, ei * 2 + 1, positionDOFStart + ei * 3);
      PGO_ALOG(jacIdx[1] >= 0);
    }

    jacIdx[2] = ES::findEntryOffset(jacobianTemplate, ei * 2 + 1, segmentDOFStart + ei * 3);
    PGO_ALOG(jacIdx[2] >= 0);

    jacobianIndices[ei] = jacIdx;
  }

  entries.clear();
  for (int ei = 0; ei < nSeg; ei++) {
    if (restPositions.size() == 0) {
      int dofStart[3] = {
        positionDOFStart + ei * 3,
        positionDOFStart + ei * 3 + 3,
        segmentDOFStart + ei * 3
      };

      for (int vi = 0; vi < 3; vi++) {
        for (int vj = 0; vj < 3; vj++) {
          for (int dofi = 0; dofi < 3; dofi++) {
            for (int dofj = 0; dofj < 3; dofj++) {
              entries.emplace_back(dofStart[vi] + dofi, dofStart[vj] + dofj, 1);
            }
          }
        }
      }
    }
    else {
      int dofStart[1] = {
        segmentDOFStart + ei * 3
      };

      for (int dofi = 0; dofi < 3; dofi++) {
        for (int dofj = 0; dofj < 3; dofj++) {
          entries.emplace_back(dofStart[0] + dofi, dofStart[0] + dofj, 1);
        }
      }
    }
  }

  lambdahTemplate.resize(nAll, nAll);
  lambdahTemplate.setFromTriplets(entries.begin(), entries.end());

  hessIndices.assign(nSeg, HessIndex::Constant(-1));
  for (int ei = 0; ei < nSeg; ei++) {
    HessIndex idx;

    if (restPositions.size() == 0) {
      int dofStart[3] = {
        positionDOFStart + ei * 3,
        positionDOFStart + ei * 3 + 3,
        segmentDOFStart + ei * 3
      };

      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          for (int k = 0; k < 3; k++) {
            for (int l = 0; l < 3; l++) {
              int row = dofStart[i] + k;
              int col = dofStart[j] + l;

              int localRow = i * 3 + k;
              int localCol = j * 3 + l;

              idx(localRow, localCol) = ES::findEntryOffset(lambdahTemplate, row, col);
              PGO_ALOG(idx(localRow, localCol) >= 0);
            }
          }
        }
      }
    }
    else {
      int dofStart[1] = {
        segmentDOFStart + ei * 3
      };

      for (int k = 0; k < 3; k++) {
        for (int l = 0; l < 3; l++) {
          int row = dofStart[0] + k;
          int col = dofStart[0] + l;

          int localRow = k;
          int localCol = l;

          idx(localRow, localCol) = ES::findEntryOffset(lambdahTemplate, row, col);
          PGO_ALOG(idx(localRow, localCol) >= 0);
        }
      }
    }

    hessIndices[ei] = idx;
  }
  hessLocks = std::vector<tbb::spin_mutex>(lambdahTemplate.nonZeros());
}

SegmentBinormalConstraintFunctions::~SegmentBinormalConstraintFunctions()
{
}

// C_i0 = dir^2 - 1 = 0
// C_i1 = (x1-x0).dir = 0
void SegmentBinormalConstraintFunctions::func(ES::ConstRefVecXd x, ES::RefVecXd g) const
{
  tbb::parallel_for(
    0, nSeg, [&](int si) {
      if (restPositions.size() == 0) {
        ES::V3d v1 = x.segment<3>(positionDOFStart + (si + 1) * 3);
        ES::V3d v0 = x.segment<3>(positionDOFStart + si * 3);
        ES::V3d d0 = x.segment<3>(segmentDOFStart + si * 3);

        g[si * 2] = d0.dot(d0) - 1;
        g[si * 2 + 1] = d0.dot(v1 - v0);
      }
      else {
        ES::V3d v1 = restPositions.segment<3>(positionDOFStart + (si + 1) * 3);
        ES::V3d v0 = restPositions.segment<3>(positionDOFStart + si * 3);
        ES::V3d d0 = x.segment<3>(segmentDOFStart + si * 3);

        g[si * 2] = d0.dot(d0) - 1;
        g[si * 2 + 1] = d0.dot(v1 - v0);
      }
    },
    tbb::static_partitioner());
}

// d C_i0 / d dir= 2 * dir
// d C_i1 / d x0= -dir
// d C_i1 / d x1= dir
// d C_i1 / d dir= x1 - x0
void SegmentBinormalConstraintFunctions::jacobian(ES::ConstRefVecXd x, ES::SpMatD &jac) const
{
  tbb::parallel_for(
    0, nSeg, [&](int si) {
      if (restPositions.size() == 0) {
        ES::V3d v1 = x.segment<3>(positionDOFStart + (si + 1) * 3);
        ES::V3d v0 = x.segment<3>(positionDOFStart + si * 3);
        ES::V3d d0 = x.segment<3>(segmentDOFStart + si * 3);

        // d C_i0
        ES::V3d ddir = 2.0 * d0;
        for (int dof = 0; dof < 3; dof++) {
          jac.valuePtr()[jacobianIndices[si][0] + dof] = ddir[dof];
        }

        // d C_i1
        ES::V3d dx0 = -d0;
        ES::V3d dx1 = d0;
        ddir = v1 - v0;
        for (int dof = 0; dof < 3; dof++) {
          jac.valuePtr()[jacobianIndices[si][1] + dof] = dx0[dof];
          jac.valuePtr()[jacobianIndices[si][1] + 3 + dof] = dx1[dof];
          jac.valuePtr()[jacobianIndices[si][2] + dof] = ddir[dof];
        }
      }
      else {
        ES::V3d v1 = restPositions.segment<3>(positionDOFStart + (si + 1) * 3);
        ES::V3d v0 = restPositions.segment<3>(positionDOFStart + si * 3);
        ES::V3d d0 = x.segment<3>(segmentDOFStart + si * 3);

        // d C_i0
        ES::V3d ddir = 2.0 * d0;
        for (int dof = 0; dof < 3; dof++) {
          jac.valuePtr()[jacobianIndices[si][0] + dof] = ddir[dof];
        }

        // d C_i1
        ddir = v1 - v0;
        for (int dof = 0; dof < 3; dof++) {
          jac.valuePtr()[jacobianIndices[si][2] + dof] = ddir[dof];
        }
      }
    },
    tbb::static_partitioner());
}

void SegmentBinormalConstraintFunctions::hessian(ES::ConstRefVecXd x, ES::ConstRefVecXd lambda, ES::SpMatD &hess) const
{
  memset(hess.valuePtr(), 0, sizeof(double) * hess.nonZeros());

  tbb::parallel_for(
    0, nSeg, [&](int si) {
      if (restPositions.size() == 0) {
        ES::M9d h1, h2;
        h1.setZero();
        h1.block<3, 3>(6, 6) = ES::M3d::Identity() * 2.0;  // d^2 C_i0 / d dir^2

        h1 *= lambda[si * 2];  // this is the lambda for the Ci0

        h2.setZero();
        h2.block<3, 3>(0, 6) = -ES::M3d::Identity();  // d^2 C_i1 / d x0 d dir= -I
        h2.block<3, 3>(6, 0) = -ES::M3d::Identity();  // d^2 C_i1 / d x0 d dir= -I

        h2.block<3, 3>(3, 6) = ES::M3d::Identity();  // d^2 C_i1 / dx1 d dir = I
        h2.block<3, 3>(6, 3) = ES::M3d::Identity();  // d^2 C_i1 / dx1 d dir = I

        // all other entries are zero
        h2 *= lambda[si * 2 + 1];

        for (int col = 0; col < 3; col++) {
          for (int dofj = 0; dofj < 3; dofj++) {
            for (int row = 0; row < 3; row++) {
              for (int dofi = 0; dofi < 3; dofi++) {
                int localRow = row * 3 + dofi;
                int localCol = col * 3 + dofj;

                ES::IDX offset = hessIndices[si](localRow, localCol);

                hessLocks[offset].lock();

                hess.valuePtr()[offset] += h1(localRow, localCol) + h2(localRow, localCol);

                hessLocks[offset].unlock();
              }
            }
          }
        }
      }
      else {
        ES::M3d h;
        h.setZero();
        h = ES::M3d::Identity() * 2.0;  // d^2 C_i0 / d dir^2
        h *= lambda[si * 2];            // this is the lambda for the Ci0
        // for the Ci1, it is a linear constraints so there is no hessian

        for (int dofj = 0; dofj < 3; dofj++) {
          for (int dofi = 0; dofi < 3; dofi++) {
            int localRow = dofi;
            int localCol = dofj;

            ES::IDX offset = hessIndices[si](localRow, localCol);

            hessLocks[offset].lock();

            hess.valuePtr()[offset] += h(localRow, localCol);

            hessLocks[offset].unlock();
          }
        }
      }
    },
    tbb::static_partitioner());
}
