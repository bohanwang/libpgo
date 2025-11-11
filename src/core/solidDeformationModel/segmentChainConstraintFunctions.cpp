/*
author: Bohan Wang
copyright to USC,MIT,NUS
*/

#include "segmentChainConstraintFunctions.h"

#include "pgoLogging.h"
#include "EigenSupport.h"

#include <tbb/parallel_for.h>
#include <tbb/partitioner.h>

using namespace pgo;
using namespace pgo::SolidDeformationModel;
namespace ES = pgo::EigenSupport;

SegmentChainConstraintFunctions::SegmentChainConstraintFunctions(int nAll, int numPoints, const double *const points, int isCircle):
  NonlinearOptimization::ConstraintFunctions(nAll)
{
  restPositions = Eigen::Map<const ES::VXd>(points, numPoints * 3);
  if (isCircle)
    nele = numPoints;
  else
    nele = numPoints - 1;
  n = numPoints;

  restLengths.resize(nele);

  for (int vi = 0; vi < nele; vi++) {
    restLengths[vi] = (restPositions.segment<3>((vi + 1) % numPoints * 3) - restPositions.segment<3>(vi * 3)).norm();
  }

  std::vector<ES::TripletD> entries;
  for (int ei = 0; ei < nele; ei++) {
    entries.emplace_back(ei, ei * 3, 1);
    entries.emplace_back(ei, ei * 3 + 1, 1);
    entries.emplace_back(ei, ei * 3 + 2, 1);

    entries.emplace_back(ei, (ei + 1) % numPoints * 3, 1);
    entries.emplace_back(ei, (ei + 1) % numPoints * 3 + 1, 1);
    entries.emplace_back(ei, (ei + 1) % numPoints * 3 + 2, 1);
  }
  jacobianTemplate.resize(nele, nAll);
  jacobianTemplate.setFromTriplets(entries.begin(), entries.end());

  jacobianIndices.assign(nele, JacIndex::Constant(-1));
  for (int ei = 0; ei < nele; ei++) {
    JacIndex jacIdx;

    jacIdx[0] = ES::findEntryOffset(jacobianTemplate, ei, ei * 3);
    PGO_ALOG(jacIdx[0] >= 0);

    jacIdx[1] = ES::findEntryOffset(jacobianTemplate, ei, (ei + 1) % numPoints * 3);
    PGO_ALOG(jacIdx[1] >= 0);

    jacobianIndices[ei] = jacIdx;
  }

  entries.clear();
  for (int ei = 0; ei < nele; ei++) {
    for (int vi = 0; vi < 2; vi++) {
      for (int vj = 0; vj < 2; vj++) {
        for (int dofi = 0; dofi < 3; dofi++) {
          for (int dofj = 0; dofj < 3; dofj++) {
            entries.emplace_back((ei + vi) % numPoints * 3 + dofi, (ei + vj) % numPoints * 3 + dofj, 1);
          }
        }
      }
    }
  }

  lambdahTemplate.resize(nAll, nAll);
  lambdahTemplate.setFromTriplets(entries.begin(), entries.end());

  hessIndices.assign(nele, HessIndex::Constant(-1));
  for (int ei = 0; ei < nele; ei++) {
    HessIndex idx;

    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++) {
        for (int k = 0; k < 3; k++) {
          for (int l = 0; l < 3; l++) {
            int row = (ei + i) % numPoints * 3 + k;
            int col = (ei + j) % numPoints * 3 + l;

            int localRow = i * 3 + k;
            int localCol = j * 3 + l;

            idx(localRow, localCol) = ES::findEntryOffset(lambdahTemplate, row, col);
            PGO_ALOG(idx(localRow, localCol) >= 0);
          }
        }
      }
    }
    hessIndices[ei] = idx;
  }
  hessLocks = std::vector<tbb::spin_mutex>(lambdahTemplate.nonZeros());
}

SegmentChainConstraintFunctions::~SegmentChainConstraintFunctions()
{
}

// C_i = (x1 - x0)^2 - l^2
void SegmentChainConstraintFunctions::func(ES::ConstRefVecXd x, ES::RefVecXd g) const
{
  tbb::parallel_for(
    0, nele, [&](int ei) {
      ES::V3d v1 = x.segment<3>((ei + 1) % n * 3);
      ES::V3d v0 = x.segment<3>(ei * 3);
      g[ei] = (v1 - v0).squaredNorm() - restLengths[ei] * restLengths[ei];
    },
    tbb::static_partitioner());
}

// d C_i = (x1 - x0)d(x1 - x0)
// be careful about the order
void SegmentChainConstraintFunctions::jacobian(ES::ConstRefVecXd x, ES::SpMatD &jac) const
{
  tbb::parallel_for(
    0, nele, [&](int ei) {
      ES::V3d v1 = x.segment<3>((ei + 1) % n * 3);
      ES::V3d v0 = x.segment<3>(ei * 3);
      ES::V3d diff = v1 - v0;

      for (int dof = 0; dof < 3; dof++) {
        jac.valuePtr()[jacobianIndices[ei][0] + dof] = -diff[dof];  // x0
        jac.valuePtr()[jacobianIndices[ei][1] + dof] = diff[dof];   // x1
      }
    },
    tbb::static_partitioner());
}

// d^2 C_i = d(x1 - x0)d(x1 - x0)
void SegmentChainConstraintFunctions::hessian(ES::ConstRefVecXd x, ES::ConstRefVecXd lambda, ES::SpMatD &hess) const
{
  memset(hess.valuePtr(), 0, sizeof(double) * hess.nonZeros());

  tbb::parallel_for(
    0, nele, [&](int ei) {
      Eigen::Matrix<double, 6, 6> h;
      h.block<3, 3>(0, 0) = ES::M3d::Identity();   // d^2 C_i / dx0^2
      h.block<3, 3>(3, 3) = ES::M3d::Identity();   // d^2 C_i / dx1^2

      h.block<3, 3>(0, 3) = -ES::M3d::Identity();  // d^2 C_i / dx0 dx1
      h.block<3, 3>(3, 0) = -ES::M3d::Identity();  // d^2 C_i / dx0 dx1

      h *= lambda[ei];

      for (int col = 0; col < 2; col++) {
        for (int dofj = 0; dofj < 3; dofj++) {
          for (int row = 0; row < 2; row++) {
            for (int dofi = 0; dofi < 3; dofi++) {
              int localRow = row * 3 + dofi;
              int localCol = col * 3 + dofj;

              ES::IDX offset = hessIndices[ei](localRow, localCol);

              hessLocks[offset].lock();

              hess.valuePtr()[offset] += h(localRow, localCol);

              hessLocks[offset].unlock();
            }
          }
        }
      }
    },
    tbb::static_partitioner());
}