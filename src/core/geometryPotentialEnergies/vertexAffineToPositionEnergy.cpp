#include "vertexAffineToPositionEnergy.h"

#include "EigenSupport.h"

#include <numeric>

using namespace pgo::PredefinedPotentialEnergies;
namespace ES = pgo::EigenSupport;

namespace pgo::PredefinedPotentialEnergies
{
struct VertexAffineToPositionEnergyBuf
{
  ES::SpMatD hess_orig, hess_orig11, hess_orig12, hess_orig22;
  ES::SpMatD H11C, CTH11C, CTH12, H21C;

  ES::SpMatI hess_orig11Map, hess_orig12Map, hess_orig22Map, H21CTransposeMap;
  ES::SpMatI CTH11CMap, CTH12Map, H21CMap, H22Map;

  ES::SymbolicMmData *mul_hc = nullptr, *mul_cthc = nullptr, *mul_cth = nullptr;

  ES::VXd grad_orig;
  ES::VXd x_orig;
};
}  // namespace pgo::PredefinedPotentialEnergies

VertexAffineToPositionEnergy::VertexAffineToPositionEnergy(int totalNumDOFs,
  const pgo::EigenSupport::VXd &restp, std::shared_ptr<const PotentialEnergy> origE):
  nAll(totalNumDOFs),
  restPosition(restp), originalEnergy(origE)
{
  buf = std::make_shared<VertexAffineToPositionEnergyBuf>();
  originalEnergy->createHessian(buf->hess_orig);
  buf->grad_orig.setZero(originalEnergy->getNumDOFs());
  buf->x_orig.setZero(originalEnergy->getNumDOFs());

  computeAToX();

  if (totalNumDOFs == (int)AToX.cols()) {
    nRestDOFs = 0;

    // Horig \in R^{3n x 3n}
    // AToX \in R^{3n x 12n}
    ES::symbolicMm(buf->hess_orig, AToX, buf->H11C, &buf->mul_hc);
    ES::symbolicMm(AToX, buf->H11C, hessTemplate, &buf->mul_cthc, 1);

    allDOFs.resize(AToX.cols());
    std::iota(allDOFs.begin(), allDOFs.end(), 0);
  }
  else {
    // have not debuged this part
    abort();

    nRestDOFs = nAll - (int)AToX.cols();

    std::vector<int> positionDOFs(AToX.rows());
    std::iota(positionDOFs.begin(), positionDOFs.end(), 0);

    allDOFs.resize(originalEnergy->getNumDOFs());
    std::iota(allDOFs.begin(), allDOFs.end(), 0);

    std::vector<int> restDOFs;
    std::set_difference(allDOFs.begin(), allDOFs.end(),
      positionDOFs.begin(), positionDOFs.end(), std::back_inserter(restDOFs));

    ES::selectRowsCols(buf->hess_orig, positionDOFs, buf->hess_orig11);
    ES::big2Small(buf->hess_orig, buf->hess_orig11, 0, 0, buf->hess_orig11Map);

    ES::selectRowsCols(buf->hess_orig, positionDOFs, restDOFs, buf->hess_orig12);
    ES::big2Small(buf->hess_orig, buf->hess_orig12, positionDOFs, restDOFs, buf->hess_orig12Map);

    ES::selectRowsCols(buf->hess_orig, restDOFs, restDOFs, buf->hess_orig22);
    ES::big2Small(buf->hess_orig, buf->hess_orig22, restDOFs, restDOFs, buf->hess_orig22Map);

    // [ h11 h12 ]
    // [ h21 h22 ]

    // [ ct h11 c,   ct h12 ]
    // [ h21 c   ,   h22    ]

    ES::symbolicMm(buf->hess_orig11, AToX, buf->H11C, &buf->mul_hc);
    ES::symbolicMm(AToX, buf->H11C, buf->CTH11C, &buf->mul_cthc, 1);

    ES::symbolicMm(AToX, buf->hess_orig12, buf->CTH12, &buf->mul_cth, 1);
    buf->H21C = buf->CTH12.transpose();
    ES::transposeMapping(buf->CTH12, buf->H21C, buf->H21CTransposeMap);

    std::vector<ES::TripletD> entries;
    for (ES::IDX rowi = 0; rowi < buf->CTH11C.rows(); rowi++) {
      for (ES::SpMatD::InnerIterator it(buf->CTH11C, rowi); it; ++it) {
        entries.emplace_back((int)it.row(), (int)it.col(), it.value());
      }
    }

    for (ES::IDX rowi = 0; rowi < buf->CTH12.rows(); rowi++) {
      for (ES::SpMatD::InnerIterator it(buf->CTH12, rowi); it; ++it) {
        entries.emplace_back((int)it.row(), (int)buf->CTH11C.cols() + (int)it.col(), it.value());
        entries.emplace_back((int)buf->CTH11C.cols() + (int)it.col(), (int)it.row(), it.value());
      }
    }

    for (ES::IDX rowi = 0; rowi < buf->hess_orig22.rows(); rowi++) {
      for (ES::SpMatD::InnerIterator it(buf->hess_orig22, rowi); it; ++it) {
        entries.emplace_back((int)buf->CTH11C.rows() + (int)it.row(), (int)buf->CTH11C.cols() + (int)it.col(), it.value());
      }
    }

    hessTemplate.resize(nAll, nAll);
    hessTemplate.setFromTriplets(entries.begin(), entries.end());

    positionDOFs.assign(AToX.cols(), 0);
    std::iota(positionDOFs.begin(), positionDOFs.end(), 0);

    allDOFs.assign(nAll, 0);
    std::iota(allDOFs.begin(), allDOFs.end(), 0);

    restDOFs.clear();
    std::set_difference(allDOFs.begin(), allDOFs.end(),
      positionDOFs.begin(), positionDOFs.end(), std::back_inserter(restDOFs));

    ES::small2Big(buf->CTH11C, hessTemplate, positionDOFs, positionDOFs, buf->CTH11CMap);
    ES::small2Big(buf->CTH12, hessTemplate, positionDOFs, restDOFs, buf->CTH12Map);
    ES::small2Big(buf->H21C, hessTemplate, restDOFs, positionDOFs, buf->H21CMap);
    ES::small2Big(buf->hess_orig22, hessTemplate, restDOFs, restDOFs, buf->H22Map);
  }
}

VertexAffineToPositionEnergy::~VertexAffineToPositionEnergy()
{
  ES::destroySymbolicMmData(buf->mul_hc);
  ES::destroySymbolicMmData(buf->mul_cthc);
  ES::destroySymbolicMmData(buf->mul_cth);
}

void VertexAffineToPositionEnergy::computeAToX()
{
  AToX.setZero();

  ES::V2i indices[] = {
    {0,  0},
    {0,  3},
    {0,  6},
    {0,  9},
    {1,  1},
    {1,  4},
    {1,  7},
    {1, 10},
    {2,  2},
    {2,  5},
    {2,  8},
    {2, 11},
  };

  std::vector<ES::TripletD> entries;
  for (int vi = 0; vi < (int)restPosition.size() / 3; vi++) {
    ES::V3d X = restPosition.segment<3>(vi * 3);
    // x = A * X
    // x = Z vec(A)
    Eigen::Matrix<double, 3, 12> Z;
    Z.setZero();

    // x0 = A00 X0 + A01 X1 + A02 X2 + A03
    Z(0, 0) = X(0), Z(0, 3) = X(1), Z(0, 6) = X(2), Z(0, 9) = 1.0;
    // x1 = A10 X0 + A11 X1 + A12 X2 + A13
    Z(1, 1) = X(0), Z(1, 4) = X(1), Z(1, 7) = X(2), Z(1, 10) = 1.0;
    // x2 = A20 X0 + A21 X1 + A22 X2 + A23
    Z(2, 2) = X(0), Z(2, 5) = X(1), Z(2, 8) = X(2), Z(2, 11) = 1.0;

    // Eigen::Matrix<double, 3, 4> A;
    // A.block<3, 3>(0, 0).setIdentity();
    // A.col(3).setZero();

    // ES::V3d x = A.block<3, 3>(0, 0) * X + A.col(3);
    // ES::V3d x1 = Z * ES::Mp<ES::V12d>(A.data());

    for (int i = 0; i < 12; i++) {
      int r = indices[i][0];
      int c = indices[i][1];
      entries.emplace_back(vi * 3 + r, vi * 12 + c, Z(r, c));
    }
  }
  AToX.resize(restPosition.size(), restPosition.size() * 4);
  AToX.setFromTriplets(entries.begin(), entries.end());
}

void VertexAffineToPositionEnergy::compute_x_orig(pgo::EigenSupport::ConstRefVecXd x, pgo::EigenSupport::RefVecXd x_orig) const
{
  // x = C A
  ES::mv(AToX, x.head(AToX.cols()), x_orig.head(AToX.rows()), 0);
  if (nRestDOFs) {
    x_orig.tail(nRestDOFs) = x.tail(nRestDOFs);
  }
}

double VertexAffineToPositionEnergy::func(pgo::EigenSupport::ConstRefVecXd x) const
{
  compute_x_orig(x, buf->x_orig);

  return originalEnergy->func(buf->x_orig);
}

void VertexAffineToPositionEnergy::gradient(pgo::EigenSupport::ConstRefVecXd x, pgo::EigenSupport::RefVecXd grad) const
{
  compute_x_orig(x, buf->x_orig);

  buf->grad_orig.setZero();
  originalEnergy->gradient(buf->x_orig, buf->grad_orig);

  ES::mv(AToX, buf->grad_orig.head(AToX.rows()), grad.head(AToX.cols()), 1);

  if (nRestDOFs) {
    grad.tail(nRestDOFs) = buf->grad_orig.tail(nRestDOFs);
  }
}

void VertexAffineToPositionEnergy::hessian(pgo::EigenSupport::ConstRefVecXd x, pgo::EigenSupport::SpMatD &hess) const
{
  compute_x_orig(x, buf->x_orig);

  std::memset(buf->hess_orig.valuePtr(), 0, sizeof(double) * buf->hess_orig.nonZeros());
  originalEnergy->hessian(buf->x_orig, buf->hess_orig);

  if (nRestDOFs == 0) {
    ES::mm(buf->hess_orig, AToX, buf->mul_hc, buf->H11C, 0);
    ES::mm(AToX, buf->H11C, buf->mul_cthc, hess, 1);
  }
  else {
    throw std::runtime_error("not done");
  }
}
