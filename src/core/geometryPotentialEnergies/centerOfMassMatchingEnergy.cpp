#include "centerOfMassMatchingEnergy.h"

#include <autodiff/reverse/var/eigen.hpp>
#include <autodiff/reverse/var.hpp>

#include <tbb/parallel_for.h>
#include <tbb/combinable.h>

#include <numeric>

using namespace pgo;
namespace ES = pgo::EigenSupport;
using namespace PredefinedPotentialEnergies;

namespace pgo::PredefinedPotentialEnergies
{
void tetMeshGeoToMatrices(const pgo::Mesh::TetMeshGeo &mesh, Eigen::MatrixXd &vtx, Eigen::MatrixXi &tet)
{
  vtx.resize(mesh.numVertices(), 3);
  tet.resize(mesh.numTets(), 4);

  for (int i = 0; i < mesh.numVertices(); i++) {
    const auto &v = mesh.pos(i);
    Eigen::Vector3d v3d(v[0], v[1], v[2]);
    vtx.row(i) = v3d;
  }

  for (int i = 0; i < mesh.numTets(); i++) {
    const auto &t = mesh.tet(i);
    tet.row(i) = Eigen::Vector4i(t[0], t[1], t[2], t[3]);
  }
}

autodiff::var dot(const autodiff::VectorXvar &a, const autodiff::VectorXvar &b)
{
  autodiff::var result = 0.0;
  for (int i = 0; i < a.size(); i++) {
    result += a(i) * b(i);
  }
  return result;
}

autodiff::Vector3var cross(const autodiff::Vector3var &a, const autodiff::Vector3var &b)
{
  autodiff::Vector3var result{ 0.0, 0.0, 0.0 };
  result(0) += a(1) * b(2) - a(2) * b(1);
  result(1) += a(2) * b(0) - a(0) * b(2);
  result(2) += a(0) * b(1) - a(1) * b(0);
  return result;
}

autodiff::var determinant(const autodiff::Vector3var &a, const autodiff::Vector3var &b, const autodiff::Vector3var &c, const autodiff::Vector3var &d)
{
  return dot(d - a, cross(b - a, c - a));
}

autodiff::Vector3var centerOfMassNumerator(const autodiff::Vector3var &a, const autodiff::Vector3var &b, const autodiff::Vector3var &c, const autodiff::Vector3var &d)
{
  // the numerator for each tet element
  autodiff::var value = determinant(a, b, c, d);
  autodiff::Vector3var vMean = (a + b + c + d) / 4.0;
  return value * vMean;
}

autodiff::var centerOfMassDenominator(const autodiff::Vector3var &a, const autodiff::Vector3var &b, const autodiff::Vector3var &c, const autodiff::Vector3var &d)
{
  // the denominator for each tet element
  autodiff::var value = determinant(a, b, c, d);
  return value;
}
}  // namespace pgo::PredefinedPotentialEnergies

CenterOfMassMatchingEnergy::CenterOfMassMatchingEnergy(const pgo::Mesh::TetMeshGeo &tetMesh, const ES::V3d &tgtCoM, const ES::M3d &projMat)
{
  ES::MXd vtxInit;
  ES::MXi tetInit;
  tetMeshGeoToMatrices(tetMesh, vtxInit, tetInit);

  m_tet = tetInit;
  m_tgtCoM = tgtCoM;
  m_projMat = projMat;
}

double CenterOfMassMatchingEnergy::func(EigenSupport::ConstRefVecXd x) const
{
  // Assuming ES::V3d and ES::MXd support reduction by using tbb::combinable
  tbb::combinable<ES::V3d> combNumerator([]() { return ES::V3d::Zero(); });
  tbb::combinable<double> combDenominator([]() { return 0.0; });

  tbb::parallel_for(tbb::blocked_range<int>(0, m_tet.rows()), [&](const tbb::blocked_range<int> &range) {
    for (int ele = range.begin(); ele < range.end(); ++ele) {
      ES::V3d a = x.segment<3>(m_tet(ele, 0) * 3);
      ES::V3d b = x.segment<3>(m_tet(ele, 1) * 3);
      ES::V3d c = x.segment<3>(m_tet(ele, 2) * 3);
      ES::V3d d = x.segment<3>(m_tet(ele, 3) * 3);

      autodiff::Vector3var aVar = a.cast<autodiff::var>();
      autodiff::Vector3var bVar = b.cast<autodiff::var>();
      autodiff::Vector3var cVar = c.cast<autodiff::var>();
      autodiff::Vector3var dVar = d.cast<autodiff::var>();

      autodiff::Vector3var valNumerator = centerOfMassNumerator(aVar, bVar, cVar, dVar);
      autodiff::var valDenominator = centerOfMassDenominator(aVar, bVar, cVar, dVar);

      combNumerator.local() += ES::V3d(autodiff::val(valNumerator(0)), autodiff::val(valNumerator(1)), autodiff::val(valNumerator(2)));
      combDenominator.local() += autodiff::val(valDenominator);
    });
    // Combine results from all threads
    ES::V3d numerator = combNumerator.combine(std::plus<>());
    double denominator = combDenominator.combine(std::plus<>());

    ES::V3d com = numerator / denominator;
    return (m_projMat * com - m_projMat * m_tgtCoM).squaredNorm();
}

void CenterOfMassMatchingEnergy::gradient(EigenSupport::ConstRefVecXd x, EigenSupport::RefVecXd grad) const
{
    // Assuming ES::V3d and ES::MXd support reduction by using tbb::combinable
    tbb::combinable<ES::V3d> combNumerator([]() { return ES::V3d::Zero(); });
    tbb::combinable<double> combDenominator([]() { return 0.0; });
    tbb::combinable<ES::MXd> combDNdxAll([&]() { return ES::MXd::Zero(3, x.size()); });
    tbb::combinable<ES::VXd> combDDdxAll([&]() { return ES::VXd::Zero(x.size()); });

    tbb::parallel_for(tbb::blocked_range<int>(0, m_tet.rows()), [&](const tbb::blocked_range<int> &range) {
      for (int ele = range.begin(); ele < range.end(); ++ele) {
        ES::V3d a = x.segment<3>(m_tet(ele, 0) * 3);
        ES::V3d b = x.segment<3>(m_tet(ele, 1) * 3);
        ES::V3d c = x.segment<3>(m_tet(ele, 2) * 3);
        ES::V3d d = x.segment<3>(m_tet(ele, 3) * 3);

        autodiff::Vector3var aVar = a.cast<autodiff::var>();
        autodiff::Vector3var bVar = b.cast<autodiff::var>();
        autodiff::Vector3var cVar = c.cast<autodiff::var>();
        autodiff::Vector3var dVar = d.cast<autodiff::var>();

        autodiff::Vector3var valNumerator = centerOfMassNumerator(aVar, bVar, cVar, dVar);
        autodiff::var valDenominator = centerOfMassDenominator(aVar, bVar, cVar, dVar);

        combNumerator.local() += ES::V3d(autodiff::val(valNumerator(0)), autodiff::val(valNumerator(1)), autodiff::val(valNumerator(2)));
        combDenominator.local() += autodiff::val(valDenominator);

        // backward
        for (int i = 0; i < 3; i++) {
          ES::V3d dNda = autodiff::gradient(valNumerator(i), aVar);
          ES::V3d dNdb = autodiff::gradient(valNumerator(i), bVar);
          ES::V3d dNdc = autodiff::gradient(valNumerator(i), cVar);
          ES::V3d dNdd = autodiff::gradient(valNumerator(i), dVar);

          for (int k = 0; k < 3; k++) {
            combDNdxAll.local()(i, m_tet(ele, 0) * 3 + k) += dNda(k);
            combDNdxAll.local()(i, m_tet(ele, 1) * 3 + k) += dNdb(k);
            combDNdxAll.local()(i, m_tet(ele, 2) * 3 + k) += dNdc(k);
            combDNdxAll.local()(i, m_tet(ele, 3) * 3 + k) += dNdd(k);
          }
        }

        ES::V3d dDda = autodiff::gradient(valDenominator, aVar);
        ES::V3d dDdb = autodiff::gradient(valDenominator, bVar);
        ES::V3d dDdc = autodiff::gradient(valDenominator, cVar);
        ES::V3d dDdd = autodiff::gradient(valDenominator, dVar);

        combDDdxAll.local().segment<3>(m_tet(ele, 0) * 3) += dDda;
        combDDdxAll.local().segment<3>(m_tet(ele, 1) * 3) += dDdb;
        combDDdxAll.local().segment<3>(m_tet(ele, 2) * 3) += dDdc;
        combDDdxAll.local().segment<3>(m_tet(ele, 3) * 3) += dDdd;
      }
    });
    // Combine results from all threads
    ES::V3d numerator = combNumerator.combine(std::plus<>());
    double denominator = combDenominator.combine(std::plus<>());
    ES::MXd dNdxAll = combDNdxAll.combine(std::plus<>());
    ES::VXd dDdxAll = combDDdxAll.combine(std::plus<>());
    ES::V3d com = numerator / denominator;

    // gradient
    ES::MXd gradP1 = 2.0 * m_projMat.transpose() * m_projMat * (com - m_tgtCoM);
    ES::MXd gradP2 = dNdxAll / denominator;
    ES::MXd gradP3 = numerator * dDdxAll.transpose();

    gradP3 /= (denominator * denominator);
    ES::MXd gradP23 = gradP2 - gradP3;
    grad.noalias() = gradP23.transpose() * gradP1;
}

void CenterOfMassMatchingEnergy::compute_com_and_energy_and_grad(ES::ConstRefVecXd x, ES::V3d &com,  double &energy, ES::RefVecXd grad) const
{
    // eng = (N / D - tgtCoM)^2
    // deng / dx = 2 * (N / D - tgtCoM) * d(N / D) / dx
    // dC / dx = 1 / D * dN / dx - N / D^2 * dD / dx
    ES::V3d numerator = ES::V3d::Zero();
    double denominator = 0.0;
    ES::MXd dNdxAll = ES::MXd::Zero(3, x.size());
    ES::VXd dDdxAll = ES::VXd::Zero(x.size());

    if (0) {
      for (int ele = 0; ele < m_tet.rows(); ele++) {
        // std::cout << "ele: " << ele << " / " << m_tet.rows() << std::endl;
        // forward
        ES::V3d a = x.segment<3>(m_tet(ele, 0) * 3);
        ES::V3d b = x.segment<3>(m_tet(ele, 1) * 3);
        ES::V3d c = x.segment<3>(m_tet(ele, 2) * 3);
        ES::V3d d = x.segment<3>(m_tet(ele, 3) * 3);

        autodiff::Vector3var aVar = a.cast<autodiff::var>();
        autodiff::Vector3var bVar = b.cast<autodiff::var>();
        autodiff::Vector3var cVar = c.cast<autodiff::var>();
        autodiff::Vector3var dVar = d.cast<autodiff::var>();

        autodiff::Vector3var valNumerator = centerOfMassNumerator(aVar, bVar, cVar, dVar);
        autodiff::var valDenominator = centerOfMassDenominator(aVar, bVar, cVar, dVar);

        numerator(0) += autodiff::val(valNumerator(0));
        numerator(1) += autodiff::val(valNumerator(1));
        numerator(2) += autodiff::val(valNumerator(2));
        denominator += autodiff::val(valDenominator);

        // backward
        // numerator
        for (int i = 0; i < 3; i++) {
          ES::V3d dNda = autodiff::gradient(valNumerator(i), aVar);
          ES::V3d dNdb = autodiff::gradient(valNumerator(i), bVar);
          ES::V3d dNdc = autodiff::gradient(valNumerator(i), cVar);
          ES::V3d dNdd = autodiff::gradient(valNumerator(i), dVar);
          for (int k = 0; k < 3; k++) {
            dNdxAll(i, m_tet(ele, 0) * 3 + k) += dNda(k);
            dNdxAll(i, m_tet(ele, 1) * 3 + k) += dNdb(k);
            dNdxAll(i, m_tet(ele, 2) * 3 + k) += dNdc(k);
            dNdxAll(i, m_tet(ele, 3) * 3 + k) += dNdd(k);
          }
        }
        // denominator
        ES::V3d dDda = autodiff::gradient(valDenominator, aVar);
        ES::V3d dDdb = autodiff::gradient(valDenominator, bVar);
        ES::V3d dDdc = autodiff::gradient(valDenominator, cVar);
        ES::V3d dDdd = autodiff::gradient(valDenominator, dVar);
        dDdxAll.segment<3>(m_tet(ele, 0) * 3) += dDda;
        dDdxAll.segment<3>(m_tet(ele, 1) * 3) += dDdb;
        dDdxAll.segment<3>(m_tet(ele, 2) * 3) += dDdc;
        dDdxAll.segment<3>(m_tet(ele, 3) * 3) += dDdd;
      }
    }
    else {
      // Assuming ES::V3d and ES::MXd support reduction by using tbb::combinable
      tbb::combinable<ES::V3d> combNumerator([]() { return ES::V3d::Zero(); });
      tbb::combinable<double> combDenominator([]() { return 0.0; });
      tbb::combinable<ES::MXd> combDNdxAll([&]() { return ES::MXd::Zero(3, x.size()); });
      tbb::combinable<ES::VXd> combDDdxAll([&]() { return ES::VXd::Zero(x.size()); });

      tbb::parallel_for(tbb::blocked_range<int>(0, m_tet.rows()), [&](const tbb::blocked_range<int> &range) {
        for (int ele = range.begin(); ele < range.end(); ++ele) {
          ES::V3d a = x.segment<3>(m_tet(ele, 0) * 3);
          ES::V3d b = x.segment<3>(m_tet(ele, 1) * 3);
          ES::V3d c = x.segment<3>(m_tet(ele, 2) * 3);
          ES::V3d d = x.segment<3>(m_tet(ele, 3) * 3);

          autodiff::Vector3var aVar = a.cast<autodiff::var>();
          autodiff::Vector3var bVar = b.cast<autodiff::var>();
          autodiff::Vector3var cVar = c.cast<autodiff::var>();
          autodiff::Vector3var dVar = d.cast<autodiff::var>();

          autodiff::Vector3var valNumerator = centerOfMassNumerator(aVar, bVar, cVar, dVar);
          autodiff::var valDenominator = centerOfMassDenominator(aVar, bVar, cVar, dVar);

          combNumerator.local() += ES::V3d(autodiff::val(valNumerator(0)), autodiff::val(valNumerator(1)), autodiff::val(valNumerator(2)));
          combDenominator.local() += autodiff::val(valDenominator);

          // backward
          for (int i = 0; i < 3; i++) {
            ES::V3d dNda = autodiff::gradient(valNumerator(i), aVar);
            ES::V3d dNdb = autodiff::gradient(valNumerator(i), bVar);
            ES::V3d dNdc = autodiff::gradient(valNumerator(i), cVar);
            ES::V3d dNdd = autodiff::gradient(valNumerator(i), dVar);

            for (int k = 0; k < 3; k++) {
              combDNdxAll.local()(i, m_tet(ele, 0) * 3 + k) += dNda(k);
              combDNdxAll.local()(i, m_tet(ele, 1) * 3 + k) += dNdb(k);
              combDNdxAll.local()(i, m_tet(ele, 2) * 3 + k) += dNdc(k);
              combDNdxAll.local()(i, m_tet(ele, 3) * 3 + k) += dNdd(k);
            }
          }

          ES::V3d dDda = autodiff::gradient(valDenominator, aVar);
          ES::V3d dDdb = autodiff::gradient(valDenominator, bVar);
          ES::V3d dDdc = autodiff::gradient(valDenominator, cVar);
          ES::V3d dDdd = autodiff::gradient(valDenominator, dVar);

          combDDdxAll.local().segment<3>(m_tet(ele, 0) * 3) += dDda;
          combDDdxAll.local().segment<3>(m_tet(ele, 1) * 3) += dDdb;
          combDDdxAll.local().segment<3>(m_tet(ele, 2) * 3) += dDdc;
          combDDdxAll.local().segment<3>(m_tet(ele, 3) * 3) += dDdd;
        }
      });
      // Combine results from all threads
      numerator = combNumerator.combine(std::plus<>());
      denominator = combDenominator.combine(std::plus<>());
      dNdxAll = combDNdxAll.combine(std::plus<>());
      dDdxAll = combDDdxAll.combine(std::plus<>());
    }

    com = numerator / denominator;
    // ES::V3d com = numerator;
    // ES::V3d proj = m_projMat * com;
    energy = (m_projMat * com - m_projMat * m_tgtCoM).squaredNorm();

    // gradient
    grad.resize(x.size());
    ES::MXd gradP1 = 2.0 * m_projMat.transpose() * m_projMat * (com - m_tgtCoM);
    ES::MXd gradP2 = dNdxAll / denominator;
    ES::MXd gradP3 = numerator * dDdxAll.transpose();
    gradP3 /= (denominator * denominator);
    ES::MXd gradP23 = gradP2 - gradP3;
    grad = gradP23.transpose() * gradP1;
    // ES::MXd gradP1 = 2.0 * m_projMat.transpose() * m_projMat * (com - m_tgtCoM);
    // grad = dNdxAll.transpose() * gradP1;

    // return energy;
}