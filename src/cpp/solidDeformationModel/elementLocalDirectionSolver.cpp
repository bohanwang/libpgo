/*
author: Bohan Wang
copyright to USC,MIT
*/

#include "elementLocalDirectionSolver.h"

#include "tetMesh.h"
#include "objMesh.h"

#include "generateLaplacian.h"

#include "sparseMatrix.h"
#include "EigenSupport.h"
#include "EigenVegaTypeConverter.h"

#include "triMeshGeo.h"
#include "triMeshNeighbor.h"
#include "triMeshPseudoNormal.h"

#include "minimizeEnergy.h"
#include "potentialEnergies.h"
#include "quadraticPotentialEnergy.h"
#include "linearConstraintFunctions.h"

#include "pgoLogging.h"

#include <iostream>
#include <vector>

using namespace pgo;
using namespace pgo::NonlinearOptimization;
using namespace pgo::SolidDeformationModel;

namespace pgo::SolidDeformationModel
{
class ElementLocalDirectionSolverData
{
public:
  TriMeshRef surfaceMesh;
  ES::VXd surfaceNormals;
  // TriMeshNeighbor surfaceMeshNeighbors;
  TriMeshPseudoNormal surfaceMeshPseudoNormal;

  ES::SpMatD DirichletBC_C;
  ES::VXd DirichletBC_d;
  ES::SpMatD Neumann_C;
  ES::VXd Neumann_d;

  const TetMesh *tetMesh = nullptr;
  ES::SpMatD tetmeshL;
  std::vector<std::vector<int>> vertexNeighborTets;

  ES::VXd x;
  ES::VXd gradx;
};
}  // namespace pgo::SolidDeformationModel

ElementLocalDirectionSolver::ElementLocalDirectionSolver()
{
  data = new ElementLocalDirectionSolverData;
}

ElementLocalDirectionSolver::~ElementLocalDirectionSolver()
{
  delete data;
}

void ElementLocalDirectionSolver::setSurfaceMesh(size_t numVertices, const Vec3d *vertices,
  size_t numTriangles, const Vec3i *triangles, const Vec3d *normals)
{
  data->surfaceMesh = TriMeshRef((int)numVertices, vertices, (int)numTriangles, triangles);
  // data->surfaceMeshNeighbors = TriMeshNeighbor(data->surfaceMesh);

  if (normals) {
    data->surfaceNormals.resize(numVertices * 3);
    for (size_t vi = 0; vi < numVertices; vi++) {
      data->surfaceNormals.segment<3>(vi * 3) = ES::V3d(normals[vi][0], normals[vi][1], normals[vi][2]);
      data->surfaceNormals.segment<3>(vi * 3).normalize();
    }
  }
  else {
    data->surfaceMeshPseudoNormal.buildPseudoNormals(data->surfaceMesh);
    ObjMesh tempMesh((int)numVertices, (const double *)vertices, (int)numTriangles, (const int *)triangles);

    // compute surface normal based on the laplacian matrix of surface mesh
    std::shared_ptr<SparseMatrix> Lsurface(generateLaplacian(&tempMesh, true));
    ES::VXd positions[3] = {
      ES::VXd::Zero(numVertices),
      ES::VXd::Zero(numVertices),
      ES::VXd::Zero(numVertices),
    };

    for (size_t vi = 0; vi < numVertices; vi++) {
      positions[0][vi] = vertices[vi][0];
      positions[1][vi] = vertices[vi][1];
      positions[2][vi] = vertices[vi][2];
    }

    ES::VXd normals[3] = {
      ES::VXd::Zero(numVertices),
      ES::VXd::Zero(numVertices),
      ES::VXd::Zero(numVertices),
    };

    Lsurface->MultiplyVector(positions[0].data(), normals[0].data());
    Lsurface->MultiplyVector(positions[1].data(), normals[1].data());
    Lsurface->MultiplyVector(positions[2].data(), normals[2].data());

    data->surfaceNormals.resize(numVertices * 3);
    for (size_t vi = 0; vi < numVertices; vi++) {
      ES::V3d n1(normals[0][vi], normals[1][vi], normals[2][vi]);
      Vec3d n2 = data->surfaceMeshPseudoNormal.vtxNormal((int)vi);

      if (n1[0] * n2[0] + n1[1] * n2[1] + n1[2] * n2[2] < 0)
        n1 *= -1.0;

      data->surfaceNormals.segment<3>(vi * 3) = ES::V3d(n2[0], n2[1], n2[2]);
      data->surfaceNormals.segment<3>(vi * 3).normalize();
    }
  }
}

void ElementLocalDirectionSolver::setTetMesh(const TetMesh *tetmesh)
{
  data->tetMesh = tetmesh;

  // compute tet mesh laplacian operator
  std::shared_ptr<SparseMatrix> L(generateLaplacian(tetmesh, 1));
  ES::CreateSparseMatrix(L.get(), data->tetmeshL);

  data->vertexNeighborTets.resize(tetmesh->getNumVertices());
  for (int ei = 0; ei < tetmesh->getNumElements(); ei++) {
    for (int j = 0; j < 4; j++) {
      data->vertexNeighborTets[tetmesh->getVertexIndex(ei, j)].push_back(ei);
    }
  }
}

void ElementLocalDirectionSolver::setDirichletBoundaryCoundition(size_t numVertices, const int *vertexIDs, const double *vertexValues)
{
  if (data->tetMesh == nullptr) {
    LGE << "Tetmesh is not set.";
    exit(1);
  }

  std::vector<ES::TripletD> entries;
  for (size_t i = 0; i < numVertices; i++) {
    entries.emplace_back((int)i, vertexIDs[i], 1.0);
  }

  data->DirichletBC_C.resize(numVertices, data->tetMesh->getNumVertices());
  data->DirichletBC_C.setFromTriplets(entries.begin(), entries.end());
  data->DirichletBC_d = Eigen::Map<const ES::VXd>(vertexValues, numVertices);
}

void ElementLocalDirectionSolver::setNeumannBoundaryCoundition(size_t numVertices, const int *vertexIDs, const double *values)
{
  // this is for the surface mesh, we don't use it
  // grad_f = (f_j - f_i) (v_i - v_k)^90 / 2A + (f_k - f_i) (v_j - v_i)^90 / 2A

  // this for the volumetric mesh, we use this
  // grad_f = (f_j - f_i) (v_i - v_k) x (v_h - v_k) / 2V
  //        + (f_k - f_i) (v_i - v_h) x (v_j - v_h) / 2V
  //        + (f_h - f_i) (v_k - v_i) x (v_j - v_i) / 2V

  // build the matrix for constraints
  std::vector<ES::TripletD> entries;
  for (size_t i = 0; i < numVertices; i++) {
    int vid = vertexIDs[i];
    Vec3d n(data->surfaceNormals.data() + vid * 3);
    const auto &vertexNeighbor = data->vertexNeighborTets[vid];

    // compute all volume
    double volAll = 0;
    for (size_t ti = 0; ti < vertexNeighbor.size(); ti++) {
      volAll += data->tetMesh->getElementVolume(vertexNeighbor[ti]);
    }

    for (size_t ti = 0; ti < vertexNeighbor.size(); ti++) {
      int vertexIndices[4] = {
        data->tetMesh->getVertexIndex(vertexNeighbor[ti], 0),
        data->tetMesh->getVertexIndex(vertexNeighbor[ti], 1),
        data->tetMesh->getVertexIndex(vertexNeighbor[ti], 2),
        data->tetMesh->getVertexIndex(vertexNeighbor[ti], 3)
      };

      Vec3d p[4] = {
        data->tetMesh->getVertex(vertexIndices[0]),
        data->tetMesh->getVertex(vertexIndices[1]),
        data->tetMesh->getVertex(vertexIndices[2]),
        data->tetMesh->getVertex(vertexIndices[3]),
      };

      Vec3d vik_x_vhk = cross(p[0] - p[2], p[3] - p[2]);
      Vec3d vih_x_vjh = cross(p[0] - p[3], p[1] - p[3]);
      Vec3d vki_x_vji = cross(p[2] - p[0], p[1] - p[0]);

      entries.emplace_back((int)i, vertexIndices[1], dot(vik_x_vhk, n) * 0.5 / volAll);
      entries.emplace_back((int)i, vertexIndices[0], -dot(vik_x_vhk, n) * 0.5 / volAll);

      entries.emplace_back((int)i, vertexIndices[2], dot(vih_x_vjh, n) * 0.5 / volAll);
      entries.emplace_back((int)i, vertexIndices[0], -dot(vih_x_vjh, n) * 0.5 / volAll);

      entries.emplace_back((int)i, vertexIndices[3], dot(vki_x_vji, n) * 0.5 / volAll);
      entries.emplace_back((int)i, vertexIndices[0], -dot(vki_x_vji, n) * 0.5 / volAll);
    }
  }

  data->Neumann_C.resize(numVertices, data->tetMesh->getNumVertices());
  data->Neumann_C.setFromTriplets(entries.begin(), entries.end());

  data->Neumann_d = Eigen::Map<const ES::VXd>(values, numVertices);
}

void ElementLocalDirectionSolver::solve(SolverMode sm, double W0, double W1, double W2)
{
  // solve equation
  // L x = 0
  // st.
  //     f(boundary) = boundary values
  //     gradf (boundary) \dot n = boundary values

  if (sm == SM_PURE_PENALITY) {
    std::vector<ES::TripletD> entries;
    std::vector<double> rhs;
    int offset = 0;
    for (ES::IDX outeri = 0; outeri < data->tetmeshL.outerSize(); outeri++) {
      for (ES::SpMatD::InnerIterator it(data->tetmeshL, outeri); it; ++it) {
        entries.emplace_back((int)it.row(), (int)it.col(), it.value() * W0);
      }
      rhs.push_back(0);
    }

    offset = data->tetMesh->getNumVertices();
    if (data->DirichletBC_C.rows() > 0) {
      for (ES::IDX outeri = 0; outeri < data->DirichletBC_C.outerSize(); outeri++) {
        for (ES::SpMatD::InnerIterator it(data->DirichletBC_C, outeri); it; ++it) {
          entries.emplace_back((int)it.row() + offset, (int)it.col(), it.value() * W1);
        }
        rhs.push_back(data->DirichletBC_d[outeri] * W1);
      }
      offset += (int)data->DirichletBC_C.rows();
    }

    if (data->Neumann_C.rows() > 0) {
      for (ES::IDX outeri = 0; outeri < data->Neumann_C.outerSize(); outeri++) {
        for (ES::SpMatD::InnerIterator it(data->Neumann_C, outeri); it; ++it) {
          entries.emplace_back((int)it.row() + offset, (int)it.col(), it.value() * W2);
        }
        rhs.push_back(data->Neumann_d[outeri] * W2);
      }
      offset += (int)data->Neumann_C.rows();
    }

    PGO_ALOG(offset == (int)rhs.size());

    ES::SpMatD A(offset, data->tetMesh->getNumVertices());
    A.setFromTriplets(entries.begin(), entries.end());

    // this solution cannot guarantee the constraints
    ES::VXd b = Eigen::Map<ES::VXd>(rhs.data(), rhs.size());
    ES::VXd ATb = ES::VXd::Zero(data->tetMesh->getNumVertices());

    ES::SpMatD ATA;
    ES::mm(A, A, ATA, 1);
    ES::mv(A, b, ATb, 1);

    LGI << b.norm() << ',' << ATb.norm();

    ES::LUSolver solver;
    solver.analyzePattern(ATA);
    solver.factorize(ATA);
    data->x = solver.solve(ATb);

    LGI << "\n"
        << "||x||=" << data->x.norm() << '\n'
        << "||Cx-d||=" << (data->Neumann_C * data->x - data->Neumann_d).norm() << '\n'
        << "||x-b||=" << (data->DirichletBC_C * data->x - data->DirichletBC_d).norm() << '\n'
        << "||Lx||=" << (data->tetmeshL * data->x).norm();
  }
  else if (sm == SM_HARD_CONSTRAINTS) {
    std::shared_ptr<QuadraticPotentialEnergy> LE = std::make_shared<QuadraticPotentialEnergy>(data->tetmeshL, 1);
    std::shared_ptr<PotentialEnergies> energyAll = std::make_shared<PotentialEnergies>(data->tetMesh->getNumVertices());
    energyAll->addPotentialEnergy(LE);
    energyAll->init();

    ES::VXd d = ES::VXd::Zero(data->Neumann_C.rows());
    std::shared_ptr<LinearConstraintFunctions> constraints = std::make_shared<LinearConstraintFunctions>(data->Neumann_C, d);

    ES::VXd x = ES::VXd::Zero(data->tetMesh->getNumVertices());
    ES::VXd xlow = ES::VXd::Ones(data->tetMesh->getNumVertices()) * -1e20;
    ES::VXd xhi = ES::VXd::Ones(data->tetMesh->getNumVertices()) * 1e20;
    ES::VXd lambda(data->Neumann_d.size()), g(data->Neumann_d.size());

    for (ES::IDX outeri = 0; outeri < data->DirichletBC_C.outerSize(); outeri++) {
      for (ES::SpMatD::InnerIterator it(data->DirichletBC_C, outeri); it; ++it) {
        xlow[it.col()] = data->DirichletBC_d[it.row()];
        xhi[it.col()] = data->DirichletBC_d[it.row()];
      }
    }

    ES::VXd clow = data->Neumann_d;
    ES::VXd chi = data->Neumann_d;

    EnergyOptimizer::minimize(x, energyAll, xlow, xhi,
      lambda, g, constraints, clow, chi,
      EnergyOptimizer::SolverType::ST_KNITRO, 300, 1e-6, 5);
    data->x = x;
    LGI << "\n"
        << "||x||=" << data->x.norm() << '\n'
        << "||Cx-d||=" << g.norm() << '\n'
        << "||x-b||=" << (data->DirichletBC_C * data->x - data->DirichletBC_d).norm() << '\n'
        << "||Lx||=" << (data->tetmeshL * x).norm();
  }
  else if (sm == SM_HARD_LAPLACIAN) {
    ES::VXd d0 = data->Neumann_d * -1.0;
    std::shared_ptr<QuadraticPotentialEnergy> CE = std::make_shared<QuadraticPotentialEnergy>(data->Neumann_C, d0, 1);
    std::shared_ptr<PotentialEnergies> energyAll = std::make_shared<PotentialEnergies>(data->tetMesh->getNumVertices());
    energyAll->addPotentialEnergy(CE);
    energyAll->init();

    ES::VXd d1 = ES::VXd::Zero(data->tetMesh->getNumVertices());
    std::shared_ptr<LinearConstraintFunctions> constraints = std::make_shared<LinearConstraintFunctions>(data->tetmeshL, d1);

    ES::VXd x = ES::VXd::Zero(data->tetMesh->getNumVertices());
    ES::VXd xlow = ES::VXd::Ones(data->tetMesh->getNumVertices()) * -1e20;
    ES::VXd xhi = ES::VXd::Ones(data->tetMesh->getNumVertices()) * 1e20;
    ES::VXd lambda(data->Neumann_d.size()), g(data->Neumann_d.size());

    for (ES::IDX outeri = 0; outeri < data->DirichletBC_C.outerSize(); outeri++) {
      for (ES::SpMatD::InnerIterator it(data->DirichletBC_C, outeri); it; ++it) {
        xlow[it.col()] = data->DirichletBC_d[it.row()];
        xhi[it.col()] = data->DirichletBC_d[it.row()];
      }
    }

    ES::VXd clow = d1;
    ES::VXd chi = d1;

    EnergyOptimizer::minimize(x, energyAll, xlow, xhi,
      lambda, g, constraints, clow, chi,
      EnergyOptimizer::SolverType::ST_KNITRO, 300, 1e-6, 5);
    data->x = x;
    LGI << "\n"
        << "||x||=" << data->x.norm() << '\n'
        << "||Cx-d||=" << g.norm() << '\n'
        << "||x-b||=" << (data->DirichletBC_C * data->x - data->DirichletBC_d).norm() << '\n'
        << "||Lx||=" << (data->tetmeshL * x).norm();
  }
  else {
    LGE << "Wrong mode.";
    exit(1);
  }

  data->gradx.resize(data->x.size() * 3);
  for (int vertexID = 0; vertexID < data->tetMesh->getNumVertices(); vertexID++) {
    const auto &vertexNeighbor = data->vertexNeighborTets[vertexID];

    double volAll = 0;
    for (size_t ti = 0; ti < vertexNeighbor.size(); ti++) {
      volAll += data->tetMesh->getElementVolume(vertexNeighbor[ti]);
    }

    Vec3d gradfAll(0.0);
    for (size_t ti = 0; ti < vertexNeighbor.size(); ti++) {
      int vertexIndices[4] = {
        data->tetMesh->getVertexIndex(vertexNeighbor[ti], 0),
        data->tetMesh->getVertexIndex(vertexNeighbor[ti], 1),
        data->tetMesh->getVertexIndex(vertexNeighbor[ti], 2),
        data->tetMesh->getVertexIndex(vertexNeighbor[ti], 3)
      };

      Vec3d p[4] = {
        data->tetMesh->getVertex(vertexIndices[0]),
        data->tetMesh->getVertex(vertexIndices[1]),
        data->tetMesh->getVertex(vertexIndices[2]),
        data->tetMesh->getVertex(vertexIndices[3]),
      };

      Vec3d vik_x_vhk = cross(p[0] - p[2], p[3] - p[2]);
      Vec3d vih_x_vjh = cross(p[0] - p[3], p[1] - p[3]);
      Vec3d vki_x_vji = cross(p[2] - p[0], p[1] - p[0]);

      Vec3d gradf = vik_x_vhk * (data->x[vertexIndices[1]] - data->x[vertexIndices[0]]) +
        vih_x_vjh * (data->x[vertexIndices[2]] - data->x[vertexIndices[0]]) +
        vki_x_vji * (data->x[vertexIndices[3]] - data->x[vertexIndices[0]]);

      gradfAll += gradf * 0.5;
    }

    data->gradx.segment<3>(vertexID * 3) = ES::V3d(gradfAll[0], gradfAll[1], gradfAll[2]) / volAll;
  }
}

const double *ElementLocalDirectionSolver::getScalarField() const
{
  return data->x.data();
}

const double *ElementLocalDirectionSolver::getGradientField() const
{
  return data->gradx.data();
}
