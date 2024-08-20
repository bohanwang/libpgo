/*
author: Bohan Wang
copyright to USC,MIT
*/

#include "elementLocalDirection.h"
#include "elementLocalDirectionSolver.h"

#include "triMeshGeo.h"
#include "triMeshNeighbor.h"
#include "triMeshPseudoNormal.h"

#include "generateLaplacian.h"
#include "sparseMatrix.h"

#include "EigenSupport.h"
#include "EigenVegaTypeConverter.h"
#include "fastInterpolationWeights.h"
#include "tetMesh.h"
#include "potentialEnergies.h"
#include "quadraticPotentialEnergy.h"
#include "linearConstraintFunctions.h"
#include "minimizeEnergy.h"
#include "objMesh.h"
#include "pgoLogging.h"

#include <tbb/parallel_for.h>

#if defined(CGAL_PMP_USE_CERES_SOLVER)
#  undef CGAL_PMP_USE_CERES_SOLVER
#endif

#include "CGALInterface.h"
#include "CGALMeshProcessing.h"
#include "tetgenInterface.h"

#include <fmt/format.h>

#include <vector>
#include <memory>
#include <filesystem>
#include <fstream>

using namespace pgo;
using namespace pgo::NonlinearOptimization;
using namespace pgo::SolidDeformationModel;

namespace pgo::SolidDeformationModel
{
class ElementLocalDirectionData
{
public:
  const TetMesh *tetMesh = nullptr;
  std::vector<TriMeshRef> muscles;
  std::vector<TriMeshGeo> muscleClosedMeshes;
  std::vector<std::vector<Vec3d>> muscleVertexNormals;

  std::vector<Vec3d> elementDirections;
  std::vector<Vec3d> elementCenters;

  std::vector<Vec3d> vertexDirections;
  std::vector<double> vertexScalarField;

  std::string workingDir;
};
}  // namespace pgo::SolidDeformationModel

ElementLocalDirection::ElementLocalDirection()
{
  data = new ElementLocalDirectionData();
}

ElementLocalDirection::~ElementLocalDirection()
{
  delete data;
}

void ElementLocalDirection::setWorkingDir(const char *filename)
{
  data->workingDir = filename;
  std::filesystem::create_directories(data->workingDir);
}

void ElementLocalDirection::addMuscle(size_t numVertices, const Vec3d *vertices, size_t numTriangles, const Vec3i *triangles, const Vec3d *normals)
{
  data->muscles.emplace_back((int)numVertices, vertices, (int)numTriangles, triangles);

  if (normals)
    data->muscleVertexNormals.emplace_back(normals, normals + numVertices);
  else
    data->muscleVertexNormals.emplace_back();
}

void ElementLocalDirection::addWrapperTetMesh(const TetMesh *tetMesh)
{
  data->tetMesh = tetMesh;
}

Vec3d ElementLocalDirection::getFiberDirection(int elementID) const
{
  return data->elementDirections[elementID];
}

Vec3d ElementLocalDirection::getVertexFiberDirection(int vertexID) const
{
  return data->vertexDirections[vertexID];
}

void ElementLocalDirection::save(const char *filenamePrefix) const
{
  std::ofstream outfile(fmt::format("{}.ele.pos.txt", filenamePrefix).c_str());
  for (size_t i = 0; i < data->elementCenters.size(); i++) {
    outfile << data->elementCenters[i][0] << ',' << data->elementCenters[i][1] << ',' << data->elementCenters[i][2] << '\n';
  }
  outfile.close();
  outfile.clear();

  outfile.open(fmt::format("{}.ele.dir.txt", filenamePrefix).c_str());
  for (size_t i = 0; i < data->elementDirections.size(); i++) {
    outfile << data->elementDirections[i][0] << ',' << data->elementDirections[i][1] << ',' << data->elementDirections[i][2] << '\n';
  }
  outfile.close();
  outfile.clear();

  outfile.open(fmt::format("{}.vtx.dir.txt", filenamePrefix).c_str());
  for (size_t i = 0; i < data->vertexDirections.size(); i++) {
    outfile << data->vertexDirections[i][0] << ',' << data->vertexDirections[i][1] << ',' << data->vertexDirections[i][2] << '\n';
  }
  outfile.close();
  outfile.clear();

  outfile.open(fmt::format("{}.vtx.scale.txt", filenamePrefix).c_str());
  for (size_t i = 0; i < data->vertexScalarField.size(); i++) {
    outfile << data->vertexScalarField[i] << '\n';
  }
  outfile.close();
  outfile.clear();

  outfile.open(fmt::format("{}.vtx.pos.txt", filenamePrefix).c_str());
  for (int i = 0; i < data->tetMesh->getNumVertices(); i++) {
    Vec3d p = data->tetMesh->getVertex(i);
    outfile << p[0] << ',' << p[1] << ',' << p[2] << '\n';
  }
  outfile.close();
}

void ElementLocalDirection::saveSurfaceVertices(const char *filenamePrefix) const
{
  std::ofstream outfile(fmt::format("{}.surface.pos.txt", filenamePrefix).c_str());
  for (int i = 0; i < data->muscleClosedMeshes[0].numVertices(); i++) {
    Vec3d p = data->muscleClosedMeshes[0].pos(i);
    outfile << p[0] << ',' << p[1] << ',' << p[2] << '\n';
  }
  outfile.close();
  outfile.clear();

  outfile.open(fmt::format("{}.surface.dir.txt", filenamePrefix).c_str());
  for (int i = 0; i < data->muscleClosedMeshes[0].numVertices(); i++) {
    outfile << data->vertexDirections[i][0] << ',' << data->vertexDirections[i][1] << ',' << data->vertexDirections[i][2] << '\n';
  }
  outfile.close();
  outfile.clear();
};

int ElementLocalDirection::load(const char *filenamePrefix) const
{
  std::string filename = fmt::format("{}.ele.pos.txt", filenamePrefix);
  std::ifstream infile(filename.c_str());
  if (!infile) {
    LGE << "Cannot open file " << filename;
    return 1;
  }

  std::string line;
  data->elementCenters.clear();
  while (std::getline(infile, line) && line.length()) {
    Vec3d p;
    sscanf(line.c_str(), "%lf,%lf,%lf", &p[0], &p[1], &p[2]);
    data->elementCenters.push_back(p);
  }
  PGO_ALOG(data->elementCenters.size() == (size_t)data->tetMesh->getNumElements());
  infile.close();
  infile.clear();

  filename = fmt::format("{}.ele.dir.txt", filenamePrefix);
  infile.open(filename.c_str());
  if (!infile) {
    LGE << "Cannot open file " << filename;
    return 1;
  }

  data->elementDirections.clear();
  while (std::getline(infile, line) && line.length()) {
    Vec3d p;
    sscanf(line.c_str(), "%lf,%lf,%lf", &p[0], &p[1], &p[2]);
    data->elementDirections.push_back(p);
  }
  PGO_ALOG(data->elementDirections.size() == (size_t)data->tetMesh->getNumElements());
  infile.close();
  infile.clear();

  filename = fmt::format("{}.vtx.dir.txt", filenamePrefix);
  infile.open(filename.c_str());
  if (!infile) {
    LGE << "Cannot open file " << filename;
    return 1;
  }

  data->vertexDirections.clear();
  while (std::getline(infile, line) && line.length()) {
    Vec3d p;
    sscanf(line.c_str(), "%lf,%lf,%lf", &p[0], &p[1], &p[2]);
    data->vertexDirections.push_back(p);
  }
  PGO_ALOG(data->vertexDirections.size() == (size_t)data->tetMesh->getNumVertices());
  infile.close();
  infile.clear();

  return 0;
}

size_t ElementLocalDirection::getMuscleClosedMeshNumVertices(int mi) const
{
  return data->muscleClosedMeshes[mi].numVertices();
}

const Vec3d *ElementLocalDirection::getMuscleClosedMeshVertices(int mi) const
{
  return data->muscleClosedMeshes[mi].positions().data();
}

size_t ElementLocalDirection::getMuscleClosedMeshNumTriangles(int mi) const
{
  return data->muscleClosedMeshes[mi].numTriangles();
}

const Vec3i *ElementLocalDirection::getMuscleClosedMeshTriangles(int mi) const
{
  return data->muscleClosedMeshes[mi].triangles().data();
}

void ElementLocalDirection::compute()
{
  std::vector<double> samplePositions;
  std::vector<double> sampleDirections;
  std::vector<double> sampleScalars;

  Vec3d previousCenter;

  for (size_t mi = 0; mi < data->muscles.size(); mi++) {
    typedef CGALInterface::KernelInexact Krnl;
    typedef CGALInterface::CGALMesh<Krnl> CM;
    typedef CM::Polyhedron Poly;

    std::vector<std::set<Poly::Vertex_handle>> holeVertices;

    ObjMesh tempMuscle(data->muscles[mi].numVertices(), data->muscles[mi].positionBuffer(), data->muscles[mi].numTriangles(), data->muscles[mi].triangleBuffer());
    Poly muscleP;
    CM::ConvertObjMesh2Polyhedron(&tempMuscle, muscleP);

    for (Poly::Halfedge_iterator hit = muscleP.halfedges_begin(); hit != muscleP.halfedges_end(); ++hit) {
      if (hit->is_border()) {
        std::vector<Poly::Facet_handle> patch_facets;
        std::vector<Poly::Vertex_handle> patch_vertices;

        Poly::Halfedge_iterator hitCur = hit;
        while (1) {
          patch_vertices.push_back(hit->vertex());
          hitCur = hitCur->next();

          if (hitCur == hit)
            break;
        }

        bool success = std::get<0>(
          CGAL::Polygon_mesh_processing::triangulate_refine_and_fair_hole(
            muscleP, hit, std::back_inserter(patch_facets), std::back_inserter(patch_vertices),
            CGAL::Polygon_mesh_processing::parameters::vertex_point_map(CGAL::get(CGAL::vertex_point, muscleP))));

        // LGI << success;

        holeVertices.emplace_back(patch_vertices.begin(), patch_vertices.end());
      }
    }

    PGO_ALOG(holeVertices.size() == 2);

    int idx = 0;
    std::map<Poly::Vertex_const_handle, int> vertexIndices;
    std::vector<std::vector<int>> boundaryVertices(holeVertices.size());
    std::vector<int> otherVertices;

    std::vector<Vec3d> vertices;
    for (Poly::Vertex_iterator itt = muscleP.vertices_begin(); itt != muscleP.vertices_end(); ++itt) {
      Poly::Point pt = itt->point();
      vertices.emplace_back(
        CGAL::to_double(pt[0]),
        CGAL::to_double(pt[1]),
        CGAL::to_double(pt[2]));

      vertexIndices.insert(make_pair(itt, idx));

      bool isBoundary = false;
      for (size_t hi = 0; hi < holeVertices.size(); hi++) {
        if (holeVertices[hi].find(itt) != holeVertices[hi].end()) {
          boundaryVertices[hi].push_back(idx);
          isBoundary = true;
          break;
        }
      }

      if (isBoundary == false) {
        otherVertices.push_back(idx);
      }

      idx++;
    }

    std::vector<Vec3i> triangles;
    for (Poly::Face_iterator itt = muscleP.facets_begin(); itt != muscleP.facets_end(); ++itt) {
      Poly::Halfedge_around_facet_circulator ctr = itt->facet_begin();
      Vec3i tri;
      int fvid = 0;
      do {
        auto idxItt = vertexIndices.find(ctr->vertex());
        if (idxItt == vertexIndices.end()) {
          LOG(ERR) << "Cannot find the vertex";
          abort();
        }
        tri[fvid++] = idxItt->second;
        ctr++;
      } while (ctr != itt->facet_begin());

      triangles.push_back(tri);
    }

    std::shared_ptr<ObjMesh> newMesh;
    std::string muscleFilename = fmt::format("{}/mm{:02d}.obj", data->workingDir, mi);
    if (std::filesystem::exists(muscleFilename)) {
      newMesh = std::make_shared<ObjMesh>(muscleFilename);
      for (int i = 0; i < (int)newMesh->getNumVertices(); i++) {
        Vec3d p = newMesh->getPosition(i) * (1.0 / outputScale);
        newMesh->setPosition(i, p);
      }
    }
    else {
      newMesh = std::make_shared<ObjMesh>(vertices, triangles);

      std::shared_ptr<ObjMesh> outputMesh = std::make_shared<ObjMesh>(*newMesh);
      for (int i = 0; i < (int)outputMesh->getNumVertices(); i++) {
        Vec3d p = outputMesh->getPosition(i) * outputScale;
        outputMesh->setPosition(i, p);
      }
      outputMesh->save(muscleFilename);
    }

    std::shared_ptr<TetMesh> muscleTetMesh = Tetgen::Tetrahedralize(newMesh.get(), "pq1.1Y");

    std::vector<Vec3d> surfaceVertices;
    std::vector<Vec3i> surfaceTriangles;
    newMesh->exportGeometry(surfaceVertices, surfaceTriangles);
    data->muscleClosedMeshes.emplace_back(surfaceVertices, surfaceTriangles);

    ElementLocalDirectionSolver le;
    le.setSurfaceMesh(surfaceVertices.size(), surfaceVertices.data(), surfaceTriangles.size(), surfaceTriangles.data());
    le.setTetMesh(muscleTetMesh.get());

    std::vector<int> bvAll;
    std::vector<double> bvValues;

    if (mi == 0) {
      bvAll.insert(bvAll.end(), boundaryVertices[0].begin(), boundaryVertices[0].end());
      for (size_t i = 0; i < boundaryVertices[0].size(); i++)
        bvValues.push_back(0.0);

      bvAll.insert(bvAll.end(), boundaryVertices[1].begin(), boundaryVertices[1].end());
      for (size_t i = 0; i < boundaryVertices[1].size(); i++)
        bvValues.push_back(1.0);

      previousCenter = asVec3d(0.0);
      for (size_t i = 0; i < boundaryVertices[0].size(); i++) {
        previousCenter += surfaceVertices[boundaryVertices[0][i]];
      }

      previousCenter *= 1.0 / (double)boundaryVertices[0].size();
    }
    else {
      Vec3d curCenter0(0.0), curCenter1(0.0);
      for (size_t i = 0; i < boundaryVertices[0].size(); i++) {
        curCenter0 += surfaceVertices[boundaryVertices[0][i]];
      }
      curCenter0 *= 1.0 / (double)boundaryVertices[0].size();

      for (size_t i = 0; i < boundaryVertices[1].size(); i++) {
        curCenter1 += surfaceVertices[boundaryVertices[1][i]];
      }
      curCenter1 *= 1.0 / (double)boundaryVertices[1].size();

      if (dot(curCenter0 - previousCenter, curCenter0 - previousCenter) <
        dot(curCenter1 - previousCenter, curCenter1 - previousCenter)) {
        bvAll.insert(bvAll.end(), boundaryVertices[0].begin(), boundaryVertices[0].end());
        for (size_t i = 0; i < boundaryVertices[0].size(); i++)
          bvValues.push_back(0.0);

        bvAll.insert(bvAll.end(), boundaryVertices[1].begin(), boundaryVertices[1].end());
        for (size_t i = 0; i < boundaryVertices[1].size(); i++)
          bvValues.push_back(1.0);
      }
      else {
        bvAll.insert(bvAll.end(), boundaryVertices[0].begin(), boundaryVertices[0].end());
        for (size_t i = 0; i < boundaryVertices[0].size(); i++)
          bvValues.push_back(1.0);

        bvAll.insert(bvAll.end(), boundaryVertices[1].begin(), boundaryVertices[1].end());
        for (size_t i = 0; i < boundaryVertices[1].size(); i++)
          bvValues.push_back(0.0);
      }
    }

    le.setDirichletBoundaryCoundition(bvAll.size(), bvAll.data(), bvValues.data());

    std::vector<double> otherValues(otherVertices.size(), 0.0);
    le.setNeumannBoundaryCoundition(otherVertices.size(), otherVertices.data(), otherValues.data());

    LGI << boundaryVertices[0].size() + boundaryVertices[1].size() + otherValues.size() << "=" << newMesh->getNumVertices();

    le.solve(ElementLocalDirectionSolver::SM_PURE_PENALITY);

    for (int i = 0; i < muscleTetMesh->getNumVertices(); i++) {
      Vec3d gradf(le.getGradientField() + i * 3);
      if (dot(gradf, gradf) > 1e-9) {
        sampleScalars.push_back(le.getScalarField()[i]);

        Vec3d p = muscleTetMesh->getVertex(i);
        samplePositions.push_back(p[0]);
        samplePositions.push_back(p[1]);
        samplePositions.push_back(p[2]);

        Vec3d d = norm(gradf);
        sampleDirections.push_back(d[0]);
        sampleDirections.push_back(d[1]);
        sampleDirections.push_back(d[2]);
      }
    }
  }

  LGI << "Sampling Done.";

  std::ofstream outfile(fmt::format("{}/sample-positions.txt", data->workingDir).c_str());
  for (size_t i = 0; i < samplePositions.size() / 3; i++) {
    outfile << samplePositions[i * 3] << ',' << samplePositions[i * 3 + 1] << ',' << samplePositions[i * 3 + 2] << '\n';
  }
  outfile.close();

  std::ofstream outfile1(fmt::format("{}/sample-dirs.txt", data->workingDir).c_str());
  for (size_t i = 0; i < sampleDirections.size() / 3; i++) {
    outfile1 << sampleDirections[i * 3] << ',' << sampleDirections[i * 3 + 1] << ',' << sampleDirections[i * 3 + 2] << '\n';
  }
  outfile1.close();

  std::ofstream outfile2(fmt::format("{}/sample-scalars.txt", data->workingDir).c_str());
  for (size_t i = 0; i < sampleDirections.size() / 3; i++) {
    outfile2 << sampleScalars[i] << '\n';
  }
  outfile2.close();

  // exit(1);
  LGI << "Re-interpolating the directions into an embedding mesh...";

  // exit(1);
  std::shared_ptr<SparseMatrix> L(generateLaplacian(data->tetMesh, 1));
  ES::SpMatD Lsp;
  ES::CreateSparseMatrix(L.get(), Lsp);

  std::shared_ptr<FastInterpolationWeights> fiw = std::make_shared<FastInterpolationWeights>(data->tetMesh);

  int *eleVertices = nullptr;
  double *eleVertexWeights = nullptr;
  fiw->generateInterpolationWeights((int)samplePositions.size() / 3, samplePositions.data(), &eleVertices, &eleVertexWeights);

  std::vector<ES::TripletD> entries;
  for (int i = 0; i < (int)samplePositions.size() / 3; i++) {
    for (int j = 0; j < 4; j++)
      entries.emplace_back(i, eleVertices[i * 4 + j], eleVertexWeights[i * 4 + j]);
  }

  ES::SpMatD W((int)samplePositions.size() / 3, data->tetMesh->getNumVertices());
  W.setFromTriplets(entries.begin(), entries.end());

  ES::VXd boundaryConditions[3] = {
    ES::VXd(samplePositions.size() / 3),
    ES::VXd(samplePositions.size() / 3),
    ES::VXd(samplePositions.size() / 3)
  };

  for (int i = 0; i < (int)samplePositions.size() / 3; i++) {
    boundaryConditions[0][i] = -sampleDirections[i * 3 + 0];
    boundaryConditions[1][i] = -sampleDirections[i * 3 + 1];
    boundaryConditions[2][i] = -sampleDirections[i * 3 + 2];
  }

  std::shared_ptr<QuadraticPotentialEnergy> LE = std::make_shared<QuadraticPotentialEnergy>(Lsp, 1);

  ES::VXd directions(data->tetMesh->getNumVertices() * 3);

  for (int i = 0; i < 3; i++) {
    std::shared_ptr<QuadraticPotentialEnergy> CE = std::make_shared<QuadraticPotentialEnergy>(W, boundaryConditions[i], 1);
    std::shared_ptr<PotentialEnergies> energyAll = std::make_shared<PotentialEnergies>((int)Lsp.rows());
    energyAll->addPotentialEnergy(LE);
    energyAll->addPotentialEnergy(CE, 100000);
    energyAll->init();

    ES::VXd x = ES::VXd::Zero(data->tetMesh->getNumVertices());
    ES::VXd xlow = ES::VXd::Ones(data->tetMesh->getNumVertices()) * -1e20;
    ES::VXd xhi = ES::VXd::Ones(data->tetMesh->getNumVertices()) * 1e20;
    ES::VXd lambda, g;

    EnergyOptimizer::minimize(x, energyAll, xlow, xhi,
      lambda, g, nullptr, ES::VXd(), ES::VXd(),
      EnergyOptimizer::SolverType::ST_KNITRO, 300, 1e-6, 5);

    LGI << "dir " << i << ": ||Wx-d||=" << (W * x + boundaryConditions[i]).norm();
    for (ES::IDX j = 0; j < x.size(); j++)
      directions[j * 3 + i] = x[j];
  }

  for (int i = 0; i < data->tetMesh->getNumVertices(); i++) {
    Vec3d dir(directions[i * 3], directions[i * 3 + 1], directions[i * 3 + 2]);
    data->vertexDirections.push_back(norm(dir));
  }

  for (int i = 0; i < data->tetMesh->getNumElements(); i++) {
    ES::V3d dir = ES::V3d::Zero();
    Vec3d p(0.0);
    for (int j = 0; j < 4; j++) {
      dir += directions.segment<3>(data->tetMesh->getVertexIndex(i, j) * 3);
      p += data->tetMesh->getVertex(i, j);
    }

    dir *= 0.25;
    p *= 0.25;

    dir.normalize();
    data->elementDirections.emplace_back(dir[0], dir[1], dir[2]);

    data->elementCenters.push_back(p);
  }
}

void ElementLocalDirection::computeNonEmbedding(size_t numZeroVertices, const int *zeroVertices,
  size_t numOneVertices, const int *oneVertices, const double weights[3])
{
  PGO_ALOG(data->muscles.size() == 1);

  data->elementCenters.resize(data->tetMesh->getNumElements());
  data->elementDirections.resize(data->tetMesh->getNumElements());
  data->vertexDirections.resize(data->tetMesh->getNumVertices());

  std::vector<int> vertexFlags(data->tetMesh->getNumVertices(), 0);

  data->muscleClosedMeshes.emplace_back(data->muscles[0].numVertices(), data->muscles[0].positions(),
    data->muscles[0].numTriangles(), data->muscles[0].triangles());

  const std::vector<Vec3d> &surfaceVertices = data->muscleClosedMeshes[0].positions();
  const std::vector<Vec3i> &surfaceTriangles = data->muscleClosedMeshes[0].triangles();

  TriMeshPseudoNormal surfaceNormals;
  surfaceNormals.buildPseudoNormals(data->muscleClosedMeshes[0].ref());

  ElementLocalDirectionSolver le;
  le.setSurfaceMesh(surfaceVertices.size(), surfaceVertices.data(),
    surfaceTriangles.size(), surfaceTriangles.data(),
    data->muscleVertexNormals[0].size() ? data->muscleVertexNormals[0].data() : nullptr);
  le.setTetMesh(data->tetMesh);

  std::vector<int> bvAll;
  std::vector<double> bvValues;

  for (size_t i = 0; i < numZeroVertices; i++) {
    bvAll.push_back(zeroVertices[i]);
    bvValues.push_back(0.0);
  }

  for (size_t i = 0; i < numOneVertices; i++) {
    bvAll.push_back(oneVertices[i]);
    bvValues.push_back(1.0);
  }

  le.setDirichletBoundaryCoundition(bvAll.size(), bvAll.data(), bvValues.data());

  std::vector<int> endVertices(bvAll);
  std::sort(endVertices.begin(), endVertices.end());

  std::vector<int> allVertices(surfaceVertices.size());
  std::iota(allVertices.begin(), allVertices.end(), 0);

  std::vector<int> otherVertices;
  std::set_difference(allVertices.begin(), allVertices.end(), endVertices.begin(), endVertices.end(), std::back_inserter(otherVertices));

  std::vector<double> otherValues(otherVertices.size(), 0.0);
  le.setNeumannBoundaryCoundition(otherVertices.size(), otherVertices.data(), otherValues.data());

  LGI << endVertices.size() + otherVertices.size() << "=" << surfaceVertices.size();

  if (weights) {
    le.solve(ElementLocalDirectionSolver::SM_PURE_PENALITY, weights[0], weights[1], weights[2]);
  }
  else {
    le.solve(ElementLocalDirectionSolver::SM_PURE_PENALITY, 1, 100, 10);
  }

  for (int i = 0; i < data->tetMesh->getNumVertices(); i++) {
    Vec3d gradf(le.getGradientField() + i * 3);
    if (dot(gradf, gradf) > 1e-9) {
      vertexFlags[i] = 1;
      data->vertexDirections[i] = norm(gradf);
    }
  }

  data->vertexScalarField.assign(le.getScalarField(), le.getScalarField() + data->tetMesh->getNumVertices());

  std::vector<std::vector<int>> vertexNeighboringElements(data->tetMesh->getNumVertices());
  for (int ele = 0; ele < data->tetMesh->getNumElements(); ele++) {
    for (int j = 0; j < 4; j++)
      vertexNeighboringElements[data->tetMesh->getVertexIndex(ele, j)].push_back(ele);
  }

  std::vector<std::vector<int>> vertexNeighboringVertices(data->tetMesh->getNumVertices());
  tbb::parallel_for(0, data->tetMesh->getNumVertices(), [&](int i) {
    if (vertexFlags[i] == 0) {
      for (int ele : vertexNeighboringElements[i]) {
        for (int j = 0; j < 4; j++)
          vertexNeighboringVertices[i].push_back(data->tetMesh->getVertexIndex(ele, j));
      }
      std::sort(vertexNeighboringVertices[i].begin(), vertexNeighboringVertices[i].end());
      auto iter = std::unique(vertexNeighboringVertices[i].begin(), vertexNeighboringVertices[i].end());
      vertexNeighboringVertices[i].erase(iter, vertexNeighboringVertices[i].end());
    }
  });

  int numTotalIter = 0;
  while (numTotalIter < 100) {
    ptrdiff_t numValidVertices = std::count_if(vertexFlags.begin(), vertexFlags.end(), [](int flag) { return flag == 1; });
    ptrdiff_t numTotalVertices = (ptrdiff_t)vertexFlags.size();
    LGI << numValidVertices << "/" << numTotalVertices;

    if (numValidVertices >= numTotalVertices)
      break;

    for (int i = 0; i < data->tetMesh->getNumVertices(); i++) {
      if (vertexFlags[i] == 0) {
        const auto &vtxNeighbors = vertexNeighboringVertices[i];
        Vec3d d(0.0);
        int validCount = 0;

        for (int vtxID : vtxNeighbors) {
          if (vertexFlags[vtxID]) {
            d += data->vertexDirections[vtxID];
            validCount++;
          }
        }

        if (validCount) {
          if (i < (int)surfaceVertices.size()) {
            Vec3d n = surfaceNormals.vtxNormal(i);
            d -= dot(n, d) * n;
          }

          data->vertexDirections[i] = norm(d);
          vertexFlags[i] = 1;
        }
      }
    }
    numTotalIter++;
  }

  for (int i = 0; i < data->tetMesh->getNumElements(); i++) {
    Vec3d dir(0.0);
    Vec3d p(0.0);
    for (int j = 0; j < 4; j++) {
      int vid = data->tetMesh->getVertexIndex(i, j);
      if (vertexFlags[vid]) {
        dir += data->vertexDirections[vid];
        p += data->tetMesh->getVertex(i, j);
      }
    }

    dir *= 0.25;
    p *= 0.25;

    if (len(dir) < 1e-8) {
      data->elementDirections[i] = Vec3d(0, 1, 0);
      data->elementCenters[i] = p;
    }
    else {
      dir.normalize();
      data->elementDirections[i] = dir;
      data->elementCenters[i] = p;
    }
  }
}

/*


    // compute tet mesh laplacian operator
    std::shared_ptr<SparseMatrix> L(generateLaplacian(muscleTetMesh.get(), 1));
    ES::SpMatD Lsp;
    ES::CreateSparseMatrix(L.get(), Lsp);

    // compute surface normal based on the laplacian matrix of surface mesh
    std::shared_ptr<SparseMatrix> Lsurface(generateLaplacian(newMesh.get(), true));
    ES::VXd positions[3] = {
      ES::VXd::Zero(vertices.size()),
      ES::VXd::Zero(vertices.size()),
      ES::VXd::Zero(vertices.size()),
    };

    for (size_t vi = 0; vi < vertices.size(); vi++) {
      positions[0][vi] = vertices[vi][0];
      positions[1][vi] = vertices[vi][1];
      positions[2][vi] = vertices[vi][2];
    }

    ES::VXd normals[3] = {
      ES::VXd::Zero(vertices.size()),
      ES::VXd::Zero(vertices.size()),
      ES::VXd::Zero(vertices.size()),
    };

    Lsurface->MultiplyVector(positions[0].data(), normals[0].data());
    Lsurface->MultiplyVector(positions[1].data(), normals[1].data());
    Lsurface->MultiplyVector(positions[2].data(), normals[2].data());

    ES::VXd surfaceNormals(vertices.size() * 3);
    for (size_t vi = 0; vi < vertices.size(); vi++) {
      surfaceNormals.segment<3>(vi * 3) = ES::V3d(normals[0][vi], normals[1][vi], normals[2][vi]);
      surfaceNormals.segment<3>(vi * 3).normalize();
    }

    TriMeshRef newMeshRef(vertices, triangles);
    TriMeshNeighbor neighbors(newMeshRef);
    TriMeshPseudoNormal newMeshPseudoNormal(newMeshRef);

    std::vector<std::vector<int>> vertexNeighborTets(muscleTetMesh->getNumVertices());
    for (int ei = 0; ei < muscleTetMesh->getNumElements(); ei++) {
      for (int j = 0; j < 4; j++) {
        vertexNeighborTets[muscleTetMesh->getVertexIndex(ei, j)].push_back(ei);
      }
    }
    // this is for the surface mesh, we don't use it
    // grad_f = (f_j - f_i) (v_i - v_k)^90 / 2A + (f_k - f_i) (v_j - v_i)^90 / 2A

    // this for the volumetric mesh, we use this
    // grad_f = (f_j - f_i) (v_i - v_k) x (v_h - v_k) / 2V
    //        + (f_k - f_i) (v_i - v_h) x (v_j - v_h) / 2V
    //        + (f_h - f_i) (v_k - v_i) x (v_j - v_i) / 2V
    int indexOrders[4][4] = {
      { 0, 1, 2, 3 },
      { 1, 2, 0, 3 },
      { 2, 0, 1, 3 },
      { 3, 1, 0, 2 }
    };
    // build the matrix for constraints
    std::vector<ES::TripletD> entries;
    for (size_t i = 0; i < otherVertices.size(); i++) {
      int vertexID = otherVertices[i];
      Vec3d n(surfaceNormals.data() + vertexID * 3);
      const auto &vertexNeighbor = vertexNeighborTets[vertexID];
      double volAll = 0;
      for (size_t ti = 0; ti < vertexNeighbor.size(); ti++) {
        volAll += muscleTetMesh->getElementVolume(vertexNeighbor[ti]);
      }

      for (size_t ti = 0; ti < vertexNeighbor.size(); ti++) {
        // find the start index, then we permute the entire vertex indices
        int startIndex = 0;
        for (; startIndex < 4; startIndex++) {
          if (muscleTetMesh->getVertexIndex(vertexNeighbor[ti], startIndex) == vertexID)
            break;
        }
        PGO_ALOG(startIndex < 4);

        int vertexIndices[4] = {
          muscleTetMesh->getVertexIndex(vertexNeighbor[ti], indexOrders[startIndex][0]),
          muscleTetMesh->getVertexIndex(vertexNeighbor[ti], indexOrders[startIndex][1]),
          muscleTetMesh->getVertexIndex(vertexNeighbor[ti], indexOrders[startIndex][2]),
          muscleTetMesh->getVertexIndex(vertexNeighbor[ti], indexOrders[startIndex][3])
        };

        Vec3d p[4] = {
          muscleTetMesh->getVertex(vertexIndices[0]),
          muscleTetMesh->getVertex(vertexIndices[1]),
          muscleTetMesh->getVertex(vertexIndices[2]),
          muscleTetMesh->getVertex(vertexIndices[3]),
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

    ES::SpMatD C(otherVertices.size(), muscleTetMesh->getNumVertices());
    C.setFromTriplets(entries.begin(), entries.end());

    std::shared_ptr<QuadraticPotentialEnergy> LE = std::make_shared<QuadraticPotentialEnergy>(Lsp, 1);
    std::shared_ptr<QuadraticPotentialEnergy> CE = std::make_shared<QuadraticPotentialEnergy>(C, 1);

    std::shared_ptr<PotentialEnergies> energyAll = std::make_shared<PotentialEnergies>(Lsp.rows());
    energyAll->addPotentialEnergy(LE);
    // energyAll->addPotentialEnergy(CE, 10000.0);
    energyAll->init();

    ES::VXd d = ES::VXd::Zero(C.rows());
    std::shared_ptr<LinearConstraintFunctions> constraints = std::make_shared<LinearConstraintFunctions>(C, d);


    ES::VXd x = ES::VXd::Zero(muscleTetMesh->getNumVertices());
    ES::VXd xlow = ES::VXd::Ones(muscleTetMesh->getNumVertices()) * -1e20;
    ES::VXd xhi = ES::VXd::Ones(muscleTetMesh->getNumVertices()) * 1e20;
    ES::VXd lambda, g;

    for (int vid : boundaryVertices[0]) {
      for (int eleID : vertexNeighborTets[vid]) {
        for (int j = 0; j < 4; j++) {
          xlow[muscleTetMesh->getVertexIndex(eleID, j)] = 0;
          xhi[muscleTetMesh->getVertexIndex(eleID, j)] = 0;
        }
      }
    }

    for (int vid : boundaryVertices[1]) {
      for (int eleID : vertexNeighborTets[vid]) {
        for (int j = 0; j < 4; j++) {
          xlow[muscleTetMesh->getVertexIndex(eleID, j)] = 1;
          xhi[muscleTetMesh->getVertexIndex(eleID, j)] = 1;
        }
      }
    }

    EnergyOptimizer::minimize(x, energyAll, xlow, xhi,
      lambda, g, constraints, d, d,
      EnergyOptimizer::ST_IPOPT, 300, 1e-6, 5);

    for (int vertexID = 0; vertexID < muscleTetMesh->getNumVertices(); vertexID++) {
      const auto &vertexNeighbor = vertexNeighborTets[vertexID];
      double volAll = 0;
      for (size_t ti = 0; ti < vertexNeighbor.size(); ti++) {
        volAll += muscleTetMesh->getElementVolume(vertexNeighbor[ti]);
      }

      Vec3d gradfAll(0.0);

      for (size_t ti = 0; ti < vertexNeighbor.size(); ti++) {
        // find the start index, then we permute the entire vertex indices
        int startIndex = 0;
        for (; startIndex < 4; startIndex++) {
          if (muscleTetMesh->getVertexIndex(vertexNeighbor[ti], startIndex) == vertexID)
            break;
        }
        PGO_ALOG(startIndex < 4);

        int vertexIndices[4] = {
          muscleTetMesh->getVertexIndex(vertexNeighbor[ti], indexOrders[startIndex][0]),
          muscleTetMesh->getVertexIndex(vertexNeighbor[ti], indexOrders[startIndex][1]),
          muscleTetMesh->getVertexIndex(vertexNeighbor[ti], indexOrders[startIndex][2]),
          muscleTetMesh->getVertexIndex(vertexNeighbor[ti], indexOrders[startIndex][3])
        };

        Vec3d p[4] = {
          muscleTetMesh->getVertex(vertexIndices[0]),
          muscleTetMesh->getVertex(vertexIndices[1]),
          muscleTetMesh->getVertex(vertexIndices[2]),
          muscleTetMesh->getVertex(vertexIndices[3]),
        };

        Vec3d vik_x_vhk = cross(p[0] - p[2], p[3] - p[2]);
        Vec3d vih_x_vjh = cross(p[0] - p[3], p[1] - p[3]);
        Vec3d vki_x_vji = cross(p[2] - p[0], p[1] - p[0]);

        Vec3d gradf = vik_x_vhk * (x[vertexIndices[1]] - x[vertexIndices[0]]) +
          vih_x_vjh * (x[vertexIndices[2]] - x[vertexIndices[0]]) +
          vki_x_vji * (x[vertexIndices[3]] - x[vertexIndices[0]]);

        gradfAll += gradf * 0.5;
      }

      if (dot(gradfAll, gradfAll) > 1e-8) {
        Vec3d p = muscleTetMesh->getVertex(vertexID);
        samplePositions.push_back(p[0]);
        samplePositions.push_back(p[1]);
        samplePositions.push_back(p[2]);

        Vec3d d = norm(gradfAll);
        sampleDirections.push_back(d[0]);
        sampleDirections.push_back(d[1]);
        sampleDirections.push_back(d[2]);

        sampleScalars.push_back(x[vertexID]);
      }
    }
*/
