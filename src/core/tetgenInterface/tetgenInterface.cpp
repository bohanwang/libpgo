#include "tetgenInterface.h"

#include "triMeshGeo.h"

#include <tetgen.h>

#include <fmt/format.h>

#include <tuple>
#include <vector>
#include <chrono>

int pgo::TetgenInterface::computeVoronoiDiagram(const std::vector<EigenSupport::V3d> &points, VoronoiDiagram &vd)
{
  namespace ES = EigenSupport;

  tetgenio tin, tout;

  // assign vertices;
  tin.numberofpoints = int(points.size());
  tin.pointlist = new double[tin.numberofpoints * 3];
  for (int i = 0; i < (int)points.size(); i++) {
    (tin.pointlist + i * 3)[0] = points[i][0];
    (tin.pointlist + i * 3)[1] = points[i][1];
    (tin.pointlist + i * 3)[2] = points[i][2];
  }

  tetgenbehavior tetb;
  tetb.voroout = 1;
  // tetb.quiet = 0;
  // tetb.verbose = 3;

  std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

  tetrahedralize(&tetb, &tin, &tout);

  std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

  std::cout << double(std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count()) * 1e-6 << "s." << std::endl;

  // points
  vd.vertexPositions.resize(3, tout.numberofvpoints);
  for (int vi = 0; vi < tout.numberofvpoints; vi++) {
    vd.vertexPositions.col(vi) = ES::Mp<ES::V3d>(tout.vpointlist + vi * 3);
  }

  // edges
  vd.edgeDirs.resize(3, tout.numberofvedges);
  vd.edges.resize(tout.numberofvedges);
  for (int i = 0; i < tout.numberofvedges; i++) {
    vd.edges[i] = std::make_tuple(tout.vedgelist[i].v1, tout.vedgelist[i].v2);
    if (tout.vedgelist[i].v2 < 0) {
      vd.edgeDirs.col(i) = ES::Mp<ES::V3d>(tout.vedgelist[i].vnormal);
    }
    else {
      vd.edgeDirs.col(i).setZero();
    }
  }

  // facets
  vd.facets.resize(tout.numberofvfacets);
  vd.facetCells.resize(tout.numberofvfacets);
  for (int i = 0; i < tout.numberofvfacets; i++) {
    vd.facets[i].resize(tout.vfacetlist[i].elist[0]);
    for (int j = 0; j < tout.vfacetlist[i].elist[0]; j++) {
      vd.facets[i][j] = tout.vfacetlist[i].elist[j + 1];
    }

    vd.facetCells[i] = std::make_tuple(tout.vfacetlist[i].c1, tout.vfacetlist[i].c2);
  }

  // cells
  vd.cells.resize(tout.numberofvcells);
  // std::cout << tout.numberofvcells << std::endl;
  for (int i = 0; i < tout.numberofvcells; i++) {
    // std::cout << "c " << i << '/' << tout.numberofvcells << std::flush;
    // std::cout << " " << (i < tout.numberofvcells) << ":" << std::flush;
    // std::cout << tout.vcelllist[i][0] << std::endl;

    vd.cells[i].resize(tout.vcelllist[i][0]);
    for (int j = 0; j < tout.vcelllist[i][0]; j++) {
      // std::cout << "  sub " << j << ": " << tout.vcelllist[i][j + 1] << ' ' << i << std::endl;
      vd.cells[i][j] = tout.vcelllist[i][j + 1];
    }
    // std::cout << "done." << std::endl;
  }

  return 0;
}

int pgo::TetgenInterface::computeTetMesh(const Mesh::TriMeshGeo &mesh, const std::string &switcher, EigenSupport::MXd &vtx, EigenSupport::MXi &tet)
{
  // assign vertices;
  tetgenio tin, tout, taddin;

  tin.numberofpoints = mesh.numVertices();
  tin.pointlist = new double[mesh.numVertices() * 3];
  for (int i = 0; i < mesh.numVertices(); i++) {
    const Vec3d &p = mesh.pos(i);
    tin.pointlist[i * 3] = p[0];
    tin.pointlist[i * 3 + 1] = p[1];
    tin.pointlist[i * 3 + 2] = p[2];
  }

  tin.numberoffacets = mesh.numTriangles();
  tin.facetlist = new tetgenio::facet[mesh.numTriangles()];

  for (int f = 0; f < mesh.numTriangles(); f++) {
    tetgenio::facet *facet = tin.facetlist + f;
    facet->numberofpolygons = 1;
    facet->polygonlist = new tetgenio::polygon[facet->numberofpolygons];
    facet->numberofholes = 0;
    facet->holelist = nullptr;

    tetgenio::polygon *poly = facet->polygonlist;
    poly->numberofvertices = 3;
    poly->vertexlist = new int[3];

    const auto &tri = mesh.tri(f);
    poly->vertexlist[0] = tri[0];
    poly->vertexlist[1] = tri[1];
    poly->vertexlist[2] = tri[2];
  }

  std::vector<char> sw;
  sw.resize(switcher.length() + 10);
  std::strcpy(sw.data(), switcher.c_str());

  try {
    std::cout << "Tetrahedralizing surface mesh...\n";
    tetrahedralize(sw.data(), &tin, &tout, nullptr);
  }
  catch (...) {
    std::cerr << "Error happened. \n";
    return 1;
  }

  if (tout.numberoftetrahedra > 0) {
    tet.resize(tout.numberoftetrahedra, 4);
    for (int i = 0; i < tout.numberoftetrahedra; i++) {
      for (int j = 0; j < 4; j++) {
        tet(i, j) = tout.tetrahedronlist[i * 4 + j];
      }
    }

    vtx.resize(tout.numberofpoints, 3);
    for (int i = 0; i < tout.numberofpoints; i++) {
      vtx(i, 0) = tout.pointlist[i * 3];
      vtx(i, 1) = tout.pointlist[i * 3 + 1];
      vtx(i, 2) = tout.pointlist[i * 3 + 2];
    }

    return 0;
  }
  else {
    std::cerr << "Error: cannot tetrahedralize the surface mesh.\n";
    return 1;
  }
}