/*
author: Bohan Wang
copyright to MIT
*/

#pragma once

#if defined(PGO_HAS_CERES)
#  define CGAL_PMP_USE_CERES_SOLVER
#endif

#include "triMeshGeo.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>

namespace pgo
{
namespace CGALInterface
{
using KernelSimple = CGAL::Simple_cartesian<double>;
using KernelInexact = CGAL::Exact_predicates_inexact_constructions_kernel;
using KernelExact = CGAL::Exact_predicates_exact_constructions_kernel;
using KernelExactWithSqrt = CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt;
using Triangle_int = CGAL::Triple<int, int, int>;

template<typename Kernel>
using Polyhedron = CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3>;

template<typename Kernel>
class PolyhedronBuilder : public CGAL::Modifier_base<typename Polyhedron<Kernel>::HalfedgeDS>
{
protected:
  const Mesh::TriMeshGeo &meshIn;

public:
  PolyhedronBuilder(const Mesh::TriMeshGeo &mesh):
    meshIn(mesh)
  {
  }

  using HalfedgeDS = typename Polyhedron<Kernel>::HalfedgeDS;
  void operator()(HalfedgeDS &hds)
  {
    typedef typename HalfedgeDS::Vertex Vertex;
    typedef typename Vertex::Point Point;

    // create a cgal incremental builder
    CGAL::Polyhedron_incremental_builder_3<HalfedgeDS> B(hds, true);
    B.begin_surface(meshIn.numVertices(), meshIn.numTriangles());

    // add the polyhedron vertices
    for (int vi = 0; vi < meshIn.numVertices(); vi++) {
      const auto &vtx = meshIn.pos(vi);
      typename HalfedgeDS::Vertex_handle vit = B.add_vertex(Point(vtx[0], vtx[1], vtx[2]));
      vit->id() = vi;
    }

    for (int fi = 0; fi < meshIn.numTriangles(); fi++) {
      const auto &f = meshIn.tri(fi);

      typename HalfedgeDS::Face_handle fit = B.begin_facet();
      fit->id() = fi;

      B.add_vertex_to_facet(f[0]);
      B.add_vertex_to_facet(f[1]);
      B.add_vertex_to_facet(f[2]);

      B.end_facet();
    }

    // finish up the surface
    B.end_surface();
  }
};

template<typename Kernel, typename VertexIterator, typename FaceIterator>
class PolyhedronBuilderNonTriangle : public CGAL::Modifier_base<typename Polyhedron<Kernel>::HalfedgeDS>
{
protected:
  VertexIterator vertexBegin;
  FaceIterator faceBegin;
  int numVertices, numFaces;

public:
  PolyhedronBuilderNonTriangle(VertexIterator vBeg, FaceIterator fBeg, int nvtx, int nface):
    vertexBegin(vBeg), faceBegin(fBeg), numVertices(nvtx), numFaces(nface)
  {
  }

  using HalfedgeDS = typename Polyhedron<Kernel>::HalfedgeDS;
  void operator()(HalfedgeDS &hds)
  {
    typedef typename HalfedgeDS::Vertex Vertex;
    typedef typename Vertex::Point Point;

    // create a cgal incremental builder
    CGAL::Polyhedron_incremental_builder_3<HalfedgeDS> B(hds, true);
    B.begin_surface(numVertices, numFaces);

    // add the polyhedron vertices
    VertexIterator vIter = vertexBegin;
    for (int vi = 0; vi < numVertices; vi++) {
      const auto &vtx = *vIter;
      typename HalfedgeDS::Vertex_handle vit = B.add_vertex(Point(vtx[0], vtx[1], vtx[2]));
      vit->id() = vi;

      ++vIter;
    }

    FaceIterator fIter = faceBegin;
    for (int fi = 0; fi < numFaces; fi++) {
      typename HalfedgeDS::Face_handle fit = B.begin_facet();
      fit->id() = fi;

      const auto &face = *fIter;
      for (int vi : face) {
        B.add_vertex_to_facet(vi);
      }

      B.end_facet();

      ++fIter;
    }

    // finish up the surface
    B.end_surface();
  }
};

template<typename Pt3>
std::array<double, 3> toDoublePt3(const Pt3 &pt)
{
  return std::array<double, 3>{ CGAL::to_double(pt.x()), CGAL::to_double(pt.y()), CGAL::to_double(pt.z()) };
}

template<typename SurfaceMesh, typename PropertyMap>
void triangleMesh2SurfaceMesh(const Mesh::TriMeshGeo &mesh, SurfaceMesh &surfaceMesh,
  PropertyMap *indexMap)
{
  std::vector<typename SurfaceMesh::Vertex_index> vertexIndices;
  int inc = 0;
  std::map<typename SurfaceMesh::Vertex_index, int> idx;
  for (const auto &pos : mesh.positions()) {
    typename SurfaceMesh::Vertex_index vi = surfaceMesh.add_vertex(typename SurfaceMesh::Point(pos[0], pos[1], pos[2]));
    vertexIndices.push_back(vi);
    idx.emplace(vi, inc);
    inc++;
  }

  for (const auto &tri : mesh.triangles()) {
    [[maybe_unused]] typename SurfaceMesh::Face_index fi = surfaceMesh.add_face(vertexIndices[tri[0]], vertexIndices[tri[1]], vertexIndices[tri[2]]);
  }

  if (indexMap) {
    *indexMap = surfaceMesh.template add_property_map<typename boost::graph_traits<SurfaceMesh>::vertex_descriptor, int>("v:idx", int(-1)).first;
    for (const auto &v : surfaceMesh.vertices()) {
      boost::put(*indexMap, v, idx[v]);
    }
  }
}

template<typename Kernel>
void surfaceMesh2TriangleMesh(const CGAL::Surface_mesh<typename Kernel::Point_3> &surfaceMesh, Mesh::TriMeshGeo &mesh)
{
  using SM = typename CGAL::Surface_mesh<typename Kernel::Point_3>;
  std::map<typename SM::Vertex_index, int> mapping;

  mesh.clear();

  int inc = 0;
  for (typename SM::Vertex_iterator vit = surfaceMesh.vertices_begin(); vit != surfaceMesh.vertices_end(); vit++) {
    const auto &pt = surfaceMesh.point(*vit);
    Vec3d ptd;
    ptd[0] = CGAL::to_double(pt.x());
    ptd[1] = CGAL::to_double(pt.y());
    ptd[2] = CGAL::to_double(pt.z());

    mesh.addPos(ptd);
    mapping.emplace(*vit, inc++);
  }

  for (typename SM::Face_iterator fit = surfaceMesh.faces_begin(); fit != surfaceMesh.faces_end(); fit++) {
    typename SM::Halfedge_index hi = surfaceMesh.halfedge(*fit);
    typename SM::Vertex_index vis = surfaceMesh.source(hi);
    typename SM::Vertex_index vit = surfaceMesh.target(hi);

    Vec3i vertexIdx;

    int idx = 0;
    while (vis != vit) {
      vertexIdx[idx++] = mapping[vit];

      hi = surfaceMesh.next(hi);
      vit = surfaceMesh.target(hi);
    }
    vertexIdx[idx] = mapping[vit];

    if (idx != 2) {
      std::cerr << "Surface mesh is not a triangle mesh!\n"
                << std::flush;
      exit(1);
    }

    mesh.addTri(vertexIdx);
  }
}

template<typename Kernel>
void triangleMesh2Polyhedron(const Mesh::TriMeshGeo &mesh, Polyhedron<Kernel> &P)
{
  // build a polyhedron from the loaded arrays
  PolyhedronBuilder<Kernel> builder(mesh);
  P.delegate(builder);
}

template<typename Kernel>
void polyhedron2TriangleMesh(const Polyhedron<Kernel> &P, Mesh::TriMeshGeo &mesh)
{
  int idx = 0;
  std::map<typename Polyhedron<Kernel>::Vertex_const_handle, int> vertexIndices;
  for (typename Polyhedron<Kernel>::Vertex_const_iterator itt = P.vertices_begin(); itt != P.vertices_end(); ++itt) {
    const typename Kernel::Point_3 &pt = itt->point();
    mesh.addPos(Vec3d{ CGAL::to_double(pt[0]),
      CGAL::to_double(pt[1]),
      CGAL::to_double(pt[2]) });

    vertexIndices.insert(make_pair(itt, idx));
    idx++;
  }

  for (typename Polyhedron<Kernel>::Face_const_iterator itt = P.facets_begin(); itt != P.facets_end(); ++itt) {
    Vec3i face;
    int inc = 0;

    typename Polyhedron<Kernel>::Halfedge_around_facet_const_circulator ctr = itt->facet_begin();
    do {
      auto idxItt = vertexIndices.find(ctr->vertex());
      if (idxItt == vertexIndices.end()) {
        std::cerr << "Cannot find the vertex" << std::endl;
        abort();
      }
      if (inc == 3) {
        std::cerr << "Not a triangle mesh" << std::endl;
        abort();
      }

      face[inc++] = idxItt->second;
      ctr++;
    } while (ctr != itt->facet_begin());
    mesh.addTri(face);
  }
}

}  // namespace CGALInterface
}  // namespace pgo