/*
author: Bohan Wang
copyright to MIT, USC
*/

#include "cgalInterface.h"
#include "cgalBasic.h"
#include "cgalTemplateUtilities.h"

#include "triMeshGeo.h"
#include "simpleSphere.h"
#include "boundingVolumeTree.h"
#include "triMeshPseudoNormal.h"
#include "initPredicates.h"
#include "EigenSupport.h"
#include "pgoLogging.h"

#if defined(_MSC_VER)
#if defined(ERROR)
#  undef ERROR
#endif
#endif

#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/Polygon_mesh_processing/angle_and_area_smoothing.h>
#include <CGAL/Polygon_mesh_processing/smooth_shape.h>
#include <CGAL/Polygon_mesh_processing/detect_features.h>
#include <CGAL/Polygon_mesh_processing/locate.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/extrude.h>
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/clip.h>

#include <CGAL/Subdivision_method_3/subdivision_hosts_3.h>

#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/IO/facets_in_complex_2_to_triangle_mesh.h>

#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>

#include <CGAL/mesh_segmentation.h>
#include <CGAL/property_map.h>

#include <CGAL/convex_hull_3.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Simple_cartesian.h>

#include <unordered_set>

namespace pgo::CGALInterface
{
template<typename SurfaceMesh>
struct User_input_is_constrained_edge_map_for_SurfaceMesh
{
  typedef bool value_type;
  typedef value_type reference;
  typedef boost::readable_property_map_tag category;
  typedef typename boost::graph_traits<SurfaceMesh>::edge_descriptor key_type;

  const SurfaceMesh *surfaceMesh;
  const std::vector<typename SurfaceMesh::Edge_index> *userSelectedEdges;

  User_input_is_constrained_edge_map_for_SurfaceMesh():
    surfaceMesh(nullptr), userSelectedEdges(nullptr)
  {
  }

  User_input_is_constrained_edge_map_for_SurfaceMesh(const SurfaceMesh &sm, const std::vector<typename SurfaceMesh::Edge_index> &usedge):
    surfaceMesh(&sm), userSelectedEdges(&usedge)
  {
  }
};
}  // namespace pgo::CGALInterface

pgo::Mesh::TriMeshGeo pgo::CGALInterface::isotropicRemeshing(const Mesh::TriMeshGeo &mesh, double targetEdgeLength, int numIter, double angleInDegree,
  std::vector<int> *constrainedVertices, std::vector<std::pair<int, int>> *constrainedEdges, std::vector<int> *faceSubset)
{
  namespace ES = pgo::EigenSupport;
  using SM = CGAL::Surface_mesh<KernelInexact::Point_3>;
  using EIFMap = boost::property_map<SM, CGAL::edge_is_feature_t>::type;

  SM surfaceMesh;
  std::vector<SM::Vertex_index> vertexIndices;
  int inc = 0;
  std::map<SM::Vertex_index, int> idxMapping;
  for (const auto &pos : mesh.positions()) {
    SM::Vertex_index vi = surfaceMesh.add_vertex(SM::Point(pos[0], pos[1], pos[2]));
    vertexIndices.push_back(vi);
    idxMapping.emplace(vi, inc);
    inc++;
  }

  std::vector<SM::Face_index> smFaceSubset;
  for (int ti = 0; const auto &tri : mesh.triangles()) {
    SM::Face_index fi = surfaceMesh.add_face(vertexIndices[tri[0]], vertexIndices[tri[1]], vertexIndices[tri[2]]);
    if (faceSubset && faceSubset->size()) {
      if (std::binary_search(faceSubset->begin(), faceSubset->end(), ti)) {
        smFaceSubset.emplace_back(fi);
      }
    }

    ++ti;
  }

  // constrain angles
  EIFMap eifMap = CGAL::get(CGAL::edge_is_feature, surfaceMesh);
  CGAL::Polygon_mesh_processing::detect_sharp_edges(surfaceMesh, angleInDegree, eifMap);

  std::vector<SM::Edge_index> hardEdges;
  for (const auto &eidx : surfaceMesh.edges()) {
    if (boost::get(eifMap, eidx)) {
      hardEdges.emplace_back(eidx);
    }
  }

  // if there is no user defined sharp features
  // split long edges
  if (constrainedEdges == nullptr && constrainedVertices == nullptr && smFaceSubset.size() == 0ull) {
    CGAL::Polygon_mesh_processing::split_long_edges(hardEdges, targetEdgeLength, surfaceMesh);
  }

  // redo checking
  CGAL::Polygon_mesh_processing::detect_sharp_edges(surfaceMesh, angleInDegree, eifMap);

  // may not work if split edge happens
  if (constrainedEdges && constrainedEdges->size()) {
    std::unordered_set<std::pair<int, int>, ES::IntPairHash, ES::IntPairEqual> edgeSet;
    for (const auto &eidx : *constrainedEdges) {
      edgeSet.emplace(eidx);
    }

    int counter = 0;
    for (const auto &eidx : surfaceMesh.edges()) {
      SM::Vertex_index v0 = surfaceMesh.vertex(eidx, 0);
      SM::Vertex_index v1 = surfaceMesh.vertex(eidx, 1);

      int i0, i1;
      auto iter = idxMapping.find(v0);
      if (iter != idxMapping.end())
        i0 = iter->second;
      else
        i0 = -1;

      iter = idxMapping.find(v1);
      if (iter != idxMapping.end())
        i1 = iter->second;
      else
        i1 = -1;

      if (i0 >= 0 && i1 >= 0) {
        auto iter1 = edgeSet.find(std::make_pair(i0, i1));
        if (iter1 != edgeSet.end()) {
          boost::put(eifMap, eidx, true);
          counter++;
        }

        iter1 = edgeSet.find(std::make_pair(i1, i0));
        if (iter1 != edgeSet.end()) {
          boost::put(eifMap, eidx, true);
          counter++;
        }
      }
    }
    std::cout << "##: " << counter << "," << edgeSet.size() << std::endl;
  }

  if (smFaceSubset.size()) {
    CGAL::Polygon_mesh_processing::isotropic_remeshing(smFaceSubset, targetEdgeLength, surfaceMesh,
      CGAL::Polygon_mesh_processing::parameters::number_of_iterations(numIter).protect_constraints(true).edge_is_constrained_map(eifMap));
  }
  else {
    CGAL::Polygon_mesh_processing::isotropic_remeshing(CGAL::faces(surfaceMesh), targetEdgeLength, surfaceMesh,
      CGAL::Polygon_mesh_processing::parameters::number_of_iterations(numIter).protect_constraints(true).edge_is_constrained_map(eifMap));
  }

  Mesh::TriMeshGeo meshOut;
  surfaceMesh2TriangleMesh<KernelInexact>(surfaceMesh, meshOut);

  return meshOut;
}

pgo::Mesh::TriMeshGeo pgo::CGALInterface::refineMesh(const Mesh::TriMeshGeo &mesh, double density)
{
  using SM = CGAL::Surface_mesh<KernelInexact::Point_3>;
  using vertex_descriptor = boost::graph_traits<SM>::vertex_descriptor;

  SM surfaceMesh;
  triangleMesh2SurfaceMesh<SM, SM::Property_map<vertex_descriptor, int>>(mesh, surfaceMesh, nullptr);

  std::vector<SM::Face_index> new_facets;
  std::vector<SM::Vertex_index> new_vertices;

  CGAL::Polygon_mesh_processing::refine(surfaceMesh, CGAL::faces(surfaceMesh),
    std::back_inserter(new_facets), std::back_inserter(new_vertices),
    CGAL::parameters::density_control_factor(density));

  Mesh::TriMeshGeo meshOut;
  surfaceMesh2TriangleMesh<KernelInexact>(surfaceMesh, meshOut);

  return meshOut;
}

pgo::Mesh::TriMeshGeo pgo::CGALInterface::refineSharpRegionOnMesh(const Mesh::TriMeshGeo &mesh, double density, double angleThreshold)
{
  using K = KernelInexact;
  using SM = CGAL::Surface_mesh<K::Point_3>;
  using vertex_descriptor = boost::graph_traits<SM>::vertex_descriptor;

  SM surfaceMesh;
  triangleMesh2SurfaceMesh<SM, SM::Property_map<vertex_descriptor, int>>(mesh, surfaceMesh, nullptr);

  std::vector<SM::Face_index> selectedFaces;
  for (SM::Face_iterator fi = surfaceMesh.faces_begin(); fi != surfaceMesh.faces_end(); ++fi) {
    SM::Halfedge_index h = surfaceMesh.halfedge(*fi);
    const SM::Point &p0 = surfaceMesh.point(surfaceMesh.target(h));
    const SM::Point &p1 = surfaceMesh.point(surfaceMesh.target(surfaceMesh.next(h)));
    const SM::Point &p2 = surfaceMesh.point(surfaceMesh.target(surfaceMesh.next(surfaceMesh.next(h))));

    double angle[3];

    K::Vector_3 e1 = p1 - p0;
    K::Vector_3 e2 = p2 - p0;
    K::FT cosAngle = CGAL::scalar_product(e1, e2) / CGAL::sqrt(e1.squared_length() * e2.squared_length());
    angle[0] = std::acos(std::clamp(CGAL::to_double(cosAngle), -1.0, 1.0));

    e1 = p0 - p1;
    e2 = p2 - p1;
    cosAngle = CGAL::scalar_product(e1, e2) / CGAL::sqrt(e1.squared_length() * e2.squared_length());
    angle[1] = std::acos(std::clamp(CGAL::to_double(cosAngle), -1.0, 1.0));

    e1 = p0 - p2;
    e2 = p1 - p2;
    cosAngle = CGAL::scalar_product(e1, e2) / CGAL::sqrt(e1.squared_length() * e2.squared_length());
    angle[2] = std::acos(std::clamp(CGAL::to_double(cosAngle), -1.0, 1.0));

    double minAngle = std::min({ angle[0], angle[1], angle[2] });
    if (minAngle < angleThreshold / 180.0 * M_PI) {
      selectedFaces.emplace_back(*fi);
    }
  }

  std::vector<SM::Face_index> new_facets;
  std::vector<SM::Vertex_index> new_vertices;

  CGAL::Polygon_mesh_processing::refine(surfaceMesh, selectedFaces,
    std::back_inserter(new_facets), std::back_inserter(new_vertices),
    CGAL::parameters::density_control_factor(density));

  Mesh::TriMeshGeo meshOut;
  surfaceMesh2TriangleMesh<KernelInexact>(surfaceMesh, meshOut);

  return meshOut;
}

pgo::Mesh::TriMeshGeo pgo::CGALInterface::smoothMesh(const Mesh::TriMeshGeo &mesh, int numIter, double angleInDegree,
  std::vector<int> *constrainedVertices, std::vector<std::pair<int, int>> *constrainedEdges)
{
  using SM = CGAL::Surface_mesh<KernelInexact::Point_3>;
  using vertex_descriptor = boost::graph_traits<SM>::vertex_descriptor;
  using edge_descriptor = boost::graph_traits<SM>::edge_descriptor;

  SM surfaceMesh;
  triangleMesh2SurfaceMesh<SM, SM::Property_map<vertex_descriptor, int>>(mesh, surfaceMesh, nullptr);

  // static_assert(std::is_same<Kernel, KernelSimple>::value);
  // Constrain edges with a dihedral angle over 60
  using EIFMap = boost::property_map<SM, CGAL::edge_is_feature_t>::type;
  using VIFMap = boost::property_map<SM, bool>::type;

  namespace PMP = CGAL::Polygon_mesh_processing;

  EIFMap eif = CGAL::get(CGAL::edge_is_feature, surfaceMesh);
  PMP::detect_sharp_edges(surfaceMesh, angleInDegree, eif);

  if (constrainedEdges) {
    for (SM::Edge_iterator eitt = surfaceMesh.edges_begin(); eitt != surfaceMesh.edges_end(); eitt++) {
      SM::Vertex_index vi0 = surfaceMesh.vertex(*eitt, 0);
      SM::Vertex_index vi1 = surfaceMesh.vertex(*eitt, 1);

      bool f1 = std::binary_search(constrainedEdges->begin(), constrainedEdges->end(), std::pair<int, int>(vi0, vi1));
      bool f2 = std::binary_search(constrainedEdges->begin(), constrainedEdges->end(), std::pair<int, int>(vi1, vi0));

      if (f1 || f2) {
        eif[*eitt] = 1;
      }
    }
  }

  auto vif = surfaceMesh.add_property_map<vertex_descriptor, bool>("v:is_feature", false).first;
  for (const auto &v : surfaceMesh.vertices()) {
    boost::put(vif, v, false);
  }

  if (constrainedVertices) {
    for (SM::Vertex_iterator it = surfaceMesh.vertices_begin(); it != surfaceMesh.vertices_end(); it++) {
      SM::Vertex_index vi0 = *it;

      bool f1 = std::binary_search(constrainedVertices->begin(), constrainedVertices->end(), int(vi0));
      if (f1) {
        vif[*it] = true;
      }
    }
  }

  int sharp_counter = 0;
  for (edge_descriptor e : CGAL::edges(surfaceMesh))
    if (boost::get(eif, e))
      ++sharp_counter;

  std::cout << sharp_counter << " sharp edges" << std::endl;
  std::cout << "Smoothing mesh... (" << numIter << " iterations)" << std::endl;

  // Smooth with both angle and area criteria + Delaunay flips
  PMP::angle_and_area_smoothing(surfaceMesh, PMP::parameters::number_of_iterations(numIter)  // # iterations
                                               .use_safety_constraints(true)                 // authorize all moves
                                               .edge_is_constrained_map(eif));

  Mesh::TriMeshGeo meshOut;
  surfaceMesh2TriangleMesh<KernelInexact>(surfaceMesh, meshOut);

  return meshOut;
}

pgo::Mesh::TriMeshGeo pgo::CGALInterface::smoothShape(const Mesh::TriMeshGeo &mesh, double time, int numIter)
{
  using SM = CGAL::Surface_mesh<KernelInexact::Point_3>;
  using vertex_descriptor = boost::graph_traits<SM>::vertex_descriptor;
  using edge_descriptor = boost::graph_traits<SM>::edge_descriptor;
  namespace PMP = CGAL::Polygon_mesh_processing;

  SM surfaceMesh;
  triangleMesh2SurfaceMesh<SM, SM::Property_map<vertex_descriptor, int>>(mesh, surfaceMesh, nullptr);

  // Smooth with both angle and area criteria + Delaunay flips
  PMP::smooth_shape(surfaceMesh, time, PMP::parameters::number_of_iterations(numIter));

  Mesh::TriMeshGeo meshOut;
  surfaceMesh2TriangleMesh<KernelInexact>(surfaceMesh, meshOut);

  return meshOut;
}

bool pgo::CGALInterface::corefineAndComputeUnion(const Mesh::TriMeshGeo &mesh1, const Mesh::TriMeshGeo &mesh2, Mesh::TriMeshGeo &unionMesh)
{
  using K = KernelExact;
  using SM = CGAL::Surface_mesh<K::Point_3>;
  using vertex_descriptor = boost::graph_traits<SM>::vertex_descriptor;
  namespace PMP = CGAL::Polygon_mesh_processing;
  SM surfMesh1;
  triangleMesh2SurfaceMesh<SM, SM::Property_map<vertex_descriptor, int>>(mesh1, surfMesh1, nullptr);
  SM surfMesh2;
  triangleMesh2SurfaceMesh<SM, SM::Property_map<vertex_descriptor, int>>(mesh2, surfMesh2, nullptr);
  SM outSurfMesh;
  bool valid_union = PMP::corefine_and_compute_union(surfMesh1, surfMesh2, outSurfMesh);
  std::cout << "CGAL valid_union: " << valid_union << std::endl;
  if (valid_union) {
    surfaceMesh2TriangleMesh<K>(outSurfMesh, unionMesh);
  }
  return valid_union;
}

pgo::Mesh::TriMeshGeo pgo::CGALInterface::triangulateHolePolyline(const std::vector<Vec3d> &polyline)
{
  using K = KernelInexact;
  std::vector<K::Point_3> points(polyline.size());
  for (int i = 0; i < (int)polyline.size(); i++) {
    points[i] = K::Point_3(polyline[i][0], polyline[i][1], polyline[i][2]);
  }

  std::vector<CGAL::Triple<int, int, int>> patch_facets;
  patch_facets.reserve(points.size() - 2);
  CGAL::Polygon_mesh_processing::triangulate_hole_polyline(points, std::back_inserter(patch_facets));

  Mesh::TriMeshGeo meshOut;
  for (int i = 0; i < (int)points.size(); i++) {
    meshOut.addPos(Vec3d(points[i][0], points[i][1], points[i][2]));
  }
  for (int i = 0; i < (int)patch_facets.size(); i++) {
    meshOut.addTri(Vec3i(patch_facets[i].first, patch_facets[i].second, patch_facets[i].third));
  }

  return meshOut;
}

pgo::Mesh::TriMeshGeo pgo::CGALInterface::triangulateRefineFairHole(const Mesh::TriMeshGeo &mesh)
{
  // using K = KernelInexact;
  using K = KernelExact;
  using SM = CGAL::Surface_mesh<K::Point_3>;
  using vertex_descriptor = boost::graph_traits<SM>::vertex_descriptor;
  using halfedge_descriptor = boost::graph_traits<SM>::halfedge_descriptor;
  using face_descriptor = boost::graph_traits<SM>::face_descriptor;
  namespace PMP = CGAL::Polygon_mesh_processing;

  SM surfMesh;
  triangleMesh2SurfaceMesh<SM, SM::Property_map<vertex_descriptor, int>>(mesh, surfMesh, nullptr);
  std::vector<halfedge_descriptor> border_cycles;
  PMP::extract_boundary_cycles(surfMesh, std::back_inserter(border_cycles));

  std::cout << "Number of border cycles: " << border_cycles.size() << std::endl;

  for (halfedge_descriptor h : border_cycles) {
    std::vector<face_descriptor> patch_facets;
    std::vector<vertex_descriptor> patch_vertices;
    bool success = std::get<0>(PMP::triangulate_refine_and_fair_hole(surfMesh, h, std::back_inserter(patch_facets),
      std::back_inserter(patch_vertices)));
    std::cout << "CGAL success: " << success << std::endl;
  }
  Mesh::TriMeshGeo meshOut;
  surfaceMesh2TriangleMesh<K>(surfMesh, meshOut);
  return meshOut;
}

pgo::Mesh::TriMeshGeo pgo::CGALInterface::clipMesh(const Mesh::TriMeshGeo &meshIn, const Vec3d &planeN, const Vec3d &planeP)
{
  // using K = KernelInexact;
  using K = KernelExact;
  using SM = CGAL::Surface_mesh<K::Point_3>;
  using vertex_descriptor = boost::graph_traits<SM>::vertex_descriptor;
  namespace PMP = CGAL::Polygon_mesh_processing;

  SM surfMesh;
  triangleMesh2SurfaceMesh<SM, SM::Property_map<vertex_descriptor, int>>(meshIn, surfMesh, nullptr);

  K::Plane_3 plane(K::Point_3(planeP[0], planeP[1], planeP[2]), K::Vector_3(planeN[0], planeN[1], planeN[2]));
  if (PMP::clip(surfMesh, plane) != true)
    throw std::domain_error("cannot clip");

  Mesh::TriMeshGeo meshOut;
  surfaceMesh2TriangleMesh<K>(surfMesh, meshOut);

  return meshOut;
}

pgo::Mesh::TriMeshGeo pgo::CGALInterface::clipMesh(const Mesh::TriMeshGeo &meshIn, const Mesh::TriMeshGeo &volMesh)
{
  // using K = KernelInexact;
  using K = KernelExact;
  using Poly = Polyhedron<K>;
  namespace PMP = CGAL::Polygon_mesh_processing;

  Poly m1, m2;
  triangleMesh2Polyhedron(meshIn, m1);
  triangleMesh2Polyhedron(volMesh, m2);

  if (PMP::clip(m1, m2) != true)
    throw std::domain_error("cannot clip");

  Mesh::TriMeshGeo meshOut;
  polyhedron2TriangleMesh(m1, meshOut);

  return meshOut;
}

void pgo::CGALInterface::corefineOnly(Mesh::TriMeshGeo &mesh1, Mesh::TriMeshGeo &mesh2, bool noModify1, bool noModify2, double edgeLengthThreshold)
{
  using K = KernelExact;
  using Poly = Polyhedron<K>;

  namespace PMP = CGAL::Polygon_mesh_processing;
  namespace SMS = CGAL::Surface_mesh_simplification;

  Poly P1, P2;
  triangleMesh2Polyhedron(mesh1, P1);
  triangleMesh2Polyhedron(mesh2, P2);

  PMP::corefine(P1, P2, CGAL::parameters::do_not_modify(noModify1), CGAL::parameters::do_not_modify(noModify2));

  if (noModify1 == false) {
    removeSmallPieces(P1, edgeLengthThreshold);
  }

  if (noModify2 == false) {
    removeSmallPieces(P2, edgeLengthThreshold);
  }

  mesh1.clear();
  polyhedron2TriangleMesh(P1, mesh1);

  mesh2.clear();
  polyhedron2TriangleMesh(P2, mesh2);
}

void pgo::CGALInterface::convexHullMesh(const std::vector<Mesh::TriMeshGeo> &meshes, Mesh::TriMeshGeo &meshOut)
{
  namespace ES = pgo::EigenSupport;
  using SM = CGAL::Surface_mesh<KernelInexact::Point_3>;

  std::vector<KernelInexact::Point_3> points;
  for (const auto &mesh : meshes) {
    for (const auto &v : mesh.positions()) {
      points.emplace_back(v[0], v[1], v[2]);
    }
  }

  SM sm;
  CGAL::convex_hull_3(points.begin(), points.end(), sm);
  surfaceMesh2TriangleMesh<KernelInexact>(sm, meshOut);
}

template<typename Poly>
class SimpleMask3
{
  typedef typename boost::graph_traits<Poly>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<Poly>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::property_map<Poly, CGAL::vertex_point_t>::type Vertex_pmap;
  typedef typename boost::property_traits<Vertex_pmap>::value_type Point;
  typedef typename boost::property_traits<Vertex_pmap>::reference Point_ref;
  Poly &pmesh;
  Vertex_pmap vpm;

public:
  SimpleMask3(Poly &pmesh):
    pmesh(pmesh), vpm(get(CGAL::vertex_point, pmesh))
  {
  }

  void edge_node(halfedge_descriptor hd, Point &pt)
  {
    Point_ref p1 = boost::get(vpm, CGAL::target(hd, pmesh));
    Point_ref p2 = boost::get(vpm, CGAL::target(CGAL::opposite(hd, pmesh), pmesh));
    pt = CGAL::midpoint(p1, p2);
  }

  void vertex_node(vertex_descriptor vd, Point &pt)
  {
    pt = boost::get(vpm, vd);
  }

  void border_node(halfedge_descriptor hd, Point &ept, Point &vpt)
  {
    Point_ref p1 = boost::get(vpm, CGAL::target(hd, pmesh));
    Point_ref p2 = boost::get(vpm, CGAL::target(CGAL::opposite(hd, pmesh), pmesh));

    ept = CGAL::midpoint(p1, p2);
    vpt = p1;
  }
};

pgo::Mesh::TriMeshGeo pgo::CGALInterface::subdivideMesh(const Mesh::TriMeshGeo &meshIn, int nIter, double smallSize)
{
  Polyhedron<KernelExact> P;
  triangleMesh2Polyhedron(meshIn, P);

  SimpleMask3<Polyhedron<KernelExact>> mask(P);

  CGAL::Subdivision_method_3::PTQ(P, mask, CGAL::parameters::number_of_iterations(nIter));

  if (smallSize > 0) {
    removeSmallPieces(P, smallSize);
  }

  Mesh::TriMeshGeo meshOut;
  polyhedron2TriangleMesh(P, meshOut);

  return meshOut;
}

pgo::Mesh::TriMeshGeo pgo::CGALInterface::triangulate(const std::vector<Vec3d> &vtx, const std::vector<std::vector<int>> &faces)
{
  Polyhedron<KernelInexact> P;
  PolyhedronBuilderNonTriangle<KernelInexact, std::vector<Vec3d>::const_iterator, std::vector<std::vector<int>>::const_iterator> builder(vtx.begin(), faces.begin(), (int)vtx.size(), (int)faces.size());

  P.delegate(builder);

  CGAL::Polygon_mesh_processing::triangulate_faces(P);

  Mesh::TriMeshGeo meshOut;
  polyhedron2TriangleMesh(P, meshOut);

  return meshOut;
}

bool pgo::CGALInterface::isSelfIntersected(const Mesh::TriMeshGeo &meshIn)
{
  Polyhedron<KernelInexact> P;
  triangleMesh2Polyhedron(meshIn, P);

  bool intersecting = CGAL::Polygon_mesh_processing::does_self_intersect<CGAL::Parallel_if_available_tag>(P, CGAL::parameters::vertex_point_map(CGAL::get(CGAL::vertex_point, P)));

  return intersecting;
}

void pgo::CGALInterface::getLargestCC(const Mesh::TriMeshGeo &meshIn, Mesh::TriMeshGeo &meshOut)
{
  Polyhedron<KernelInexact> P;
  triangleMesh2Polyhedron(meshIn, P);

  CGAL::Polygon_mesh_processing::keep_largest_connected_components(P, 1ull);

  meshOut.clear();
  polyhedron2TriangleMesh(P, meshOut);
}

bool pgo::CGALInterface::isManifold(const Mesh::TriMeshGeo &meshIn)
{
  std::vector<std::vector<size_t>> faces;
  for (int ti = 0; ti < meshIn.numTriangles(); ti++) {
    faces.emplace_back();
    faces.back().resize(3);

    faces.back()[0] = meshIn.tri(ti)[0];
    faces.back()[1] = meshIn.tri(ti)[1];
    faces.back()[2] = meshIn.tri(ti)[2];
  }

  bool ret = CGAL::Polygon_mesh_processing::is_polygon_soup_a_polygon_mesh(faces);
  return ret;
}

void pgo::CGALInterface::segmentMesh(const Mesh::TriMeshGeo &meshIn, int nClusters, std::vector<int> &classID)
{
  using Poly = Polyhedron<KernelInexact>;
  using face_descriptor = boost::graph_traits<Poly>::face_descriptor;
  // create a property-map for SDF values
  using Facet_double_map = std::map<face_descriptor, double>;

  Polyhedron<KernelInexact> P;
  triangleMesh2Polyhedron(meshIn, P);

  Facet_double_map internal_sdf_map;
  boost::associative_property_map<Facet_double_map> sdf_property_map(internal_sdf_map);

  // compute SDF values using default parameters for number of rays, and cone angle
  CGAL::sdf_values(P, sdf_property_map);

  // create a property-map for segment-ids
  using Facet_int_map = std::map<face_descriptor, std::size_t>;

  Facet_int_map internal_segment_map;
  boost::associative_property_map<Facet_int_map> segment_property_map(internal_segment_map);

  // segment the mesh using default parameters for number of levels, and smoothing lambda
  // Any other scalar values can be used instead of using SDF values computed using the CGAL function
  std::size_t number_of_segments = CGAL::segmentation_from_sdf_values(P, sdf_property_map, segment_property_map);
  std::cout << "Number of segments: " << number_of_segments << std::endl;

  // print segment-ids
  // for (face_descriptor f : CGAL::faces(P)) {
  //  // ids are between [0, number_of_segments -1]
  //  std::cout << segment_property_map[f] << " ";
  //}

  // std::cout << std::endl;

  std::size_t number_of_clusters = std::min((std::size_t)nClusters, number_of_segments);
  const double smoothing_lambda = 0.3;  // importance of surface features, suggested to be in-between [0,1]

  // Note that we can use the same SDF values (sdf_property_map) over and over again for segmentation.
  // This feature is relevant for segmenting the mesh several times with different parameters.
  CGAL::segmentation_from_sdf_values(P, sdf_property_map, segment_property_map, number_of_clusters, smoothing_lambda);

  // print segment-ids
  classID.resize(meshIn.numTriangles());
  for (face_descriptor f : CGAL::faces(P)) {
    int fid = f->id();
    PGO_ALOG(fid < meshIn.numTriangles());
    classID[fid] = segment_property_map[f];
  }

  std::cout << "Done." << std::endl;
}