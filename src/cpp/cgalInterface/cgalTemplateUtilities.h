/*
author: Bohan Wang
copyright to MIT, USC
*/

#pragma once

#include "cgalBasic.h"

#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_cost.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_placement.h>
#include <CGAL/Surface_mesh_simplification/Edge_collapse_visitor_base.h>

namespace pgo
{
namespace CGALInterface
{
template<typename Kernel>
struct MyVisitor : CGAL::Surface_mesh_simplification::Edge_collapse_visitor_base<Polyhedron<Kernel>>
{
  using Profile = CGAL::Surface_mesh_simplification::Edge_profile<Polyhedron<Kernel>>;
  using vertex_descriptor = boost::graph_traits<Polyhedron<Kernel>>::vertex_descriptor;
  MyVisitor() {}

  // Called during the collecting phase for each edge collected.
  void OnCollected(const Profile &, const boost::optional<typename Kernel::FT> &) {}

  // Called during the processing phase for each edge selected.
  // If cost is absent the edge won't be collapsed.
  void OnSelected(const Profile &,
    boost::optional<typename Kernel::FT> cost, std::size_t initial, std::size_t current)
  {
    std::cerr << *cost << '\n'
              << std::flush;
  }

  // Called during the processing phase for each edge being collapsed.
  // If placement is absent the edge is left uncollapsed.
  void OnCollapsing(const Profile &, boost::optional<typename Kernel::Point_3> placement)
  {
  }

  // Called for each edge which failed the so called link-condition,
  // that is, which cannot be collapsed because doing so would
  // turn the surface mesh into a non-manifold.
  void OnNonCollapsable(const Profile &)
  {
    std::cout << "Link\n"
              << std::flush;
  }

  // Called after each edge has been collapsed
  void OnCollapsed(const Profile &, vertex_descriptor)
  {
  }
};

template<typename Kernel>
void removeSmallPieces(Polyhedron<Kernel> &p, double edgeLengthThreshold)
{
  using Poly = Polyhedron<Kernel>;
  namespace SMS = CGAL::Surface_mesh_simplification;

  typename Kernel::FT eps(edgeLengthThreshold);
  typename Kernel::FT eps2 = eps * eps;
  // My_visitor vis;

  constexpr int debug = 1;

  while (1) {
    int r = SMS::edge_collapse(p,
      CGAL::Surface_mesh_simplification::Edge_length_stop_predicate<typename Kernel::FT>(edgeLengthThreshold),
      CGAL::parameters::get_cost(SMS::Edge_length_cost<Poly>())
        .get_placement(SMS::Midpoint_placement<Poly>()));

    std::cout << "#collapsed edges: " << r << std::endl;
    break;

    if constexpr (debug) {
      if (r) {
        std::cout << "#vtx" << p.size_of_vertices() << std::endl;
        std::cout << "#tri" << p.size_of_facets() << std::endl;
        std::ofstream("aa.off") << p;
      }
    }

    bool isModified = false;
    for (auto it : p.facet_handles()) {
      typename Poly::Halfedge_handle h[3] = {
        it->halfedge(),
      };
      h[1] = h[0]->next();
      h[2] = h[1]->next();

      for (int j = 0; j < 3; j++) {
        typename Kernel::Segment_3 seg(h[(j + 1) % 3]->vertex()->point(), h[(j + 2) % 3]->vertex()->point());
        typename Kernel::FT d2 = CGAL::squared_distance(h[j]->vertex()->point(), seg);
        if (d2 > eps2) {
          continue;
        }

        std::cout << "Find one bad triangle: " << CGAL::to_double(d2) << std::endl;

        // (t (p1 - p0) + p0 - x)T (p1 - p0)
        // t (p1- p0)^T(p1 - p0
        typename Kernel::Vector_3 dir = h[(j + 2) % 3]->vertex()->point() - h[(j + 1) % 3]->vertex()->point();
        typename Kernel::Vector_3 diff = h[(j + 1) % 3]->vertex()->point() - h[j]->vertex()->point();

        typename Kernel::FT a = dir.squared_length();
        typename Kernel::FT b = CGAL::scalar_product(dir, diff);
        typename Kernel::FT t = -b / a;

        if (t <= 0 || t >= 1) {
          std::cout << "t not in [0, 1]; t=" << CGAL::to_double(t) << std::endl;
          continue;
        }

        if constexpr (0) {
          typename Kernel::FT edgeLength[3] = {
            CGAL::squared_distance(h[0]->vertex()->point(), h[0]->opposite()->vertex()->point()),
            CGAL::squared_distance(h[1]->vertex()->point(), h[1]->opposite()->vertex()->point()),
            CGAL::squared_distance(h[2]->vertex()->point(), h[2]->opposite()->vertex()->point()),
          };

          if (edgeLength[0] < edgeLength[1] && edgeLength[0] < edgeLength[2]) {
            CGAL::Euler::collapse_edge(CGAL::edge(h[0], p), p);
          }
          else if (edgeLength[1] < edgeLength[0] && edgeLength[1] < edgeLength[2]) {
            CGAL::Euler::collapse_edge(CGAL::edge(h[1], p), p);
          }
          else {
            CGAL::Euler::collapse_edge(CGAL::edge(h[2], p), p);
          }

          std::cout << "After splitting: is triangle mesh: " << p.is_pure_triangle() << std::endl;

          isModified = true;
          break;
        }
        else {
          // create a point in the middle of the segment
          typename Poly::Halfedge_handle hnew = p.split_edge(h[(j + 2) % 3]);
          typename Poly::Halfedge_handle hnew_oppo = hnew->next()->opposite();

          // std::cout << hnew->vertex()->point() << std::endl;
          // std::cout << hnew->next()->vertex()->point() << std::endl;
          // std::cout << hnew->prev()->vertex()->point() << std::endl;
          // std::cout << "==\n";

          if constexpr (debug) {
            Poly pp;
            pp.make_triangle(h[0]->vertex()->point(), h[1]->vertex()->point(), h[2]->vertex()->point());
            std::ofstream("aa.off") << pp;
          }

          typename Kernel::Point_3 pnew = hnew->vertex()->point() + dir * t;
          hnew->vertex()->point() = pnew;

          // std::cout << hnew->vertex()->point() << std::endl;
          // std::cout << hnew->next()->vertex()->point() << std::endl;
          // std::cout << hnew->prev()->vertex()->point() << std::endl;
          typename Poly::Halfedge_handle h1 = p.split_facet(hnew, h[j]);

          // Poly pp1;
          // pp1.make_triangle(h1->vertex()->point(), h1->next()->vertex()->point(), h1->next()->next()->vertex()->point());
          // std::ofstream("aa.off") << pp1;

          typename Poly::Halfedge_handle g = hnew_oppo->next()->next();
          typename Poly::Halfedge_handle h2 = p.split_facet(hnew_oppo, g);

          CGAL::Euler::collapse_edge(CGAL::edge(h1, p), p);

          std::cout << "After splitting: is triangle mesh: " << p.is_pure_triangle() << std::endl;

          isModified = true;
          break;
        }
      }

      if (isModified)
        break;
    }

    if constexpr (debug) {
      if (isModified) {
        std::cout << "#vtx" << p.size_of_vertices() << std::endl;
        std::cout << "#tri" << p.size_of_facets() << std::endl;

        std::ofstream("aa.off") << p;

        // std::cout << p.is_pure_triangle() << std::endl;
        // double minArea = 1e100;

        // for (auto it : p.facet_handles()) {
        //   Poly::Halfedge_handle h[3] = {
        //     it->halfedge(),
        //   };
        //   h[1] = h[0]->next();
        //   h[2] = h[1]->next();

        //  Kernel::Triangle_3 tri(h[0]->vertex()->point(), h[1]->vertex()->point(), h[2]->vertex()->point());
        //  minArea = std::min(minArea, std::sqrt(CGAL::to_double(tri.squared_area())));
        //}
        // std::cout << minArea;
        // }
      }
    }

    if (isModified == false)
      break;
  }
}

}  // namespace CGALUtilities
}  // namespace pgo
