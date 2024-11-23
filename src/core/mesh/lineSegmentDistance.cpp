#include "geometryQuery.h"

#if defined(PGO_USE_CGAL_PREDICATES)
#include <CGAL/squared_distance_3.h>
#include <CGAL/Segment_3.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

double pgo::Mesh::minimalDistance2OfTwoLineSegments(const Vec3d &p0, const Vec3d &p1, const Vec3d &q0, const Vec3d &q1)
{
  using KernelExact = CGAL::Exact_predicates_exact_constructions_kernel;
  using Seg = KernelExact::Segment_3;

  Seg s1(KernelExact::Point_3(p0[0], p0[1], p0[2]), KernelExact::Point_3(p1[0], p1[1], p1[2]));
  Seg s2(KernelExact::Point_3(q0[0], q0[1], q0[2]), KernelExact::Point_3(q1[0], q1[1], q1[2]));

  KernelExact::FT dist2 = CGAL::squared_distance(s1, s2);

  return CGAL::to_double(dist2);
}
#else
#include <iostream>

double pgo::Mesh::minimalDistance2OfTwoLineSegments(const Vec3d &p0, const Vec3d &p1, const Vec3d &q0, const Vec3d &q1)
{
  std::cerr << "unsupported routine." << std::endl;
  return 0;
}
#endif