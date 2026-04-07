/*
copyright to Bohan Wang
*/

// Full implementation of Codimensional IPC collision handling.
// =============================================================================

#include "CIPC.h"
#include "CIPC_autogen.h"
#include "CIPC_autogen_ll.h"

#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>
#include <tbb/blocked_range.h>
#include <tbb/enumerable_thread_specific.h>

#include <iostream>
#include <cassert>
#include <limits>
#include <numeric>
#include <cmath>
#include <functional>
#include <fstream>
#include <atomic>
#include <cstdint>

namespace pgo
{
namespace Contact
{
namespace CIPC
{
static constexpr double kEps = 1e-20;  // numerical guard

// helper: pack an edge key for adjacency lookup
static int64_t packPair(int a, int b)
{
  return (int64_t)a * 10000000LL + (int64_t)b;
}

// =========================================================================
//  Distance type classification
// =========================================================================
namespace distance
{

PTDistType classifyPT(const V3d &p, const V3d &t0,
  const V3d &t1, const V3d &t2)
{
  // Compute barycentric coordinates of the projection of p onto the
  // triangle plane, then decide which feature is closest.
  const V3d e0 = t1 - t0;
  const V3d e1 = t2 - t0;
  const V3d v = p - t0;

  double d00 = e0.dot(e0);
  double d01 = e0.dot(e1);
  double d11 = e1.dot(e1);
  double d20 = v.dot(e0);
  double d21 = v.dot(e1);
  double denom = d00 * d11 - d01 * d01;

  // Degenerate triangle guard
  if (std::abs(denom) < kEps) {
    // fall back to closest vertex
    double d0 = (p - t0).squaredNorm();
    double d1 = (p - t1).squaredNorm();
    double d2 = (p - t2).squaredNorm();
    if (d0 <= d1 && d0 <= d2)
      return PTDistType::PP_PT0;
    if (d1 <= d2)
      return PTDistType::PP_PT1;
    return PTDistType::PP_PT2;
  }

  double s = (d11 * d20 - d01 * d21) / denom;
  double t = (d00 * d21 - d01 * d20) / denom;

  // Inside triangle
  if (s >= 0 && t >= 0 && (s + t) <= 1.0)
    return PTDistType::PT;

  // Project onto each edge and pick the closest
  // Edge t0-t1:  parameterize by u in [0,1]
  auto edgeParam = [](const V3d &pp, const V3d &a, const V3d &b) -> double {
    V3d ab = b - a;
    double len2 = ab.squaredNorm();
    if (len2 < kEps)
      return 0.0;
    return std::clamp(ab.dot(pp - a) / len2, 0.0, 1.0);
  };

  auto edgeDist2 = [&](const V3d &pp, const V3d &a, const V3d &b) -> double {
    double u = edgeParam(pp, a, b);
    return (pp - a - u * (b - a)).squaredNorm();
  };

  double dE01 = edgeDist2(p, t0, t1);
  double dE12 = edgeDist2(p, t1, t2);
  double dE20 = edgeDist2(p, t2, t0);

  double u01 = edgeParam(p, t0, t1);
  double u12 = edgeParam(p, t1, t2);
  double u20 = edgeParam(p, t2, t0);

  // find minimum among edges
  int bestEdge = 0;
  double bestD = dE01;
  if (dE12 < bestD) {
    bestD = dE12;
    bestEdge = 1;
  }
  if (dE20 < bestD) {
    bestD = dE20;
    bestEdge = 2;
  }

  if (bestEdge == 0) {
    if (u01 <= 0.0)
      return PTDistType::PP_PT0;
    else if (u01 >= 1.0)
      return PTDistType::PP_PT1;
    else
      return PTDistType::PE_PT0T1;
  }
  else if (bestEdge == 1) {
    if (u12 <= 0.0)
      return PTDistType::PP_PT1;
    else if (u12 >= 1.0)
      return PTDistType::PP_PT2;
    else
      return PTDistType::PE_PT1T2;
  }
  else {
    if (u20 <= 0.0)
      return PTDistType::PP_PT2;
    else if (u20 >= 1.0)
      return PTDistType::PP_PT0;
    else
      return PTDistType::PE_PT2T0;
  }
}

EEDistType classifyEE(const V3d &ea0, const V3d &ea1,
  const V3d &eb0, const V3d &eb1)
{
  // Based on Codim-IPC reference (DISTANCE_TYPE.h) which itself is based on
  // http://geomalgorithms.com/a07-_distance.html with a near-parallel guard.
  V3d u = ea1 - ea0;
  V3d v = eb1 - eb0;
  V3d w = ea0 - eb0;

  double a = u.squaredNorm();
  double b = u.dot(v);
  double c = v.squaredNorm();
  double d = u.dot(w);
  double e = v.dot(w);
  double D = a * c - b * b;  // always >= 0

  double tD = D;
  double sN, tN;

  int defaultCase = 8;  // EE

  sN = b * e - c * d;
  if (sN <= 0.0) {
    tN = e;
    tD = c;
    defaultCase = 4;  // PE_Ea0_Eb
  }
  else if (sN >= D) {
    tN = e + b;
    tD = c;
    defaultCase = 5;  // PE_Ea1_Eb
  }
  else {
    tN = a * e - b * d;
    // Near-parallel or coplanar edges: avoid degenerate EE distance.
    // If the interior closest point falls on the line-line segment but
    // the cross product is nearly zero, force to a PE case instead.
    if (tN > 0.0 && tN < tD &&
      (u.cross(v).dot(w) == 0.0 ||
        u.cross(v).squaredNorm() < 1.0e-20 * a * c)) {
      if (sN < D / 2) {
        tN = e;
        tD = c;
        defaultCase = 4;  // PE_Ea0_Eb
      }
      else {
        tN = e + b;
        tD = c;
        defaultCase = 5;  // PE_Ea1_Eb
      }
    }
  }

  if (tN <= 0.0) {
    if (-d <= 0.0)
      return EEDistType::PP_Ea0Eb0;
    else if (-d >= a)
      return EEDistType::PP_Ea1Eb0;
    else
      return EEDistType::PE_Eb0_Ea;
  }
  else if (tN >= tD) {
    if ((-d + b) <= 0.0)
      return EEDistType::PP_Ea0Eb1;
    else if ((-d + b) >= a)
      return EEDistType::PP_Ea1Eb1;
    else
      return EEDistType::PE_Eb1_Ea;
  }

  switch (defaultCase) {
  case 4: return EEDistType::PE_Ea0_Eb;
  case 5: return EEDistType::PE_Ea1_Eb;
  default: return EEDistType::EE;
  }
}

// =========================================================================
//  PP squared distance   x = [a; b]  (6 DOF)
// =========================================================================
double ppSqDist(const V3d &a, const V3d &b)
{
  return (a - b).squaredNorm();
}

V6d ppSqDistGrad(const V3d &a, const V3d &b)
{
  V3d d = a - b;
  V6d g;
  g.head<3>() = 2.0 * d;
  g.tail<3>() = -2.0 * d;
  return g;
}

M6d ppSqDistHess(const V3d & /*a*/, const V3d & /*b*/)
{
  M6d H;
  H.setZero();
  H.block<3, 3>(0, 0) = 2.0 * M3d::Identity();
  H.block<3, 3>(0, 3) = -2.0 * M3d::Identity();
  H.block<3, 3>(3, 0) = -2.0 * M3d::Identity();
  H.block<3, 3>(3, 3) = 2.0 * M3d::Identity();
  return H;
}

// =========================================================================
//  PE squared distance   x = [p; e0; e1]  (9 DOF)
// =========================================================================
double peSqDist(const V3d &p, const V3d &e0, const V3d &e1)
{
  V3d e = e1 - e0;
  V3d v = p - e0;
  double len2 = e.squaredNorm();
  if (len2 < kEps)
    return v.squaredNorm();
  double t = v.dot(e) / len2;
  V3d r = v - t * e;
  return r.squaredNorm();
}

// Gradient via analytical derivation:
//   s = alpha - beta^2/gamma
//   alpha = ||p-e0||^2,  beta = (p-e0).(e1-e0),  gamma = ||e1-e0||^2
V9d peSqDistGrad(const V3d &p, const V3d &e0, const V3d &e1)
{
  V3d a = p - e0;   // 3
  V3d b = e1 - e0;  // 3
  double gamma = b.squaredNorm();
  if (gamma < kEps) {
    // degenerate edge: fall back to PP(p, e0)
    V9d g;
    g.setZero();
    g.head<3>() = 2.0 * a;
    g.segment<3>(3) = -2.0 * a;
    return g;
  }
  double beta = a.dot(b);
  double t = beta / gamma;
  V3d r = a - t * b;  // residual

  // ds/dp   = 2 r
  // ds/de0  = -2(1-t) r   +  ... actually, let me re-derive cleanly.
  //
  //   s = ||a||^2 - (a.b)^2 / ||b||^2
  //
  //   ds/dp  = 2a - 2*(a.b)*b / ||b||^2  = 2*(a - t*b) = 2*r
  //   ds/de0 = -2a - d/de0 [ (a.b)^2 / ||b||^2 ]
  //
  //   d/de0 [(a.b)^2 / ||b||^2]
  //     = [2*(a.b)*(a-b)*||b||^2 - (a.b)^2*(-2*b)] / ||b||^4
  //     = 2*beta*[ (a-b)*gamma + beta*b ] / gamma^2
  //     = 2*t*[ (a-b) + t*b ]
  //     = 2*t*[ a - b + t*b ]
  //     = 2*t*[ a - (1-t)*b ]
  //
  //   ds/de0 = -2*a - 2*t*(a - (1-t)*b)
  //          = -2*a*(1+t) + 2*t*(1-t)*b
  //
  // Similarly:
  //   ds/de1 = -d/de1 [ (a.b)^2 / ||b||^2 ]
  //     d/de1: d(beta)/de1 = a,  d(gamma)/de1 = 2b
  //     = [2*beta*a*gamma - beta^2*2*b] / gamma^2
  //     = 2*t*a - 2*t^2*b
  //     = 2*t*(a - t*b)
  //     = 2*t*r
  //   ds/de1 = -2*t*r
  //
  //   But ds/de0 should be: let me verify via r:
  //   ds/de0 = -2(1+t)a + 2t(1-t)b
  //          = -2a -2ta + 2tb - 2t^2 b
  //          = -2(a + t(a - b + tb))
  //   Hmm let me just use  ds/de0 = -(ds/dp + ds/de1)  since s depends
  //   only on the relative geometry and shifting all vertices by the same
  //   amount doesn't change s.
  //
  //   ds/dp + ds/de0 + ds/de1 = 0  (translation invariance)
  //   ds/de0 = -ds/dp - ds/de1 = -2r -(-2t*r) = -2r + 2tr = -2(1-t)r

  V9d g;
  g.segment<3>(0) = 2.0 * r;
  g.segment<3>(6) = -2.0 * t * r;                        // ds/de1
  g.segment<3>(3) = -g.segment<3>(0) - g.segment<3>(6);  // translation invariance
  return g;
}

M9d peSqDistHess(const V3d &p, const V3d &e0, const V3d &e1)
{
  // Analytical Hessian via auto-generated symbolic derivatives
  double H_data[81];
  pgo::Contact::CIPC::autogen::point_line_distance_hessian_3D(
    p[0], p[1], p[2],
    e0[0], e0[1], e0[2],
    e1[0], e1[1], e1[2],
    H_data);
  // IPC toolkit stores column-major (Eigen default)
  return Eigen::Map<M9d>(H_data);
}

// =========================================================================
//  PT face squared distance   x = [p; t0; t1; t2]  (12 DOF)
// =========================================================================
//  d^2 = [(p - t0) . n]^2 / ||n||^2      n = (t1-t0) x (t2-t0)
//
double ptSqDist(const V3d &p, const V3d &t0,
  const V3d &t1, const V3d &t2)
{
  V3d n = (t1 - t0).cross(t2 - t0);
  double nn = n.squaredNorm();
  if (nn < kEps)
    return (p - t0).squaredNorm();  // degenerate
  double dn = (p - t0).dot(n);
  return dn * dn / nn;
}

V12d ptSqDistGrad(const V3d &p, const V3d &t0,
  const V3d &t1, const V3d &t2)
{
  // Analytical gradient via auto-generated symbolic derivatives
  double g_data[12];
  pgo::Contact::CIPC::autogen::point_plane_distance_gradient(
    p[0], p[1], p[2],
    t0[0], t0[1], t0[2],
    t1[0], t1[1], t1[2],
    t2[0], t2[1], t2[2],
    g_data);
  return Eigen::Map<V12d>(g_data);
}

M12d ptSqDistHess(const V3d &p, const V3d &t0,
  const V3d &t1, const V3d &t2)
{
  double H_data[144];
  pgo::Contact::CIPC::autogen::point_plane_distance_hessian(
    p[0], p[1], p[2],
    t0[0], t0[1], t0[2],
    t1[0], t1[1], t1[2],
    t2[0], t2[1], t2[2],
    H_data);
  return Eigen::Map<M12d>(H_data);
}

// =========================================================================
//  EE interior squared distance   x = [ea0; ea1; eb0; eb1]  (12 DOF)
// =========================================================================
double eeSqDist(const V3d &ea0, const V3d &ea1,
  const V3d &eb0, const V3d &eb1)
{
  V3d da = ea1 - ea0;
  V3d db = eb1 - eb0;
  V3d dg = ea0 - eb0;
  double a = da.dot(da);
  double b_val = da.dot(db);
  double c = db.dot(db);
  double d = da.dot(dg);
  double e = db.dot(dg);
  double denom = a * c - b_val * b_val;

  double s, t;
  if (std::abs(denom) < kEps) {
    s = 0.0;
    t = (c > kEps) ? e / c : 0.0;
  }
  else {
    s = (b_val * e - c * d) / denom;
    t = (a * e - b_val * d) / denom;
  }
  // Note: for interior EE type, s,t should be in (0,1).
  // We don't clamp here as the classification already ensures interior.
  V3d diff = dg + s * da - t * db;
  return diff.squaredNorm();
}

V12d eeSqDistGrad(const V3d &ea0, const V3d &ea1,
  const V3d &eb0, const V3d &eb1)
{
  double g_data[12];
  pgo::Contact::CIPC::autogen::line_line_distance_gradient(
    ea0[0], ea0[1], ea0[2],
    ea1[0], ea1[1], ea1[2],
    eb0[0], eb0[1], eb0[2],
    eb1[0], eb1[1], eb1[2],
    g_data);
  return Eigen::Map<V12d>(g_data);
}

M12d eeSqDistHess(const V3d &ea0, const V3d &ea1,
  const V3d &eb0, const V3d &eb1)
{
  double H_data[144];
  pgo::Contact::CIPC::autogen::line_line_distance_hessian(
    ea0[0], ea0[1], ea0[2],
    ea1[0], ea1[1], ea1[2],
    eb0[0], eb0[1], eb0[2],
    eb1[0], eb1[1], eb1[2],
    H_data);
  return Eigen::Map<M12d>(H_data);
}

// =========================================================================
//  EE Mollifier  (Codim-IPC §4.2  —  smoothly disables EE barrier when
//  edges are near-parallel so the PT barrier takes over)
// =========================================================================
// m(e_a, e_b) uses the squared cross-product magnitude relative to a
// rest-pose threshold eps_x.
//
//   c = ||ea x eb||^2 / (||ea||^2 * ||eb||^2)       in [0, 1]
//   m(c) = { 0                           if c < eps_x
//          { (-c/eps_x + 2) * c/eps_x    otherwise      (C1 blend)
//   For c >= eps_x  m is 1 essentially  (we use a smooth ramp)
//
// In practice the mollifier is:
//   m(x) = (-x + 2*eps_x)^2 * x / eps_x^3    for x in [0, eps_x)
//   m(x) = 1                                    for x >= eps_x
// where x = || (ea1-ea0) x (eb1-eb0) ||^2

double eeMollifier(const V3d &ea0, const V3d &ea1,
  const V3d &eb0, const V3d &eb1,
  double eps_x)
{
  V3d ea = ea1 - ea0;
  V3d eb = eb1 - eb0;
  double x = ea.cross(eb).squaredNorm();
  if (x >= eps_x)
    return 1.0;
  double r = x / eps_x;
  // m = (-r + 2) * r  =  r*(2 - r)      [C1 at r=1: m(1)=1, m'(1)=0]
  return r * (2.0 - r);
}

// Analytical mollifier gradient:
//   m(x) = r*(2-r),  r = x/eps_x,  x = ||ea x eb||^2
//   dm/dx = (2/eps_x)*(1 - x/eps_x)
//   dm/dvertices = dm/dx * dx/dvertices
V12d eeMollifierGrad(const V3d &ea0, const V3d &ea1,
  const V3d &eb0, const V3d &eb1,
  double eps_x)
{
  V3d ea = ea1 - ea0;
  V3d eb = eb1 - eb0;
  double x = ea.cross(eb).squaredNorm();
  if (x >= eps_x)
    return V12d::Zero();

  // dm/dx
  double one_div_eps = 1.0 / eps_x;
  double dmoll_dx = 2.0 * one_div_eps * (1.0 - one_div_eps * x);

  // dx/dvertices (analytical)
  double gx_data[12];
  pgo::Contact::CIPC::autogen::edge_edge_cross_squarednorm_gradient(
    ea0[0], ea0[1], ea0[2], ea1[0], ea1[1], ea1[2],
    eb0[0], eb0[1], eb0[2], eb1[0], eb1[1], eb1[2],
    gx_data);

  return dmoll_dx * Eigen::Map<V12d>(gx_data);
}

// Analytical mollifier Hessian:
//   H_m = (d2m/dx2) * gx * gx^T  +  (dm/dx) * Hx
//   d2m/dx2 = -2 / eps_x^2
M12d eeMollifierHess(const V3d &ea0, const V3d &ea1,
  const V3d &eb0, const V3d &eb1,
  double eps_x)
{
  V3d ea = ea1 - ea0;
  V3d eb = eb1 - eb0;
  double x = ea.cross(eb).squaredNorm();
  if (x >= eps_x)
    return M12d::Zero();

  double one_div_eps = 1.0 / eps_x;
  double dmoll_dx = 2.0 * one_div_eps * (1.0 - one_div_eps * x);
  double d2moll_dx2 = -2.0 * one_div_eps * one_div_eps;

  double gx_data[12];
  pgo::Contact::CIPC::autogen::edge_edge_cross_squarednorm_gradient(
    ea0[0], ea0[1], ea0[2], ea1[0], ea1[1], ea1[2],
    eb0[0], eb0[1], eb0[2], eb1[0], eb1[1], eb1[2],
    gx_data);
  V12d gx = Eigen::Map<V12d>(gx_data);

  double Hx_data[144];
  pgo::Contact::CIPC::autogen::edge_edge_cross_squarednorm_hessian(
    ea0[0], ea0[1], ea0[2], ea1[0], ea1[1], ea1[2],
    eb0[0], eb0[1], eb0[2], eb1[0], eb1[1], eb1[2],
    Hx_data);
  M12d Hx = Eigen::Map<M12d>(Hx_data);

  return dmoll_dx * Hx + d2moll_dx2 * gx * gx.transpose();
}

// =========================================================================
//  Unified dispatchers  —  classify distance type, call the right sub-
//  routine, and embed into the 12-DOF canonical ordering.
// =========================================================================

// Helper: embed a 6-DOF (PP) result into 12-DOF at given slot indices
static void embedPP(int slotA, int slotB,
  const V6d &g6, const M6d &H6,
  V12d &g12, M12d &H12)
{
  g12.setZero();
  H12.setZero();
  g12.segment<3>(3 * slotA) = g6.head<3>();
  g12.segment<3>(3 * slotB) = g6.tail<3>();
  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 2; ++j) {
      int ri = (i == 0 ? slotA : slotB);
      int ci = (j == 0 ? slotA : slotB);
      H12.block<3, 3>(3 * ri, 3 * ci) = H6.block<3, 3>(3 * i, 3 * j);
    }
}

// Helper: embed a 9-DOF (PE) result into 12-DOF
//  peSlots[3] = which of {0,1,2,3} the 3 PE vertices map to
static void embedPE(const int peSlots[3],
  const V9d &g9, const M9d &H9,
  V12d &g12, M12d &H12)
{
  g12.setZero();
  H12.setZero();
  for (int i = 0; i < 3; ++i)
    g12.segment<3>(3 * peSlots[i]) = g9.segment<3>(3 * i);

  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      H12.block<3, 3>(3 * peSlots[i], 3 * peSlots[j]) =
        H9.block<3, 3>(3 * i, 3 * j);
}

// ----- PT dispatcher -----
double computePTSqDist(const V3d &p, const V3d &t0,
  const V3d &t1, const V3d &t2)
{
  auto tp = classifyPT(p, t0, t1, t2);
  switch (tp) {
  case PTDistType::PP_PT0: return ppSqDist(p, t0);
  case PTDistType::PP_PT1: return ppSqDist(p, t1);
  case PTDistType::PP_PT2: return ppSqDist(p, t2);
  case PTDistType::PE_PT0T1: return peSqDist(p, t0, t1);
  case PTDistType::PE_PT1T2: return peSqDist(p, t1, t2);
  case PTDistType::PE_PT2T0: return peSqDist(p, t2, t0);
  case PTDistType::PT: return ptSqDist(p, t0, t1, t2);
  }
  return 0.0;
}

V12d computePTSqDistGrad(const V3d &p, const V3d &t0,
  const V3d &t1, const V3d &t2)
{
  V12d g12;
  g12.setZero();
  M12d dummy;
  dummy.setZero();
  auto tp = classifyPT(p, t0, t1, t2);
  switch (tp) {
  case PTDistType::PP_PT0: {
    V6d g6 = ppSqDistGrad(p, t0);
    M6d H6;
    H6.setZero();
    embedPP(0, 1, g6, H6, g12, dummy);
    break;
  }
  case PTDistType::PP_PT1: {
    V6d g6 = ppSqDistGrad(p, t1);
    M6d H6;
    H6.setZero();
    embedPP(0, 2, g6, H6, g12, dummy);
    break;
  }
  case PTDistType::PP_PT2: {
    V6d g6 = ppSqDistGrad(p, t2);
    M6d H6;
    H6.setZero();
    embedPP(0, 3, g6, H6, g12, dummy);
    break;
  }
  case PTDistType::PE_PT0T1: {
    V9d g9 = peSqDistGrad(p, t0, t1);
    M9d H9;
    H9.setZero();
    int slots[3] = { 0, 1, 2 };
    embedPE(slots, g9, H9, g12, dummy);
    break;
  }
  case PTDistType::PE_PT1T2: {
    V9d g9 = peSqDistGrad(p, t1, t2);
    M9d H9;
    H9.setZero();
    int slots[3] = { 0, 2, 3 };
    embedPE(slots, g9, H9, g12, dummy);
    break;
  }
  case PTDistType::PE_PT2T0: {
    V9d g9 = peSqDistGrad(p, t2, t0);
    M9d H9;
    H9.setZero();
    int slots[3] = { 0, 3, 1 };
    embedPE(slots, g9, H9, g12, dummy);
    break;
  }
  case PTDistType::PT:
    g12 = ptSqDistGrad(p, t0, t1, t2);
    break;
  }
  return g12;
}

M12d computePTSqDistHess(const V3d &p, const V3d &t0,
  const V3d &t1, const V3d &t2)
{
  V12d dummyG;
  dummyG.setZero();
  M12d H12;
  H12.setZero();
  auto tp = classifyPT(p, t0, t1, t2);
  switch (tp) {
  case PTDistType::PP_PT0: {
    V6d g6;
    g6.setZero();
    M6d H6 = ppSqDistHess(p, t0);
    embedPP(0, 1, g6, H6, dummyG, H12);
    break;
  }
  case PTDistType::PP_PT1: {
    V6d g6;
    g6.setZero();
    M6d H6 = ppSqDistHess(p, t1);
    embedPP(0, 2, g6, H6, dummyG, H12);
    break;
  }
  case PTDistType::PP_PT2: {
    V6d g6;
    g6.setZero();
    M6d H6 = ppSqDistHess(p, t2);
    embedPP(0, 3, g6, H6, dummyG, H12);
    break;
  }
  case PTDistType::PE_PT0T1: {
    V9d g9;
    g9.setZero();
    M9d H9 = peSqDistHess(p, t0, t1);
    int slots[3] = { 0, 1, 2 };
    embedPE(slots, g9, H9, dummyG, H12);
    break;
  }
  case PTDistType::PE_PT1T2: {
    V9d g9;
    g9.setZero();
    M9d H9 = peSqDistHess(p, t1, t2);
    int slots[3] = { 0, 2, 3 };
    embedPE(slots, g9, H9, dummyG, H12);
    break;
  }
  case PTDistType::PE_PT2T0: {
    V9d g9;
    g9.setZero();
    M9d H9 = peSqDistHess(p, t2, t0);
    int slots[3] = { 0, 3, 1 };
    embedPE(slots, g9, H9, dummyG, H12);
    break;
  }
  case PTDistType::PT:
    H12 = ptSqDistHess(p, t0, t1, t2);
    break;
  }
  return H12;
}

// ----- EE dispatcher -----
double computeEESqDist(const V3d &ea0, const V3d &ea1,
  const V3d &eb0, const V3d &eb1)
{
  auto tp = classifyEE(ea0, ea1, eb0, eb1);
  switch (tp) {
  case EEDistType::PP_Ea0Eb0: return ppSqDist(ea0, eb0);
  case EEDistType::PP_Ea0Eb1: return ppSqDist(ea0, eb1);
  case EEDistType::PP_Ea1Eb0: return ppSqDist(ea1, eb0);
  case EEDistType::PP_Ea1Eb1: return ppSqDist(ea1, eb1);
  case EEDistType::PE_Ea0_Eb: return peSqDist(ea0, eb0, eb1);
  case EEDistType::PE_Ea1_Eb: return peSqDist(ea1, eb0, eb1);
  case EEDistType::PE_Eb0_Ea: return peSqDist(eb0, ea0, ea1);
  case EEDistType::PE_Eb1_Ea: return peSqDist(eb1, ea0, ea1);
  case EEDistType::EE: return eeSqDist(ea0, ea1, eb0, eb1);
  }
  return 0.0;
}

V12d computeEESqDistGrad(const V3d &ea0, const V3d &ea1,
  const V3d &eb0, const V3d &eb1)
{
  V12d g12;
  g12.setZero();
  M12d dummy;
  dummy.setZero();
  auto tp = classifyEE(ea0, ea1, eb0, eb1);
  switch (tp) {
  case EEDistType::PP_Ea0Eb0: {
    embedPP(0, 2, ppSqDistGrad(ea0, eb0), M6d::Zero(), g12, dummy);
    break;
  }
  case EEDistType::PP_Ea0Eb1: {
    embedPP(0, 3, ppSqDistGrad(ea0, eb1), M6d::Zero(), g12, dummy);
    break;
  }
  case EEDistType::PP_Ea1Eb0: {
    embedPP(1, 2, ppSqDistGrad(ea1, eb0), M6d::Zero(), g12, dummy);
    break;
  }
  case EEDistType::PP_Ea1Eb1: {
    embedPP(1, 3, ppSqDistGrad(ea1, eb1), M6d::Zero(), g12, dummy);
    break;
  }
  case EEDistType::PE_Ea0_Eb: {
    // PE with [ea0, eb0, eb1] -> slots [0, 2, 3]
    int s[3] = { 0, 2, 3 };
    embedPE(s, peSqDistGrad(ea0, eb0, eb1), M9d::Zero(), g12, dummy);
    break;
  }
  case EEDistType::PE_Ea1_Eb: {
    int s[3] = { 1, 2, 3 };
    embedPE(s, peSqDistGrad(ea1, eb0, eb1), M9d::Zero(), g12, dummy);
    break;
  }
  case EEDistType::PE_Eb0_Ea: {
    // PE with [eb0, ea0, ea1] -> slots [2, 0, 1]
    int s[3] = { 2, 0, 1 };
    embedPE(s, peSqDistGrad(eb0, ea0, ea1), M9d::Zero(), g12, dummy);
    break;
  }
  case EEDistType::PE_Eb1_Ea: {
    int s[3] = { 3, 0, 1 };
    embedPE(s, peSqDistGrad(eb1, ea0, ea1), M9d::Zero(), g12, dummy);
    break;
  }
  case EEDistType::EE:
    g12 = eeSqDistGrad(ea0, ea1, eb0, eb1);
    break;
  }
  return g12;
}

M12d computeEESqDistHess(const V3d &ea0, const V3d &ea1,
  const V3d &eb0, const V3d &eb1)
{
  V12d dummyG;
  dummyG.setZero();
  M12d H12;
  H12.setZero();
  auto tp = classifyEE(ea0, ea1, eb0, eb1);
  switch (tp) {
  case EEDistType::PP_Ea0Eb0:
    embedPP(0, 2, V6d::Zero(), ppSqDistHess(ea0, eb0), dummyG, H12);
    break;
  case EEDistType::PP_Ea0Eb1:
    embedPP(0, 3, V6d::Zero(), ppSqDistHess(ea0, eb1), dummyG, H12);
    break;
  case EEDistType::PP_Ea1Eb0:
    embedPP(1, 2, V6d::Zero(), ppSqDistHess(ea1, eb0), dummyG, H12);
    break;
  case EEDistType::PP_Ea1Eb1:
    embedPP(1, 3, V6d::Zero(), ppSqDistHess(ea1, eb1), dummyG, H12);
    break;
  case EEDistType::PE_Ea0_Eb: {
    int s[3] = { 0, 2, 3 };
    embedPE(s, V9d::Zero(), peSqDistHess(ea0, eb0, eb1), dummyG, H12);
    break;
  }
  case EEDistType::PE_Ea1_Eb: {
    int s[3] = { 1, 2, 3 };
    embedPE(s, V9d::Zero(), peSqDistHess(ea1, eb0, eb1), dummyG, H12);
    break;
  }
  case EEDistType::PE_Eb0_Ea: {
    int s[3] = { 2, 0, 1 };
    embedPE(s, V9d::Zero(), peSqDistHess(eb0, ea0, ea1), dummyG, H12);
    break;
  }
  case EEDistType::PE_Eb1_Ea: {
    int s[3] = { 3, 0, 1 };
    embedPE(s, V9d::Zero(), peSqDistHess(eb1, ea0, ea1), dummyG, H12);
    break;
  }
  case EEDistType::EE:
    H12 = eeSqDistHess(ea0, ea1, eb0, eb1);
    break;
  }
  return H12;
}

}  // namespace distance

// =========================================================================
//  Barrier function
// =========================================================================
namespace barrier
{

// Barrier on squared distance (matching Codim-IPC elastic formulation):
//   b(s, shat) = -(s/shat - 1)^2 * ln(s/shat)    for 0 < s < shat
// where s = d^2 (squared distance), shat = dhat^2
double b(double s, double shat)
{
  if (s >= shat || s <= 0.0)
    return 0.0;
  double r = s / shat;
  double rm1 = r - 1.0;
  return -rm1 * rm1 * std::log(r);
}

// db/ds = -(1/shat) * (r - 1) * (2 ln(r) + (r-1)/r)
double dbds(double s, double shat)
{
  if (s >= shat || s <= 0.0)
    return 0.0;
  double r = s / shat;
  double rm1 = r - 1.0;
  return -(rm1 * (2.0 * std::log(r) + rm1 / r)) / shat;
}

// d2b/ds2 = -(1/shat^2) * (2 ln(r) + 4(r-1)/r - (r-1)^2/r^2)
double d2bds2(double s, double shat)
{
  if (s >= shat || s <= 0.0)
    return 0.0;
  double r = s / shat;
  double rm1 = r - 1.0;
  double invr = 1.0 / r;
  double shat2 = shat * shat;
  return -(2.0 * std::log(r) + 4.0 * rm1 * invr - rm1 * rm1 * invr * invr) / shat2;
}

}  // namespace barrier

// =========================================================================
//  CCD — Additive Continuous Collision Detection (ACCD)
//
//  Based on Codim-IPC [Li et al. 2021]. Uses conservative advancement:
//  iteratively steps forward by a safe amount based on current distance
//  and maximum displacement magnitude. Much more robust than cubic-solver
//  CCD in floating point.
// =========================================================================
namespace ccd
{

static constexpr long ACCD_MAX_ITER = 1000000;
static constexpr double ACCD_CONSERVATIVE_RESCALING = 0.1;

// Core ACCD algorithm: given stacked positions x, displacement dx,
// a distance function, and max displacement magnitude, find TOI.
// Returns true if collision found, false otherwise.
static bool additiveCCD(
  V12d x,  // mutable copy of stacked positions
  const V12d &dx,
  const std::function<double(const V12d &)> &distanceSquared,
  double maxDispMag,
  double &toi,
  double minDistance,
  double tmax)
{
  const double eta = ACCD_CONSERVATIVE_RESCALING;
  const double minDistSq = minDistance * minDistance;

  double d_sq = distanceSquared(x);
  double d = std::sqrt(d_sq);
  if (d <= minDistance) {
    toi = 0.0;
    return true;
  }

  double d_func = d_sq - minDistSq;
  // gap = η * (d² - ξ²) / (d + ξ)
  const double gap = eta * d_func / (d + minDistance);

  toi = 0.0;
  for (long i = 0; i < ACCD_MAX_ITER; ++i) {
    // Conservative lower bound on step:
    //   t_l = (1 - η) * (d² - ξ²) / ((d + ξ) * l_p)
    double toi_lower = (1.0 - eta) * d_func / ((d + minDistance) * maxDispMag);

    x += toi_lower * dx;

    d = std::sqrt(d_sq = distanceSquared(x));
    d_func = d_sq - minDistSq;

    if (d_func <= 0.0) {
      // Shouldn't happen with conservative stepping, but handle it
      break;
    }

    if (toi > 0.0 && d_func / (d + minDistance) < gap) {
      break;  // distance has decreased below gap threshold
    }

    toi += toi_lower;
    if (toi > tmax) {
      return false;  // no collision within tmax
    }
  }

  return true;
}

// Helper: subtract mean displacement for translation invariance
template<typename... Args>
static void subtractMean(Args &...args)
{
  constexpr int n = sizeof...(args);
  V3d mean = V3d::Zero();
  for (const V3d &v : { args... })
    mean += v;
  mean /= n;
  for (V3d *v : { &args... })
    *v -= mean;
}

double pointTriangleCCD(const V3d &p, const V3d &t0,
  const V3d &t1, const V3d &t2,
  const V3d &dp, const V3d &dt0,
  const V3d &dt1, const V3d &dt2,
  double thickness, double tMax)
{
  // Check initial distance
  double initDist = distance::computePTSqDist(p, t0, t1, t2);
  if (initDist <= thickness * thickness) {
    return 0.0;
  }

  // Compute displacements and subtract mean for translation invariance
  V3d ddp = dp;
  V3d ddt0 = dt0;
  V3d ddt1 = dt1;
  V3d ddt2 = dt2;
  subtractMean(ddp, ddt0, ddt1, ddt2);

  // Max displacement magnitude
  double maxDisp = ddp.norm() + std::sqrt(std::max({ ddt0.squaredNorm(), ddt1.squaredNorm(), ddt2.squaredNorm() }));
  if (maxDisp < 1e-30)
    return tMax;

  // Stack positions and displacements
  V12d x, dx;
  x.segment<3>(0) = p;
  x.segment<3>(3) = t0;
  x.segment<3>(6) = t1;
  x.segment<3>(9) = t2;

  dx.segment<3>(0) = ddp;
  dx.segment<3>(3) = ddt0;
  dx.segment<3>(6) = ddt1;
  dx.segment<3>(9) = ddt2;

  auto distSq = [](const V12d &xx) -> double {
    return distance::computePTSqDist(
      xx.segment<3>(0), xx.segment<3>(3),
      xx.segment<3>(6), xx.segment<3>(9));
  };

  double toi;
  if (additiveCCD(x, dx, distSq, maxDisp, toi, thickness, tMax)) {
    return toi;
  }
  return tMax;
}

double edgeEdgeCCD(const V3d &ea0, const V3d &ea1,
  const V3d &eb0, const V3d &eb1,
  const V3d &dea0, const V3d &dea1,
  const V3d &deb0, const V3d &deb1,
  double thickness, double tMax)
{
  // Check initial distance
  double initDist = distance::computeEESqDist(ea0, ea1, eb0, eb1);
  if (initDist <= thickness * thickness) {
    return 0.0;
  }

  // Compute displacements and subtract mean
  V3d ddea0 = dea0, ddea1 = dea1;
  V3d ddeb0 = deb0, ddeb1 = deb1;
  subtractMean(ddea0, ddea1, ddeb0, ddeb1);

  double maxDisp =
    std::sqrt(std::max(ddea0.squaredNorm(), ddea1.squaredNorm())) + std::sqrt(std::max(ddeb0.squaredNorm(), ddeb1.squaredNorm()));
  if (maxDisp < 1e-30)
    return tMax;

  double minDistSq = thickness * thickness;
  V12d x, dx;
  x.segment<3>(0) = ea0;
  x.segment<3>(3) = ea1;
  x.segment<3>(6) = eb0;
  x.segment<3>(9) = eb1;
  dx.segment<3>(0) = ddea0;
  dx.segment<3>(3) = ddea1;
  dx.segment<3>(6) = ddeb0;
  dx.segment<3>(9) = ddeb1;

  auto distSq = [minDistSq](const V12d &xx) -> double {
    double d_sq = distance::computeEESqDist(
      xx.segment<3>(0), xx.segment<3>(3),
      xx.segment<3>(6), xx.segment<3>(9));

    if (d_sq - minDistSq <= 0) {
      // Near-parallel edges: fall back to minimum PP distance
      const auto &a0 = xx.segment<3>(0);
      const auto &a1 = xx.segment<3>(3);
      const auto &b0 = xx.segment<3>(6);
      const auto &b1 = xx.segment<3>(9);
      d_sq = std::min({ (a0 - b0).squaredNorm(), (a0 - b1).squaredNorm(),
        (a1 - b0).squaredNorm(), (a1 - b1).squaredNorm() });
    }
    return d_sq;
  };

  double toi;
  if (additiveCCD(x, dx, distSq, maxDisp, toi, thickness, tMax)) {
    return toi;
  }
  return tMax;
}

}  // namespace ccd

// =========================================================================
//  PSD projection of 12x12 symmetric matrix
// =========================================================================
M12d projectToPSD(const M12d &H)
{  
  Eigen::SelfAdjointEigenSolver<M12d> es(H);
  const auto &evals = es.eigenvalues();

  // Early exit: eigenvalues are sorted ascending, so if the smallest is >= 0
  // the matrix is already PSD.
  if (evals(0) >= 0.0)
    return H;

  // Clamp negative eigenvalues to zero
  Eigen::DiagonalMatrix<double, 12> D(evals);
  for (int i = 0; i < 12; ++i) {
    if (D.diagonal()(i) < 0.0)
      D.diagonal()(i) = 0.0;
    else
      break;  // remaining eigenvalues are >= 0
  }

  return es.eigenvectors() * D * es.eigenvectors().transpose();
}

// =========================================================================
//  CollisionIPC  —  mesh setup
// =========================================================================

void CIPCSolver::setMesh(const MXd &V, const MXi &F)
{
  numVerts_ = (int)V.rows();
  triangles_.resize(F.rows());
  for (int i = 0; i < (int)F.rows(); ++i)
    triangles_[i] = { F(i, 0), F(i, 1), F(i, 2) };

  buildEdges();
  buildAdjacency();
  buildAreaWeights(V);
}

void CIPCSolver::setMesh(int numVerts,
  const std::vector<std::array<int, 3>> &triangles)
{
  numVerts_ = numVerts;
  triangles_ = triangles;
  buildEdges();
  buildAdjacency();
  // Area weights not computed — caller must provide positions via the other overload
  vertexArea_.assign(numVerts_, 1.0);
  triArea_.assign(triangles_.size(), 1.0);
  edgeLength_.assign(edges_.size(), 1.0);
}

void CIPCSolver::buildEdges()
{
  std::set<std::pair<int, int>> edgeSet;
  for (auto &tri : triangles_) {
    for (int j = 0; j < 3; ++j) {
      int a = tri[j], b = tri[(j + 1) % 3];
      if (a > b)
        std::swap(a, b);
      edgeSet.insert({ a, b });
    }
  }
  edges_.clear();
  edges_.reserve(edgeSet.size());
  for (auto &e : edgeSet)
    edges_.push_back({ e.first, e.second });
}

void CIPCSolver::buildAdjacency()
{
  vertexTriAdj_.clear();
  for (int fi = 0; fi < (int)triangles_.size(); ++fi) {
    for (int j = 0; j < 3; ++j)
      vertexTriAdj_.insert(packPair(triangles_[fi][j], fi));
  }

  edgeVertAdj_.clear();
  for (int ei = 0; ei < (int)edges_.size(); ++ei) {
    edgeVertAdj_.insert(packPair(ei, edges_[ei][0]));
    edgeVertAdj_.insert(packPair(ei, edges_[ei][1]));
  }
}

void CIPCSolver::buildAreaWeights(const MXd &V)
{
  int nTri = (int)triangles_.size();
  int nEdge = (int)edges_.size();

  // Triangle areas
  triArea_.resize(nTri);
  for (int fi = 0; fi < nTri; ++fi) {
    V3d v0 = V.row(triangles_[fi][0]).transpose();
    V3d v1 = V.row(triangles_[fi][1]).transpose();
    V3d v2 = V.row(triangles_[fi][2]).transpose();
    triArea_[fi] = 0.5 * (v1 - v0).cross(v2 - v0).norm();
  }

  // Lumped vertex area: 1/3 of sum of adjacent triangle areas
  vertexArea_.assign(numVerts_, 0.0);
  for (int fi = 0; fi < nTri; ++fi) {
    double a3 = triArea_[fi] / 3.0;
    for (int j = 0; j < 3; ++j)
      vertexArea_[triangles_[fi][j]] += a3;
  }

  // Edge rest lengths
  edgeLength_.resize(nEdge);
  for (int ei = 0; ei < nEdge; ++ei) {
    V3d e0 = V.row(edges_[ei][0]).transpose();
    V3d e1 = V.row(edges_[ei][1]).transpose();
    edgeLength_[ei] = (e1 - e0).norm();
  }
}

// =========================================================================
//  Spatial hash grid for O(n) broad-phase collision detection
// =========================================================================

namespace
{

struct AABB
{
  V3d lo, hi;
  void init(const V3d &v, double pad)
  {
    lo = v.array() - pad;
    hi = v.array() + pad;
  }
  void expand(const V3d &v, double pad)
  {
    lo = (v.array() - pad).cwiseMin(lo.array());
    hi = (v.array() + pad).cwiseMax(hi.array());
  }
  void expand(const V3d &v)
  {
    lo = lo.cwiseMin(v);
    hi = hi.cwiseMax(v);
  }
  bool overlaps(const AABB &o) const
  {
    return (lo.array() <= o.hi.array()).all() &&
      (o.lo.array() <= hi.array()).all();
  }
};

// Spatial hash: maps cell keys to lists of primitive indices.
// Uses insert-then-query pattern: insert one primitive type, query with another.
struct SpatialHash
{
  double cellSize;
  std::unordered_map<int64_t, std::vector<int>> cells;

  SpatialHash(int capacity)
  {
    cells.reserve(capacity);
  }

  void clear() { cells.clear(); }

  static int64_t hashCoord(int ix, int iy, int iz)
  {
    // Large primes for spatial hashing (from Teschner et al. 2003)
    constexpr int64_t p1 = 73856093LL;
    constexpr int64_t p2 = 19349663LL;
    constexpr int64_t p3 = 83492791LL;
    return ((int64_t)ix * p1) ^ ((int64_t)iy * p2) ^ ((int64_t)iz * p3);
  }

  void toGrid(const V3d &p, int &ix, int &iy, int &iz) const
  {
    ix = (int)std::floor(p.x() / cellSize);
    iy = (int)std::floor(p.y() / cellSize);
    iz = (int)std::floor(p.z() / cellSize);
  }

  // Insert a primitive's AABB into all overlapping cells
  void insert(const AABB &box, int index)
  {
    int lo_ix, lo_iy, lo_iz, hi_ix, hi_iy, hi_iz;
    toGrid(box.lo, lo_ix, lo_iy, lo_iz);
    toGrid(box.hi, hi_ix, hi_iy, hi_iz);

    for (int iz = lo_iz; iz <= hi_iz; ++iz)
      for (int iy = lo_iy; iy <= hi_iy; ++iy)
        for (int ix = lo_ix; ix <= hi_ix; ++ix)
          cells[hashCoord(ix, iy, iz)].push_back(index);
  }

  // Query: collect all indices stored in cells that overlap the given AABB.
  // Uses a visited flag array to deduplicate and avoid self-matches.
  void query(const AABB &box, int selfIdx,
    std::vector<int> &visited_stamp, int stamp,
    std::vector<int> &result) const
  {
    int lo_ix, lo_iy, lo_iz, hi_ix, hi_iy, hi_iz;
    toGrid(box.lo, lo_ix, lo_iy, lo_iz);
    toGrid(box.hi, hi_ix, hi_iy, hi_iz);

    for (int iz = lo_iz; iz <= hi_iz; ++iz)
      for (int iy = lo_iy; iy <= hi_iy; ++iy)
        for (int ix = lo_ix; ix <= hi_ix; ++ix) {
          auto it = cells.find(hashCoord(ix, iy, iz));
          if (it == cells.end())
            continue;
          for (int idx : it->second) {
            if (idx == selfIdx)
              continue;
            if (visited_stamp[idx] == stamp)
              continue;
            visited_stamp[idx] = stamp;
            result.push_back(idx);
          }
        }
  }
};

}  // anonymous namespace

// =========================================================================
//  Broad phase: find candidate PT and EE pairs using spatial hashing
//  Uses insert-then-query: insert one type, query with the other.
// =========================================================================
void CIPCSolver::findCollisionPairs(const VXd &x)
{
  ptPairs_.clear();
  eePairs_.clear();

  const double inflate = dhat;

  auto getV = [&](int i) -> V3d {
    return x.segment<3>(3 * i);
  };

  int nTri = (int)triangles_.size();
  int nEdge = (int)edges_.size();

  // --- Build AABBs (parallel) ---
  std::vector<AABB> vertBox(numVerts_);
  tbb::parallel_for(tbb::blocked_range<int>(0, numVerts_),
    [&](const tbb::blocked_range<int> &r) {
      for (int vi = r.begin(); vi < r.end(); ++vi)
        vertBox[vi].init(getV(vi), inflate);
    });

  std::vector<AABB> triBox(nTri);
  tbb::parallel_for(tbb::blocked_range<int>(0, nTri),
    [&](const tbb::blocked_range<int> &r) {
      for (int fi = r.begin(); fi < r.end(); ++fi) {
        auto &tri = triangles_[fi];
        triBox[fi].init(getV(tri[0]), inflate);
        triBox[fi].expand(getV(tri[1]), inflate);
        triBox[fi].expand(getV(tri[2]), inflate);
      }
    });

  std::vector<AABB> edgeBox(nEdge);
  tbb::parallel_for(tbb::blocked_range<int>(0, nEdge),
    [&](const tbb::blocked_range<int> &r) {
      for (int ei = r.begin(); ei < r.end(); ++ei) {
        edgeBox[ei].init(getV(edges_[ei][0]), inflate);
        edgeBox[ei].expand(getV(edges_[ei][1]), inflate);
      }
    });

  double avgBoxDiag = 0;
  for (const auto &aabb : triBox) {
    avgBoxDiag += (aabb.hi - aabb.lo).norm();
  }

  double cellSize = nTri > 0 ? std::max(avgBoxDiag / nTri, 1e-6) : std::max(1e-6, dhat);
  double dhat2 = dhat * dhat;

  // --- PT pairs: insert triangles (serial), query with vertices (parallel) ---
  {
    SpatialHash triHash(nTri);
    triHash.cellSize = cellSize;
    for (int fi = 0; fi < nTri; ++fi)
      triHash.insert(triBox[fi], fi);

    // Thread-local buffers: each thread gets its own visited array,
    // candidates vector, and result pairs. The visited stamp dedup is
    // per-query (prevents same triangle appearing from multiple cells),
    // not cross-query, so thread-local is correct.
    tbb::enumerable_thread_specific<std::vector<int>> tls_visited(
      [nTri]() { return std::vector<int>(nTri, 0); });
    tbb::enumerable_thread_specific<std::vector<int>> tls_candidates;
    tbb::enumerable_thread_specific<std::vector<PTPair>> tls_pairs;

    tbb::parallel_for(tbb::blocked_range<int>(0, numVerts_),
      [&](const tbb::blocked_range<int> &range) {
        auto &visited = tls_visited.local();
        auto &candidates = tls_candidates.local();
        auto &localPairs = tls_pairs.local();

        for (int vi = range.begin(); vi < range.end(); ++vi) {
          candidates.clear();
          triHash.query(vertBox[vi], -1, visited, vi + 1, candidates);

          for (int fi : candidates) {
            auto &tri = triangles_[fi];
            if (vi == tri[0] || vi == tri[1] || vi == tri[2])
              continue;
            if (!vertBox[vi].overlaps(triBox[fi]))
              continue;

            V3d vp = getV(vi);
            V3d vt0 = getV(tri[0]), vt1 = getV(tri[1]), vt2 = getV(tri[2]);
            double d2 = distance::computePTSqDist(vp, vt0, vt1, vt2);
            if (d2 < dhat2)
              localPairs.push_back({ vi, tri[0], tri[1], tri[2],
                vertexArea_[vi] * triArea_[fi] });
          }
        }
      });

    for (auto &lp : tls_pairs)
      ptPairs_.insert(ptPairs_.end(), lp.begin(), lp.end());
  }

  // --- EE pairs: insert edges (serial), query with edges (parallel) ---
  {
    SpatialHash edgeHash(nEdge);
    edgeHash.cellSize = cellSize;
    for (int ei = 0; ei < nEdge; ++ei)
      edgeHash.insert(edgeBox[ei], ei);

    tbb::enumerable_thread_specific<std::vector<int>> tls_visited(
      [nEdge]() { return std::vector<int>(nEdge, 0); });
    tbb::enumerable_thread_specific<std::vector<int>> tls_candidates;
    tbb::enumerable_thread_specific<std::vector<EEPair>> tls_pairs;

    tbb::parallel_for(tbb::blocked_range<int>(0, nEdge),
      [&](const tbb::blocked_range<int> &range) {
        auto &visited = tls_visited.local();
        auto &candidates = tls_candidates.local();
        auto &localPairs = tls_pairs.local();

        for (int ei = range.begin(); ei < range.end(); ++ei) {
          candidates.clear();
          edgeHash.query(edgeBox[ei], ei, visited, ei + 1, candidates);

          int a0 = edges_[ei][0], a1 = edges_[ei][1];
          for (int ej : candidates) {
            if (ej <= ei)
              continue;

            int b0 = edges_[ej][0], b1 = edges_[ej][1];
            if (a0 == b0 || a0 == b1 || a1 == b0 || a1 == b1)
              continue;
            if (!edgeBox[ei].overlaps(edgeBox[ej]))
              continue;

            V3d va0 = getV(a0), va1 = getV(a1);
            V3d vb0 = getV(b0), vb1 = getV(b1);
            double d2 = distance::computeEESqDist(va0, va1, vb0, vb1);
            if (d2 < dhat2)
              localPairs.push_back({ a0, a1, b0, b1,
                edgeLength_[ei] * edgeLength_[ej] });
          }
        }
      });

    for (auto &lp : tls_pairs)
      eePairs_.insert(eePairs_.end(), lp.begin(), lp.end());
  }
}

// =========================================================================
//  1)  Maximum step size  (CCD-based line search with spatial hashing)
// =========================================================================
double CIPCSolver::computeMaxStepSize(const VXd &x, const VXd &dx,
  double slackness) const
{
  auto getV = [&](int i) -> V3d {
    return x.segment<3>(3 * i);
  };
  auto getdV = [&](int i) -> V3d {
    return dx.segment<3>(3 * i);
  };

  int nTri = (int)triangles_.size();
  int nEdge = (int)edges_.size();

  // Build swept-volume AABBs (parallel)
  std::vector<AABB> vertBox(numVerts_);
  tbb::parallel_for(tbb::blocked_range<int>(0, numVerts_),
    [&](const tbb::blocked_range<int> &r) {
      for (int vi = r.begin(); vi < r.end(); ++vi) {
        V3d p0 = getV(vi), p1 = p0 + getdV(vi);
        vertBox[vi].init(p0, 0.0);
        vertBox[vi].expand(p1);
      }
    });

  std::vector<AABB> triBox(nTri);
  tbb::parallel_for(tbb::blocked_range<int>(0, nTri),
    [&](const tbb::blocked_range<int> &r) {
      for (int fi = r.begin(); fi < r.end(); ++fi) {
        auto &tri = triangles_[fi];
        V3d v0 = getV(tri[0]), v1 = getV(tri[1]), v2 = getV(tri[2]);
        V3d d0 = getdV(tri[0]), d1 = getdV(tri[1]), d2 = getdV(tri[2]);
        triBox[fi].init(v0, 0.0);
        triBox[fi].expand(v1);
        triBox[fi].expand(v2);
        triBox[fi].expand(v0 + d0);
        triBox[fi].expand(v1 + d1);
        triBox[fi].expand(v2 + d2);
      }
    });

  std::vector<AABB> edgeBox(nEdge);
  tbb::parallel_for(tbb::blocked_range<int>(0, nEdge),
    [&](const tbb::blocked_range<int> &r) {
      for (int ei = r.begin(); ei < r.end(); ++ei) {
        V3d a0 = getV(edges_[ei][0]), a1 = getV(edges_[ei][1]);
        V3d da0 = getdV(edges_[ei][0]), da1 = getdV(edges_[ei][1]);
        edgeBox[ei].init(a0, 0.0);
        edgeBox[ei].expand(a1);
        edgeBox[ei].expand(a0 + da0);
        edgeBox[ei].expand(a1 + da1);
      }
    });

  // Cell size: average swept-AABB diagonal (parallel reduction)
  double avgBoxDiag = tbb::parallel_reduce(
    tbb::blocked_range<int>(0, nTri), 0.0,
    [&](const tbb::blocked_range<int> &r, double sum) {
      for (int fi = r.begin(); fi < r.end(); ++fi)
        sum += (triBox[fi].hi - triBox[fi].lo).norm();
      return sum;
    },
    std::plus<double>());
  double cellSize = nTri > 0 ? std::max(avgBoxDiag / nTri, 1e-6) : std::max(1e-6, dhat);

  double alpha = 1.0;

  // --- PT CCD: insert triangles (serial), query with vertices (parallel reduce) ---
  {
    SpatialHash triHash(nTri);
    triHash.cellSize = cellSize;
    for (int fi = 0; fi < nTri; ++fi)
      triHash.insert(triBox[fi], fi);

    tbb::enumerable_thread_specific<std::vector<int>> tls_visited(
      [nTri]() { return std::vector<int>(nTri, 0); });
    tbb::enumerable_thread_specific<std::vector<int>> tls_candidates;

    // {
    //   double alpha = 1.0;
    //   auto &visited = tls_visited.local();
    //   auto &candidates = tls_candidates.local();

    //   for (int vi = 0; vi < numVerts_; ++vi) {
    //     candidates.clear();
    //     triHash.query(vertBox[vi], -1, visited, vi + 1, candidates);

    //     for (int fi : candidates) {
    //       auto &tri = triangles_[fi];
    //       if (vi == tri[0] || vi == tri[1] || vi == tri[2])
    //         continue;
    //       if (!vertBox[vi].overlaps(triBox[fi]))
    //         continue;

    //       V3d p = getV(vi), dp = getdV(vi);
    //       V3d t0 = getV(tri[0]), dt0 = getdV(tri[0]);
    //       V3d t1 = getV(tri[1]), dt1 = getdV(tri[1]);
    //       V3d t2 = getV(tri[2]), dt2 = getdV(tri[2]);

    //       double toi = ccd::pointTriangleCCD(p, t0, t1, t2,
    //         dp, dt0, dt1, dt2,
    //         0.0, alpha);

    //       if (toi < 1e-3) {
    //         toi = ccd::pointTriangleCCD(p, t0, t1, t2,
    //           dp, dt0, dt1, dt2,
    //           0.0, 1.0);

    //         V3d q = p + dp;
    //         V3d s0 = t0 + dt0;
    //         V3d s1 = t1 + dt1;
    //         V3d s2 = t2 + dt2;

    //         std::cout << p[0] << ',' << p[1] << ',' << p[2] << std::endl;
    //         std::cout << q[0] << ',' << q[1] << ',' << q[2] << std::endl;

    //         std::cout << t0[0] << ',' << t0[1] << ',' << t0[2] << std::endl;
    //         std::cout << t1[0] << ',' << t1[1] << ',' << t1[2] << std::endl;
    //         std::cout << t2[0] << ',' << t2[1] << ',' << t2[2] << std::endl;

    //         std::cout << s0[0] << ',' << s0[1] << ',' << s0[2] << std::endl;
    //         std::cout << s1[0] << ',' << s1[1] << ',' << s1[2] << std::endl;
    //         std::cout << s2[0] << ',' << s2[1] << ',' << s2[2] << std::endl;

    //         std::ofstream outfile("/home/user/code/libpgo-private/data/deb.obj");
    //         V3d e(1e-5, 1e-5, 1e-5);

    //         // Triangle: p, p+dp, p+e
    //         auto wv = [&](const V3d &v) {
    //           outfile << "v " << v[0] << " " << v[1] << " " << v[2] << "\n";
    //         };
    //         wv(p);       // 1
    //         wv(p + dp);  // 2
    //         wv(p + e);   // 3
    //         outfile << "f 1 2 3\n";

    //         // Prism: base t0,t1,t2 / top t0+dt0,t1+dt1,t2+dt2
    //         wv(t0);                    // 4
    //         wv(t1);                    // 5
    //         wv(t2);                    // 6
    //         wv(t0 + dt0);              // 7
    //         wv(t1 + dt1);              // 8
    //         wv(t2 + dt2);              // 9
    //         outfile << "f 4 5 6\n";    // bottom
    //         outfile << "f 7 9 8\n";    // top (flipped for outward normal)
    //         outfile << "f 4 5 8 7\n";  // side 01
    //         outfile << "f 5 6 9 8\n";  // side 12
    //         outfile << "f 6 4 7 9\n";  // side 20
    //       }

    //       if (toi < alpha)
    //         alpha = toi * slackness;
    //     }
    //   }
    //   std::cout << "!" << alpha << std::endl;
    // }

    alpha = tbb::parallel_reduce(
      tbb::blocked_range<int>(0, numVerts_), 1.0,
      [&](const tbb::blocked_range<int> &range, double localAlpha) {
        auto &visited = tls_visited.local();
        auto &candidates = tls_candidates.local();

        for (int vi = range.begin(); vi < range.end(); ++vi) {
          candidates.clear();
          triHash.query(vertBox[vi], -1, visited, vi + 1, candidates);

          for (int fi : candidates) {
            auto &tri = triangles_[fi];
            if (vi == tri[0] || vi == tri[1] || vi == tri[2])
              continue;
            if (!vertBox[vi].overlaps(triBox[fi]))
              continue;

            V3d p = getV(vi), dp = getdV(vi);
            V3d t0 = getV(tri[0]), dt0 = getdV(tri[0]);
            V3d t1 = getV(tri[1]), dt1 = getdV(tri[1]);
            V3d t2 = getV(tri[2]), dt2 = getdV(tri[2]);

            double toi = ccd::pointTriangleCCD(p, t0, t1, t2,
              dp, dt0, dt1, dt2,
              0.0, localAlpha);
            if (toi < localAlpha)
              localAlpha = toi * slackness;
          }
        }
        return localAlpha;
      },
      [](double a, double b) { return std::min(a, b); });
  }

  // --- EE CCD: insert edges (serial), query with edges (parallel reduce) ---
  {
    SpatialHash edgeHash(nEdge);
    edgeHash.cellSize = cellSize;
    for (int ei = 0; ei < nEdge; ++ei)
      edgeHash.insert(edgeBox[ei], ei);

    tbb::enumerable_thread_specific<std::vector<int>> tls_visited(
      [nEdge]() { return std::vector<int>(nEdge, 0); });
    tbb::enumerable_thread_specific<std::vector<int>> tls_candidates;

    alpha = tbb::parallel_reduce(
      tbb::blocked_range<int>(0, nEdge), alpha,
      [&](const tbb::blocked_range<int> &range, double localAlpha) {
        auto &visited = tls_visited.local();
        auto &candidates = tls_candidates.local();

        for (int ei = range.begin(); ei < range.end(); ++ei) {
          candidates.clear();
          edgeHash.query(edgeBox[ei], ei, visited, ei + 1, candidates);

          int a0 = edges_[ei][0], a1 = edges_[ei][1];
          for (int ej : candidates) {
            if (ej <= ei)
              continue;

            int b0 = edges_[ej][0], b1 = edges_[ej][1];
            if (a0 == b0 || a0 == b1 || a1 == b0 || a1 == b1)
              continue;
            if (!edgeBox[ei].overlaps(edgeBox[ej]))
              continue;

            V3d va0 = getV(a0), va1 = getV(a1);
            V3d vb0 = getV(b0), vb1 = getV(b1);
            V3d da0 = getdV(a0), da1 = getdV(a1);
            V3d db0 = getdV(b0), db1 = getdV(b1);

            double toi = ccd::edgeEdgeCCD(va0, va1, vb0, vb1,
              da0, da1, db0, db1,
              0.0, localAlpha);
            if (toi < localAlpha)
              localAlpha = toi * slackness;
          }
        }
        return localAlpha;
      },
      [](double a, double b) { return std::min(a, b); });
  }

  return std::max(alpha, 1e-12);
}

// =========================================================================
//  2)  Energy
// =========================================================================
double CIPCSolver::computeEnergy(const VXd &x)
{
  double ee_eps = eps_ee;

  double dhat2 = dhat * dhat;

  // PT pairs
  double ptEnergy = tbb::parallel_reduce(
    tbb::blocked_range<int>(0, (int)ptPairs_.size()), 0.0,
    [&](const tbb::blocked_range<int> &range, double localE) {
      for (int i = range.begin(); i < range.end(); ++i) {
        auto &pair = ptPairs_[i];
        V3d p = vtx(x, pair.p);
        V3d t0 = vtx(x, pair.t0);
        V3d t1 = vtx(x, pair.t1);
        V3d t2 = vtx(x, pair.t2);
        double d2 = distance::computePTSqDist(p, t0, t1, t2);
        if (d2 < dhat2 && d2 > 0.0)
          localE += pair.weight * kappa * barrier::b(d2, dhat2);
      }
      return localE;
    },
    std::plus<double>());

  // EE pairs
  double eeEnergy = tbb::parallel_reduce(
    tbb::blocked_range<int>(0, (int)eePairs_.size()), 0.0,
    [&](const tbb::blocked_range<int> &range, double localE) {
      for (int i = range.begin(); i < range.end(); ++i) {
        auto &pair = eePairs_[i];
        V3d ea0 = vtx(x, pair.ea0);
        V3d ea1 = vtx(x, pair.ea1);
        V3d eb0 = vtx(x, pair.eb0);
        V3d eb1 = vtx(x, pair.eb1);
        double d2 = distance::computeEESqDist(ea0, ea1, eb0, eb1);
        if (d2 < dhat2 && d2 > 0.0) {
          double m = 1.0;
          if (ee_eps > 0.0)
            m = distance::eeMollifier(ea0, ea1, eb0, eb1, ee_eps);
          localE += pair.weight * kappa * m * barrier::b(d2, dhat2);
        }
      }
      return localE;
    },
    std::plus<double>());

  // Floor penalty: 0.5 * kappa * (z - h)^2  when z < h
  double floorEnergy = 0.0;
  if (useFloor) {
    floorEnergy = tbb::parallel_reduce(
      tbb::blocked_range<int>(0, numVerts_), 0.0,
      [&](const tbb::blocked_range<int> &range, double localE) {
        for (int vi = range.begin(); vi < range.end(); ++vi) {
          double dz = x[3 * vi + 2] - floorHeight;
          if (dz < 0.0)
            localE += 0.5 * floorKappa * dz * dz;
        }
        return localE;
      },
      std::plus<double>());
  }

  return ptEnergy + eeEnergy + floorEnergy;
}

// =========================================================================
//  2)  Gradient
// =========================================================================
void CIPCSolver::computeGradient(const VXd &x, VXd &grad)
{
  int n = 3 * numVerts_;
  if (grad.size() != n)
    grad.setZero(n);

  // Atomic scatter: accumulate a 12-vector into global gradient
  auto scatter = [&](const V12d &local, const int idx[4]) {
    double *gdata = grad.data();
    for (int i = 0; i < 4; ++i)
      if (idx[i] >= 0)
        for (int d = 0; d < 3; ++d)
          std::atomic_ref<double>(gdata[3 * idx[i] + d])
            .fetch_add(local[3 * i + d], std::memory_order_relaxed);
  };

  double dhat2 = dhat * dhat;

  // PT pairs
  tbb::parallel_for(
    tbb::blocked_range<int>(0, (int)ptPairs_.size()),
    [&](const tbb::blocked_range<int> &range) {
      for (int i = range.begin(); i < range.end(); ++i) {
        auto &pair = ptPairs_[i];
        V3d p = vtx(x, pair.p);
        V3d t0 = vtx(x, pair.t0);
        V3d t1 = vtx(x, pair.t1);
        V3d t2 = vtx(x, pair.t2);

        double d2 = distance::computePTSqDist(p, t0, t1, t2);
        if (d2 >= dhat2 || d2 <= 0.0)
          continue;

        V12d gd2 = distance::computePTSqDistGrad(p, t0, t1, t2);
        double coeff = pair.weight * kappa * barrier::dbds(d2, dhat2);
        V12d gE = coeff * gd2;

        int idx[4] = { pair.p, pair.t0, pair.t1, pair.t2 };
        scatter(gE, idx);
      }
    });

  // EE pairs
  double ee_eps = eps_ee;
  tbb::parallel_for(
    tbb::blocked_range<int>(0, (int)eePairs_.size()),
    [&](const tbb::blocked_range<int> &range) {
      for (int i = range.begin(); i < range.end(); ++i) {
        auto &pair = eePairs_[i];
        V3d ea0 = vtx(x, pair.ea0);
        V3d ea1 = vtx(x, pair.ea1);
        V3d eb0 = vtx(x, pair.eb0);
        V3d eb1 = vtx(x, pair.eb1);

        double d2 = distance::computeEESqDist(ea0, ea1, eb0, eb1);
        if (d2 >= dhat2 || d2 <= 0.0)
          continue;

        double wk = pair.weight * kappa;
        V12d gd2 = distance::computeEESqDistGrad(ea0, ea1, eb0, eb1);
        double b_val = barrier::b(d2, dhat2);
        double dbv = barrier::dbds(d2, dhat2);

        V12d gE;
        if (ee_eps > 0.0) {
          double m = distance::eeMollifier(ea0, ea1, eb0, eb1, ee_eps);
          V12d gm = distance::eeMollifierGrad(ea0, ea1, eb0, eb1, ee_eps);
          gE = wk * (gm * b_val + m * dbv * gd2);
        }
        else {
          gE = wk * dbv * gd2;
        }

        int idx[4] = { pair.ea0, pair.ea1, pair.eb0, pair.eb1 };
        scatter(gE, idx);
      }
    });

  // Floor penalty gradient: floorKappa * (z - h) when z < h
  if (useFloor) {
    tbb::parallel_for(
      tbb::blocked_range<int>(0, numVerts_),
      [&](const tbb::blocked_range<int> &range) {
        double *gdata = grad.data();
        for (int vi = range.begin(); vi < range.end(); ++vi) {
          double dz = x[3 * vi + 2] - floorHeight;
          if (dz < 0.0)
            gdata[3 * vi + 2] += floorKappa * dz;
        }
      });
  }
}

// =========================================================================
//  2)  Sparse Hessian
// =========================================================================
void CIPCSolver::computeHessian(const VXd &x, SpMatD &hess)
{
  int n = 3 * numVerts_;
  int nPT = (int)ptPairs_.size();
  int nEE = (int)eePairs_.size();
  int totalPairs = nPT + nEE;

  // Preallocate exactly 144 triplet slots per pair (12x12 block).
  // Each pair writes to its own contiguous slice — no conflicts.
  std::vector<TripletD> triplets(144 * totalPairs, TripletD(0, 0, 0.0));

  // Write 144 triplets from a 12x12 local hessian into a preallocated slice.
  auto scatterH = [&](int pairIdx, const M12d &localH, const int idx[4]) {
    TripletD *base = triplets.data() + 144 * pairIdx;
    int k = 0;
    for (int i = 0; i < 4; ++i) {
      int ri = (idx[i] >= 0) ? 3 * idx[i] : 0;
      bool vi = (idx[i] >= 0);
      for (int j = 0; j < 4; ++j) {
        int cj = (idx[j] >= 0) ? 3 * idx[j] : 0;
        bool valid = vi && (idx[j] >= 0);
        for (int di = 0; di < 3; ++di)
          for (int dj = 0; dj < 3; ++dj, ++k)
            base[k] = TripletD(ri + di, cj + dj,
              valid ? localH(3 * i + di, 3 * j + dj) : 0.0);
      }
    }
  };

  double dhat2 = dhat * dhat;

  // PT pairs
  tbb::parallel_for(
    tbb::blocked_range<int>(0, nPT),
    [&](const tbb::blocked_range<int> &range) {
      for (int i = range.begin(); i < range.end(); ++i) {
        auto &pair = ptPairs_[i];
        V3d p = vtx(x, pair.p);
        V3d t0 = vtx(x, pair.t0);
        V3d t1 = vtx(x, pair.t1);
        V3d t2 = vtx(x, pair.t2);

        double d2 = distance::computePTSqDist(p, t0, t1, t2);
        if (d2 >= dhat2 || d2 <= 0.0)
          continue;  // slot already zeroed

        V12d gd2 = distance::computePTSqDistGrad(p, t0, t1, t2);
        M12d Hd2 = distance::computePTSqDistHess(p, t0, t1, t2);
        double wk = pair.weight * kappa;
        double gp = barrier::dbds(d2, dhat2);
        double gpp = barrier::d2bds2(d2, dhat2);
        M12d localH = wk * (gpp * gd2 * gd2.transpose() + gp * Hd2);
        localH = projectToPSD(localH);

        int idx[4] = { pair.p, pair.t0, pair.t1, pair.t2 };
        scatterH(i, localH, idx);
      }
    });

  // EE pairs
  double ee_eps = eps_ee;
  tbb::parallel_for(
    tbb::blocked_range<int>(0, nEE),
    [&](const tbb::blocked_range<int> &range) {
      for (int i = range.begin(); i < range.end(); ++i) {
        auto &pair = eePairs_[i];
        V3d ea0 = vtx(x, pair.ea0);
        V3d ea1 = vtx(x, pair.ea1);
        V3d eb0 = vtx(x, pair.eb0);
        V3d eb1 = vtx(x, pair.eb1);

        double d2 = distance::computeEESqDist(ea0, ea1, eb0, eb1);
        if (d2 >= dhat2 || d2 <= 0.0)
          continue;  // slot already zeroed

        double wk = pair.weight * kappa;
        V12d gd2 = distance::computeEESqDistGrad(ea0, ea1, eb0, eb1);
        M12d Hd2 = distance::computeEESqDistHess(ea0, ea1, eb0, eb1);
        double gp = barrier::dbds(d2, dhat2);
        double gpp = barrier::d2bds2(d2, dhat2);

        M12d localH;
        if (ee_eps > 0.0) {
          double m = distance::eeMollifier(ea0, ea1, eb0, eb1, ee_eps);
          V12d gm = distance::eeMollifierGrad(ea0, ea1, eb0, eb1, ee_eps);
          M12d Hm = distance::eeMollifierHess(ea0, ea1, eb0, eb1, ee_eps);
          double bv = barrier::b(d2, dhat2);
          double dbv = barrier::dbds(d2, dhat2);
          V12d gb = dbv * gd2;
          localH = wk * (bv * Hm + gm * gb.transpose() + gb * gm.transpose() + m * (gpp * gd2 * gd2.transpose() + gp * Hd2));
        }
        else {
          localH = wk * (gpp * gd2 * gd2.transpose() + gp * Hd2);
        }
        localH = projectToPSD(localH);

        int idx[4] = { pair.ea0, pair.ea1, pair.eb0, pair.eb1 };
        scatterH(nPT + i, localH, idx);
      }
    });

  // Floor penalty hessian: floorKappa at diagonal (3*vi+2, 3*vi+2) when z < h
  if (useFloor) {
    tbb::enumerable_thread_specific<std::vector<TripletD>> tls_floorTrip;
    tbb::parallel_for(
      tbb::blocked_range<int>(0, numVerts_),
      [&](const tbb::blocked_range<int> &range) {
        auto &local = tls_floorTrip.local();
        for (int vi = range.begin(); vi < range.end(); ++vi) {
          double dz = x[3 * vi + 2] - floorHeight;
          if (dz < 0.0) {
            int row = 3 * vi + 2;
            local.emplace_back(row, row, floorKappa);
          }
        }
      });
    for (auto &local : tls_floorTrip)
      triplets.insert(triplets.end(), local.begin(), local.end());
  }

  hess.resize(n, n);
  hess.setFromTriplets(triplets.begin(), triplets.end());
}

// =========================================================================
//  Combined computation (single broad-phase pass)
// =========================================================================
void CIPCSolver::computeAll(const VXd &x,
  double &energy, VXd &grad, SpMatD &hess)
{
  findCollisionPairs(x);  // one broad-phase pass

  int n = 3 * numVerts_;
  int nPT = (int)ptPairs_.size();
  int nEE = (int)eePairs_.size();
  int totalPairs = nPT + nEE;

  energy = 0.0;
  grad.setZero(n);
  std::vector<TripletD> triplets(144 * totalPairs, TripletD(0, 0, 0.0));

  // Atomic scatter for gradient
  auto scatterG = [&](const V12d &local, const int idx[4]) {
    double *gdata = grad.data();
    for (int i = 0; i < 4; ++i)
      if (idx[i] >= 0)
        for (int d = 0; d < 3; ++d)
          std::atomic_ref<double>(gdata[3 * idx[i] + d])
            .fetch_add(local[3 * i + d], std::memory_order_relaxed);
  };

  // Preallocated scatter for hessian
  auto scatterH = [&](int pairIdx, const M12d &localH, const int idx[4]) {
    TripletD *base = triplets.data() + 144 * pairIdx;
    int k = 0;
    for (int i = 0; i < 4; ++i) {
      int ri = (idx[i] >= 0) ? 3 * idx[i] : 0;
      bool vi = (idx[i] >= 0);
      for (int j = 0; j < 4; ++j) {
        int cj = (idx[j] >= 0) ? 3 * idx[j] : 0;
        bool valid = vi && (idx[j] >= 0);
        for (int di = 0; di < 3; ++di)
          for (int dj = 0; dj < 3; ++dj, ++k)
            base[k] = TripletD(ri + di, cj + dj,
              valid ? localH(3 * i + di, 3 * j + dj) : 0.0);
      }
    }
  };

  double ee_eps = eps_ee;
  double dhat2 = dhat * dhat;

  // ---- PT pairs ----
  double ptEnergy = tbb::parallel_reduce(
    tbb::blocked_range<int>(0, nPT), 0.0,
    [&](const tbb::blocked_range<int> &range, double localE) {
      for (int i = range.begin(); i < range.end(); ++i) {
        auto &pair = ptPairs_[i];
        V3d p = vtx(x, pair.p);
        V3d t0 = vtx(x, pair.t0);
        V3d t1 = vtx(x, pair.t1);
        V3d t2 = vtx(x, pair.t2);

        double d2 = distance::computePTSqDist(p, t0, t1, t2);
        if (d2 >= dhat2 || d2 <= 0.0)
          continue;  // triplet slot already zeroed

        double wk = pair.weight * kappa;

        // energy
        localE += wk * barrier::b(d2, dhat2);

        // gradient
        V12d gd2 = distance::computePTSqDistGrad(p, t0, t1, t2);
        double gp = barrier::dbds(d2, dhat2);
        V12d gE = wk * gp * gd2;
        int idx[4] = { pair.p, pair.t0, pair.t1, pair.t2 };
        scatterG(gE, idx);

        // hessian
        M12d Hd2 = distance::computePTSqDistHess(p, t0, t1, t2);
        double gpp = barrier::d2bds2(d2, dhat2);
        M12d localH = wk * (gpp * gd2 * gd2.transpose() + gp * Hd2);
        localH = projectToPSD(localH);
        scatterH(i, localH, idx);
      }
      return localE;
    },
    std::plus<double>());

  // ---- EE pairs ----
  double eeEnergy = tbb::parallel_reduce(
    tbb::blocked_range<int>(0, nEE), 0.0,
    [&](const tbb::blocked_range<int> &range, double localE) {
      for (int i = range.begin(); i < range.end(); ++i) {
        auto &pair = eePairs_[i];
        V3d ea0 = vtx(x, pair.ea0);
        V3d ea1 = vtx(x, pair.ea1);
        V3d eb0 = vtx(x, pair.eb0);
        V3d eb1 = vtx(x, pair.eb1);

        double d2 = distance::computeEESqDist(ea0, ea1, eb0, eb1);
        if (d2 >= dhat2 || d2 <= 0.0)
          continue;  // triplet slot already zeroed

        double wk = pair.weight * kappa;
        double m = 1.0;
        V12d gm;
        gm.setZero();
        M12d Hm;
        Hm.setZero();
        if (ee_eps > 0.0) {
          m = distance::eeMollifier(ea0, ea1, eb0, eb1, ee_eps);
          gm = distance::eeMollifierGrad(ea0, ea1, eb0, eb1, ee_eps);
          Hm = distance::eeMollifierHess(ea0, ea1, eb0, eb1, ee_eps);
        }

        double bv = barrier::b(d2, dhat2);

        // energy
        localE += wk * m * bv;

        // gradient
        V12d gd2 = distance::computeEESqDistGrad(ea0, ea1, eb0, eb1);
        double dbv = barrier::dbds(d2, dhat2);
        V12d gE = wk * (gm * bv + m * dbv * gd2);
        int idx[4] = { pair.ea0, pair.ea1, pair.eb0, pair.eb1 };
        scatterG(gE, idx);

        // hessian
        M12d Hd2 = distance::computeEESqDistHess(ea0, ea1, eb0, eb1);
        double gp = barrier::dbds(d2, dhat2);
        double gpp = barrier::d2bds2(d2, dhat2);
        V12d gb = dbv * gd2;
        M12d localH;
        if (ee_eps > 0.0) {
          localH = wk * (bv * Hm + gm * gb.transpose() + gb * gm.transpose() + m * (gpp * gd2 * gd2.transpose() + gp * Hd2));
        }
        else {
          localH = wk * (gpp * gd2 * gd2.transpose() + gp * Hd2);
        }
        localH = projectToPSD(localH);
        scatterH(nPT + i, localH, idx);
      }
      return localE;
    },
    std::plus<double>());

  energy = ptEnergy + eeEnergy;

  // Floor penalty (energy + gradient + hessian): 0.5 * kappa * (z-h)^2 when z < h
  if (useFloor) {
    tbb::enumerable_thread_specific<double> tls_floorE(0.0);
    tbb::enumerable_thread_specific<std::vector<TripletD>> tls_floorTrip;
    tbb::parallel_for(
      tbb::blocked_range<int>(0, numVerts_),
      [&](const tbb::blocked_range<int> &range) {
        double &localE = tls_floorE.local();
        auto &localTrip = tls_floorTrip.local();
        double *gdata = grad.data();
        for (int vi = range.begin(); vi < range.end(); ++vi) {
          double dz = x[3 * vi + 2] - floorHeight;
          if (dz >= 0.0)
            continue;

          // energy
          localE += 0.5 * floorKappa * dz * dz;

          // gradient
          std::atomic_ref<double>(gdata[3 * vi + 2])
            .fetch_add(floorKappa * dz, std::memory_order_relaxed);

          // hessian
          int row = 3 * vi + 2;
          localTrip.emplace_back(row, row, floorKappa);
        }
      });

    for (auto &e : tls_floorE)
      energy += e;
    for (auto &local : tls_floorTrip)
      triplets.insert(triplets.end(), local.begin(), local.end());
  }

  hess.resize(n, n);
  hess.setFromTriplets(triplets.begin(), triplets.end());
}


}  // namespace CIPC
}  // namespace Contact
}  // namespace pgo
