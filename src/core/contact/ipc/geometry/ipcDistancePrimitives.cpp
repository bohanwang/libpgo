/*
copyright to Bohan Wang
*/

#include "ipcDistancePrimitives.h"

#include "ipc/geometry/generated/CIPC_autogen.h"
#include "ipc/geometry/generated/CIPC_autogen_ll.h"

#include <algorithm>
#include <cmath>

namespace pgo {
namespace Contact {
namespace IPC {
static constexpr double kEps = 1e-20;  // numerical guard

// =========================================================================
//  Distance type classification
// =========================================================================
namespace distance {

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
  pgo::Contact::IPC::autogen::point_line_distance_hessian_3D(
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
  pgo::Contact::IPC::autogen::point_plane_distance_gradient(
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
  pgo::Contact::IPC::autogen::point_plane_distance_hessian(
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
  pgo::Contact::IPC::autogen::line_line_distance_gradient(
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
  pgo::Contact::IPC::autogen::line_line_distance_hessian(
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
  pgo::Contact::IPC::autogen::edge_edge_cross_squarednorm_gradient(
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
  pgo::Contact::IPC::autogen::edge_edge_cross_squarednorm_gradient(
    ea0[0], ea0[1], ea0[2], ea1[0], ea1[1], ea1[2],
    eb0[0], eb0[1], eb0[2], eb1[0], eb1[1], eb1[2],
    gx_data);
  V12d gx = Eigen::Map<V12d>(gx_data);

  double Hx_data[144];
  pgo::Contact::IPC::autogen::edge_edge_cross_squarednorm_hessian(
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

// ----- Combined dispatchers (classify once) -----

PTDistAll computePTSqDistAll(const V3d &p, const V3d &t0,
  const V3d &t1, const V3d &t2)
{
  PTDistAll result;
  auto tp = classifyPT(p, t0, t1, t2);
  switch (tp) {
  case PTDistType::PP_PT0: {
    V6d g6 = ppSqDistGrad(p, t0);
    M6d H6 = ppSqDistHess(p, t0);
    result.d2 = ppSqDist(p, t0);
    embedPP(0, 1, g6, H6, result.grad, result.hess);
    break;
  }
  case PTDistType::PP_PT1: {
    V6d g6 = ppSqDistGrad(p, t1);
    M6d H6 = ppSqDistHess(p, t1);
    result.d2 = ppSqDist(p, t1);
    embedPP(0, 2, g6, H6, result.grad, result.hess);
    break;
  }
  case PTDistType::PP_PT2: {
    V6d g6 = ppSqDistGrad(p, t2);
    M6d H6 = ppSqDistHess(p, t2);
    result.d2 = ppSqDist(p, t2);
    embedPP(0, 3, g6, H6, result.grad, result.hess);
    break;
  }
  case PTDistType::PE_PT0T1: {
    V9d g9 = peSqDistGrad(p, t0, t1);
    M9d H9 = peSqDistHess(p, t0, t1);
    result.d2 = peSqDist(p, t0, t1);
    int slots[3] = { 0, 1, 2 };
    embedPE(slots, g9, H9, result.grad, result.hess);
    break;
  }
  case PTDistType::PE_PT1T2: {
    V9d g9 = peSqDistGrad(p, t1, t2);
    M9d H9 = peSqDistHess(p, t1, t2);
    result.d2 = peSqDist(p, t1, t2);
    int slots[3] = { 0, 2, 3 };
    embedPE(slots, g9, H9, result.grad, result.hess);
    break;
  }
  case PTDistType::PE_PT2T0: {
    V9d g9 = peSqDistGrad(p, t2, t0);
    M9d H9 = peSqDistHess(p, t2, t0);
    result.d2 = peSqDist(p, t2, t0);
    int slots[3] = { 0, 3, 1 };
    embedPE(slots, g9, H9, result.grad, result.hess);
    break;
  }
  case PTDistType::PT:
    result.d2 = ptSqDist(p, t0, t1, t2);
    result.grad = ptSqDistGrad(p, t0, t1, t2);
    result.hess = ptSqDistHess(p, t0, t1, t2);
    break;
  }
  return result;
}

EEDistAll computeEESqDistAll(const V3d &ea0, const V3d &ea1,
  const V3d &eb0, const V3d &eb1)
{
  EEDistAll result;
  auto tp = classifyEE(ea0, ea1, eb0, eb1);
  switch (tp) {
  case EEDistType::PP_Ea0Eb0: {
    V6d g6 = ppSqDistGrad(ea0, eb0);
    M6d H6 = ppSqDistHess(ea0, eb0);
    result.d2 = ppSqDist(ea0, eb0);
    embedPP(0, 2, g6, H6, result.grad, result.hess);
    break;
  }
  case EEDistType::PP_Ea0Eb1: {
    V6d g6 = ppSqDistGrad(ea0, eb1);
    M6d H6 = ppSqDistHess(ea0, eb1);
    result.d2 = ppSqDist(ea0, eb1);
    embedPP(0, 3, g6, H6, result.grad, result.hess);
    break;
  }
  case EEDistType::PP_Ea1Eb0: {
    V6d g6 = ppSqDistGrad(ea1, eb0);
    M6d H6 = ppSqDistHess(ea1, eb0);
    result.d2 = ppSqDist(ea1, eb0);
    embedPP(1, 2, g6, H6, result.grad, result.hess);
    break;
  }
  case EEDistType::PP_Ea1Eb1: {
    V6d g6 = ppSqDistGrad(ea1, eb1);
    M6d H6 = ppSqDistHess(ea1, eb1);
    result.d2 = ppSqDist(ea1, eb1);
    embedPP(1, 3, g6, H6, result.grad, result.hess);
    break;
  }
  case EEDistType::PE_Ea0_Eb: {
    V9d g9 = peSqDistGrad(ea0, eb0, eb1);
    M9d H9 = peSqDistHess(ea0, eb0, eb1);
    result.d2 = peSqDist(ea0, eb0, eb1);
    int s[3] = { 0, 2, 3 };
    embedPE(s, g9, H9, result.grad, result.hess);
    break;
  }
  case EEDistType::PE_Ea1_Eb: {
    V9d g9 = peSqDistGrad(ea1, eb0, eb1);
    M9d H9 = peSqDistHess(ea1, eb0, eb1);
    result.d2 = peSqDist(ea1, eb0, eb1);
    int s[3] = { 1, 2, 3 };
    embedPE(s, g9, H9, result.grad, result.hess);
    break;
  }
  case EEDistType::PE_Eb0_Ea: {
    V9d g9 = peSqDistGrad(eb0, ea0, ea1);
    M9d H9 = peSqDistHess(eb0, ea0, ea1);
    result.d2 = peSqDist(eb0, ea0, ea1);
    int s[3] = { 2, 0, 1 };
    embedPE(s, g9, H9, result.grad, result.hess);
    break;
  }
  case EEDistType::PE_Eb1_Ea: {
    V9d g9 = peSqDistGrad(eb1, ea0, ea1);
    M9d H9 = peSqDistHess(eb1, ea0, ea1);
    result.d2 = peSqDist(eb1, ea0, ea1);
    int s[3] = { 3, 0, 1 };
    embedPE(s, g9, H9, result.grad, result.hess);
    break;
  }
  case EEDistType::EE:
    result.d2 = eeSqDist(ea0, ea1, eb0, eb1);
    result.grad = eeSqDistGrad(ea0, ea1, eb0, eb1);
    result.hess = eeSqDistHess(ea0, ea1, eb0, eb1);
    break;
  }
  return result;
}

}  // namespace distance

}  // namespace IPC
}  // namespace Contact
}  // namespace pgo
