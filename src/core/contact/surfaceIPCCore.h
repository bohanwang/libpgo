/*
copyright to Bohan Wang
*/

#pragma once

#include "EigenDef.h"
#include "potentialEnergy.h"

#include <vector>
#include <array>
#include <cmath>
#include <algorithm>
#include <set>
#include <unordered_set>
#include <unordered_map>

namespace pgo
{
namespace Contact
{
namespace CIPC
{

using namespace pgo::EigenSupport;

// -------------------------------------------------------------------------
// Distance-type classification
// -------------------------------------------------------------------------
// Point-Triangle: which feature of the triangle is closest to the point?
enum class PTDistType
{
  PP_PT0,    // point  <-> triangle vertex 0
  PP_PT1,    // point  <-> triangle vertex 1
  PP_PT2,    // point  <-> triangle vertex 2
  PE_PT0T1,  // point  <-> edge t0-t1
  PE_PT1T2,  // point  <-> edge t1-t2
  PE_PT2T0,  // point  <-> edge t2-t0
  PT         // point  <-> triangle interior (face)
};

// Edge-Edge: which sub-features are closest?
enum class EEDistType
{
  PP_Ea0Eb0,
  PP_Ea0Eb1,
  PP_Ea1Eb0,
  PP_Ea1Eb1,
  PE_Ea0_Eb,  // vertex ea0 to edge eb
  PE_Ea1_Eb,  // vertex ea1 to edge eb
  PE_Eb0_Ea,  // vertex eb0 to edge ea
  PE_Eb1_Ea,  // vertex eb1 to edge ea
  EE          // interior edge-edge
};

// -------------------------------------------------------------------------
// Collision pair descriptors
// -------------------------------------------------------------------------
struct PTPair
{
  int p, t0, t1, t2;  // vertex indices
  double weight;      // area(point) * area(triangle)
};

struct EEPair
{
  int ea0, ea1, eb0, eb1;  // vertex indices
  double weight;           // length(edgeA) * length(edgeB)
};

// =========================================================================
// Distance primitives  (value, gradient, hessian of SQUARED distance)
// =========================================================================
namespace distance
{

// --- classification -------------------------------------------------------
PTDistType classifyPT(const V3d &p, const V3d &t0,
  const V3d &t1, const V3d &t2);

EEDistType classifyEE(const V3d &ea0, const V3d &ea1,
  const V3d &eb0, const V3d &eb1);

// --- Point-Point (PP)  stacked as x = [a; b]  (6 DOF) -------------------
double ppSqDist(const V3d &a, const V3d &b);
V6d ppSqDistGrad(const V3d &a, const V3d &b);
M6d ppSqDistHess(const V3d &a, const V3d &b);

// --- Point-Edge  (PE)  stacked as x = [p; e0; e1]  (9 DOF) -------------
double peSqDist(const V3d &p, const V3d &e0, const V3d &e1);
V9d peSqDistGrad(const V3d &p, const V3d &e0, const V3d &e1);
M9d peSqDistHess(const V3d &p, const V3d &e0, const V3d &e1);

// --- Point-Triangle face  x = [p; t0; t1; t2]  (12 DOF) ----------------
double ptSqDist(const V3d &p, const V3d &t0,
  const V3d &t1, const V3d &t2);
V12d ptSqDistGrad(const V3d &p, const V3d &t0,
  const V3d &t1, const V3d &t2);
M12d ptSqDistHess(const V3d &p, const V3d &t0,
  const V3d &t1, const V3d &t2);

// --- Edge-Edge interior  x = [ea0; ea1; eb0; eb1]  (12 DOF) -------------
double eeSqDist(const V3d &ea0, const V3d &ea1,
  const V3d &eb0, const V3d &eb1);
V12d eeSqDistGrad(const V3d &ea0, const V3d &ea1,
  const V3d &eb0, const V3d &eb1);
M12d eeSqDistHess(const V3d &ea0, const V3d &ea1,
  const V3d &eb0, const V3d &eb1);

// --- Mollified EE (for near-parallel edges in codim-IPC) -----------------
double eeMollifier(const V3d &ea0, const V3d &ea1,
  const V3d &eb0, const V3d &eb1,
  double eps_x);
V12d eeMollifierGrad(const V3d &ea0, const V3d &ea1,
  const V3d &eb0, const V3d &eb1,
  double eps_x);
M12d eeMollifierHess(const V3d &ea0, const V3d &ea1,
  const V3d &eb0, const V3d &eb1,
  double eps_x);

// --- Unified dispatchers (classify then call sub-routine) -----------------
//   They write into 12-dim vectors/matrices with the canonical ordering:
//   PT:  [p; t0; t1; t2]      EE: [ea0; ea1; eb0; eb1]
double computePTSqDist(const V3d &p, const V3d &t0,
  const V3d &t1, const V3d &t2);
V12d computePTSqDistGrad(const V3d &p, const V3d &t0,
  const V3d &t1, const V3d &t2);
M12d computePTSqDistHess(const V3d &p, const V3d &t0,
  const V3d &t1, const V3d &t2);

double computeEESqDist(const V3d &ea0, const V3d &ea1,
  const V3d &eb0, const V3d &eb1);
V12d computeEESqDistGrad(const V3d &ea0, const V3d &ea1,
  const V3d &eb0, const V3d &eb1);
M12d computeEESqDistHess(const V3d &ea0, const V3d &ea1,
  const V3d &eb0, const V3d &eb1);

}  // namespace distance

// =========================================================================
// Barrier function   b(d, dhat)
// =========================================================================
namespace barrier
{

// Barrier on squared distance (Codim-IPC elastic formulation):
//   b(s, shat) = -(s/shat - 1)^2 * ln(s/shat)   if 0 < s < shat,  else 0
// where s = d^2, shat = dhat^2
double b(double s, double shat);
double dbds(double s, double shat);    // first derivative  w.r.t. s
double d2bds2(double s, double shat);  // second derivative w.r.t. s

}  // namespace barrier

// =========================================================================
// Continuous Collision Detection (CCD)
// =========================================================================
namespace ccd
{

// Returns the earliest time of impact in (0, tMax].
// If no collision, returns tMax.
double pointTriangleCCD(const V3d &p, const V3d &t0,
  const V3d &t1, const V3d &t2,
  const V3d &dp, const V3d &dt0,
  const V3d &dt1, const V3d &dt2,
  double thickness = 0.0,
  double tMax = 1.0);

double edgeEdgeCCD(const V3d &ea0, const V3d &ea1,
  const V3d &eb0, const V3d &eb1,
  const V3d &dea0, const V3d &dea1,
  const V3d &deb0, const V3d &deb1,
  double thickness = 0.0,
  double tMax = 1.0);

}  // namespace ccd

// =========================================================================
// Utility: project a symmetric matrix to PSD  (clamp negative eigenvalues)
// =========================================================================
M12d projectToPSD(const M12d &H);

// =========================================================================
//  Surface-space IPC core
// =========================================================================
class SurfaceIPCCore
{
public:
  struct Parameters
  {
    double dhat = 1e-1;
    double kappa = 0.1;
    double eps_ee = 0.0;
    double slackness = 1.0;
  };

  SurfaceIPCCore() = default;
  explicit SurfaceIPCCore(const Parameters &params) { setParameters(params); }

  void setParameters(const Parameters &params);
  Parameters getParameters() const;

  void setMesh(const MXd &V, const MXi &F);

  double computeEnergy(EigenSupport::ConstRefVecXd x_surf) const;
  void computeGradient(EigenSupport::ConstRefVecXd x_surf, EigenSupport::RefVecXd g_surf) const;
  void computeHessian(EigenSupport::ConstRefVecXd x_surf, EigenSupport::SpMatD &H_surf) const;
  void computeAll(EigenSupport::ConstRefVecXd x_surf, double &energy, VXd &g_surf, SpMatD &H_surf) const;
  double computeMaxStepSize(EigenSupport::ConstRefVecXd x_surf, EigenSupport::ConstRefVecXd dx_surf) const;

  const std::vector<PTPair> &getPTPairs() const { return ptPairs_; }
  const std::vector<EEPair> &getEEPairs() const { return eePairs_; }

  int getNumSurfaceVertices() const { return numVerts_; }
  int getNumSurfaceDOFs() const { return 3 * numVerts_; }

private:
  void findCollisionPairs(const VXd &positions) const;
  void buildEdges();
  void buildAdjacency();
  void buildAreaWeights(const MXd &V);

  static V3d vtx(const VXd &x, int i)
  {
    return x.segment<3>(3 * i);
  }

  double dhat = 1e-1;
  double kappa = 0.1;
  double eps_ee = 0.0;
  double slackness = 1.0;
  int numVerts_ = 0;
  std::vector<std::array<int, 3>> triangles_;
  std::vector<std::array<int, 2>> edges_;
  std::vector<int> allDOFs_;
  std::unordered_set<long long> vertexTriAdj_;
  std::unordered_set<long long> edgeVertAdj_;
  std::vector<double> vertexArea_;
  std::vector<double> triArea_;
  std::vector<double> edgeLength_;
  mutable std::vector<PTPair> ptPairs_;
  mutable std::vector<EEPair> eePairs_;
};

}  // namespace CIPC
}  // namespace Contact
}  // namespace pgo
