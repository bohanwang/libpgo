/*
author: Bohan Wang
copyright to USC
*/

#include "CCDKernel.h"

#include "predicates.h"
#include "triangle.h"

#include "safeCCDWrapper.h"
#include "exactCCDWrapper.h"

#include <iostream>
#include <algorithm>

using namespace pgo::CCDKernel;
using namespace pgo::CCDKernel::TriangleCCD;

namespace ES = pgo::EigenSupport;

namespace pgo
{
namespace CCDKernel
{
namespace TriangleCCD
{
bool CCDTestVega(int numVertices, const double *pos_t0, const double *pos_t1, int numTriangles, const ES::V3i *triangles, int triA, int triB, CCDData *ccdDataPtr = nullptr);
bool CCDTest3rdParty(int numVertices, const double *pos_t0, const double *pos_t1, int numTriangles, const ES::V3i *triangles, int triA, int triB, CCDData *ccdDataPtr = nullptr);
}  // namespace TriangleCCD
}  // namespace CCDKernel
}  // namespace pgo

bool pgo::CCDKernel::TriangleCCD::CCDTest(int numVertices, const double *pos_t0, const double *pos_t1, int numTriangles, const ES::V3i *triangles,
  int triA, int triB, CCD_METHOD m, CCDData *ccdDataPtr)
{
  if (m == CCDM_VEGA) {
    return pgo::CCDKernel::TriangleCCD::CCDTestVega(numVertices, pos_t0, pos_t1, numTriangles, triangles, triA, triB, ccdDataPtr);
  }
  else if (m == CCDM_3RDPARTY) {
    return pgo::CCDKernel::TriangleCCD::CCDTest3rdParty(numVertices, pos_t0, pos_t1, numTriangles, triangles, triA, triB, ccdDataPtr);
  }
  else {
    std::cerr << "Unknown method: " << m << std::endl;
    return false;
  }
}

bool pgo::CCDKernel::TriangleCCD::computeCCDPlane(int numVertices, const double *pos_t0, const double *pos_t1, int numTriangles, const ES::V3i *triangles,
  const CCDData &ccdData, bool doSign, double planePoint[3], double planeNormal[3], double &signA, double &signB)
{
  constexpr double eps = 1e-10;

  const double *vtxA0[3] = {
    pos_t0 + triangles[ccdData.triA][0] * 3,
    pos_t0 + triangles[ccdData.triA][1] * 3,
    pos_t0 + triangles[ccdData.triA][2] * 3
  };

  const double *vtxA1[3] = {
    pos_t1 + triangles[ccdData.triA][0] * 3,
    pos_t1 + triangles[ccdData.triA][1] * 3,
    pos_t1 + triangles[ccdData.triA][2] * 3
  };

  const double *vtxB0[3] = {
    pos_t0 + triangles[ccdData.triB][0] * 3,
    pos_t0 + triangles[ccdData.triB][1] * 3,
    pos_t0 + triangles[ccdData.triB][2] * 3
  };

  const double *vtxB1[3] = {
    pos_t1 + triangles[ccdData.triB][0] * 3,
    pos_t1 + triangles[ccdData.triB][1] * 3,
    pos_t1 + triangles[ccdData.triB][2] * 3
  };

  ES::V3d contactPoint;

  if (ccdData.ccdCase == CCDC_EDGEA_EDGEB) {
    int vid0 = ccdData.eid[0];
    int vid1 = (ccdData.eid[0] + 1) % 3;

    int vid2 = ccdData.eid[1];
    int vid3 = (ccdData.eid[1] + 1) % 3;

    ES::V3d edgeA0_t0 = asVec3d(vtxA0[vid0]);
    ES::V3d edgeA1_t0 = asVec3d(vtxA0[vid1]);
    ES::V3d edgeA0_t1 = asVec3d(vtxA1[vid0]);
    ES::V3d edgeA1_t1 = asVec3d(vtxA1[vid1]);
    ES::V3d edgeA0_t = edgeA0_t0 * (1 - ccdData.t) + edgeA0_t1 * ccdData.t;
    ES::V3d edgeA1_t = edgeA1_t0 * (1 - ccdData.t) + edgeA1_t1 * ccdData.t;

    ES::V3d edgeB0_t0 = asVec3d(vtxB0[vid2]);
    ES::V3d edgeB1_t0 = asVec3d(vtxB0[vid3]);
    ES::V3d edgeB0_t1 = asVec3d(vtxB1[vid2]);
    ES::V3d edgeB1_t1 = asVec3d(vtxB1[vid3]);
    ES::V3d edgeB0_t = edgeB0_t0 * (1 - ccdData.t) + edgeB0_t1 * ccdData.t;
    ES::V3d edgeB1_t = edgeB1_t0 * (1 - ccdData.t) + edgeB1_t1 * ccdData.t;

    ES::V3d ea = edgeA1_t - edgeA0_t;
    ES::V3d eb = edgeB1_t - edgeB0_t;

    ES::V3d pn(1, 0, 0);
    // it these two edge are in parallel,
    // we ignore it, because it should be captured by the vtx-triangle case
    if (std::abs(ea.dot(eb)) > 1.0 - eps) {
      return false;
    }
    else {
      pn = ea.cross(eb);
    }
    pn.normalize();
    (ES::Mp<ES::V3d>(planeNormal)) = pn;
    // edgeA0_t.convertToArray(planePoint);

    contactPoint = edgeA0_t * ccdData.weights[0] + edgeA1_t * ccdData.weights[1];
    (ES::Mp<ES::V3d>(planePoint)) = contactPoint;
  }
  else if (ccdData.ccdCase == CCDC_VTXA_TRIB) {
    ES::V3d vtxA_t0 = asVec3d(vtxA0[ccdData.vid]);
    ES::V3d vtxA_t1 = asVec3d(vtxA1[ccdData.vid]);
    ES::V3d vtxA_t = vtxA_t0 * (1 - ccdData.t) + vtxA_t1 * ccdData.t;

    ES::V3d triB0_t0 = asVec3d(vtxB0[0]);
    ES::V3d triB1_t0 = asVec3d(vtxB0[1]);
    ES::V3d triB2_t0 = asVec3d(vtxB0[2]);
    ES::V3d triB0_t1 = asVec3d(vtxB1[0]);
    ES::V3d triB1_t1 = asVec3d(vtxB1[1]);
    ES::V3d triB2_t1 = asVec3d(vtxB1[2]);
    ES::V3d triB0_t = triB0_t0 * (1 - ccdData.t) + triB0_t1 * ccdData.t;
    ES::V3d triB1_t = triB1_t0 * (1 - ccdData.t) + triB1_t1 * ccdData.t;
    ES::V3d triB2_t = triB2_t0 * (1 - ccdData.t) + triB2_t1 * ccdData.t;

    ES::V3d ea = triB1_t - triB0_t;
    ES::V3d eb = triB2_t - triB0_t;
    ES::V3d pn(1, 0, 0);

    // if we encounter a degenerated triangle,
    // the contact should be captured by the neighboring triangles.
    if (std::abs(ea.dot(eb)) > 1 - eps) {
      return false;
    }
    else
      pn = ea.cross(eb);

    pn.normalize();
    (ES::Mp<ES::V3d>(planeNormal)) = pn;
    // triB0_t.convertToArray(planePoint);

    contactPoint = vtxA_t;
    (ES::Mp<ES::V3d>(planePoint)) = contactPoint;
  }
  else if (ccdData.ccdCase == CCDC_VTXB_TRIA) {
    ES::V3d vtxB_t0 = asVec3d(vtxB0[ccdData.vid]);
    ES::V3d vtxB_t1 = asVec3d(vtxB1[ccdData.vid]);
    ES::V3d vtxB_t = vtxB_t0 * (1 - ccdData.t) + vtxB_t1 * ccdData.t;

    ES::V3d triA0_t0 = asVec3d(vtxA0[0]);
    ES::V3d triA1_t0 = asVec3d(vtxA0[1]);
    ES::V3d triA2_t0 = asVec3d(vtxA0[2]);
    ES::V3d triA0_t1 = asVec3d(vtxA1[0]);
    ES::V3d triA1_t1 = asVec3d(vtxA1[1]);
    ES::V3d triA2_t1 = asVec3d(vtxA1[2]);
    ES::V3d triA0_t = triA0_t0 * (1 - ccdData.t) + triA0_t1 * ccdData.t;
    ES::V3d triA1_t = triA1_t0 * (1 - ccdData.t) + triA1_t1 * ccdData.t;
    ES::V3d triA2_t = triA2_t0 * (1 - ccdData.t) + triA2_t1 * ccdData.t;

    ES::V3d ea = triA1_t - triA0_t;
    ES::V3d eb = triA2_t - triA0_t;
    ES::V3d pn(1, 0, 0);

    // if we encounter a degenerated triangle,
    // the contact should be captured by the neighboring triangles.
    if (std::abs(ea.dot(eb)) > 1 - eps) {
      return false;
    }
    else
      pn = ea.cross(eb);
    pn.normalize();
    (ES::Mp<ES::V3d>(planeNormal)) = pn;
    // triA0_t.convertToArray(planePoint);

    contactPoint = vtxB_t;
    (ES::Mp<ES::V3d>(planePoint)) = contactPoint;
  }
  else if (ccdData.ccdCase == CCDC_COLLIDED) {
    return false;
  }
  else {
    std::cerr << "Un-addressed case" << std::endl;
    abort();
  }

  if (doSign == false)
    return true;

  computeSign(numVertices, pos_t0, pos_t1, numTriangles, triangles, ccdData, planePoint, planeNormal, signA, signB);

  if (signA * signB > 0) {
    ES::V3d triA0_t0 = asVec3d(vtxA0[0]);
    ES::V3d triA1_t0 = asVec3d(vtxA0[1]);
    ES::V3d triA2_t0 = asVec3d(vtxA0[2]);
    ES::V3d triA0_t1 = asVec3d(vtxA1[0]);
    ES::V3d triA1_t1 = asVec3d(vtxA1[1]);
    ES::V3d triA2_t1 = asVec3d(vtxA1[2]);
    ES::V3d triA0_t = triA0_t0 * (1 - ccdData.t) + triA0_t1 * ccdData.t;
    ES::V3d triA1_t = triA1_t0 * (1 - ccdData.t) + triA1_t1 * ccdData.t;
    ES::V3d triA2_t = triA2_t0 * (1 - ccdData.t) + triA2_t1 * ccdData.t;

    ES::V3d triB0_t0 = asVec3d(vtxB0[0]);
    ES::V3d triB1_t0 = asVec3d(vtxB0[1]);
    ES::V3d triB2_t0 = asVec3d(vtxB0[2]);
    ES::V3d triB0_t1 = asVec3d(vtxB1[0]);
    ES::V3d triB1_t1 = asVec3d(vtxB1[1]);
    ES::V3d triB2_t1 = asVec3d(vtxB1[2]);
    ES::V3d triB0_t = triB0_t0 * (1 - ccdData.t) + triB0_t1 * ccdData.t;
    ES::V3d triB1_t = triB1_t0 * (1 - ccdData.t) + triB1_t1 * ccdData.t;
    ES::V3d triB2_t = triB2_t0 * (1 - ccdData.t) + triB2_t1 * ccdData.t;

    ES::V3d eA0 = triA1_t - triA0_t;
    ES::V3d eA1 = triA2_t - triA0_t;

    ES::V3d eB0 = triB1_t - triB0_t;
    ES::V3d eB1 = triB2_t - triB0_t;

    ES::V3d nA = eA0.cross(eA1);
    nA.normalize();

    ES::V3d nB = eB0.cross(eB1);
    nB.normalize();

#if 0
    ES::V3d pn;
    if (dot(nA, nB) > 0) {
      // std::cerr << "Two faces are collided but has the same orientation!" << std::endl;
      pn = (nA - nB) * 0.5;
    }
    else {
      pn = (nA - nB) * 0.5;
    }
    pn.normalize();

    pn.convertToArray(planeNormal);
    contactPoint.convertToArray(planePoint);

    signA = -1;
    signB = 1;
#else
    ES::V3d eA12 = triA2_t - triA1_t;
    ES::V3d eB12 = triB2_t - triB1_t;

    ES::V3d axes[] = {
      nA,
      nB,
      (nA - nB) * 0.5,
      eA0.cross(eB0),
      eA0.cross(eB1),
      eA0.cross(eB12),

      eA1.cross(eB0),
      eA1.cross(eB1),
      eA1.cross(eB12),

      eA12.cross(eB0),
      eA12.cross(eB1),
      eA12.cross(eB12)
    };

    size_t numAxes = sizeof(axes) / sizeof(ES::V3d);
    size_t axisId = 0;
    double minOverlap = 1e100;
    double minProj = 1e100;

    for (size_t ai = 0; ai < numAxes; ai++) {
      axes[ai].normalize();

      double projA[3] = {
        triA0_t.dot(axes[ai]),
        triA1_t.dot(axes[ai]),
        triA2_t.dot(axes[ai])
      };

      double projB[3] = {
        triB0_t.dot(axes[ai]),
        triB1_t.dot(axes[ai]),
        triB2_t.dot(axes[ai])
      };

      double minA = std::min({ projA[0], projA[1], projA[2] });
      double maxA = std::max({ projA[0], projA[1], projA[2] });

      double minB = std::min({ projB[0], projB[1], projB[2] });
      double maxB = std::max({ projB[0], projB[1], projB[2] });

      double overlap = 0;
      if (minB > maxA || minA > maxB)
        overlap = -1;
      else {
        double minO = std::max(minA, minB);
        double maxO = std::min(maxA, maxB);
        overlap = std::abs(maxO - minO);
      }

      // find the axis with the smallest overlap,
      // -1 means no overlap
      // in this case, we find the one with the minimal projection area.
      if (overlap < minOverlap) {
        minOverlap = overlap;
        axisId = ai;
      }
      else if (overlap == minOverlap) {
        double proj = (maxA - minA) + (maxB - minB);
        if (proj < minProj) {
          axisId = ai;
          minProj = proj;
        }
      }
    }

    // after we found the axis,
    // we check the sign
    double dir = (triA0_t0 - contactPoint).dot(axes[axisId]);
    ES::V3d pn = axes[axisId];
    if (dir > 0)
      pn *= -1;

    (ES::Mp<ES::V3d>(planeNormal)) = pn;
    (ES::Mp<ES::V3d>(planePoint)) = contactPoint;

    signA = -1;
    signB = 1;
#endif
  }

  return true;
}

bool pgo::CCDKernel::TriangleCCD::computeDCDPlane(int, const double *pos_t0, int, const ES::V3i *triangles,
  CCDData &ccdData, double planePoint[3], double planeNormal[3])
{
  // constexpr double eps = 1e-10;

  const double *vtxA0[3] = {
    pos_t0 + triangles[ccdData.triA][0] * 3,
    pos_t0 + triangles[ccdData.triA][1] * 3,
    pos_t0 + triangles[ccdData.triA][2] * 3
  };

  const double *vtxB0[3] = {
    pos_t0 + triangles[ccdData.triB][0] * 3,
    pos_t0 + triangles[ccdData.triB][1] * 3,
    pos_t0 + triangles[ccdData.triB][2] * 3
  };

  ES::V3d vtxA[3] = { asVec3d(vtxA0[0]), asVec3d(vtxA0[1]), asVec3d(vtxA0[2]) };
  ES::V3d vtxB[3] = { asVec3d(vtxB0[0]), asVec3d(vtxB0[1]), asVec3d(vtxB0[2]) };

  ES::V3d nA = (vtxA[1] - vtxA[0]).cross(vtxA[2] - vtxA[0]);
  if (nA.dot(nA) < 1e-40) {
    return false;
  }
  nA.normalize();

  ES::V3d nB = (vtxB[1] - vtxB[0]).cross(vtxB[2] - vtxB[0]);
  if (nB.dot(nB) < 1e-40) {
    return false;
  }
  nB.normalize();

  Mesh::TriangleWithCollisionInfo infoA(vtxA[0], vtxA[1], vtxA[2]), infoB(vtxB[0], vtxB[1], vtxB[2]);
  double closestDistSq = 1e100;

  ccdData.t = 0;

  // vtx A vs tri B
  for (int vi = 0; vi < 3; vi++) {
    // we only check the penetrated vertices
    if ((vtxA[vi] - vtxB[0]).dot(nB) > 0)
      continue;

    int feature;
    double w[3];
    double distSq = infoB.distanceToPoint2(vtxA[vi], &feature, w, w + 1, w + 2);
    if (distSq < closestDistSq) {
      ccdData.vid = vi;
      ccdData.ccdCase = CCDKernel::TriangleCCD::CCDC_VTXA_TRIB;
      ccdData.weights[0] = w[0];
      ccdData.weights[1] = w[1];
      ccdData.weights[2] = w[2];

      closestDistSq = distSq;
    }
  }

  // vtx B vs tri A
  for (int vi = 0; vi < 3; vi++) {
    // we only check the penetrated vertices
    if ((vtxB[vi] - vtxA[0]).dot(nA) > 0)
      continue;

    int feature;
    double w[3];
    double distSq = infoA.distanceToPoint2(vtxB[vi], &feature, w, w + 1, w + 2);
    if (distSq < closestDistSq) {
      ccdData.vid = vi;
      ccdData.ccdCase = CCDKernel::TriangleCCD::CCDC_VTXB_TRIA;
      ccdData.weights[0] = w[0];
      ccdData.weights[1] = w[1];
      ccdData.weights[2] = w[2];

      closestDistSq = distSq;
    }
  }

  // edge vs edge
  // for (int ei = 0; ei < 3; ei++) {
  //  ES::V3d e0(vtxA[ei]), e1(vtxA[(ei + 1) % 3]);

  //  for (int ej = 0; ej < 3; ej++) {
  //    ES::V3d e2(vtxB[ej]), e3(vtxB[(ej + 1) % 3]);

  //    double distSq()

  //  }
  //}

  // std::cout << "min dist: " << sqrt(closestDistSq) << std::endl;
  if (closestDistSq > 1.e-3 * 1.e-3) {
    std::cout << "Detect a bad min dist: " << sqrt(closestDistSq) << std::endl;
    return false;
  }

  if (ccdData.ccdCase == CCDC_VTXA_TRIB) {
    (ES::Mp<ES::V3d>(planeNormal)) = nB;

    ES::V3d contactPoint = vtxB[0] * ccdData.weights[0] + vtxB[1] * ccdData.weights[1] + vtxB[2] * ccdData.weights[2];
    (ES::Mp<ES::V3d>(planePoint)) = contactPoint;
  }
  else if (ccdData.ccdCase == CCDC_VTXB_TRIA) {
    (ES::Mp<ES::V3d>(planeNormal)) = nA;
    ES::V3d contactPoint = vtxA[0] * ccdData.weights[0] + vtxA[1] * ccdData.weights[1] + vtxA[2] * ccdData.weights[2];
    (ES::Mp<ES::V3d>(planePoint)) = contactPoint;
  }

  return true;
}

void pgo::CCDKernel::TriangleCCD::computeSign(int, const double *pos_t0, const double *pos_t1, int, const ES::V3i *triangles,
  const CCDData &ccdData, const double planePoint[3], const double planeNormal[3], double &signA, double &signB)
{
  const double *vtxA0[3] = {
    pos_t0 + triangles[ccdData.triA][0] * 3,
    pos_t0 + triangles[ccdData.triA][1] * 3,
    pos_t0 + triangles[ccdData.triA][2] * 3
  };

  (void)pos_t1;
  // const double *vtxA1[3] = {
  //   pos_t1 + triangles[ccdData.triA][0] * 3,
  //   pos_t1 + triangles[ccdData.triA][1] * 3,
  //   pos_t1 + triangles[ccdData.triA][2] * 3
  // };

  const double *vtxB0[3] = {
    pos_t0 + triangles[ccdData.triB][0] * 3,
    pos_t0 + triangles[ccdData.triB][1] * 3,
    pos_t0 + triangles[ccdData.triB][2] * 3
  };

  // const double *vtxB1[3] = {
  //   pos_t1 + triangles[ccdData.triB][0] * 3,
  //   pos_t1 + triangles[ccdData.triB][1] * 3,
  //   pos_t1 + triangles[ccdData.triB][2] * 3
  // };

  ES::V3d p = asVec3d(planePoint), pn = asVec3d(planeNormal);

  if (ccdData.ccdCase == CCDC_EDGEA_EDGEB) {
    int vid0 = (ccdData.eid[0] + 2) % 3;
    int vid1 = (ccdData.eid[1] + 2) % 3;

    ES::V3d otherA_t0 = asVec3d(vtxA0[vid0]);
    // ES::V3d otherA_t1(vtxA1[vid0]);
    // ES::V3d otherA_t = otherA_t0 * (1 - ccdData.t) + otherA_t1 * ccdData.t;
    // signA = dot(pn, otherA_t - p);
    signA = pn.dot(otherA_t0 - p);

    ES::V3d otherB_t0 = asVec3d(vtxB0[vid1]);
    // ES::V3d otherB_t1(vtxB1[vid1]);
    // ES::V3d otherB_t = otherB_t0 * (1 - ccdData.t) + otherB_t1 * ccdData.t;
    // signB = dot(pn, otherB_t - p);
    signB = pn.dot(otherB_t0 - p);

    // ALOG(signA * signB <= 0);
  }
  else if (ccdData.ccdCase == CCDC_VTXA_TRIB) {
    ES::V3d vtxA_t0 = asVec3d(vtxA0[(ccdData.vid + 1) % 3]);
    // ES::V3d vtxA_t1(vtxA1[(ccdData.vid + 1) % 3]);
    // ES::V3d vtxA_t = vtxA_t0 * (1 - ccdData.t) + vtxA_t1 * ccdData.t;

    ES::V3d vtxB_t0 = asVec3d(vtxB0[0]);

    // signA = dot(pn, vtxA_t - p);
    signA = pn.dot(vtxA_t0 - p);
    signB = pn.dot(vtxB_t0 - p);
    // signB = -signA;
  }
  else if (ccdData.ccdCase == CCDC_VTXB_TRIA) {
    ES::V3d vtxB_t0 = asVec3d(vtxB0[(ccdData.vid + 1) % 3]);
    // ES::V3d vtxB_t1(vtxB1[(ccdData.vid + 1) % 3]);
    // ES::V3d vtxB_t = vtxB_t0 * (1 - ccdData.t) + vtxB_t1 * ccdData.t;

    ES::V3d vtxA_t0 = asVec3d(vtxA0[0]);

    // signB = dot(pn, vtxB_t - p);
    signB = pn.dot(vtxB_t0 - p);
    signA = pn.dot(vtxA_t0 - p);
    // signA = -signB;
  }
  else if (ccdData.ccdCase == CCDC_COLLIDED) {
  }
  else {
    std::cerr << "Unaddressed case" << std::endl;
  }
}

void pgo::CCDKernel::TriangleCCD::dumpCCD(int, const double *pos_t0, const double *pos_t1, int, const ES::V3i *triangles,
  const CCDData &ccdData, std::vector<ES::V3d> &outVertices, std::vector<ES::V3i> &outTriangles)
{
  const double *vtxA0[3] = {
    pos_t0 + triangles[ccdData.triA][0] * 3,
    pos_t0 + triangles[ccdData.triA][1] * 3,
    pos_t0 + triangles[ccdData.triA][2] * 3
  };

  const double *vtxA1[3] = {
    pos_t1 + triangles[ccdData.triA][0] * 3,
    pos_t1 + triangles[ccdData.triA][1] * 3,
    pos_t1 + triangles[ccdData.triA][2] * 3
  };

  const double *vtxB0[3] = {
    pos_t0 + triangles[ccdData.triB][0] * 3,
    pos_t0 + triangles[ccdData.triB][1] * 3,
    pos_t0 + triangles[ccdData.triB][2] * 3
  };

  const double *vtxB1[3] = {
    pos_t1 + triangles[ccdData.triB][0] * 3,
    pos_t1 + triangles[ccdData.triB][1] * 3,
    pos_t1 + triangles[ccdData.triB][2] * 3
  };

  int numVtx = static_cast<int>(outVertices.size());
  // bot
  outVertices.emplace_back(vtxA0[0]);
  outVertices.emplace_back(vtxA0[1]);
  outVertices.emplace_back(vtxA0[2]);
  outTriangles.emplace_back(numVtx, numVtx + 1, numVtx + 2);
  numVtx += 3;

  // top
  outVertices.emplace_back(vtxA1[0]);
  outVertices.emplace_back(vtxA1[1]);
  outVertices.emplace_back(vtxA1[2]);
  outTriangles.emplace_back(numVtx, numVtx + 1, numVtx + 2);
  numVtx += 3;

  ES::V3d triA0_t0 = asVec3d(vtxA0[0]);
  ES::V3d triA1_t0 = asVec3d(vtxA0[1]);
  ES::V3d triA2_t0 = asVec3d(vtxA0[2]);
  ES::V3d triA0_t1 = asVec3d(vtxA1[0]);
  ES::V3d triA1_t1 = asVec3d(vtxA1[1]);
  ES::V3d triA2_t1 = asVec3d(vtxA1[2]);
  ES::V3d triA0_t = triA0_t0 * (1 - ccdData.t) + triA0_t1 * ccdData.t;
  ES::V3d triA1_t = triA1_t0 * (1 - ccdData.t) + triA1_t1 * ccdData.t;
  ES::V3d triA2_t = triA2_t0 * (1 - ccdData.t) + triA2_t1 * ccdData.t;

  // mid
  outVertices.emplace_back(triA0_t);
  outVertices.emplace_back(triA1_t);
  outVertices.emplace_back(triA2_t);
  outTriangles.emplace_back(numVtx, numVtx + 1, numVtx + 2);
  numVtx += 3;

  // bot
  outVertices.emplace_back(vtxB0[0]);
  outVertices.emplace_back(vtxB0[1]);
  outVertices.emplace_back(vtxB0[2]);
  outTriangles.emplace_back(numVtx, numVtx + 1, numVtx + 2);
  numVtx += 3;

  // top
  outVertices.emplace_back(vtxB1[0]);
  outVertices.emplace_back(vtxB1[1]);
  outVertices.emplace_back(vtxB1[2]);
  outTriangles.emplace_back(numVtx, numVtx + 1, numVtx + 2);
  numVtx += 3;

  ES::V3d triB0_t0 = asVec3d(vtxB0[0]);
  ES::V3d triB1_t0 = asVec3d(vtxB0[1]);
  ES::V3d triB2_t0 = asVec3d(vtxB0[2]);
  ES::V3d triB0_t1 = asVec3d(vtxB1[0]);
  ES::V3d triB1_t1 = asVec3d(vtxB1[1]);
  ES::V3d triB2_t1 = asVec3d(vtxB1[2]);
  ES::V3d triB0_t = triB0_t0 * (1 - ccdData.t) + triB0_t1 * ccdData.t;
  ES::V3d triB1_t = triB1_t0 * (1 - ccdData.t) + triB1_t1 * ccdData.t;
  ES::V3d triB2_t = triB2_t0 * (1 - ccdData.t) + triB2_t1 * ccdData.t;

  // mid
  outVertices.emplace_back(triB0_t);
  outVertices.emplace_back(triB1_t);
  outVertices.emplace_back(triB2_t);
  outTriangles.emplace_back(numVtx, numVtx + 1, numVtx + 2);
  numVtx += 3;
}

bool pgo::CCDKernel::TriangleCCD::CCDTest3rdParty(int, const double *pos_t0, const double *pos_t1, int, const ES::V3i *triangles, int triA, int triB, CCDData *ccdDataPtr)
{
  CCDData ccdData;
  ccdData.triA = triA;
  ccdData.triB = triB;
  ccdData.t = 1000;
  ccdData.ccdCase = CCDC_NONE;

  const double *vtxA0[3] = {
    pos_t0 + triangles[ccdData.triA][0] * 3,
    pos_t0 + triangles[ccdData.triA][1] * 3,
    pos_t0 + triangles[ccdData.triA][2] * 3
  };

  const double *vtxA1[3] = {
    pos_t1 + triangles[ccdData.triA][0] * 3,
    pos_t1 + triangles[ccdData.triA][1] * 3,
    pos_t1 + triangles[ccdData.triA][2] * 3
  };

  const double *vtxB0[3] = {
    pos_t0 + triangles[ccdData.triB][0] * 3,
    pos_t0 + triangles[ccdData.triB][1] * 3,
    pos_t0 + triangles[ccdData.triB][2] * 3
  };

  const double *vtxB1[3] = {
    pos_t1 + triangles[ccdData.triB][0] * 3,
    pos_t1 + triangles[ccdData.triB][1] * 3,
    pos_t1 + triangles[ccdData.triB][2] * 3
  };

  // if it is already collided
  if (Mesh::intersectTriTri(vtxA0[0], vtxA0[1], vtxA0[2], vtxB0[0], vtxB0[1], vtxB0[2])) {
    ccdData.ccdCase = CCDC_COLLIDED;
    ccdData.t = 0;

    if (ccdDataPtr)
      *ccdDataPtr = ccdData;

    return true;
  }

  for (int i = 0; i < 3; i++) {
    const double *vtx0 = vtxA0[i];
    const double *vtx1 = vtxA1[i];

    double thist, thisb2, thisb3;
    if (SafeCCDWrapper::VertexTriangleCCD(vtx0, vtx1, vtxB0[0], vtxB1[0], vtxB0[1], vtxB1[1], vtxB0[2], vtxB1[2], thist, &thisb2, &thisb3) &&
      ExactCCDWrapper::VertexTriangleCCD(vtx0, vtx1, vtxB0[0], vtxB1[0], vtxB0[1], vtxB1[1], vtxB0[2], vtxB1[2])) {
      if (thist < ccdData.t) {
        ccdData.t = thist;

        ccdData.ccdCase = CCDC_VTXA_TRIB;
        ccdData.weights[0] = 1 - thisb2 - thisb3;
        ccdData.weights[1] = thisb2;
        ccdData.weights[2] = thisb3;

        ccdData.vid = i;
      }
    }
  }

  for (int i = 0; i < 3; i++) {
    const double *vtx0 = vtxB0[i];
    const double *vtx1 = vtxB1[i];

    double thist, thisb2, thisb3;
    if (SafeCCDWrapper::VertexTriangleCCD(vtx0, vtx1, vtxA0[0], vtxA1[0], vtxA0[1], vtxA1[1], vtxA0[2], vtxA1[2], thist, &thisb2, &thisb3) &&
      ExactCCDWrapper::VertexTriangleCCD(vtx0, vtx1, vtxA0[0], vtxA1[0], vtxA0[1], vtxA1[1], vtxA0[2], vtxA1[2])) {
      if (thist < ccdData.t) {
        ccdData.t = thist;

        ccdData.ccdCase = CCDC_VTXB_TRIA;
        ccdData.weights[0] = 1 - thisb2 - thisb3;
        ccdData.weights[1] = thisb2;
        ccdData.weights[2] = thisb3;

        ccdData.vid = i;
      }
    }
  }

  for (int i = 0; i < 3; i++) {
    int v0 = i;
    int v1 = (i + 1) % 3;

    for (int j = 0; j < 3; j++) {
      int v2 = j;
      int v3 = (j + 1) % 3;

      double thist, thisr, thiss;
      if (SafeCCDWrapper::EdgeEdgeCCD(vtxA0[v0], vtxA1[v0], vtxA0[v1], vtxA1[v1], vtxB0[v2], vtxB1[v2], vtxB0[v3], vtxB1[v3], thist, &thisr, &thiss) &&
        ExactCCDWrapper::EdgeEdgeCCD(vtxA0[v0], vtxA1[v0], vtxA0[v1], vtxA1[v1], vtxB0[v2], vtxB1[v2], vtxB0[v3], vtxB1[v3])) {
        if (thist < ccdData.t) {
          ccdData.t = thist;

          ccdData.ccdCase = CCDC_EDGEA_EDGEB;
          ccdData.weights[0] = 1 - thisr;
          ccdData.weights[1] = thisr;

          ccdData.weights[2] = 1 - thiss;
          ccdData.weights[3] = thiss;

          ccdData.eid[0] = v0;
          ccdData.eid[1] = v2;
        }
      }
    }
  }

  if (ccdDataPtr)
    *ccdDataPtr = ccdData;

  if (ccdData.ccdCase != CCDC_NONE)
    return true;
  else
    return false;
}

bool pgo::CCDKernel::TriangleCCD::CCDTestVega(int, const double *pos_t0, const double *pos_t1, int, const ES::V3i *triangles, int triA, int triB, CCDData *ccdDataPtr)
{
#if 0
  CCDData ccdData;
  ccdData.triA = triA;
  ccdData.triB = triB;
  ccdData.t = 1000;
  ccdData.ccdCase = CCDC_NONE;

  const double *vtxA0[3] = {
    pos_t0 + triangles[ccdData.triA][0] * 3,
    pos_t0 + triangles[ccdData.triA][1] * 3,
    pos_t0 + triangles[ccdData.triA][2] * 3
  };

  const double *vtxA1[3] = {
    pos_t1 + triangles[ccdData.triA][0] * 3,
    pos_t1 + triangles[ccdData.triA][1] * 3,
    pos_t1 + triangles[ccdData.triA][2] * 3
  };

  const double *vtxB0[3] = {
    pos_t0 + triangles[ccdData.triB][0] * 3,
    pos_t0 + triangles[ccdData.triB][1] * 3,
    pos_t0 + triangles[ccdData.triB][2] * 3
  };

  const double *vtxB1[3] = {
    pos_t1 + triangles[ccdData.triB][0] * 3,
    pos_t1 + triangles[ccdData.triB][1] * 3,
    pos_t1 + triangles[ccdData.triB][2] * 3
  };

  if (triangleVsTriangleOverlap(vtxA0[0], vtxA0[1], vtxA0[2], vtxB0[0], vtxB0[1], vtxB0[2])) {
    ccdData.t = 0;
    ccdData.ccdCase = CCDC_COLLIDED;

    if (ccdDataPtr)
      *ccdDataPtr = ccdData;

    return true;
  }

  const double *triangleAPos[3] = { vtxA0[0], vtxA0[1], vtxA0[2] };
  const double *triangleBPos[3] = { vtxB0[0], vtxB0[1], vtxB0[2] };

  double p0Motion[3] = { vtxA1[0][0] - vtxA0[0][0], vtxA1[0][1] - vtxA0[0][1], vtxA1[0][2] - vtxA0[0][2] };
  double p1Motion[3] = { vtxA1[1][0] - vtxA0[1][0], vtxA1[1][1] - vtxA0[1][1], vtxA1[1][2] - vtxA0[1][2] };
  double p2Motion[3] = { vtxA1[2][0] - vtxA0[2][0], vtxA1[2][1] - vtxA0[2][1], vtxA1[2][2] - vtxA0[2][2] };
  double q0Motion[3] = { vtxB1[0][0] - vtxB0[0][0], vtxB1[0][1] - vtxB0[0][1], vtxB1[0][2] - vtxB0[0][2] };
  double q1Motion[3] = { vtxB1[1][0] - vtxB0[1][0], vtxB1[1][1] - vtxB0[1][1], vtxB1[1][2] - vtxB0[1][2] };
  double q2Motion[3] = { vtxB1[2][0] - vtxB0[2][0], vtxB1[2][1] - vtxB0[2][1], vtxB1[2][2] - vtxB0[2][2] };

  double *triangleAMotion[3] = { p0Motion, p1Motion, p2Motion };
  double *triangleBMotion[3] = { q0Motion, q1Motion, q2Motion };

  //printf("====\n");
  //printf("Triangle A:\n");
  //printf("%G %G %G\n", triangleAPos[0][0], triangleAPos[0][1], triangleAPos[0][2]);
  //printf("%G %G %G\n", triangleAPos[1][0], triangleAPos[1][1], triangleAPos[1][2]);
  //printf("%G %G %G\n", triangleAPos[2][0], triangleAPos[2][1], triangleAPos[2][2]);
  //printf("Triangle B:\n");
  //printf("%G %G %G\n", triangleBPos[0][0], triangleBPos[0][1], triangleBPos[0][2]);
  //printf("%G %G %G\n", triangleBPos[1][0], triangleBPos[1][1], triangleBPos[1][2]);
  //printf("%G %G %G\n", triangleBPos[2][0], triangleBPos[2][1], triangleBPos[2][2]);

  double radius = 1.1;

  double comA[3];
  double comB[3];
  for (int i = 0; i < 3; i++) {
    comA[i] = 1.0 / 3 * (triangleAPos[0][i] + triangleAPos[1][i] + triangleAPos[2][i]);
    comB[i] = 1.0 / 3 * (triangleBPos[0][i] + triangleBPos[1][i] + triangleBPos[2][i]);
  }

  ES::V3d comAv = asVec3d(comA);
  ES::V3d comBv = asVec3d(comB);
  ES::V3d APosv[3] = { asVec3d(triangleAPos[0]), asVec3d(triangleAPos[1]), asVec3d(triangleAPos[2]) };
  ES::V3d BPosv[3] = { asVec3d(triangleBPos[0]), asVec3d(triangleBPos[1]), asVec3d(triangleBPos[2]) };

  double r0 = 0;
  for (int i = 0; i < 3; i++) {
    double l = len(APosv[i] - comAv);
    if (l > r0)
      r0 = l;
  }

  double r1 = 0;
  for (int i = 0; i < 3; i++) {
    double l = len(BPosv[i] - comBv);
    if (l > r1)
      r1 = l;
  }

  double triDist = len(comAv - comBv) - r0 - r1;
  if (triDist < 0)
    triDist = 0;

  double modeMax = 0;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      ES::V3d modeA(triangleAMotion[i][0], triangleAMotion[i][1], triangleAMotion[i][2]);
      ES::V3d modeB(triangleBMotion[j][0], triangleBMotion[j][1], triangleBMotion[j][2]);
      double l = len(modeA - modeB);
      if (l > modeMax)
        modeMax = l;
    }

  if (modeMax < triDist)
    return 0;

  // check edge vs edge
  int mask[3][2] = { { 0, 1 }, { 1, 2 }, { 2, 0 } };
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      double Up0[3] = { triangleAMotion[mask[i][0]][0],
        triangleAMotion[mask[i][0]][1],
        triangleAMotion[mask[i][0]][2] };

      double Up1[3] = { triangleAMotion[mask[i][1]][0],
        triangleAMotion[mask[i][1]][1],
        triangleAMotion[mask[i][1]][2] };

      double Us0[3] = { triangleBMotion[mask[j][0]][0],
        triangleBMotion[mask[j][0]][1],
        triangleBMotion[mask[j][0]][2] };

      double Us1[3] = { triangleBMotion[mask[j][1]][0],
        triangleBMotion[mask[j][1]][1],
        triangleBMotion[mask[j][1]][2] };

      double edgeRadius = TriTriCCD::EdgevsEdgeCCD(triangleAPos[mask[i][0]], triangleAPos[mask[i][1]], triangleBPos[mask[j][0]], triangleBPos[mask[j][1]], Up0, Up1, Us0, Us1, radius);
      //printf("e %d %d %G\n", i, j, edgeRadius);

      if (edgeRadius < radius) {
        radius = edgeRadius;
        //printf("e %d %d %G\n", i, j, radius);

        ccdData.t = radius;
        ccdData.eid[0] = mask[i][0];
        ccdData.eid[1] = mask[j][0];
        ccdData.ccdCase = CCDC_EDGEA_EDGEB;
      }
    }

  // vtx vs face
  for (int i = 0; i < 3; i++) {
    double Up0[3] = { triangleAMotion[i][0],
      triangleAMotion[i][1],
      triangleAMotion[i][2] };

    double Us0[3] = { triangleBMotion[0][0],
      triangleBMotion[0][1],
      triangleBMotion[0][2] };

    double Us1[3] = { triangleBMotion[1][0],
      triangleBMotion[1][1],
      triangleBMotion[1][2] };

    double Us2[3] = { triangleBMotion[2][0],
      triangleBMotion[2][1],
      triangleBMotion[2][2] };

    //printf("=!!!===\n");
    //printf("Triangle A:\n");
    //printf("%G %G %G\n", triangleAPos[0][0], triangleAPos[0][1], triangleAPos[0][2]);
    //printf("%G %G %G\n", triangleAPos[1][0], triangleAPos[1][1], triangleAPos[1][2]);
    //printf("%G %G %G\n", triangleAPos[2][0], triangleAPos[2][1], triangleAPos[2][2]);
    //printf("Triangle B:\n");
    //printf("%G %G %G\n", triangleBPos[0][0], triangleBPos[0][1], triangleBPos[0][2]);
    //printf("%G %G %G\n", triangleBPos[1][0], triangleBPos[1][1], triangleBPos[1][2]);
    //printf("%G %G %G\n", triangleBPos[2][0], triangleBPos[2][1], triangleBPos[2][2]);

    double vtxFaceRadius = TriTriCCD::VtxvsFaceCCD(triangleAPos[i], triangleBPos[0], triangleBPos[1], triangleBPos[2], Up0, Us0, Us1, Us2, radius);

    //printf("f1 %d %G\n", i, radius);

    if (vtxFaceRadius < radius) {
      radius = vtxFaceRadius;
      //printf("f1 %d %G\n", i, radius);

      ccdData.t = radius;
      ccdData.vid = i;
      ccdData.ccdCase = CCDC_VTXA_TRIB;
    }
  }

  //printf("=!!!===\n");
  //printf("Triangle A:\n");
  //printf("%G %G %G\n", triangleAPos[0][0], triangleAPos[0][1], triangleAPos[0][2]);
  //printf("%G %G %G\n", triangleAPos[1][0], triangleAPos[1][1], triangleAPos[1][2]);
  //printf("%G %G %G\n", triangleAPos[2][0], triangleAPos[2][1], triangleAPos[2][2]);
  //printf("Triangle B:\n");
  //printf("%G %G %G\n", triangleBPos[0][0], triangleBPos[0][1], triangleBPos[0][2]);
  //printf("%G %G %G\n", triangleBPos[1][0], triangleBPos[1][1], triangleBPos[1][2]);
  //printf("%G %G %G\n", triangleBPos[2][0], triangleBPos[2][1], triangleBPos[2][2]);

  // vtx vs face
  for (int i = 0; i < 3; i++) {
    double Up0[3] = { triangleBMotion[i][0],
      triangleBMotion[i][1],
      triangleBMotion[i][2] };

    double Us0[3] = { triangleAMotion[0][0],
      triangleAMotion[0][1],
      triangleAMotion[0][2] };

    double Us1[3] = { triangleAMotion[1][0],
      triangleAMotion[1][1],
      triangleAMotion[1][2] };

    double Us2[3] = { triangleAMotion[2][0],
      triangleAMotion[2][1],
      triangleAMotion[2][2] };

    //printf("i==%d\n", i);
    double vtxFaceRadius = TriTriCCD::VtxvsFaceCCD(triangleBPos[i], triangleAPos[0], triangleAPos[1], triangleAPos[2], Up0, Us0, Us1, Us2, radius);

    //printf("f2 %d %d %G\n", i, modeIndex, radius);

    if (vtxFaceRadius < radius) {
      radius = vtxFaceRadius;
      //printf("f2 %d %G\n", i, radius);

      ccdData.t = radius;
      ccdData.vid = i;
      ccdData.ccdCase = CCDC_VTXB_TRIA;
    }
  }

  if (ccdDataPtr)
    *ccdDataPtr = ccdData;

  if (ccdData.ccdCase != CCDC_NONE)
    return true;
  else
    return false;
#else
  return true;
#endif
}

// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2019
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)
// Compute the closest points on the line segments P(s) = (1-s)*P0 + s*P1 and
// Q(t) = (1-t)*Q0 + t*Q1 for 0 <= s <= 1 and 0 <= t <= 1.  The algorithm
// relies on exact rational arithmetic.

double pgo::CCDKernel::distanceSegmentSegmentExact(const double P0[3], const double P1[3], const double Q0[3], const double Q1[3], double w[2])
{
  ES::V3d P1mP0 = ES::Mp<const ES::V3d>(P1) - ES::Mp<const ES::V3d>(P0);
  ES::V3d Q1mQ0 = ES::Mp<const ES::V3d>(Q1) - ES::Mp<const ES::V3d>(Q0);
  ES::V3d P0mQ0 = ES::Mp<const ES::V3d>(P0) - ES::Mp<const ES::V3d>(Q0);
  double a = P1mP0.dot(P1mP0);
  double b = P1mP0.dot(Q1mQ0);
  double c = Q1mQ0.dot(Q1mQ0);
  double d = P1mP0.dot(P0mQ0);
  double e = Q1mQ0.dot(P0mQ0);
  constexpr double zero = (double)0;
  constexpr double one = (double)1;
  double det = a * c - b * b;
  double s, t, nd, bmd, bte, ctd, bpe, ate, btd;

  if (det > zero) {
    bte = b * e;
    ctd = c * d;
    if (bte <= ctd)  // s <= 0
    {
      s = zero;
      if (e <= zero)  // t <= 0
      {
        // region 6
        t = zero;
        nd = -d;
        if (nd >= a) {
          s = one;
        }
        else if (nd > zero) {
          s = nd / a;
        }
        // else: s is already zero
      }
      else if (e < c)  // 0 < t < 1
      {
        // region 5
        t = e / c;
      }
      else  // t >= 1
      {
        // region 4
        t = one;
        bmd = b - d;
        if (bmd >= a) {
          s = one;
        }
        else if (bmd > zero) {
          s = bmd / a;
        }
        // else:  s is already zero
      }
    }
    else  // s > 0
    {
      s = bte - ctd;
      if (s >= det)  // s >= 1
      {
        // s = 1
        s = one;
        bpe = b + e;
        if (bpe <= zero)  // t <= 0
        {
          // region 8
          t = zero;
          nd = -d;
          if (nd <= zero) {
            s = zero;
          }
          else if (nd < a) {
            s = nd / a;
          }
          // else: s is already one
        }
        else if (bpe < c)  // 0 < t < 1
        {
          // region 1
          t = bpe / c;
        }
        else  // t >= 1
        {
          // region 2
          t = one;
          bmd = b - d;
          if (bmd <= zero) {
            s = zero;
          }
          else if (bmd < a) {
            s = bmd / a;
          }
          // else:  s is already one
        }
      }
      else  // 0 < s < 1
      {
        ate = a * e;
        btd = b * d;
        if (ate <= btd)  // t <= 0
        {
          // region 7
          t = zero;
          nd = -d;
          if (nd <= zero) {
            s = zero;
          }
          else if (nd >= a) {
            s = one;
          }
          else {
            s = nd / a;
          }
        }
        else  // t > 0
        {
          t = ate - btd;
          if (t >= det)  // t >= 1
          {
            // region 3
            t = one;
            bmd = b - d;
            if (bmd <= zero) {
              s = zero;
            }
            else if (bmd >= a) {
              s = one;
            }
            else {
              s = bmd / a;
            }
          }
          else  // 0 < t < 1
          {
            // region 0
            s /= det;
            t /= det;
          }
        }
      }
    }
  }
  else {
    // The segments are parallel.  The quadratic factors to R(s,t) =
    // a*(s-(b/a)*t)^2 + 2*d*(s - (b/a)*t) + f, where a*c = b^2,
    // e = b*d/a, f = |P0-Q0|^2, and b is not zero.  R is constant along
    // lines of the form s-(b/a)*t = k, and the minimum of R occurs on the
    // line a*s - b*t + d = 0.  This line must intersect both the s-axis
    // and the t-axis because 'a' and 'b' are not zero.  Because of
    // parallelism, the line is also represented by -b*s + c*t - e = 0.
    //
    // The code determines an edge of the domain [0,1]^2 that intersects
    // the minimum line, or if none of the edges intersect, it determines
    // the closest corner to the minimum line.  The conditionals are
    // designed to test first for intersection with the t-axis (s = 0)
    // using -b*s + c*t - e = 0 and then with the s-axis (t = 0) using
    // a*s - b*t + d = 0.

    // When s = 0, solve c*t - e = 0 (t = e/c).
    if (e <= zero)  // t <= 0
    {
      // Now solve a*s - b*t + d = 0 for t = 0 (s = -d/a).
      t = zero;
      nd = -d;
      if (nd <= zero)  // s <= 0
      {
        // region 6
        s = zero;
      }
      else if (nd >= a)  // s >= 1
      {
        // region 8
        s = one;
      }
      else  // 0 < s < 1
      {
        // region 7
        s = nd / a;
      }
    }
    else if (e >= c)  // t >= 1
    {
      // Now solve a*s - b*t + d = 0 for t = 1 (s = (b-d)/a).
      t = one;
      bmd = b - d;
      if (bmd <= zero)  // s <= 0
      {
        // region 4
        s = zero;
      }
      else if (bmd >= a)  // s >= 1
      {
        // region 2
        s = one;
      }
      else  // 0 < s < 1
      {
        // region 3
        s = bmd / a;
      }
    }
    else  // 0 < t < 1
    {
      // The point (0,e/c) is on the line and domain, so we have one
      // point at which R is a minimum.
      s = zero;
      t = e / c;
    }
  }

  ES::V3d c0 = ES::Mp<const ES::V3d>(P0) + s * P1mP0;
  ES::V3d c1 = ES::Mp<const ES::V3d>(Q0) + t * Q1mQ0;
  ES::V3d diff = c1 - c0;

  if (w) {
    w[0] = s;
    w[1] = t;
  }

  return diff.squaredNorm();
}
