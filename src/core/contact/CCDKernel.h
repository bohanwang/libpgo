/*
author: Bohan Wang
copyright to USC
*/

#pragma once

#include "EigenSupport.h"

#include <vector>

namespace pgo
{
namespace CCDKernel
{
namespace TriangleCCD
{
enum CCD_CASE
{
  CCDC_NONE,
  CCDC_VTXA_TRIB,
  CCDC_VTXB_TRIA,
  CCDC_EDGEA_EDGEB,
  CCDC_COLLIDED
};

enum CCD_METHOD
{
  CCDM_VEGA,
  CCDM_3RDPARTY
};

struct CCDData
{
  int triA, triB;
  CCD_CASE ccdCase;
  double t;
  int vid;
  int eid[2];
  double weights[4];
};

bool CCDTest(int numVertices, const double *pos_t0, const double *pos_t1, int numTriangles, const EigenSupport::V3i *triangles,
  int triA, int triB, CCD_METHOD m, CCDData *ccdDataPtr = nullptr);

bool computeCCDPlane(int numVertices, const double *pos_t0, const double *pos_t1, int numTriangles, const EigenSupport::V3i *triangles,
  const CCDData &ccdData, bool doSign, double planePoint[3], double planeNormal[3], double &signA, double &signB);

bool computeDCDPlane(int numVertices, const double *pos_t0, int numTriangles, const EigenSupport::V3i *triangles,
  CCDData &ccdData, double planePoint[3], double planeNormal[3]);

void computeSign(int numVertices, const double *pos_t0, const double *pos_t1, int numTriangles, const EigenSupport::V3i *triangles,
  const CCDData &ccdData, const double planePoint[3], const double planeNormal[3], double &signA, double &signB);

void dumpCCD(int numVertices, const double *pos_t0, const double *pos_t1, int numTriangles, const EigenSupport::V3i *triangles,
  const CCDData &ccdData, std::vector<EigenSupport::V3d> &outVertices, std::vector<EigenSupport::V3i> &outTriangles);
}  // namespace TriangleCCD

double distanceSegmentSegmentExact(const double P0[3], const double P1[3], const double Q0[3], const double Q1[3], double w[2]);
}  // namespace CCDKernel
}  // namespace pgo
