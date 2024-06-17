#include "safeCCDWrapper.h"

#include "SAFE_CCD.h"


#include <cmath>

bool SafeCCDWrapper::VertexTriangleCCD(
    const double x0_0[3], const double x0_1[3], 
    const double x1_0[3], const double x1_1[3],
    const double x2_0[3], const double x2_1[3], 
    const double x3_0[3], const double x3_1[3], 
    double &t, double *b2, double *b3)
{
  SAFE_CCD<double> ccd;
  double B = computeB(x0_0, x0_1, x1_0, x1_1, x2_0, x2_1, x3_0, x3_1);
  ccd.Set_Coefficients(B);

  return ccd.Vertex_Triangle_CCD(
    const_cast<double*>(x0_0), const_cast<double*>(x0_1), 
    const_cast<double*>(x1_0), const_cast<double*>(x1_1), 
    const_cast<double*>(x2_0), const_cast<double*>(x2_1), 
    const_cast<double*>(x3_0), const_cast<double*>(x3_1),
    t, b2, b3);
}

bool SafeCCDWrapper::EdgeEdgeCCD(
    const double x0_0[], const double x0_1[], 
    const double x1_0[], const double x1_1[], 
    const double x2_0[], const double x2_1[], 
    const double x3_0[], const double x3_1[], 
    double &t, double *r, double *s)
{
  SAFE_CCD<double> ccd;

  double B = computeB(x0_0, x0_1, x1_0, x1_1, x2_0, x2_1, x3_0, x3_1);
  ccd.Set_Coefficients(B);
  return ccd.Edge_Edge_CCD(
    const_cast<double*>(x0_0), const_cast<double*>(x0_1), 
    const_cast<double*>(x1_0), const_cast<double*>(x1_1), 
    const_cast<double*>(x2_0), const_cast<double*>(x2_1),
    const_cast<double*>(x3_0), const_cast<double*>(x3_1),
    t, r, s);
}

double SafeCCDWrapper::computeB(
  const double x0_0[3], const double x0_1[3], 
  const double x1_0[3], const double x1_1[3], 
  const double x2_0[3], const double x2_1[3],
  const double x3_0[3], const double x3_1[3])
{

  double pos[6][3];
  SUB(x0_0, x1_0, pos[0]);
  SUB(x0_0, x2_0, pos[1]);
  SUB(x0_0, x3_0, pos[2]);
  SUB(x1_0, x2_0, pos[3]);
  SUB(x1_0, x3_0, pos[4]);
  SUB(x2_0, x3_0, pos[5]);

  double E = 0;
  for (int i = 0; i < 6; i++) {
    E = std::fmax(E, std::fabs(pos[i][0]));
    E = std::fmax(E, std::fabs(pos[i][1]));
    E = std::fmax(E, std::fabs(pos[i][2]));
  }

  double V = 0;
  SUB(x0_0, x0_1, pos[0]);
  SUB(x1_0, x1_1, pos[1]);
  SUB(x2_0, x2_1, pos[2]);
  SUB(x3_0, x3_1, pos[3]);
  for (int i = 0; i < 4; i++) {
    V = std::fmax(V, std::fabs(pos[i][0]));
    V = std::fmax(V, std::fabs(pos[i][1]));
    V = std::fmax(V, std::fabs(pos[i][2]));
  }

  return (E + 2 * V) * (1.0 + 0x1.p-53) * (1.0 + 0x1.p-53);
}
