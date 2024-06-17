#pragma once

class SafeCCDWrapper
{
public:
  static bool VertexTriangleCCD(
    const double x0_0[3], const double x0_1[3], 
    const double x1_0[3], const double x1_1[3],
    const double x2_0[3], const double x2_1[3], 
    const double x3_0[3], const double x3_1[3], 
    double &t, double *b2 = nullptr, double *b3 = nullptr);

  static bool EdgeEdgeCCD(
    const double x0_0[3], const double x0_1[3], 
    const double x1_0[3], const double x1_1[3], 
    const double x2_0[3], const double x2_1[3], 
    const double x3_0[3], const double x3_1[3], 
    double &t, double *r = nullptr, double *s = nullptr);

protected:
  static double computeB(
    const double x0_0[3], const double x0_1[3],
    const double x1_0[3], const double x1_1[3],
    const double x2_0[3], const double x2_1[3],
    const double x3_0[3], const double x3_1[3]);
};