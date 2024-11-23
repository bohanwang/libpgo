#pragma once

class ExactCCDWrapper
{
public:
  static bool VertexTriangleCCD(
    const double x0_0[3], const double x0_1[3], 
    const double x1_0[3], const double x1_1[3],
    const double x2_0[3], const double x2_1[3], 
    const double x3_0[3], const double x3_1[3]);

  static bool EdgeEdgeCCD(
    const double x0_0[], const double x0_1[], 
    const double x1_0[], const double x1_1[], 
    const double x2_0[], const double x2_1[], 
    const double x3_0[], const double x3_1[]);
};