#include "exactCCDWrapper.h"

#include "rootparitycollisiontest.h"

bool ExactCCDWrapper::VertexTriangleCCD(
    const double x0_0[3], const double x0_1[3], 
    const double x1_0[3], const double x1_1[3],
    const double x2_0[3], const double x2_1[3], 
    const double x3_0[3], const double x3_1[3])
{
  Vec3d x0old(x0_0);
  Vec3d x1old(x1_0);
  Vec3d x2old(x2_0);
  Vec3d x3old(x3_0);


  Vec3d x0new(x0_1);
  Vec3d x1new(x1_1);
  Vec3d x2new(x2_1);
  Vec3d x3new(x3_1);

  rootparity::RootParityCollisionTest kernel(x0old, x1old, x2old, x3old, x0new, x1new, x2new, x3new, false);

  return kernel.run_test();
}

bool ExactCCDWrapper::EdgeEdgeCCD(
    const double x0_0[], const double x0_1[], 
    const double x1_0[], const double x1_1[], 
    const double x2_0[], const double x2_1[], 
    const double x3_0[], const double x3_1[])
{
  Vec3d x0old(x0_0);
  Vec3d x1old(x1_0);
  Vec3d x2old(x2_0);
  Vec3d x3old(x3_0);


  Vec3d x0new(x0_1);
  Vec3d x1new(x1_1);
  Vec3d x2new(x2_1);
  Vec3d x3new(x3_1);

  rootparity::RootParityCollisionTest kernel(x0old, x1old, x2old, x3old, x0new, x1new, x2new, x3new, true);

  return kernel.run_test();
}