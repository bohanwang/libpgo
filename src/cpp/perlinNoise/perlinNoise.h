#pragma once

/*
  Code author: Yili Zhao

  This is an implementation of Ken Perlin's "Improved Noise".

  Reference:
  Perlin, Ken. 2002. "Improving Noise." ACM Transactions on Graphics (Proceedings of SIGGRAPH 2002) 21(3), pp. 681--682.
*/

namespace pgo
{
namespace PerlinNoise
{
class PerlinNoise
{
public:
  PerlinNoise();

  double noise1D(double x);
  double noise2D(double x[2]);
  double noise3D(double x[3]);
  double noise4D(double x[4]);

protected:
  static constexpr int TWO_FIVE_SIX = 256;
  static constexpr int TWO_FIVE_FIVE = 255;

  int permutation[TWO_FIVE_SIX];
  void generatePermutationTable(void);
  int getPermutation(int pos);
  double computeFadeCurves(double t);
  double linearInterpolation(double t, double a, double b);
  double gradient1D(int hash, double x);
  double gradient2D(int hash, double x, double y);
  double gradient3D(int hash, double x, double y, double z);
  double gradient4D(int hash, double x, double y, double z, double t);

  // the following gradient directions are constant (they could be made static)
  double g2[4][2];   // gradients for 2d noise (2 * 2 groups)
  double g3[16][3];  // gradients for 3d noise (3 * 4 groups, however, to avoid dividing by 12, we add an extra (1,1,0), (-1,1,0), (0,-1,1) and (0,-1,-1). Please check the reference paper.)
  double g4[32][4];  // gradients for 4d noise (4 * 8 groups)
};
}  // namespace PerlinNoise
}  // namespace pgo
