/*
  Code author: Yili Zhao

  This is an implementation of Ken Perlin's "Improved Noise".

  Reference:
  Perlin, Ken. 2002. "Improving Noise." ACM Transactions on Graphics (Proceedings of SIGGRAPH 2002) 21(3), pp. 681--682. 
*/

#include "perlinNoise.h"

#include <cstdlib>
#include <cstdio>
#include <cstring>

using namespace pgo::PerlinNoise;

PerlinNoise::PerlinNoise()
{
  generatePermutationTable();

  // init g2
  {
    int g[4][2] = 
    {
     { 1,  0}, {-1,  0}, { 0,  1}, { 0, -1},
    };

    for(int i=0; i<4; i++)
    {
      for(int j=0; j<2; j++)
        g2[i][j] = g[i][j];
    }
  }

  // init g3
  {
    int g[16][3] =
    {
      { 1, 1,  0}, {-1, 1,  0}, { 1,-1,  0}, {-1,-1,  0},
      { 1, 0,  1}, {-1, 0,  1}, { 1, 0, -1}, {-1, 0, -1}, 
      { 0, 1,  1}, { 0,-1,  1}, { 0, 1, -1}, { 0,-1, -1},
      { 1, 1,  0}, {-1, 1,  0}, { 0,-1,  1}, { 0,-1, -1},
    };

    for(int i=0; i<16; i++)
    {
      for(int j=0; j<3; j++)
        g3[i][j] = g[i][j];
    }
  }

  // init g4
  {
    int g[32][4] =
    { 
      { 1, 1,  1,  0}, {-1,  1,  1,  0}, { 1, -1,  1,  0}, {-1, -1,   1,  0},            
      { 1, 1, -1,  0}, {-1,  1, -1,  0}, { 1, -1, -1,  0}, {-1, -1,  -1,  0},            
      { 1, 1,  0,  1}, {-1,  1,  0,  1}, { 1, -1,  0,  1}, {-1, -1,   0,  1}, 
      { 1, 1,  0, -1}, {-1,  1,  0, -1}, { 1, -1,  0, -1}, {-1, -1,   0, -1}, 
      { 1, 0,  1,  1}, {-1,  0,  1,  1}, { 1, 0,  -1,  1}, {-1,  0,  -1,  1},
      { 1, 0,  1, -1}, {-1,  0,  1, -1}, { 1, 0,  -1, -1}, {-1,  0,  -1, -1},
      { 0, 1,  1,  1}, { 0, -1,  1,  1}, { 0, 1,  -1,  1}, { 0, -1,  -1,  1},
      { 0, 1,  1, -1}, { 0, -1,  1, -1}, { 0, 1,  -1, -1}, { 0, -1,  -1, -1},
    };

    for(int i=0; i<32; i++)
    {
      for(int j=0; j<4; j++)
        g4[i][j] = g[i][j];
    }
  }
}

void PerlinNoise::generatePermutationTable(void)
{
  for(int i=0; i<TWO_FIVE_SIX; i++)
    permutation[i] = i;
  for(int i=0; i<TWO_FIVE_SIX; i++)
  {
    int temp = permutation[i];
    int j = rand() % TWO_FIVE_SIX;
    permutation[i] = permutation[j];
    permutation[j] = temp;
  }
}

double PerlinNoise::computeFadeCurves(double t)
{
  return t * t * t * (t * (t * 6.0 - 15.0) + 10.0); // new curve (good, see "Improving Noise")
  //	return t * t * (3 - 2 * t); // old curve (suboptimal since it does not have non-zero 2nd derivative at either t=0 or t=1) 
}

double PerlinNoise::linearInterpolation(double t, double a, double b)
{
  return (a + t * (b - a));
}

int PerlinNoise::getPermutation(int pos)
{
  return (permutation[pos & TWO_FIVE_FIVE]);
}

double PerlinNoise::gradient1D(int hash, double x)
{
  if (hash % 2 == 0)
    return x;
  else
    return -x;
}

double PerlinNoise::noise1D(double x)
{
  // find the line segment containing the 1D point
  int lower;
  if (x > 0.0)
    lower = (int)x;
  else
    lower = (int)(x - 1);

  double relative = x - lower;
  double u = computeFadeCurves(relative);

  return linearInterpolation(u, gradient1D(getPermutation(lower), relative), gradient1D(getPermutation(lower+1), relative-1));
}

double PerlinNoise::gradient2D(int hash, double x, double y)
{
  int index = hash & 3;
  return (g2[index][0] * x + g2[index][1] * y);
}

double PerlinNoise::noise2D(double x[2])
{
  enum {X=0, Y, NUM_DIM};
  
  int lower[NUM_DIM];
  double relative[NUM_DIM];  // relative value of p in the grid
  double u[NUM_DIM];
  
  for(int dim=0; dim<NUM_DIM; dim++)
  {
    if (x[dim] > 0.0)
      lower[dim] = (int)(x[dim]);
    else
      lower[dim] = (int)(x[dim] - 1);

    relative[dim] = x[dim] - lower[dim];
    u[dim] = computeFadeCurves(relative[dim]);
  }

  int A = getPermutation(lower[X]) + lower[Y];
  int B = getPermutation(lower[X] + 1) + lower[Y];

  // add blended results from four grid corners
  return linearInterpolation(u[Y], linearInterpolation(u[X], gradient2D(getPermutation(A),   relative[X],   relative[Y]),  
                                                             gradient2D(getPermutation(B),   relative[X]-1, relative[Y])),
                                   linearInterpolation(u[X], gradient2D(getPermutation(A+1), relative[X],   relative[Y]-1),
                                                             gradient2D(getPermutation(B+1), relative[X]-1, relative[Y]-1)));
}

double PerlinNoise::gradient3D(int hash, double x, double y, double z)
{
  int index = hash & 15;
  return (g3[index][0] * x + g3[index][1] * y + g3[index][2] * z);
}

double PerlinNoise::noise3D(double x[3])
{
  enum {X=0, Y, Z, NUM_DIM};
  int lower[NUM_DIM];
  double relative[NUM_DIM];  // relative value of x in the cube
  double u[NUM_DIM];

  for(int dim=0; dim<NUM_DIM; dim++)
  {
    if (x[dim] > 0.0)
      lower[dim] = (int)(x[dim]);
    else
      lower[dim] = (int)(x[dim] - 1);

    relative[dim] = x[dim] - lower[dim];
    u[dim] = computeFadeCurves(relative[dim]);
  }

  // hash coordinates of the 8 corners of the cube
  int A = getPermutation(lower[X]) + lower[Y];
  int B = getPermutation(lower[X] + 1) + lower[Y];

  int AA = getPermutation(A) + lower[Z];
  int AB = getPermutation(A+1) + lower[Z];
  int BA = getPermutation(B) + lower[Z];
  int BB = getPermutation(B+1) + lower[Z];

  return linearInterpolation(u[Z], linearInterpolation(u[Y], linearInterpolation(u[X], gradient3D(getPermutation(AA),   relative[X],   relative[Y],   relative[Z]),
                                                                                       gradient3D(getPermutation(BA),   relative[X]-1, relative[Y],   relative[Z])),
                                                             linearInterpolation(u[X], gradient3D(getPermutation(AB),   relative[X],   relative[Y]-1, relative[Z]),  
                                                                                       gradient3D(getPermutation(BB),   relative[X]-1, relative[Y]-1, relative[Z]))),
                                   linearInterpolation(u[Y], linearInterpolation(u[X], gradient3D(getPermutation(AA+1), relative[X],   relative[Y],   relative[Z]-1), 
                                                                                       gradient3D(getPermutation(BA+1), relative[X]-1, relative[Y]  , relative[Z]-1)),
                                                             linearInterpolation(u[X], gradient3D(getPermutation(AB+1), relative[X],   relative[Y]-1, relative[Z]-1),
                                                                                       gradient3D(getPermutation(BB+1), relative[X]-1, relative[Y]-1, relative[Z]-1))));
}

double PerlinNoise::gradient4D(int hash, double x, double y, double z, double t)
{
  int index = hash & 31;
  return (g4[index][0] * x + g4[index][1] * y + g4[index][2] * z + g4[index][3] * t);
}

double PerlinNoise::noise4D(double x[4])
{
  enum {X, Y, Z, T, NUM_DIM};
  int lower[NUM_DIM];
  double relative[NUM_DIM];  // relative value of x in the hyper-cube
  double u[NUM_DIM];

  for(int dim=0; dim<NUM_DIM; dim++)
  {
    if (x[dim] > 0.0)
      lower[dim] = (int)(x[dim]);
    else
      lower[dim] = (int)(x[dim] - 1);

    relative[dim] = x[dim] - lower[dim];
    u[dim] = computeFadeCurves(relative[dim]);
  }

  // hash coordinates of the 16 corners of the cube

  int A = getPermutation(lower[X])     + lower[Y];
  int B = getPermutation(lower[X] + 1) + lower[Y];

  int AA = getPermutation(A)   + lower[Z];
  int AB = getPermutation(A+1) + lower[Z];
  int BA = getPermutation(B)   + lower[Z];
  int BB = getPermutation(B+1) + lower[Z];

  int AAA = getPermutation(AA)   + lower[T];
  int AAB = getPermutation(AA+1) + lower[T];
  int ABA = getPermutation(AB)   + lower[T];
  int ABB = getPermutation(AB+1) + lower[T];
  int BAA = getPermutation(BA)   + lower[T];
  int BAB = getPermutation(BA+1) + lower[T];
  int BBA = getPermutation(BB)   + lower[T];
  int BBB = getPermutation(BB+1) + lower[T];

  double d0 = linearInterpolation(u[X], gradient4D(getPermutation(AAA),   relative[0],   relative[1],   relative[2],   relative[3]),
                                        gradient4D(getPermutation(BAA),   relative[0]-1, relative[1],   relative[2],   relative[3]));

  double d1 = linearInterpolation(u[X], gradient4D(getPermutation(ABA),   relative[0],   relative[1]-1, relative[2],   relative[3]),
                                        gradient4D(getPermutation(BBA),   relative[0]-1, relative[1]-1, relative[2],   relative[3]));

  double d2 = linearInterpolation(u[X], gradient4D(getPermutation(AAB),   relative[0],   relative[1],   relative[2]-1, relative[3]),
                                        gradient4D(getPermutation(BAB),   relative[0]-1, relative[1],   relative[2]-1, relative[3]));
    
  double d3 = linearInterpolation(u[X], gradient4D(getPermutation(ABB),   relative[0],   relative[1]-1, relative[2]-1, relative[3]),
                                        gradient4D(getPermutation(BBB),   relative[0]-1, relative[1]-1, relative[2]-1, relative[3]));

  double d4 = linearInterpolation(u[X], gradient4D(getPermutation(AAA+1), relative[0],   relative[1],   relative[2],   relative[3]-1),
                                        gradient4D(getPermutation(BAA+1), relative[0]-1, relative[1],   relative[2],   relative[3]-1));
    
  double d5 = linearInterpolation(u[X], gradient4D(getPermutation(ABA+1), relative[0],   relative[1]-1, relative[2],   relative[3]-1),
                                        gradient4D(getPermutation(BBA+1), relative[0]-1, relative[1]-1, relative[2],   relative[3]-1));

  double d6 = linearInterpolation(u[X], gradient4D(getPermutation(AAB+1), relative[0],   relative[1],   relative[2]-1, relative[3]-1),
                                        gradient4D(getPermutation(BAB+1), relative[0]-1, relative[1],   relative[2]-1, relative[3]-1));
      
  double d7 = linearInterpolation(u[X], gradient4D(getPermutation(ABB+1), relative[0],   relative[1]-1, relative[2]-1, relative[3]-1),
                                        gradient4D(getPermutation(BBB+1), relative[0]-1, relative[1]-1, relative[2]-1, relative[3]-1));


  return linearInterpolation(u[T], linearInterpolation(u[Z], linearInterpolation(u[Y], d0, d1),
                                                             linearInterpolation(u[Y], d2, d3)),
                                   linearInterpolation(u[Z], linearInterpolation(u[Y], d4, d5),
                                                             linearInterpolation(u[Y], d6, d7)));
}

