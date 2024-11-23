#pragma once

// Perlin wind
// Code authors: Jernej Barbic, Somya Sharma, USC, 2012

#include "perlinNoise.h"

namespace pgo
{
namespace PerlinNoise
{

class PerlinWind
{
public:
  // direction is 3D vector
  // randomMagnitude is a scalar
  // the wind will be a sum of a constant wind with the given direction and magnitude, and a randomized spacetime 4D perlin noise
  PerlinWind(double *direction, double randomMagnitude);
  virtual ~PerlinWind();

  void SetDirection(double *direction);  // direction is 3D vector
  void SetRandomMagnitude(double randomMagnitude);

  // x,y,z,t,persistence must lie on the interval [0,1]
  void GetPerlinWind(double x, double y, double z, double t, int numFrequencies, double persistence, double *value);

protected:
  double randomMagnitude;
  double direction[3];
  double GetPerlinNoise1D(double x, int numFrequencies, double persistence);
  double GetPerlinNoise3D(double *v, int dof, int numFrequencies, double persistence);
  double GetPerlinNoise4D(double *v, int dof, int numFrequencies, double persistence);
  PerlinNoise *perlinNoise[4];
};
}  // namespace PerlinNoise
}  // namespace pgo