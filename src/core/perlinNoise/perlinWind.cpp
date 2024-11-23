#include "perlinWind.h"

#include <cstdio>
#include <iostream>
#include <cmath>

using namespace pgo::PerlinNoise;

PerlinWind::PerlinWind(double * direction_, double randomMagnitude_)
{
  direction[0] = direction_[0];
  direction[1] = direction_[1];
  direction[2] = direction_[2];
  randomMagnitude = randomMagnitude_;

  for(int i=0; i<4; i++)
    perlinNoise[i] = new PerlinNoise();
}

PerlinWind::~PerlinWind()
{
  for(int i=0; i<4; i++)
    delete(perlinNoise[i]);
}

void PerlinWind::SetDirection(double * direction_)
{
  direction[0] = direction_[0];
  direction[1] = direction_[1];
  direction[2] = direction_[2];
}

void PerlinWind::SetRandomMagnitude(double randomMagnitude_)
{
  randomMagnitude = randomMagnitude_;
}

double PerlinWind::GetPerlinNoise1D(double x, int numFrequencies, double persistence)
{
  double noise = 0.0;
  double frequency = 1.0;
  double amplitude = persistence;

  for(int i=0; i<numFrequencies - 1; i++)
  {
    noise += perlinNoise[3]->noise1D(x * frequency) * amplitude;
    frequency *= 2;
    amplitude *= persistence;
  }

  return noise;
}

double PerlinWind::GetPerlinNoise3D(double * v, int dof, int numFrequencies, double persistence)
{
  double noise = 0.0;
  double frequency = 1.0;
  double amplitude = persistence;
  double vf[3] = {0, 0, 0};
  for(int i=0; i<numFrequencies - 1; i++)
  {
    vf[0] = v[0] * frequency;
    vf[1] = v[1] * frequency;
    vf[2] = v[2] * frequency;
    noise += perlinNoise[dof]->noise3D(vf) * amplitude;
    frequency *= 2;
    amplitude *= persistence;
  }

  return noise;
}

double PerlinWind::GetPerlinNoise4D(double * v, int dof, int numFrequencies, double persistence)
{
  double noise = 0.0;
  double frequency = 1.0;
  double amplitude = persistence;
  double vf[4] = {0, 0, 0, 0};
  for(int i=0; i<numFrequencies - 1; i++)
  {
    vf[0] = v[0] * frequency;
    vf[1] = v[1] * frequency;
    vf[2] = v[2] * frequency;
    vf[3] = v[3] * frequency;
    noise += perlinNoise[dof]->noise4D(vf) * amplitude;
    frequency *= 2;
    amplitude *= persistence;
  }

  return noise;
}

void PerlinWind::GetPerlinWind(double x, double y, double z, double t, int numFrequencies, double persistence, double * value)
{
  double v[4] = {x, y, z, t};
  value[0] = direction[0] + randomMagnitude * GetPerlinNoise4D(v, 0, numFrequencies, persistence);
  value[1] = direction[1] + randomMagnitude * GetPerlinNoise4D(v, 1, numFrequencies, persistence);
  value[2] = direction[2] + randomMagnitude * GetPerlinNoise4D(v, 2, numFrequencies, persistence);
}

