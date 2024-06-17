#include "pgo_c.h"

int main(int argc, char* argv[])
{
  pgo_init();

  pgoTetMeshGeoStructHandle tetmesh = pgo_create_tetmeshgeo_from_file(argv[1]);
  pgoSmoothRSEnergyStructHandle energy = pgo_create_smooth_rs_energy(tetmesh, 0.0, 1.0);

  return 0;
}