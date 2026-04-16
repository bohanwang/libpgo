#include <gtest/gtest.h>

#include "pgo_c.h"

namespace
{
constexpr const char *kTorusVegPath = LIBPGO_TEST_TORUS_VEG;
}

TEST(PgoCTest, LoadsTetMeshFromExampleFile)
{
  pgo_init();

  auto tetmesh = pgo_create_tetmeshgeo_from_file(const_cast<char *>(kTorusVegPath));
  ASSERT_NE(tetmesh, nullptr);
  EXPECT_GT(pgo_tetmeshgeo_get_num_vertices(tetmesh), 0);
  EXPECT_GT(pgo_tetmeshgeo_get_num_tets(tetmesh), 0);

  pgo_destroy_tetmeshgeo(tetmesh);
}

TEST(PgoCTest, CreatesSmoothRSEnergyWhenSupported)
{
  pgo_init();

  auto tetmesh = pgo_create_tetmeshgeo_from_file(const_cast<char *>(kTorusVegPath));
  ASSERT_NE(tetmesh, nullptr);

  auto energy = pgo_create_smooth_rs_energy(tetmesh, 0.0, 1.0);
#if defined(PGO_HAS_MKL)
  ASSERT_NE(energy, nullptr);
  pgo_destroy_smooth_rs_energy(energy);
#else
  EXPECT_EQ(energy, nullptr);
#endif

  pgo_destroy_tetmeshgeo(tetmesh);
}
