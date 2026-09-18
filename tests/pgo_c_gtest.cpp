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

TEST(PgoCTest, UpdatesTetMeshVerticesWithoutChangingSourceMesh)
{
  pgo_init();

  double vertices[] = {
    0.0, 0.0, 0.0,
    1.0, 0.0, 0.0,
    0.0, 1.0, 0.0,
    0.0, 0.0, 1.0,
  };
  int elements[] = { 0, 1, 2, 3 };
  double translatedVertices[] = {
    1.0, 2.0, 3.0,
    2.0, 2.0, 3.0,
    1.0, 3.0, 3.0,
    1.0, 2.0, 4.0,
  };

  auto tetmesh = pgo_create_tetmesh(4, vertices, 1, elements, 1.0e5, 0.45, 1000.0);
  ASSERT_NE(tetmesh, nullptr);

  auto updatedTetmesh = pgo_tetmesh_update_vertices(tetmesh, translatedVertices);
  ASSERT_NE(updatedTetmesh, nullptr);
  ASSERT_NE(updatedTetmesh, tetmesh);

  double sourceVertices[12] = {};
  double updatedVertices[12] = {};
  pgo_tetmesh_get_vertices(tetmesh, sourceVertices);
  pgo_tetmesh_get_vertices(updatedTetmesh, updatedVertices);

  for (int i = 0; i < 12; ++i) {
    EXPECT_DOUBLE_EQ(sourceVertices[i], vertices[i]);
    EXPECT_DOUBLE_EQ(updatedVertices[i], translatedVertices[i]);
  }

  int updatedElements[4] = {};
  pgo_tetmesh_get_elements(updatedTetmesh, updatedElements);
  for (int i = 0; i < 4; ++i)
    EXPECT_EQ(updatedElements[i], elements[i]);

  pgo_destroy_tetmesh(updatedTetmesh);
  pgo_destroy_tetmesh(tetmesh);
}
