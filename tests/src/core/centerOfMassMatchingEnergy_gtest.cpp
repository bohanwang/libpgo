#include <gtest/gtest.h>

#include "centerOfMassMatchingEnergy.h"
#include "pgoLogging.h"
#include "potentialEnergies.h"
#include "tetMeshGeo.h"

#include <cmath>
#include <memory>
#include <vector>

namespace
{
namespace ES = pgo::EigenSupport;
}

TEST(CenterOfMassMatchingEnergyGTest, ExposesAllVertexDOFsSoAggregationMatchesDirectEvaluation)
{
  pgo::Logging::init();
  const double vertices[] = { 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1 };
  const int tets[] = { 0, 1, 2, 3, 1, 2, 3, 4 };
  const pgo::Mesh::TetMeshGeo mesh(5, vertices, 2, tets);

  const auto energy = std::make_shared<pgo::PredefinedPotentialEnergies::CenterOfMassMatchingEnergy>(
    mesh, ES::V3d(0.4, 0.3, 0.2), ES::M3d::Identity());
  EXPECT_EQ(energy->getNumDOFs(), 15);
  std::vector<int> dofs;
  energy->getDOFs(dofs);
  ASSERT_EQ(dofs.size(), 15u);
  for (int i = 0; i < 15; ++i)
    EXPECT_EQ(dofs[i], i);
  ES::SpMatD hessian;
  energy->createHessian(hessian);
  EXPECT_EQ(hessian.rows(), 15);

  ES::VXd x(15);
  for (int i = 0; i < 15; ++i)
    x[i] = vertices[i] + 0.05 * std::sin(1.0 + i);
  const double direct = energy->func(x);
  ES::VXd directGradient = ES::VXd::Zero(15);
  energy->gradient(x, directGradient);
  EXPECT_TRUE(std::isfinite(direct));
  EXPECT_GT(directGradient.norm(), 0.0);

  pgo::NonlinearOptimization::PotentialEnergies aggregate(15);
  aggregate.addPotentialEnergy(energy);
  aggregate.init();
  EXPECT_DOUBLE_EQ(aggregate.func(x), direct);
  ES::VXd aggregateGradient = ES::VXd::Zero(15);
  aggregate.gradient(x, aggregateGradient);
  EXPECT_LT((aggregateGradient - directGradient).norm(), 1e-12 * std::max(1.0, directGradient.norm()));
}
