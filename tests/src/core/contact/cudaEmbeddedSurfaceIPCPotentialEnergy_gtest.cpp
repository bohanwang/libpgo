#include "ipc/cuda/cudaEmbeddedSurfaceIPCPotentialEnergy.h"
#include "ipc/cuda/cudaIPCGeometry.h"
#include "ipc/cuda/cudaSurfaceIPCBroadPhase.h"
#include "embeddedSurfaceIPCPotentialEnergy.h"
#include "cudaDeviceBuffer.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <array>
#include <limits>
#include <random>
#include <set>

bool cudaIPCIntervalCertificateForTest(const double *positions, const double *motion,
  const double *error, int kind, double firstTime);

namespace
{
namespace ES = pgo::EigenSupport;
using pgo::Contact::CIPC::CudaEmbeddedSurfaceIPCPotentialEnergy;
using pgo::Contact::CIPC::EmbeddedSurfaceIPCPotentialEnergy;
using pgo::Contact::CIPC::SurfaceIPCCore;
using pgo::CudaUtilities::DeviceBuffer;

class CudaEmbeddedSurfaceIPCPotentialEnergyGTest : public ::testing::Test
{
protected:
  void SetUp() override
  {
    int count = 0;
    if (cudaGetDeviceCount(&count) != cudaSuccess || count == 0)
      GTEST_SKIP() << "A CUDA device is required.";
  }
};

struct Mesh
{
  ES::MXd vertices;
  ES::MXi faces;
  ES::SpMatD map;
  std::vector<uint8_t> mask;
};

Mesh makeMesh(bool coupled, bool obstacle, double scale = 1.0)
{
  Mesh mesh;
  mesh.vertices.resize(6, 3);
  mesh.vertices << 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0,
    0.0, 1.0, 0.0,
    0.12, 0.08, 0.045,
    0.86, 0.13, 0.055,
    0.09, 0.82, 0.05;
  mesh.vertices *= scale;
  mesh.faces.resize(2, 3);
  mesh.faces << 0, 1, 2, 3, 4, 5;
  mesh.mask.assign(6, uint8_t{ 1 });
  if (obstacle)
    std::fill(mesh.mask.begin(), mesh.mask.begin() + 3, uint8_t{ 0 });
  const int n = obstacle ? 9 : 18;
  mesh.map.resize(18, n);
  std::vector<ES::TripletD> entries;
  for (int row = obstacle ? 9 : 0; row < 18; ++row) {
    const int col = obstacle ? row - 9 : row;
    if (coupled) {
      entries.emplace_back(row, col, 0.7);
      entries.emplace_back(row, (col + 3) % n, 0.3);
    }
    else
      entries.emplace_back(row, col, 1.0);
  }
  mesh.map.setFromTriplets(entries.begin(), entries.end());
  return mesh;
}

void expectNear(const ES::VXd &a, const ES::VXd &b, double tolerance)
{
  EXPECT_LE((a - b).norm(), tolerance * std::max(1e-12, b.norm()));
}

void compare(const Mesh &mesh, double epsEE, bool project, double scale)
{
  SurfaceIPCCore::Parameters params;
  params.dhat = 0.1 * scale;
  params.kappa = 2.3;
  params.eps_ee = epsEE;
  params.projectHessianToPSD = project;
  EmbeddedSurfaceIPCPotentialEnergy cpu(mesh.vertices, mesh.faces, mesh.map, mesh.mask, params);
  CudaEmbeddedSurfaceIPCPotentialEnergy gpu(mesh.vertices, mesh.faces, mesh.map, mesh.mask, params);
  const int n = gpu.getNumDOFs();
  ES::VXd x = ES::VXd::Zero(n), direction(n), baseline(n), result(n), gradient(n);
  for (int i = 0; i < n; ++i) {
    x[i] = scale * 0.0001 * std::sin(0.23 * i);
    direction[i] = std::cos(0.7 * i + 0.4);
    baseline[i] = 0.013 * std::sin(0.3 * i);
  }
  DeviceBuffer<double> dx, dv, out;
  dx.upload(x.data(), n);
  dv.upload(direction.data(), n);
  const double referenceEnergy = cpu.func(x);
  EXPECT_NEAR(gpu.func(dx.data()), referenceEnergy, 2e-10 * std::max(1e-14, std::abs(referenceEnergy)));
  EXPECT_GT(gpu.getNumActivePairs(), 0u);
  cpu.gradient(x, gradient);
  out.upload(baseline.data(), n);
  gpu.addGradient(dx.data(), out.data());
  out.download(result.data(), n);
  expectNear(result, gradient + baseline, 2e-9);

  ES::SpMatD hessian;
  cpu.hessianDirect(x, hessian);
  gpu.prepareHessian(dx.data());
  out.upload(baseline.data(), n);
  gpu.addHessianVector(dv.data(), out.data());
  out.download(result.data(), n);
  expectNear(result, hessian * direction + baseline, 5e-8);
  out.upload(baseline.data(), n);
  gpu.addDiagonal(out.data());
  out.download(result.data(), n);
  const ES::VXd diagonal = hessian.diagonal();
  expectNear(result, diagonal + baseline, 5e-8);

  // Line search may evaluate a different x and CCD without invalidating the
  // cached Hessian used by a subsequent linear solve.
  const ES::VXd trial = x + 1e-5 * scale * direction;
  dx.upload(trial.data(), n);
  gpu.func(dx.data());
  gpu.computeMaxStepSize(dx.data(), dv.data());
  out.zero();
  gpu.addHessianVector(dv.data(), out.data());
  out.download(result.data(), n);
  expectNear(result, hessian * direction, 5e-8);

  const double h = 1e-7 * scale;
  ES::VXd perturbed = x + h * direction;
  dx.upload(perturbed.data(), n);
  const double plus = gpu.func(dx.data());
  perturbed = x - h * direction;
  dx.upload(perturbed.data(), n);
  const double minus = gpu.func(dx.data());
  EXPECT_NEAR((plus - minus) / (2 * h), gradient.dot(direction),
    2e-6 * std::max(1e-12, std::abs(gradient.dot(direction))));
}

void compareBroadPhase(const ES::MXd &vertices, const ES::MXi &faces,
  const std::vector<uint8_t> &mask, const ES::MXd *motion, double padding)
{
  using namespace pgo::Contact::CIPC;
  std::vector<int> triangles, edges;
  std::set<std::array<int, 2>> uniqueEdges;
  for (int i = 0; i < faces.rows(); ++i) {
    for (int j = 0; j < 3; ++j) {
      triangles.push_back(faces(i, j));
      const int a = faces(i, j), b = faces(i, (j + 1) % 3);
      uniqueEdges.insert({ std::min(a, b), std::max(a, b) });
    }
  }
  for (const auto &edge : uniqueEdges)
    edges.insert(edges.end(), edge.begin(), edge.end());
  std::vector<double> positions(3 * vertices.rows()), motions(positions.size());
  for (int i = 0; i < vertices.rows(); ++i)
    for (int j = 0; j < 3; ++j) {
      positions[3 * i + j] = vertices(i, j);
      motions[3 * i + j] = motion ? (*motion)(i, j) : 0.0;
    }
  DeviceBuffer<int> deviceTriangles, deviceEdges;
  DeviceBuffer<uint8_t> deviceMask;
  DeviceBuffer<double> devicePositions, deviceMotions;
  deviceTriangles.upload(triangles.data(), triangles.size());
  deviceEdges.upload(edges.data(), edges.size());
  deviceMask.upload(mask.data(), mask.size());
  devicePositions.upload(positions.data(), positions.size());
  deviceMotions.upload(motions.data(), motions.size());
  CudaSurfaceIPCBroadPhase gpu({ deviceTriangles.data(), deviceEdges.data(), deviceMask.data(),
    static_cast<int>(vertices.rows()), static_cast<int>(faces.rows()), static_cast<int>(edges.size() / 2) });

  // Independent exhaustive oracle, with the original per-pair bound formula.
  auto overlaps = [&](const std::array<int, 4> &indices, int split) {
    bool moving = false;
    for (int i : indices)
      moving |= mask[i] != 0;
    if (!moving)
      return false;
    for (int axis = 0; axis < 3; ++axis) {
      double lo[2] = { INFINITY, INFINITY }, hi[2] = { -INFINITY, -INFINITY };
      for (int i = 0; i < 4; ++i) {
        const double p = positions[3 * indices[i] + axis];
        const double q = motion ? p + motions[3 * indices[i] + axis] : p;
        const int side = i < split ? 0 : 1;
        lo[side] = std::min(lo[side], std::min(p, q));
        hi[side] = std::max(hi[side], std::max(p, q));
      }
      if (std::nextafter(lo[0] - padding, -INFINITY) > std::nextafter(hi[1] + padding, INFINITY) ||
        std::nextafter(lo[1] - padding, -INFINITY) > std::nextafter(hi[0] + padding, INFINITY))
        return false;
    }
    return true;
  };
  std::set<std::array<int, 3>> expected;
  for (int point = 0; point < vertices.rows(); ++point)
    for (int triangle = 0; triangle < faces.rows(); ++triangle) {
      std::array<int, 4> indices{ point, faces(triangle, 0), faces(triangle, 1), faces(triangle, 2) };
      if (point != indices[1] && point != indices[2] && point != indices[3] && overlaps(indices, 1))
        expected.insert({ 0, point, triangle });
    }
  const int edgeCount = static_cast<int>(edges.size() / 2);
  for (int a = 0; a < edgeCount; ++a)
    for (int b = a + 1; b < edgeCount; ++b) {
      std::array<int, 4> indices{ edges[2 * a], edges[2 * a + 1], edges[2 * b], edges[2 * b + 1] };
      if (indices[0] != indices[2] && indices[0] != indices[3] && indices[1] != indices[2] && indices[1] != indices[3] && overlaps(indices, 2))
        expected.insert({ 1, a, b });
    }
  // Rebuilding the same object also exercises retained allocations and resets.
  for (int repeat = 0; repeat < 2; ++repeat) {
    gpu.detect(devicePositions.data(), motion ? deviceMotions.data() : nullptr, padding);
    std::vector<CudaSurfaceIPCCandidate> candidates(gpu.getNumCandidates());
    if (!candidates.empty())
      pgo::CudaUtilities::checkCuda(cudaMemcpy(candidates.data(), gpu.getDeviceCandidates(),
                                      sizeof(CudaSurfaceIPCCandidate) * candidates.size(), cudaMemcpyDeviceToHost),
        "test candidate download");
    std::set<std::array<int, 3>> actual;
    for (const auto &candidate : candidates)
      actual.insert({ candidate.kind, candidate.first, candidate.second });
    EXPECT_EQ(actual.size(), candidates.size()) << "Duplicate candidates";
    EXPECT_EQ(actual, expected);
  }
}
}  // namespace

TEST_F(CudaEmbeddedSurfaceIPCPotentialEnergyGTest, MatchesCpuBarrierEmbeddingAndProjectedHessian)
{
  for (bool obstacle : { false, true })
    for (bool coupled : { false, true })
      for (bool project : { false, true })
        compare(makeMesh(coupled, obstacle), 0.0, project, 1.0);
}

TEST_F(CudaEmbeddedSurfaceIPCPotentialEnergyGTest, MatchesCpuMollifierAndSmallScale)
{
  compare(makeMesh(true, true), 0.1, true, 1.0);
  compare(makeMesh(false, true, 1e-5), 0.0, true, 1e-5);
}

TEST_F(CudaEmbeddedSurfaceIPCPotentialEnergyGTest, SweptCcdPreventsCrossingFixedTriangle)
{
  const Mesh mesh = makeMesh(false, true);
  SurfaceIPCCore::Parameters params;
  params.dhat = 0.1;
  CudaEmbeddedSurfaceIPCPotentialEnergy gpu(mesh.vertices, mesh.faces, mesh.map, mesh.mask, params);
  ES::VXd x = ES::VXd::Zero(9), motion = ES::VXd::Zero(9);
  DeviceBuffer<double> dx, dm;
  dx.upload(x.data(), 9);
  dm.upload(motion.data(), 9);
  EXPECT_DOUBLE_EQ(gpu.computeMaxStepSize(dx.data(), dm.data()), 1.0);
  for (int v = 0; v < 3; ++v)
    motion[3 * v + 2] = -0.2;
  dm.upload(motion.data(), 9);
  const double alpha = gpu.computeMaxStepSize(dx.data(), dm.data());
  EXPECT_GT(alpha, 0.0);
  EXPECT_LT(alpha, 0.045 / 0.2);
  EmbeddedSurfaceIPCPotentialEnergy cpu(mesh.vertices, mesh.faces, mesh.map, mesh.mask, params);
  EXPECT_NO_THROW(cpu.validateCollisionFreeState(alpha * motion));
}

TEST_F(CudaEmbeddedSurfaceIPCPotentialEnergyGTest, RejectsUnsupportedFrictionAndExternalMapRows)
{
  Mesh mesh = makeMesh(false, true);
  SurfaceIPCCore::Parameters params;
  params.frictionCoeff = 0.1;
  EXPECT_THROW(CudaEmbeddedSurfaceIPCPotentialEnergy(mesh.vertices, mesh.faces, mesh.map, mesh.mask, params), std::invalid_argument);
  params.frictionCoeff = 0.0;
  mesh.map.coeffRef(0, 0) = 1.0;
  EXPECT_THROW(CudaEmbeddedSurfaceIPCPotentialEnergy(mesh.vertices, mesh.faces, mesh.map, mesh.mask, params), std::invalid_argument);
}

TEST_F(CudaEmbeddedSurfaceIPCPotentialEnergyGTest, GrowsCandidateBufferWithoutDroppingPairs)
{
  constexpr int sheets = 20;
  ES::MXd vertices(3 * sheets, 3);
  ES::MXi faces(sheets, 3);
  for (int i = 0; i < sheets; ++i) {
    vertices.row(3 * i) << 0.0, 0.0, 0.001 * i;
    vertices.row(3 * i + 1) << 1.0, 0.0, 0.001 * i;
    vertices.row(3 * i + 2) << 0.0, 1.0, 0.001 * i;
    faces.row(i) << 3 * i, 3 * i + 1, 3 * i + 2;
  }
  ES::SpMatD map(9 * sheets, 9 * sheets);
  map.setIdentity();
  SurfaceIPCCore::Parameters params;
  EmbeddedSurfaceIPCPotentialEnergy cpu(vertices, faces, map, params);
  CudaEmbeddedSurfaceIPCPotentialEnergy gpu(vertices, faces, map, params);
  const ES::VXd x = ES::VXd::Zero(map.cols());
  DeviceBuffer<double> deviceX;
  deviceX.upload(x.data(), x.size());
  const double energy = cpu.func(x);
  EXPECT_NEAR(gpu.func(deviceX.data()), energy, 1e-9 * std::abs(energy));
  EXPECT_GT(gpu.getNumActivePairs(), 1024u);
}

TEST_F(CudaEmbeddedSurfaceIPCPotentialEnergyGTest, NearParallelEdgesAndRigidTranslationHaveSafeCcd)
{
  Mesh mesh = makeMesh(false, false);
  mesh.vertices.row(3) << 0.0, 0.0, 0.05;
  mesh.vertices.row(4) << 1.0, 1e-10, 0.05;
  mesh.vertices.row(5) << 0.0, 1.0, 0.05;
  CudaEmbeddedSurfaceIPCPotentialEnergy gpu(mesh.vertices, mesh.faces, mesh.map);
  ES::VXd x = ES::VXd::Zero(18), motion(18);
  for (int vertex = 0; vertex < 6; ++vertex)
    motion.segment<3>(3 * vertex) << 0.1, 0.2, -0.2;
  DeviceBuffer<double> deviceX, deviceMotion;
  deviceX.upload(x.data(), x.size());
  deviceMotion.upload(motion.data(), motion.size());
  EXPECT_DOUBLE_EQ(gpu.computeMaxStepSize(deviceX.data(), deviceMotion.data()), 1.0);
  motion.head(9).setZero();
  motion.tail(9).setZero();
  for (int vertex = 3; vertex < 6; ++vertex)
    motion[3 * vertex + 2] = -0.2;
  deviceMotion.upload(motion.data(), motion.size());
  const double alpha = gpu.computeMaxStepSize(deviceX.data(), deviceMotion.data());
  EXPECT_GT(alpha, 0.0);
  EXPECT_LT(alpha, 0.25);
}

TEST_F(CudaEmbeddedSurfaceIPCPotentialEnergyGTest, CcdUsesGeometricDistanceForCoplanarInteriorEdges)
{
  // The triangles cross only through edge interiors: no vertex lies inside the
  // other triangle. The barrier's coplanar PE convention gives positive edge
  // distances here, so CCD must use an independent segment-distance routine.
  for (double scale : { 1.0, 1e-12, 1e-32 }) {
    ES::MXd vertices(6, 3);
    vertices << -1, 0, 0, 1, 0, 0, -1, -0.1, 0,
      0, -1, 0, 0, 1, 0, 0.1, 1, 0;
    vertices *= scale;
    ES::MXi faces(2, 3);
    faces << 0, 1, 2, 3, 4, 5;
    ES::SpMatD map(18, 18);
    map.setIdentity();
    CudaEmbeddedSurfaceIPCPotentialEnergy gpu(vertices, faces, map);
    ES::VXd x = ES::VXd::Zero(18), motion = ES::VXd::Zero(18);
    DeviceBuffer<double> deviceX, deviceMotion;
    deviceX.upload(x.data(), 18);
    deviceMotion.upload(motion.data(), 18);
    EXPECT_DOUBLE_EQ(gpu.computeMaxStepSize(deviceX.data(), deviceMotion.data()), 0.0);
    for (int vertex = 3; vertex < 6; ++vertex)
      motion[3 * vertex + 2] = 0.1 * scale;
    deviceMotion.upload(motion.data(), 18);
    EXPECT_DOUBLE_EQ(gpu.computeMaxStepSize(deviceX.data(), deviceMotion.data()), 0.0);

    for (int vertex = 3; vertex < 6; ++vertex) {
      x[3 * vertex + 2] = 0.1 * scale;
      motion[3 * vertex + 2] = -0.2 * scale;
    }
    deviceX.upload(x.data(), 18);
    deviceMotion.upload(motion.data(), 18);
    const double alpha = gpu.computeMaxStepSize(deviceX.data(), deviceMotion.data());
    EXPECT_GT(alpha, 0.0);
    EXPECT_LT(alpha, 0.5);
  }
}

TEST_F(CudaEmbeddedSurfaceIPCPotentialEnergyGTest, PartialCcdRetainsInitialGapForNewtonConditioning)
{
  Mesh mesh = makeMesh(false, false);
  mesh.vertices.bottomRows(3).col(2).setConstant(0.05);
  CudaEmbeddedSurfaceIPCPotentialEnergy gpu(mesh.vertices, mesh.faces, mesh.map);
  ES::VXd x = ES::VXd::Zero(18), motion = ES::VXd::Zero(18);
  for (int vertex = 3; vertex < 6; ++vertex)
    motion[3 * vertex + 2] = -0.2;
  DeviceBuffer<double> deviceX, deviceMotion;
  deviceX.upload(x.data(), x.size());
  deviceMotion.upload(motion.data(), motion.size());
  const double alpha = gpu.computeMaxStepSize(deviceX.data(), deviceMotion.data());
  // The first collision is at .25. A partial query must retain 10% of the
  // initial gap instead of returning the nearly singular final CCD iterate.
  EXPECT_GT(alpha, 0.20);
  EXPECT_LE(alpha, 0.225000000001);
  EXPECT_GE(0.05 - 0.2 * alpha, 0.005 - 1e-12);
  // An entirely clear motion still receives the complete certificate.
  motion *= 0.1;
  deviceMotion.upload(motion.data(), motion.size());
  EXPECT_DOUBLE_EQ(gpu.computeMaxStepSize(deviceX.data(), deviceMotion.data()), 1.0);
}

TEST_F(CudaEmbeddedSurfaceIPCPotentialEnergyGTest, OriginalChordCertifiesLongTangentialNearMiss)
{
  Mesh mesh = makeMesh(false, false);
  mesh.vertices.row(3) << -2, 0.1, 1e-8;
  mesh.vertices.row(4) << -1.5, 0.1, 1e-8;
  mesh.vertices.row(5) << -2, 0.6, 1e-8;
  CudaEmbeddedSurfaceIPCPotentialEnergy gpu(mesh.vertices, mesh.faces, mesh.map);
  ES::VXd x = ES::VXd::Zero(18), motion = ES::VXd::Zero(18);
  for (int vertex = 3; vertex < 6; ++vertex)
    motion[3 * vertex] = 4.0;
  DeviceBuffer<double> deviceX, deviceMotion;
  deviceX.upload(x.data(), x.size());
  deviceMotion.upload(motion.data(), motion.size());
  // A speed-only lower bound would need many more than 1000 iterations near
  // this tiny constant gap. A separating axis across the original whole
  // chord proves the complete near miss without subdivision or rebasing.
  EXPECT_DOUBLE_EQ(gpu.computeMaxStepSize(deviceX.data(), deviceMotion.data()), 1.0);
}

TEST_F(CudaEmbeddedSurfaceIPCPotentialEnergyGTest, IntervalFallbackCertifiesDragonAdvancementExhaustion)
{
  // Actual PT [2593,1894,1950,1964] from cubic dragon frame 410. The clear
  // near-grazing chord exhausted 1000 additive steps at t=.998948, with a
  // remaining gap 2.85e-8. It is not a collision or mapping-error uncertainty.
  const double positions[12] = {
    -1.3902367855748896, -0.08987065473271055, 0.12782760605539445,
    -1.3919291709570423, -0.08665921233977701, 0.12753054256679588,
    -1.3890982715932283, -0.09206302614812611, 0.13222600107054389,
    -1.3890513553960604, -0.09209218027811888, 0.12419292794348799
  };
  const double motion[12] = {
    0.0009322397330936909, -0.001132779925460271, -0.0008704657583936122,
    0.0009029940722515738, -0.0011118090485401824, -0.0009274181412626267,
    0.0009183970362710596, -0.0011170059262299774, -0.0008721057685721573,
    0.0009632825813015974, -0.0011342655988609454, -0.0009047861266704289
  };
  const double error[12] = {
    5.383058763828114e-15, 1.3417369674235698e-15, 5.88095360143924e-16,
    5.385811096929352e-15, 8.827288708451997e-16, 1.0285444592246692e-15,
    5.3752488321162655e-15, 8.563999121358938e-16, 1.0301592814113028e-15,
    5.3752622718681095e-15, 8.785625122625168e-16, 1.0044769501073484e-15
  };
  EXPECT_TRUE(cudaIPCIntervalCertificateForTest(positions, motion, error, 0, 0.0));
  EXPECT_TRUE(cudaIPCIntervalCertificateForTest(positions, motion, error, 0, 0.99894801678691));
  ES::MXd vertices(4, 3);
  for (int i = 0; i < 4; ++i)
    vertices.row(i) = Eigen::Map<const ES::V3d>(positions + 3 * i).transpose();
  ES::MXi faces(1, 3);
  faces << 1, 2, 3;
  ES::SpMatD map(12, 12);
  map.setIdentity();
  CudaEmbeddedSurfaceIPCPotentialEnergy gpu(vertices, faces, map);
  DeviceBuffer<double> deviceX(12), deviceMotion;
  deviceX.zero();
  deviceMotion.upload(motion, 12);
  EXPECT_DOUBLE_EQ(gpu.computeMaxStepSize(deviceX.data(), deviceMotion.data()), 1.0);

  // A nearby genuinely crossing chord must still fail the very same fallback.
  ES::VXd crossing = Eigen::Map<const ES::VXd>(motion, 12);
  const ES::V3d edge0 = (vertices.row(2) - vertices.row(1)).transpose();
  const ES::V3d edge1 = (vertices.row(3) - vertices.row(1)).transpose();
  const ES::V3d normal = edge0.cross(edge1).normalized();
  const double signedGap = normal.dot((vertices.row(0) - vertices.row(1)).transpose());
  crossing.head<3>() -= std::copysign(1e-6, signedGap) * normal;
  EXPECT_FALSE(cudaIPCIntervalCertificateForTest(positions, crossing.data(), error, 0, 0.0));
  deviceMotion.upload(crossing.data(), crossing.size());
  EXPECT_LT(gpu.computeMaxStepSize(deviceX.data(), deviceMotion.data()), 1.0);
}

TEST_F(CudaEmbeddedSurfaceIPCPotentialEnergyGTest, IntervalFallbackRejectsCrossingsAndUncertainSeparation)
{
  const double error[12] = {};
  const double pt[12] = { 0.25, 0.25, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0 };
  double motion[12] = { 0, 0, -2 };
  // Absolute endpoint distances are both positive. Signed support hulls
  // must nevertheless detect the through-plane interval at t=.5.
  EXPECT_FALSE(cudaIPCIntervalCertificateForTest(pt, motion, error, 0, 0.0));
  motion[2] = 0;
  motion[0] = 2;
  EXPECT_TRUE(cudaIPCIntervalCertificateForTest(pt, motion, error, 0, 0.0));
  double uncertain[12];
  std::fill(std::begin(uncertain), std::end(uncertain), 2.0);
  // Fixed error intervals overlap at the first endpoint. Subdivision cannot
  // remove this uncertainty and must terminate without a false certificate.
  EXPECT_FALSE(cudaIPCIntervalCertificateForTest(pt, motion, uncertain, 0, 0.0));
  const double ee[12] = { -1, 0, 0, 1, 0, 0, 0, -1, 1, 0, 1, 1 };
  std::fill(std::begin(motion), std::end(motion), 0.0);
  motion[8] = motion[11] = -2;
  EXPECT_FALSE(cudaIPCIntervalCertificateForTest(ee, motion, error, 1, 0.0));
  motion[8] = motion[11] = 0;
  motion[6] = motion[9] = 4;
  EXPECT_TRUE(cudaIPCIntervalCertificateForTest(ee, motion, error, 1, 0.0));
}

TEST_F(CudaEmbeddedSurfaceIPCPotentialEnergyGTest, RepeatedCcdCannotCrossDragonPointTriangleRoot)
{
  // Extracted from the original dragon-tet frame 431. This vertex enters and
  // leaves a triangle at t=.11627438567556836 and .9304462759850638, so both
  // endpoints are clear. Rounded positive distances previously let repeated
  // conservative-advancement suffixes jump across both roots and certify one.
  ES::MXd vertices(6, 3), velocity(6, 3);
  vertices.row(0) << -1.267687175222354, -0.27383300068135347, 0.15184045565443113;
  vertices.row(1) = vertices.row(0) + ES::V3d(1e-6, 0, 0).transpose();
  vertices.row(2) = vertices.row(0) + ES::V3d(0, 1e-6, 0).transpose();
  vertices.row(3) << -1.272781006433854, -0.27206091024380585, 0.1576676778371424;
  vertices.row(4) << -1.2550081601772165, -0.2813765309687974, 0.15556259714647316;
  vertices.row(5) << -1.2675543945013559, -0.2719720561437912, 0.14050790887572215;
  velocity.row(0) << 0.0002529060242464354, -0.002683482965247508, -0.0013030811841920065;
  velocity.row(1) = velocity.row(2) = velocity.row(0);
  velocity.row(3) << -0.0009267287516814271, -0.0018603615414542474, -0.0015532419446279144;
  velocity.row(4) << -0.0014158529183845303, -0.0020562010885734616, -0.0017674687776073061;
  velocity.row(5) << -0.0010244743031080361, -0.0017888450058192173, -0.0016116669923125215;
  ES::MXi faces(2, 3);
  faces << 0, 1, 2, 3, 4, 5;
  ES::SpMatD mapping(18, 18);
  mapping.setIdentity();
  for (double restOffset : { 0.0, 1e6 }) {
    ES::MXd rest = vertices.array() + restOffset;
    CudaEmbeddedSurfaceIPCPotentialEnergy gpu(rest, faces, mapping);
    ES::VXd start = ES::VXd::Constant(18, -restOffset), delta(18);
    for (int j = 0; j < 6; ++j)
      delta.segment<3>(3 * j) = velocity.row(j).transpose();
    const ES::VXd endpoint = start + delta;
    DeviceBuffer<double> deviceX, deviceMotion;
    double progress = 0.0;
    bool certified = false;
    for (int iteration = 0; iteration < 128; ++iteration) {
      ES::VXd x = start + progress * delta, motion = endpoint - x;
      deviceX.upload(x.data(), x.size());
      deviceMotion.upload(motion.data(), motion.size());
      const double bound = gpu.computeMaxStepSize(deviceX.data(), deviceMotion.data());
      ASSERT_TRUE(std::isfinite(bound));
      ASSERT_GE(bound, 0.0);
      if (bound == 1.0) {
        certified = true;
        break;
      }
      const double next = std::nextafter(std::fma(bound, 1.0 - progress, progress), progress);
      if (!(next > progress))
        break;
      progress = next;
    }
    EXPECT_FALSE(certified);
    EXPECT_LT(progress, 0.11627438567556836);
  }
}

TEST_F(CudaEmbeddedSurfaceIPCPotentialEnergyGTest, CcdHandlesCollapsedMovingPrimitives)
{
  const Mesh mesh = makeMesh(false, false);
  CudaEmbeddedSurfaceIPCPotentialEnergy gpu(mesh.vertices, mesh.faces, mesh.map);
  for (bool collapseToPoint : { false, true }) {
    ES::VXd x = ES::VXd::Zero(18), motion = ES::VXd::Zero(18);
    // Collapse one triangle to an edge, then to a point. A vertex of the other
    // triangle moves through that primitive at a known positive time.
    x.segment<3>(6) = mesh.vertices.row(0).transpose() - mesh.vertices.row(2).transpose();
    if (collapseToPoint)
      x.segment<3>(3) = mesh.vertices.row(0).transpose() - mesh.vertices.row(1).transpose();
    for (int vertex = 3; vertex < 6; ++vertex) {
      x.segment<3>(3 * vertex) = mesh.vertices.row(0).transpose() + ES::V3d(0, 0, 0.05) - mesh.vertices.row(vertex).transpose();
      motion[3 * vertex + 2] = -0.2;
    }
    DeviceBuffer<double> deviceX, deviceMotion;
    deviceX.upload(x.data(), 18);
    deviceMotion.upload(motion.data(), 18);
    const double alpha = gpu.computeMaxStepSize(deviceX.data(), deviceMotion.data());
    EXPECT_GT(alpha, 0.0);
    EXPECT_LT(alpha, 0.25);
  }
}

TEST_F(CudaEmbeddedSurfaceIPCPotentialEnergyGTest, BvhCandidatesMatchExhaustiveStaticAndSweptBounds)
{
  constexpr int count = 57, triangleCount = 44;
  ES::MXd vertices(count, 3), motion(count, 3);
  ES::MXi faces(triangleCount, 3);
  std::mt19937 random(81723);
  std::uniform_real_distribution<double> uniform(-1.0, 1.0);
  for (int i = 0; i < count; ++i)
    for (int j = 0; j < 3; ++j) {
      vertices(i, j) = uniform(random);
      motion(i, j) = 0.7 * uniform(random);
    }
  for (int i = 0; i < triangleCount; ++i)
    faces.row(i) << i, (i + 7) % count, (i + 19) % count;
  for (int ownership = 0; ownership < 3; ++ownership) {
    std::vector<uint8_t> mask(count);
    for (int i = 0; i < count; ++i)
      mask[i] = ownership == 0 ? 1 : ownership == 1 ? i % 3 == 0 :
                                                      0;
    compareBroadPhase(vertices, faces, mask, nullptr, 0.03);
    compareBroadPhase(vertices, faces, mask, &motion, 0.0);
  }
  // Broad-phase decisions remain exact at translated and tiny scales.
  vertices.array() += 1e8;
  compareBroadPhase(vertices, faces, std::vector<uint8_t>(count, 1), &motion, 1e-8);
  vertices.array() -= 1e8;
  vertices *= 1e-10;
  motion *= 1e-10;
  compareBroadPhase(vertices, faces, std::vector<uint8_t>(count, 1), &motion, 1e-12);
}

TEST_F(CudaEmbeddedSurfaceIPCPotentialEnergyGTest, BvhEmptySingletonAndCoincidentBoundsAreComplete)
{
  ES::MXd vertices = ES::MXd::Zero(4, 3);
  ES::MXi faces(0, 3);
  compareBroadPhase(vertices, faces, std::vector<uint8_t>(4, 1), nullptr, 0.0);
  vertices << 0, 0, 0, 1, 0, 0, 0, 1, 0, 0.2, 0.2, 0;
  faces.resize(1, 3);
  faces << 0, 1, 2;
  compareBroadPhase(vertices, faces, std::vector<uint8_t>(4, 1), nullptr, 0.0);
  // All Morton codes tie. Unique primitive IDs must still yield a valid tree.
  constexpr int triangles = 30;
  vertices = ES::MXd::Zero(3 * triangles, 3);
  faces.resize(triangles, 3);
  for (int i = 0; i < triangles; ++i)
    faces.row(i) << 3 * i, 3 * i + 1, 3 * i + 2;
  compareBroadPhase(vertices, faces, std::vector<uint8_t>(vertices.rows(), 1), nullptr, 0.0);
}

TEST_F(CudaEmbeddedSurfaceIPCPotentialEnergyGTest, SweptBoundsIncludeSurfaceMappingUncertainty)
{
  Mesh mesh = makeMesh(false, true);
  for (int vertex = 3; vertex < 6; ++vertex)
    mesh.vertices(vertex, 2) = 1e6;
  CudaEmbeddedSurfaceIPCPotentialEnergy gpu(mesh.vertices, mesh.faces, mesh.map, mesh.mask);
  DeviceBuffer<double> deviceX, deviceMotion(9);
  deviceMotion.zero();
  for (double gap : { 1e-10, 1e-6 }) {
    ES::VXd x = ES::VXd::Zero(9);
    for (int vertex = 0; vertex < 3; ++vertex)
      x[3 * vertex + 2] = -1e6 + gap;
    deviceX.upload(x.data(), x.size());
    // The small gap lies inside the rest+W*u arithmetic uncertainty. Its
    // nominal AABBs are disjoint, so the broad phase must expand them too.
    EXPECT_DOUBLE_EQ(gpu.computeMaxStepSize(deviceX.data(), deviceMotion.data()), gap < 1e-9 ? 0.0 : 1.0);
  }
}

TEST_F(CudaEmbeddedSurfaceIPCPotentialEnergyGTest, NonfiniteSurfaceStateFailsClosed)
{
  const Mesh mesh = makeMesh(false, true);
  CudaEmbeddedSurfaceIPCPotentialEnergy gpu(mesh.vertices, mesh.faces, mesh.map, mesh.mask);
  ES::VXd x = ES::VXd::Zero(9);
  x[0] = std::numeric_limits<double>::quiet_NaN();
  DeviceBuffer<double> deviceX, gradient(9), motion(9);
  deviceX.upload(x.data(), x.size());
  gradient.zero();
  motion.zero();
  EXPECT_TRUE(std::isinf(gpu.func(deviceX.data())));
  EXPECT_THROW(gpu.addGradient(deviceX.data(), gradient.data()), std::domain_error);
  EXPECT_THROW(gpu.computeMaxStepSize(deviceX.data(), motion.data()), std::domain_error);
  // Squared primitive lengths overflow at these finite coordinates. Scaled
  // support axes must still reject an actual touching point in the fixed
  // triangle, independently of the two distant endpoints.
  x.setZero();
  x[3] = 1e200;
  x[6] = -1e200;
  x[2] = -mesh.vertices(3, 2);  // The first moving vertex touches the floor.
  deviceX.upload(x.data(), x.size());
  EXPECT_DOUBLE_EQ(gpu.computeMaxStepSize(deviceX.data(), motion.data()), 0.0);
}

TEST_F(CudaEmbeddedSurfaceIPCPotentialEnergyGTest, NearParallelInteriorContactUsesCpuFeatureConvention)
{
  for (double angle : { 0.5e-6, 0.9e-6, 1.1e-6, 2e-6 }) {
    Mesh mesh = makeMesh(false, false);
    mesh.vertices << -0.5, 0, 0, 0.5, 0, 0, -0.2, -0.8, 0,
      -0.5, -0.5 * angle, 0.05, 0.5, 0.5 * angle, 0.05, -0.2, 0.8, 0.05;
    SurfaceIPCCore::Parameters params;
    params.dhat = 0.1;
    params.kappa = 2.3;
    EmbeddedSurfaceIPCPotentialEnergy cpu(mesh.vertices, mesh.faces, mesh.map, params);
    CudaEmbeddedSurfaceIPCPotentialEnergy gpu(mesh.vertices, mesh.faces, mesh.map, params);
    ES::VXd x = ES::VXd::Zero(18), expected(18), actual(18);
    DeviceBuffer<double> deviceX, gradient(18);
    deviceX.upload(x.data(), x.size());
    gradient.zero();
    const double energy = cpu.func(x);
    EXPECT_NEAR(gpu.func(deviceX.data()), energy, 1e-8 * energy);
    cpu.gradient(x, expected);
    gpu.addGradient(deviceX.data(), gradient.data());
    gradient.download(actual.data(), actual.size());
    expectNear(actual, expected, 1e-7);
  }
}

TEST_F(CudaEmbeddedSurfaceIPCPotentialEnergyGTest, EnergyChangesResolveSmallStepsAndChangingContactSupport)
{
  for (double mollifier : { 0.0, 0.1 }) {
    const Mesh mesh = makeMesh(true, true);
    SurfaceIPCCore::Parameters params;
    params.dhat = 0.1;
    params.kappa = 2.3;
    params.eps_ee = mollifier;
    CudaEmbeddedSurfaceIPCPotentialEnergy gpu(mesh.vertices, mesh.faces, mesh.map, mesh.mask, params);
    ES::VXd base = ES::VXd::Zero(9), trial(9), direction(9), gradient(9), before(9), after(9);
    for (int i = 0; i < 9; ++i)
      direction[i] = std::sin(0.3 * i + 0.8);
    DeviceBuffer<double> deviceBase, deviceTrial, deviceDirection, result(9);
    deviceBase.upload(base.data(), base.size());
    deviceDirection.upload(direction.data(), direction.size());
    result.zero();
    gpu.addGradient(deviceBase.data(), result.data());
    result.download(gradient.data(), gradient.size());
    gpu.prepareHessian(deviceBase.data());
    result.zero();
    gpu.addHessianVector(deviceDirection.data(), result.data());
    result.download(before.data(), before.size());
    EXPECT_DOUBLE_EQ(gpu.energyChange(deviceBase.data(), deviceBase.data()).value, 0.0);
    for (double step : { 1e-7, 1e-11 }) {
      trial = step * direction;
      deviceTrial.upload(trial.data(), trial.size());
      const auto change = gpu.energyChange(deviceBase.data(), deviceTrial.data());
      const auto reverse = gpu.energyChange(deviceTrial.data(), deviceBase.data());
      const double predicted = step * gradient.dot(direction);
      EXPECT_NEAR(change.value, predicted, 5e-5 * std::abs(predicted));
      EXPECT_NEAR(change.value, -reverse.value, 1e-9 * std::abs(change.value));
      EXPECT_GE(change.roundoffScale, std::abs(change.value));
    }
    // Contact may be absent at the base and present only at the trial (or the
    // reverse). A base-only active set would silently lose this energy change.
    trial.setZero();
    for (int vertex = 0; vertex < 3; ++vertex)
      trial[3 * vertex + 2] = 0.2;
    deviceTrial.upload(trial.data(), trial.size());
    const double activeEnergy = gpu.func(deviceBase.data());
    EXPECT_DOUBLE_EQ(gpu.func(deviceTrial.data()), 0.0);
    EXPECT_NEAR(gpu.energyChange(deviceTrial.data(), deviceBase.data()).value, activeEnergy, 2e-12 * activeEnergy);
    EXPECT_NEAR(gpu.energyChange(deviceBase.data(), deviceTrial.data()).value, -activeEnergy, 2e-12 * activeEnergy);
    result.zero();
    gpu.addHessianVector(deviceDirection.data(), result.data());
    result.download(after.data(), after.size());
    expectNear(after, before, 1e-12);
  }
}

TEST_F(CudaEmbeddedSurfaceIPCPotentialEnergyGTest, EnergyChangesPreserveSubUlpWorldMotion)
{
  Mesh mesh = makeMesh(true, true);
  mesh.vertices.array() += 1024.0;
  SurfaceIPCCore::Parameters params;
  params.dhat = 0.1;
  params.kappa = 2.3;
  CudaEmbeddedSurfaceIPCPotentialEnergy gpu(mesh.vertices, mesh.faces, mesh.map, mesh.mask, params);
  ES::VXd base = ES::VXd::Zero(9), trial(9), gradient(9);
  for (int i = 0; i < 9; ++i)
    trial[i] = 1e-14 * std::sin(0.3 * i + 0.8);
  DeviceBuffer<double> deviceBase, deviceTrial, result(9);
  deviceBase.upload(base.data(), 9);
  deviceTrial.upload(trial.data(), 9);
  result.zero();
  gpu.addGradient(deviceBase.data(), result.data());
  result.download(gradient.data(), 9);
  const double predicted = gradient.dot(trial);
  EXPECT_NE(predicted, 0.0);
  EXPECT_NEAR(gpu.energyChange(deviceBase.data(), deviceTrial.data()).value, predicted, 1e-10 * std::abs(predicted));
  EXPECT_DOUBLE_EQ(gpu.energyChange(deviceBase.data(), deviceBase.data()).value, 0.0);
  EXPECT_DOUBLE_EQ(gpu.energyChange(deviceTrial.data(), deviceTrial.data()).value, 0.0);
}

TEST(CudaIPCGeometryGTest, NearParallelEnvelopeDerivativesRemainConsistent)
{
  using V3 = Eigen::Matrix<long double, 3, 1>;
  using V12 = Eigen::Matrix<long double, 12, 1>;
  using M12 = Eigen::Matrix<long double, 12, 12>;
  for (long double angle : { 0.2L, 1e-4L, 1e-7L }) {
    V3 points[4] = { V3(0, 0, 0), V3(1, 0, 0), V3(0.3L, -0.4L * angle, -0.02L), V3(1.3L, 0.6L * angle, -0.02L) };
    long double distance;
    V12 gradient;
    M12 hessian, numerical;
    pgo::Contact::CIPC::CudaIPCGeometry::interiorEdgeEdge(points[0], points[1], points[2], points[3], distance, &gradient, &hessian);
    EXPECT_NEAR(static_cast<double>(distance), 0.0004, 1e-18);
    const long double weights[4] = { 0.3L, 0.7L, -0.6L, -0.4L };
    for (int i = 0; i < 4; ++i)
      EXPECT_NEAR(static_cast<double>(gradient[3 * i + 2]), static_cast<double>(0.04L * weights[i]), 1e-15);
    const long double step = std::min(1e-7L, angle * angle * 1e-3L);
    for (int column = 0; column < 12; ++column) {
      V3 plus[4], minus[4];
      for (int i = 0; i < 4; ++i)
        plus[i] = minus[i] = points[i];
      plus[column / 3][column % 3] += step;
      minus[column / 3][column % 3] -= step;
      V12 gp, gm;
      pgo::Contact::CIPC::CudaIPCGeometry::interiorEdgeEdge(plus[0], plus[1], plus[2], plus[3], distance, &gp, static_cast<M12 *>(nullptr));
      pgo::Contact::CIPC::CudaIPCGeometry::interiorEdgeEdge(minus[0], minus[1], minus[2], minus[3], distance, &gm, static_cast<M12 *>(nullptr));
      numerical.col(column) = (gp - gm) / (2 * step);
    }
    EXPECT_LE(static_cast<double>((numerical - hessian).norm() / hessian.norm()), 2e-4);
    EXPECT_LE(static_cast<double>((hessian - hessian.transpose()).norm()), 1e-12);
  }
}
