/*
copyright to Bohan Wang
*/

#include "cudaSurfaceIPCFrictionPotentialEnergy.h"
#include "cudaSurfaceIPCBroadPhase.h"
#include "cudaIPCGeometry.h"
#include "cudaDeviceBuffer.h"
#include "ipc/topology/surfaceIPCTopology.h"

#include <libuipc/distance/distance_flagged.h>
#include <libuipc/distance/edge_edge_mollifier.h>

#include <climits>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace pgo::Contact::CIPC
{
namespace
{
using CudaUtilities::DeviceBuffer;
using CudaUtilities::checkCuda;
using Vec3 = Eigen::Matrix<double, 3, 1>;
using Vec12 = Eigen::Matrix<double, 12, 1>;
using Mat3 = Eigen::Matrix<double, 3, 3>;
namespace distance = uipc::backend::cuda::distance;
constexpr int blockSize = 128;
int blocks(int n)
{
  return n / blockSize + (n % blockSize != 0);
}
struct Map
{
  const int *offsets, *columns;
  const double *weights;
  int rows;
};
struct Contact
{
  int vertices[4];
  double coefficients[4];
  double projector[9];
  double force;
};
__device__ Vec3 position(const double *x, int v)
{
  return Eigen::Map<const Vec3>(x + 3 * v);
}

// Map a state or a represented displacement difference. Computing differences
// before adding rest coordinates preserves small frictional slips.
__global__ void mapKernel(Map map, const double *x, const double *reference,
  const double *rest, double *out, int *invalid)
{
  const int row = blockIdx.x * blockDim.x + threadIdx.x;
  if (row >= map.rows)
    return;
  double value = rest ? rest[row] : 0.0;
  for (int k = map.offsets[row]; k < map.offsets[row + 1]; ++k) {
    const int col = map.columns[k];
    value += map.weights[k] * (reference ? x[col] - reference[col] : x[col]);
  }
  out[row] = value;
  if (invalid && !isfinite(value))
    atomicExch(invalid, 1);
}
__global__ void transposeKernel(Map map, const double *x, double *out)
{
  const int row = blockIdx.x * blockDim.x + threadIdx.x;
  if (row >= map.rows)
    return;
  for (int k = map.offsets[row]; k < map.offsets[row + 1]; ++k)
    atomicAdd(out + map.columns[k], map.weights[k] * x[row]);
}

__global__ void freezeKernel(const CudaSurfaceIPCCandidate *pairs, int count,
  const int *triangles, const int *edges, const double *vertexArea,
  const double *triangleArea, const double *edgeLength, const double *x,
  double dhat2, double kappa, double epsEE, double mu,
  Contact *contacts, int *contactCount, int *invalid)
{
  const int index = blockIdx.x * blockDim.x + threadIdx.x;
  if (index >= count)
    return;
  const auto pair = pairs[index];
  Contact contact;
  double weight;
  if (pair.kind == 0) {
    contact.vertices[0] = pair.first;
    for (int i = 0; i < 3; ++i)
      contact.vertices[i + 1] = triangles[3 * pair.second + i];
    weight = vertexArea[pair.first] * triangleArea[pair.second];
  }
  else {
    for (int i = 0; i < 2; ++i) {
      contact.vertices[i] = edges[2 * pair.first + i];
      contact.vertices[i + 2] = edges[2 * pair.second + i];
    }
    weight = edgeLength[pair.first] * edgeLength[pair.second];
  }
  const Vec3 a = position(x, contact.vertices[0]), b = position(x, contact.vertices[1]);
  const Vec3 c = position(x, contact.vertices[2]), d = position(x, contact.vertices[3]);
  double d2;
  Vec12 gradient;
  Vec3 separation;
  if (pair.kind == 0) {
    const auto flag = distance::point_triangle_distance_flag(a, b, c, d);
    distance::point_triangle_distance2(flag, a, b, c, d, d2);
    if (d2 >= dhat2)
      return;
    distance::point_triangle_distance2_gradient(flag, a, b, c, d, gradient);
    separation = gradient.head<3>() / 2.0;
  }
  else {
    const auto flag = distance::edge_edge_distance_flag(a, b, c, d);
    if (flag.sum() == 4)
      CudaIPCGeometry::interiorEdgeEdge(a, b, c, d, d2, &gradient, static_cast<Eigen::Matrix<double, 12, 12> *>(nullptr));
    else {
      distance::edge_edge_distance2(flag, a, b, c, d, d2);
      distance::edge_edge_distance2_gradient(flag, a, b, c, d, gradient);
    }
    if (d2 >= dhat2)
      return;
    separation = (gradient.head<3>() + gradient.segment<3>(3)) / 2.0;
    if (epsEE > 0.0) {
      double mollifier;
      distance::edge_edge_mollifier(a, b, c, d, epsEE, mollifier);
      weight *= mollifier;
    }
  }
  const double separation2 = separation.squaredNorm();
  bool finiteGradient = true;
  for (int i = 0; i < 12; ++i)
    finiteGradient = finiteGradient && isfinite(gradient[i]);
  if (!isfinite(d2) || d2 <= 0.0 || !isfinite(separation2) || separation2 <= 0.0 || !finiteGradient) {
    atomicExch(invalid, 1);
    return;
  }
  for (int i = 0; i < 4; ++i) {
    contact.coefficients[i] = gradient.segment<3>(3 * i).dot(separation) / (2.0 * separation2);
    if (!isfinite(contact.coefficients[i])) {
      atomicExch(invalid, 1);
      return;
    }
  }
  const Vec3 normal = separation / sqrt(separation2);
  Eigen::Map<Mat3>(contact.projector) = Mat3::Identity() - normal * normal.transpose();
  const double r = d2 / dhat2, rm1 = r - 1.0;
  const double derivative = -(rm1 * (2.0 * log(r) + rm1 / r)) / dhat2;
  contact.force = mu * (-2.0 * sqrt(d2) * weight * kappa * derivative);
  bool finiteProjector = true;
  for (int i = 0; i < 9; ++i)
    finiteProjector = finiteProjector && isfinite(contact.projector[i]);
  if (!isfinite(contact.force) || contact.force < 0.0 || !finiteProjector) {
    atomicExch(invalid, 1);
    return;
  }
  if (contact.force > 0.0)
    contacts[atomicAdd(contactCount, 1)] = contact;
}
__device__ Vec3 slip(const Contact &contact, const double *motion)
{
  Vec3 relative = Vec3::Zero();
  for (int i = 0; i < 4; ++i)
    relative += contact.coefficients[i] * position(motion, contact.vertices[i]);
  return Eigen::Map<const Mat3>(contact.projector) * relative;
}

template<int Mode>
__global__ void evaluateKernel(const Contact *contacts, int count, const double *motion,
  double epsilon, double *energy, double *gradient, double *hessians)
{
  const int index = blockIdx.x * blockDim.x + threadIdx.x;
  if (index >= count)
    return;
  const Contact contact = contacts[index];
  const Vec3 r = slip(contact, motion);
  const double length = r.norm(), t = length / epsilon;
  if constexpr (Mode == 0) {
    const double value = length < epsilon ? epsilon * (t * t * (1.0 - t / 3.0) + 1.0 / 3.0) : length;
    atomicAdd(energy, contact.force * value);
  }
  if constexpr (Mode == 1) {
    const double factor = length < epsilon ? (2.0 - t) / epsilon : 1.0 / length;
    const Vec3 force = contact.force * factor * (Eigen::Map<const Mat3>(contact.projector) * r);
    for (int i = 0; i < 4; ++i)
      for (int axis = 0; axis < 3; ++axis)
        atomicAdd(gradient + 3 * contact.vertices[i] + axis, contact.coefficients[i] * force[axis]);
  }
  if constexpr (Mode == 2) {
    Mat3 h;
    if (length == 0.0)
      h = (2.0 / epsilon) * Mat3::Identity();
    else {
      const Vec3 direction = r / length;
      const Mat3 outer = direction * direction.transpose();
      h = length < epsilon ? Mat3(((2.0 - t) * Mat3::Identity() - t * outer) / epsilon) : Mat3((Mat3::Identity() - outer) / length);
    }
    Eigen::Map<Mat3>(hessians + 9 * index) = contact.force * Eigen::Map<const Mat3>(contact.projector) * h * Eigen::Map<const Mat3>(contact.projector);
  }
}

__global__ void changeKernel(const Contact *contacts, int count, const double *baseMotion,
  const double *difference, double epsilon, double *result)
{
  const int index = blockIdx.x * blockDim.x + threadIdx.x;
  if (index >= count)
    return;
  const Contact contact = contacts[index];
  const Vec3 r = slip(contact, baseMotion), dr = slip(contact, difference);
  if (dr.isZero(0.0))
    return;
  const double l0 = r.norm(), l1 = (r + dr).norm();
  const double squaredChange = (2.0 * r + dr).dot(dr);
  const double dl = l0 + l1 > 0.0 ? squaredChange / (l0 + l1) : 0.0;
  // Preserve uncancelled dot-product terms in the error scale when the
  // tangential increment is nearly orthogonal to the existing slip.
  const double dlScale = l0 + l1 > 0.0 ?
    (2.0 * r.cwiseAbs() + dr.cwiseAbs()).dot(dr.cwiseAbs()) / (l0 + l1) :
    0.0;
  double value, scale;
  if (l0 < epsilon && l1 < epsilon) {
    const double term1 = (l0 + l1) * dl / epsilon;
    const double term2 = (l0 * l0 + l0 * l1 + l1 * l1) * dl / (3.0 * epsilon * epsilon);
    value = term1 - term2;
    scale = ((l0 + l1) / epsilon +
              (l0 * l0 + l0 * l1 + l1 * l1) / (3.0 * epsilon * epsilon)) *
      dlScale;
  }
  else {
    // f0(l)=l-(l-epsilon)^3/(3 epsilon^2) on the regularized side.
    const double q0 = fmin(0.0, l0 - epsilon), q1 = fmin(0.0, l1 - epsilon);
    const double term = (q1 * q1 * q1 - q0 * q0 * q0) / (3.0 * epsilon * epsilon);
    value = dl - term;
    scale = dlScale + fabs(term);
  }
  atomicAdd(result, contact.force * value);
  atomicAdd(result + 1, contact.force * scale);
}
__global__ void hessianVectorKernel(const Contact *contacts, int count, const double *hessians,
  const double *v, double *out)
{
  const int index = blockIdx.x * blockDim.x + threadIdx.x;
  if (index >= count)
    return;
  const Contact contact = contacts[index];
  Vec3 relative = Vec3::Zero();
  for (int i = 0; i < 4; ++i)
    relative += contact.coefficients[i] * position(v, contact.vertices[i]);
  const Vec3 product = Eigen::Map<const Mat3>(hessians + 9 * index) * relative;
  for (int i = 0; i < 4; ++i)
    for (int a = 0; a < 3; ++a)
      atomicAdd(out + 3 * contact.vertices[i] + a, contact.coefficients[i] * product[a]);
}
__global__ void diagonalKernel(const Contact *contacts, int count, const double *hessians, Map map, double *out)
{
  const int index = blockIdx.x * blockDim.x + threadIdx.x;
  if (index >= count)
    return;
  const Contact contact = contacts[index];
  const Eigen::Map<const Mat3> h(hessians + 9 * index);
  // Keep all W^T H W cross terms when different surface coordinates share a DOF.
  for (int i = 0; i < 12; ++i) {
    const int rowI = 3 * contact.vertices[i / 3] + i % 3;
    for (int ki = map.offsets[rowI]; ki < map.offsets[rowI + 1]; ++ki) {
      double sum = 0.0;
      for (int j = 0; j < 12; ++j) {
        const int rowJ = 3 * contact.vertices[j / 3] + j % 3;
        for (int kj = map.offsets[rowJ]; kj < map.offsets[rowJ + 1]; ++kj)
          if (map.columns[ki] == map.columns[kj])
            sum += map.weights[kj] * contact.coefficients[j / 3] * h(i % 3, j % 3);
      }
      atomicAdd(out + map.columns[ki], map.weights[ki] * contact.coefficients[i / 3] * sum);
    }
  }
}
}  // namespace

class CudaSurfaceIPCFrictionPotentialEnergy::Impl
{
public:
  SurfaceIPCCore::Parameters parameters;
  SurfaceIPCTopology topology;
  EigenSupport::MXd hostRest;
  int numDOFs = 0, numContacts = 0;
  double epsilon = 0.0;
  bool ready = false, hessianReady = false;
  DeviceBuffer<int> triangles, edges, offsets, columns, count, invalid;
  DeviceBuffer<uint8_t> mask;
  DeviceBuffer<double> rest, weights, vertexArea, triangleArea, edgeLength;
  DeviceBuffer<double> reference, lagged, motion, difference, scratch, scalar, hessians;
  DeviceBuffer<Contact> contacts;
  std::unique_ptr<CudaSurfaceIPCBroadPhase> broadPhase;
  Map map() const { return { offsets.data(), columns.data(), weights.data(), 3 * topology.numVerts }; }
  void requireReady() const
  {
    if (!ready)
      throw std::logic_error("CUDA IPC friction must be updated before evaluation.");
  }
  void mapMotion(const double *x, const double *ref, double *out)
  {
    invalid.zero();
    mapKernel<<<blocks(map().rows), blockSize>>>(map(), x, ref, nullptr, out, invalid.data());
    checkCuda(cudaGetLastError(), "CUDA friction displacement mapping");
  }
  bool valid() const
  {
    int bad;
    invalid.download(&bad, 1);
    return bad == 0;
  }
  void uploadRestMeasures()
  {
    std::vector<double> values(3 * topology.numVerts);
    for (int v = 0; v < topology.numVerts; ++v)
      for (int a = 0; a < 3; ++a)
        values[3 * v + a] = hostRest(v, a);
    rest.upload(values.data(), values.size());
    vertexArea.upload(topology.vertexArea.data(), topology.vertexArea.size());
    triangleArea.upload(topology.triArea.data(), topology.triArea.size());
    edgeLength.upload(topology.edgeLength.data(), topology.edgeLength.size());
  }
};

CudaSurfaceIPCFrictionPotentialEnergy::CudaSurfaceIPCFrictionPotentialEnergy(
  const EigenSupport::MXd &vertices, const EigenSupport::MXi &faces,
  const EigenSupport::SpMatD &mapping, const std::vector<uint8_t> &mask,
  const SurfaceIPCCore::Parameters &parameters):
  impl_(std::make_unique<Impl>())
{
  SurfaceIPCCore validator(parameters);
  validator.setMesh(vertices, faces, mask);
  if (mapping.rows() != 3 * vertices.rows() || mapping.cols() <= 0 ||
    mapping.rows() > INT_MAX || mapping.cols() > INT_MAX || mapping.nonZeros() > INT_MAX)
    throw std::invalid_argument("CUDA friction displacement map has invalid dimensions.");
  auto &data = *impl_;
  data.parameters = parameters;
  data.hostRest = vertices;
  data.topology.setMesh(vertices, faces, mask);
  data.numDOFs = static_cast<int>(mapping.cols());
  Eigen::SparseMatrix<double, Eigen::RowMajor, int> rowMap = mapping;
  rowMap.makeCompressed();
  for (int row = 0; row < rowMap.rows(); ++row)
    for (decltype(rowMap)::InnerIterator entry(rowMap, row); entry; ++entry) {
      if (!std::isfinite(entry.value()) || (!mask[row / 3] && entry.value() != 0.0))
        throw std::invalid_argument("CUDA friction requires finite map entries and zero external map rows.");
    }
  data.offsets.upload(rowMap.outerIndexPtr(), rowMap.rows() + 1);
  data.columns.upload(rowMap.innerIndexPtr(), rowMap.nonZeros());
  data.weights.upload(rowMap.valuePtr(), rowMap.nonZeros());
  std::vector<int> indices;
  for (const auto &face : data.topology.triangles)
    indices.insert(indices.end(), face.begin(), face.end());
  data.triangles.upload(indices.data(), indices.size());
  indices.clear();
  for (const auto &edge : data.topology.edges)
    indices.insert(indices.end(), edge.begin(), edge.end());
  data.edges.upload(indices.data(), indices.size());
  data.mask.upload(mask.data(), mask.size());
  data.uploadRestMeasures();
  data.reference.resize(data.numDOFs);
  data.lagged.resize(mapping.rows());
  data.motion.resize(mapping.rows());
  data.difference.resize(mapping.rows());
  data.scratch.resize(mapping.rows());
  data.scalar.resize(2);
  data.count.resize(1);
  data.invalid.resize(1);
  data.broadPhase = std::make_unique<CudaSurfaceIPCBroadPhase>(CudaSurfaceIPCTopologyView{
    data.triangles.data(), data.edges.data(), data.mask.data(), data.topology.numVerts,
    static_cast<int>(data.topology.triangles.size()), static_cast<int>(data.topology.edges.size()) });
}
CudaSurfaceIPCFrictionPotentialEnergy::~CudaSurfaceIPCFrictionPotentialEnergy() = default;
int CudaSurfaceIPCFrictionPotentialEnergy::getNumDOFs() const
{
  return impl_->numDOFs;
}
std::size_t CudaSurfaceIPCFrictionPotentialEnergy::getNumFrictionPairs() const
{
  return impl_->numContacts;
}
void CudaSurfaceIPCFrictionPotentialEnergy::clear()
{
  impl_->ready = impl_->hessianReady = false;
  impl_->numContacts = 0;
}
void CudaSurfaceIPCFrictionPotentialEnergy::updateFriction(const double *reference, const double *lagged, double timestep)
{
  clear();
  auto &data = *impl_;
  const double epsilon = data.parameters.frictionEpsV * timestep;
  if (!std::isfinite(timestep) || timestep <= 0.0 || !std::isfinite(epsilon) || epsilon <= 0.0)
    throw std::invalid_argument("CUDA friction requires finite positive timestep and epsv*timestep.");
  data.epsilon = epsilon;
  checkCuda(cudaMemcpyAsync(data.reference.data(), reference, sizeof(double) * data.numDOFs, cudaMemcpyDeviceToDevice), "CUDA friction reference copy");
  data.mapMotion(reference, reference, data.motion.data());
  if (!data.valid())
    throw std::invalid_argument("CUDA friction requires a finite reference state.");
  data.invalid.zero();
  mapKernel<<<blocks(data.map().rows), blockSize>>>(data.map(), lagged, nullptr, data.rest.data(), data.lagged.data(), data.invalid.data());
  checkCuda(cudaGetLastError(), "CUDA friction lagged mapping");
  if (!data.valid())
    throw std::invalid_argument("CUDA friction requires a finite lagged state.");
  if (data.parameters.frictionCoeff != 0.0) {
    data.broadPhase->detect(data.lagged.data(), nullptr, data.parameters.dhat);
    const int candidates = static_cast<int>(data.broadPhase->getNumCandidates());
    data.contacts.resize(candidates);
    data.count.zero();
    data.invalid.zero();
    if (candidates) {
      freezeKernel<<<blocks(candidates), blockSize>>>(data.broadPhase->getDeviceCandidates(), candidates,
        data.triangles.data(), data.edges.data(), data.vertexArea.data(), data.triangleArea.data(), data.edgeLength.data(),
        data.lagged.data(), data.parameters.dhat * data.parameters.dhat, data.parameters.kappa, data.parameters.eps_ee,
        data.parameters.frictionCoeff, data.contacts.data(), data.count.data(), data.invalid.data());
      checkCuda(cudaGetLastError(), "CUDA friction lagged contact preparation");
    }
    if (!data.valid())
      throw std::runtime_error("CUDA friction requires finite strictly separated lagged contacts.");
    data.count.download(&data.numContacts, 1);
  }
  data.ready = true;
}
void CudaSurfaceIPCFrictionPotentialEnergy::updateExternalSurface(int first, const EigenSupport::MXd &positions)
{
  auto &data = *impl_;
  if (positions.cols() != 3 || first < 0 || positions.rows() > data.topology.numVerts - first || !positions.allFinite())
    throw std::invalid_argument("CUDA friction external update requires finite positions in range.");
  for (int i = 0; i < positions.rows(); ++i)
    if (data.topology.isVertexDeformable(first + i))
      throw std::invalid_argument("CUDA friction external update addresses a deformable vertex.");
  clear();
  data.hostRest.middleRows(first, positions.rows()) = positions;
  EigenSupport::VXd flat(3 * data.topology.numVerts);
  for (int v = 0; v < data.topology.numVerts; ++v)
    flat.segment<3>(3 * v) = data.hostRest.row(v).transpose();
  data.topology.updateMeasures(flat);
  data.uploadRestMeasures();
}
double CudaSurfaceIPCFrictionPotentialEnergy::func(const double *x)
{
  auto &data = *impl_;
  data.requireReady();
  data.mapMotion(x, data.reference.data(), data.motion.data());
  if (!data.valid())
    return std::numeric_limits<double>::infinity();
  data.scalar.zero();
  if (data.numContacts)
    evaluateKernel<0><<<blocks(data.numContacts), blockSize>>>(data.contacts.data(), data.numContacts, data.motion.data(), data.epsilon, data.scalar.data(), nullptr, nullptr);
  checkCuda(cudaGetLastError(), "CUDA friction energy");
  double result;
  data.scalar.download(&result, 1);
  return result;
}
NonlinearOptimization::CudaPotentialEnergy::EnergyChange CudaSurfaceIPCFrictionPotentialEnergy::energyChange(const double *base, const double *trial)
{
  auto &data = *impl_;
  data.requireReady();
  data.mapMotion(base, data.reference.data(), data.motion.data());
  if (!data.valid())
    return { std::numeric_limits<double>::infinity(), 0.0 };
  data.mapMotion(trial, base, data.difference.data());
  if (!data.valid())
    return { std::numeric_limits<double>::infinity(), 0.0 };
  data.scalar.zero();
  if (data.numContacts)
    changeKernel<<<blocks(data.numContacts), blockSize>>>(data.contacts.data(), data.numContacts, data.motion.data(), data.difference.data(), data.epsilon, data.scalar.data());
  checkCuda(cudaGetLastError(), "CUDA friction local energy change");
  double result[2];
  data.scalar.download(result, 2);
  return { result[0], result[1] };
}
void CudaSurfaceIPCFrictionPotentialEnergy::addGradient(const double *x, double *out)
{
  auto &data = *impl_;
  data.requireReady();
  data.mapMotion(x, data.reference.data(), data.motion.data());
  if (!data.valid())
    throw std::invalid_argument("CUDA friction gradient requires a finite state.");
  data.scratch.zero();
  if (data.numContacts)
    evaluateKernel<1><<<blocks(data.numContacts), blockSize>>>(data.contacts.data(), data.numContacts, data.motion.data(), data.epsilon, nullptr, data.scratch.data(), nullptr);
  transposeKernel<<<blocks(data.map().rows), blockSize>>>(data.map(), data.scratch.data(), out);
  checkCuda(cudaGetLastError(), "CUDA friction mapped gradient");
}
void CudaSurfaceIPCFrictionPotentialEnergy::prepareHessian(const double *x)
{
  auto &data = *impl_;
  data.requireReady();
  data.hessianReady = false;
  data.mapMotion(x, data.reference.data(), data.motion.data());
  if (!data.valid())
    throw std::invalid_argument("CUDA friction Hessian requires a finite state.");
  data.hessians.resize(9 * static_cast<std::size_t>(data.numContacts));
  if (data.numContacts)
    evaluateKernel<2><<<blocks(data.numContacts), blockSize>>>(data.contacts.data(), data.numContacts, data.motion.data(), data.epsilon, nullptr, nullptr, data.hessians.data());
  checkCuda(cudaGetLastError(), "CUDA friction Hessian preparation");
  data.hessianReady = true;
}
void CudaSurfaceIPCFrictionPotentialEnergy::addHessianVector(const double *v, double *out) const
{
  auto &data = *impl_;
  if (!data.hessianReady)
    throw std::logic_error("CUDA friction Hessian was not prepared.");
  if (!data.numContacts)
    return;
  data.mapMotion(v, nullptr, data.difference.data());
  data.scratch.zero();
  hessianVectorKernel<<<blocks(data.numContacts), blockSize>>>(data.contacts.data(), data.numContacts, data.hessians.data(), data.difference.data(), data.scratch.data());
  transposeKernel<<<blocks(data.map().rows), blockSize>>>(data.map(), data.scratch.data(), out);
  checkCuda(cudaGetLastError(), "CUDA friction mapped Hessian-vector product");
}
void CudaSurfaceIPCFrictionPotentialEnergy::addDiagonal(double *out) const
{
  const auto &data = *impl_;
  if (!data.hessianReady)
    throw std::logic_error("CUDA friction Hessian was not prepared.");
  if (data.numContacts)
    diagonalKernel<<<blocks(data.numContacts), blockSize>>>(data.contacts.data(), data.numContacts, data.hessians.data(), data.map(), out);
  checkCuda(cudaGetLastError(), "CUDA friction mapped Hessian diagonal");
}
}  // namespace pgo::Contact::CIPC
