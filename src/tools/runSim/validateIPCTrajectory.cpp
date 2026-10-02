#include "configFileJSON.h"
#include "hexJacobianValidation.h"
#include "triangleAreaValidation.h"
#include "embeddedSurfaceIPCPotentialEnergy.h"
#include "EigenSupport.h"
#include "initPredicates.h"
#include "ipc/broadPhase/surfaceIPCSelfBroadPhase.h"
#include "ipc/core/surfaceIPCMaxStep.h"
#include "ipc/geometry/ipcDistancePrimitives.h"
#include "ipc/topology/surfaceIPCTopology.h"
#include "pgoLogging.h"
#include "runIPCSimSetup.h"
#include "simulationMesh.h"

#include <argparse/argparse.hpp>
#include <fmt/format.h>
#include <tbb/parallel_for.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace
{
namespace ES = pgo::EigenSupport;
namespace IPC = pgo::Contact::CIPC;
namespace fs = std::filesystem;
using Json = nlohmann::json;

bool readState(const fs::path &path, int rows, ES::MXd &state)
{
  // Validate the header before allocating: the generic Eigen reader resizes
  // directly from file dimensions and only warns on an unexpected scalar size.
  std::ifstream stream(path, std::ios::binary);
  std::array<int, 3> header{};
  if (!stream.read(reinterpret_cast<char *>(header.data()), sizeof(header)) ||
    header[0] != rows || header[1] != 3 || header[2] != sizeof(double))
    return false;
  state.resize(rows, 3);
  return bool(stream.read(reinterpret_cast<char *>(state.data()), sizeof(double) * state.size())) && state.allFinite();
}

struct Tet
{
  std::array<int, 4> vertices;
  ES::M3d inverseRestEdges;
};

std::vector<Tet> prepareTets(const pgo::SolidDeformationModel::SimulationMesh &mesh)
{
  if (mesh.getElementType() != pgo::SolidDeformationModel::SimulationMeshType::TET || mesh.getNumElementVertices() != 4)
    throw std::invalid_argument("Trajectory determinant validation currently requires a tetrahedral simulation mesh.");
  std::vector<Tet> tets(mesh.getNumElements());
  for (int e = 0; e < mesh.getNumElements(); ++e) {
    auto &tet = tets[e];
    ES::V3d origin;
    mesh.getVertex(e, 0, origin.data());
    ES::M3d edges;
    for (int i = 0; i < 4; ++i) {
      tet.vertices[i] = mesh.getVertexIndex(e, i);
      if (i) {
        ES::V3d vertex;
        mesh.getVertex(e, i, vertex.data());
        edges.col(i - 1) = vertex - origin;
      }
    }
    if (!std::isfinite(edges.determinant()) || edges.determinant() == 0)
      throw std::invalid_argument("Rest tetrahedron has a zero or nonfinite determinant.");
    tet.inverseRestEdges = edges.inverse();
  }
  return tets;
}

ES::M3d displacementGradient(const Tet &tet, ES::ConstRefVecXd u)
{
  ES::M3d edges;
  for (int i = 0; i < 3; ++i)
    edges.col(i) = u.segment<3>(3 * tet.vertices[i + 1]) - u.segment<3>(3 * tet.vertices[0]);
  return edges * tet.inverseRestEdges;
}

// det(F + t D) is a cubic. Its minimum on [0,1] occurs at an endpoint
// or a real root of its quadratic derivative. Coefficients use determinant
// multilinearity, avoiding subtraction of almost equal endpoint determinants.
double minimumSegmentDeterminant(const ES::M3d &F, const ES::M3d &D)
{
  std::array<double, 4> c{};
  for (int mask = 0; mask < 8; ++mask) {
    ES::M3d mixed;
    int degree = 0;
    for (int column = 0; column < 3; ++column) {
      const bool moving = (mask & (1 << column)) != 0;
      mixed.col(column) = moving ? D.col(column) : F.col(column);
      degree += moving;
    }
    c[degree] += mixed.determinant();
  }
  auto evaluate = [&](double t) {
    return ((c[3] * t + c[2]) * t + c[1]) * t + c[0];
  };
  double minimum = std::min(c[0], evaluate(1.0));
  auto consider = [&](double t) {
    if (std::isfinite(t) && t > 0 && t < 1)
      minimum = std::min(minimum, evaluate(t));
  };
  const double a = 3 * c[3], b = 2 * c[2], d = c[1];
  if (a == 0) {
    if (b != 0)
      consider(-d / b);
  }
  else {
    const double discriminant = b * b - 4 * a * d;
    if (discriminant >= 0 && std::isfinite(discriminant)) {
      const double q = -0.5 * (b + std::copysign(std::sqrt(discriminant), b));
      consider(q / a);
      if (q != 0)
        consider(d / q);
    }
  }
  return minimum;
}

Json checkDeterminants(const std::vector<Tet> &tets, ES::ConstRefVecXd u, const ES::VXd *previous)
{
  double minimum = std::numeric_limits<double>::infinity();
  double sweptMinimum = minimum;
  int invalid = 0, sweptInvalid = 0, worst = -1;
  for (std::size_t e = 0; e < tets.size(); ++e) {
    const ES::M3d F = ES::M3d::Identity() + displacementGradient(tets[e], u);
    const double determinant = F.determinant();
    if (!std::isfinite(determinant) || determinant <= 0)
      ++invalid;
    if (determinant < minimum) {
      minimum = determinant;
      worst = static_cast<int>(e);
    }
    if (previous) {
      const ES::M3d F0 = ES::M3d::Identity() + displacementGradient(tets[e], *previous);
      const double segmentMinimum = minimumSegmentDeterminant(F0, F - F0);
      sweptMinimum = std::min(sweptMinimum, segmentMinimum);
      if (!std::isfinite(segmentMinimum) || segmentMinimum <= 0)
        ++sweptInvalid;
    }
  }
  Json result = { { "passed", invalid == 0 }, { "minimum", minimum },
    { "worst_tetrahedron", worst }, { "nonpositive_or_nonfinite_count", invalid } };
  if (previous) {
    result["linear_segment_minimum"] = sweptMinimum;
    result["linear_segment_nonpositive_or_nonfinite_count"] = sweptInvalid;
  }
  return result;
}

namespace Hex = pgo::RunSim::HexValidation;
struct Hexahedron
{
  std::array<int, 8> vertices;
  Hex::Positions rest;
  std::array<Hex::Interval, 8> gaussRest, cornerRest;
  std::array<double, 8> gaussOrientation;
  double orientation = 1, restUpper = 0;
  bool wholeRestCertified = false;
};

constexpr double gaussDeterminantFloor = 1e-8;

bool finiteInterval(Hex::Interval value)
{
  return std::isfinite(value.lo) && std::isfinite(value.hi) && value.lo <= value.hi;
}

// Keep invalid arithmetic sticky: min/max of products must never hide NaNs.
Hex::Interval intervalOperation(Hex::Interval a, Hex::Interval b, char operation)
{
  if (!finiteInterval(a) || !finiteInterval(b))
    return { -INFINITY, INFINITY };
  const auto value = operation == '+' ? a + b : operation == '-' ? a - b :
                                                                   a * b;
  return finiteInterval(value) ? value : Hex::Interval(-INFINITY, INFINITY);
}

std::array<Hex::Interval, 3> gaussCoordinates(int point)
{
  const double offset = 0.5 / std::sqrt(3.0);
  const double coordinate[2] = { 0.5 - offset, 0.5 + offset };
  return { coordinate[point / 4], coordinate[(point / 2) % 2], coordinate[point % 2] };
}

Hex::Interval intervalDeterminant(const Hex::Matrix &matrix, double orientation)
{
  auto add = [](auto a, auto b) {
    return intervalOperation(a, b, '+');
  };
  auto sub = [](auto a, auto b) {
    return intervalOperation(a, b, '-');
  };
  auto mul = [](auto a, auto b) {
    return intervalOperation(a, b, '*');
  };
  const auto value = add(sub(mul(matrix[0][0], sub(mul(matrix[1][1], matrix[2][2]), mul(matrix[1][2], matrix[2][1]))),
                           mul(matrix[0][1], sub(mul(matrix[1][0], matrix[2][2]), mul(matrix[1][2], matrix[2][0])))),
    mul(matrix[0][2], sub(mul(matrix[1][0], matrix[2][1]), mul(matrix[1][1], matrix[2][0]))));
  return mul(value, Hex::Interval(orientation));
}

Hex::Interval parameterJacobian(const Hex::Positions &positions,
  const std::array<Hex::Interval, 3> &parameter, double orientation)
{
  Hex::Matrix matrix{};
  auto add = [](auto a, auto b) {
    return intervalOperation(a, b, '+');
  };
  auto sub = [](auto a, auto b) {
    return intervalOperation(a, b, '-');
  };
  auto mul = [](auto a, auto b) {
    return intervalOperation(a, b, '*');
  };
  for (int axis = 0; axis < 3; ++axis)
    for (int bits = 0; bits < 8; ++bits) {
      if (bits & (1 << axis))
        continue;
      Hex::Interval weight(1);
      for (int other = 0; other < 3; ++other)
        if (other != axis)
          weight = mul(weight, (bits & (1 << other)) ? parameter[other] : sub(Hex::Interval(1), parameter[other]));
      for (int row = 0; row < 3; ++row)
        matrix[row][axis] = add(matrix[row][axis], mul(weight, sub(positions[Hex::bitToVertex[bits | (1 << axis)]][row], positions[Hex::bitToVertex[bits]][row])));
    }
  return intervalDeterminant(matrix, orientation);
}

Hex::Interval gaussJacobian(const Hex::Positions &positions, int point, double orientation)
{
  return parameterJacobian(positions, gaussCoordinates(point), orientation);
}

std::vector<int> negativeHexCorners(const Hexahedron &element, const Hex::Positions &positions)
{
  std::vector<int> result;
  for (int corner = 0; corner < 8; ++corner) {
    Hex::Matrix matrix{};
    for (int axis = 0; axis < 3; ++axis)
      for (int row = 0; row < 3; ++row)
        matrix[row][axis] = intervalOperation(positions[Hex::bitToVertex[corner | (1 << axis)]][row],
          positions[Hex::bitToVertex[corner & ~(1 << axis)]][row], '-');
    const auto rest = element.cornerRest[corner];
    const auto current = intervalDeterminant(matrix, 1);
    if (finiteInterval(rest) && finiteInterval(current) &&
      ((rest.lo > 0 && current.hi < 0) || (rest.hi < 0 && current.lo > 0)))
      result.push_back(corner);
  }
  return result;
}

// A fixed Gauss-point Jacobian is affine in time, so its determinant is cubic.
// Subdivide interval Bernstein coefficients; shared dyadic split points cover
// the entire original chord. Exhausted work or uncertain signs fail closed.
Hex::Certificate certifyGaussCoefficients(const std::array<Hex::Interval, 4> &coefficients)
{
  struct Node
  {
    std::array<Hex::Interval, 4> coefficients;
    int depth;
  };
  std::array<Node, 33> stack;
  stack[0] = { coefficients, 0 };
  int pending = 1, work = 0;
  double lower = INFINITY, upper = -INFINITY;
  while (pending) {
    const auto node = stack[--pending];
    double lo = INFINITY, hi = -INFINITY;
    for (const auto value : node.coefficients) {
      if (!finiteInterval(value))
        return {};
      lo = std::min(lo, value.lo);
      hi = std::max(hi, value.hi);
    }
    if (lo > 0) {
      lower = std::min(lower, lo);
      upper = std::max(upper, hi);
      continue;
    }
    if (++work > 4096 || node.depth == 32 || hi <= 0 || node.coefficients.front().hi <= 0 || node.coefficients.back().hi <= 0)
      return { false, hi <= 0 || node.coefficients.front().hi <= 0 || node.coefficients.back().hi <= 0, lo, hi, work };
    auto temporary = node.coefficients;
    Node left{ {}, node.depth + 1 }, right{ {}, node.depth + 1 };
    left.coefficients[0] = temporary[0];
    right.coefficients[3] = temporary[3];
    for (int level = 1; level < 4; ++level) {
      for (int i = 0; i < 4 - level; ++i)
        temporary[i] = intervalOperation(Hex::Interval(0.5), intervalOperation(temporary[i], temporary[i + 1], '+'), '*');
      left.coefficients[level] = temporary[0];
      right.coefficients[3 - level] = temporary[3 - level];
    }
    stack[pending++] = right;
    stack[pending++] = left;
  }
  return { true, false, lower, upper, work };
}

Hex::Certificate certifyGaussPoint(const Hexahedron &element, int point,
  const Hex::Positions &positions, const Hex::Positions *previous)
{
  const auto threshold = intervalOperation(Hex::Interval(gaussDeterminantFloor), element.gaussRest[point], '*');
  auto margin = [&](const Hex::Positions &x) {
    return intervalOperation(gaussJacobian(x, point, element.gaussOrientation[point]), threshold, '-');
  };
  if (!previous) {
    const auto value = margin(positions);
    return { finiteInterval(value) && value.lo > 0, !finiteInterval(value) || value.hi <= 0, value.lo, value.hi, 0 };
  }
  std::array<Hex::Interval, 4> coefficients;
  for (int t = 0; t < 4; ++t) {
    auto interpolated = t == 0 ? *previous : positions;
    if (t == 1 || t == 2) {
      const auto time = Hex::fraction(t, 3);
      for (int vertex = 0; vertex < 8; ++vertex)
        for (int axis = 0; axis < 3; ++axis)
          interpolated[vertex][axis] = intervalOperation(
            intervalOperation((*previous)[vertex][axis], intervalOperation(Hex::Interval(1), time, '-'), '*'),
            intervalOperation(positions[vertex][axis], time, '*'), '+');
    }
    coefficients[t] = margin(interpolated);
  }
  const auto a = intervalOperation(intervalOperation(intervalOperation(Hex::Interval(27), coefficients[1], '*'),
                                     intervalOperation(Hex::Interval(8), coefficients[0], '*'), '-'),
    coefficients[3], '-');
  const auto b = intervalOperation(intervalOperation(intervalOperation(Hex::Interval(27), coefficients[2], '*'),
                                     coefficients[0], '-'),
    intervalOperation(Hex::Interval(8), coefficients[3], '*'), '-');
  coefficients[1] = intervalOperation(intervalOperation(intervalOperation(Hex::Interval(2), a, '*'), b, '-'), Hex::fraction(1, 18), '*');
  coefficients[2] = intervalOperation(intervalOperation(intervalOperation(Hex::Interval(2), b, '*'), a, '-'), Hex::fraction(1, 18), '*');
  return certifyGaussCoefficients(coefficients);
}

std::vector<Hexahedron> prepareHexahedra(const pgo::SolidDeformationModel::SimulationMesh &mesh, const std::string &criterion)
{
  if (mesh.getElementType() != pgo::SolidDeformationModel::SimulationMeshType::CUBIC || mesh.getNumElementVertices() != 8)
    throw std::invalid_argument("Expected eight-node cubic elements.");
  std::vector<Hexahedron> elements(mesh.getNumElements());
  for (int e = 0; e < mesh.getNumElements(); ++e) {
    auto &element = elements[e];
    ES::V3d origin;
    mesh.getVertex(e, 0, origin.data());
    for (int v = 0; v < 8; ++v) {
      element.vertices[v] = mesh.getVertexIndex(e, v);
      ES::V3d vertex;
      mesh.getVertex(e, v, vertex.data());
      for (int axis = 0; axis < 3; ++axis)
        element.rest[v][axis] = v == 0 ? Hex::Interval(0) : Hex::Interval(vertex[axis]) - Hex::Interval(origin[axis]);
    }
    const auto center = Hex::jacobian(element.rest, { 0.5, 0.5, 0.5 });
    element.orientation = center.hi < 0 ? -1 : 1;
    const auto certificate = Hex::certify(element.rest, nullptr, element.orientation);
    element.wholeRestCertified = certificate.certified;
    if (!certificate.certified && criterion == "whole-cell")
      throw std::invalid_argument("Could not certify a nonzero, constant-sign rest hexahedron Jacobian throughout its parameter domain.");
    element.restUpper = certificate.upper;
    for (int point = 0; point < 8; ++point) {
      const auto raw = gaussJacobian(element.rest, point, 1);
      element.gaussOrientation[point] = raw.hi < 0 ? -1 : 1;
      element.gaussRest[point] = intervalOperation(raw, Hex::Interval(element.gaussOrientation[point]), '*');
      if (!finiteInterval(element.gaussRest[point]) || element.gaussRest[point].lo <= 0)
        throw std::invalid_argument("Could not certify a positive rest Jacobian at a cubic Gauss point.");
      const std::array<Hex::Interval, 3> corner{ double(point & 1), double((point >> 1) & 1), double((point >> 2) & 1) };
      element.cornerRest[point] = parameterJacobian(element.rest, corner, 1);
    }
  }
  return elements;
}

Hex::Positions deformedHexahedron(const Hexahedron &element, ES::ConstRefVecXd u)
{
  auto positions = element.rest;
  for (int v = 1; v < 8; ++v)
    for (int axis = 0; axis < 3; ++axis)
      positions[v][axis] = positions[v][axis] + (Hex::Interval(u[3 * element.vertices[v] + axis]) - Hex::Interval(u[3 * element.vertices[0] + axis]));
  return positions;
}

Json checkWholeHexDeterminants(const std::vector<Hexahedron> &elements, ES::ConstRefVecXd u, const ES::VXd *previous)
{
  struct Result
  {
    Hex::Certificate endpoint, segment;
    double minimum = 0, sweptMinimum = 0;
    std::vector<int> negativeCorners;
    bool endpointInverted = false, segmentInverted = false;
  };
  std::vector<Result> results(elements.size());
  tbb::parallel_for(std::size_t(0), elements.size(), [&](std::size_t e) {
    const auto positions = deformedHexahedron(elements[e], u);
    auto &result = results[e];
    result.negativeCorners = negativeHexCorners(elements[e], positions);
    result.endpointInverted = !result.negativeCorners.empty();
    if (elements[e].wholeRestCertified) {
      result.endpoint = Hex::certify(positions, nullptr, elements[e].orientation);
      result.minimum = Hex::down(result.endpoint.lower / elements[e].restUpper);
      // A negative sample or an entirely negative Bernstein sub-box is a
      // strict witness. An inconclusive/zero certificate is not an inversion.
      result.endpointInverted = result.endpointInverted ||
        (std::isfinite(result.endpoint.upper) && result.endpoint.upper < 0);
    }
    if (previous) {
      const auto before = deformedHexahedron(elements[e], *previous);
      result.segmentInverted = result.endpointInverted || !negativeHexCorners(elements[e], before).empty();
      if (elements[e].wholeRestCertified) {
        result.segment = Hex::certify(before, &positions, elements[e].orientation);
        result.sweptMinimum = Hex::down(result.segment.lower / elements[e].restUpper);
        result.segmentInverted = result.segmentInverted ||
          (std::isfinite(result.segment.upper) && result.segment.upper < 0);
      }
    }
  });
  double minimum = std::numeric_limits<double>::infinity(), sweptMinimum = minimum;
  int failed = 0, sweptFailed = 0, negative = 0, sweptNegative = 0, worst = -1, invalidRest = 0;
  std::vector<int> inverted, sweptInverted, unresolved, sweptUnresolved, negativeCornerElements;
  Json cornerWitnesses = Json::array();
  for (std::size_t e = 0; e < results.size(); ++e) {
    invalidRest += !elements[e].wholeRestCertified;
    const auto &result = results[e];
    failed += !result.endpoint.certified;
    negative += result.endpoint.nonpositiveSample;
    if (result.endpointInverted)
      inverted.push_back(static_cast<int>(e));
    else if (!result.endpoint.certified)
      unresolved.push_back(static_cast<int>(e));
    if (!result.negativeCorners.empty()) {
      negativeCornerElements.push_back(static_cast<int>(e));
      cornerWitnesses.push_back({ { "element", e }, { "corners", result.negativeCorners } });
    }
    if (result.minimum < minimum) {
      minimum = result.minimum;
      worst = static_cast<int>(e);
    }
    if (previous) {
      sweptFailed += !result.segment.certified;
      sweptNegative += result.segment.nonpositiveSample;
      if (result.segmentInverted)
        sweptInverted.push_back(static_cast<int>(e));
      else if (!result.segment.certified)
        sweptUnresolved.push_back(static_cast<int>(e));
      sweptMinimum = std::min(sweptMinimum, result.sweptMinimum);
    }
  }
  Json result = { { "passed", failed == 0 }, { "minimum", minimum }, { "worst_hexahedron", worst },
    { "minimum_is_certified_lower_bound", failed == 0 },
    { "method", "outward-rounded Jacobian interval bound, with tensor Bernstein subdivision fallback over the entire [0,1]^3 parameter domain and optional linear time interval" },
    { "nonpositive_or_nonfinite_count", negative }, { "uncertified_count", failed },
    { "uncertified_rest_count", invalidRest } };
  result["confirmed_inverted_element_ids"] = inverted;
  result["confirmed_inverted_element_count"] = inverted.size();
  result["negative_corner_element_ids"] = negativeCornerElements;
  result["negative_corner_element_count"] = negativeCornerElements.size();
  result["negative_corner_witnesses"] = std::move(cornerWitnesses);
  result["corner_indexing"] = "x + 2*y + 4*z, with x,y,z in {0,1}";
  result["uncertified_without_negative_witness_element_ids"] = unresolved;
  result["uncertified_without_negative_witness_element_count"] = unresolved.size();
  result["inversion_witness_scope"] = "strictly negative relative Jacobian at any of eight corners, or a strictly negative whole-cell certificate sample/sub-box with certified rest orientation; conservative witness count may miss unproven interior inversions";
  if (invalidRest) {
    result["minimum"] = nullptr;
    result["normalization_scope"] = "whole-cell relative determinant is undefined where a constant-sign nonzero rest Jacobian cannot be certified";
  }
  if (previous) {
    result["linear_segment_minimum"] = sweptMinimum;
    result["linear_segment_minimum_is_certified_lower_bound"] = sweptFailed == 0;
    // A failed certificate is conservatively rejected, even when a negative
    // sample was not found. Expose both counts rather than claim an inversion.
    result["linear_segment_nonpositive_or_nonfinite_count"] = sweptNegative;
    result["linear_segment_uncertified_count"] = sweptFailed;
    result["linear_segment_passed"] = sweptFailed == 0;
    result["linear_segment_confirmed_inverted_element_ids"] = sweptInverted;
    result["linear_segment_confirmed_inverted_element_count"] = sweptInverted.size();
    result["linear_segment_uncertified_without_negative_witness_element_ids"] = sweptUnresolved;
    result["linear_segment_uncertified_without_negative_witness_element_count"] = sweptUnresolved.size();
    if (invalidRest)
      result["linear_segment_minimum"] = nullptr;
  }
  return result;
}

Json checkHexDeterminants(const std::vector<Hexahedron> &elements, ES::ConstRefVecXd u,
  const ES::VXd *previous, const std::string &criterion)
{
  Json whole = checkWholeHexDeterminants(elements, u, previous);
  whole["criterion"] = "whole-cell";
  whole["required_for_validation"] = criterion == "whole-cell";
  Json result = whole;
  if (criterion == "gauss") {
    struct Result
    {
      bool endpointPassed = true, segmentPassed = true;
      int endpointNonpositive = 0, endpointBelowFloor = 0, segmentBelowFloor = 0;
      double minimum = INFINITY, sweptMinimum = INFINITY;
      int worstPoint = -1;
    };
    std::vector<Result> results(elements.size());
    tbb::parallel_for(std::size_t(0), elements.size(), [&](std::size_t e) {
      const auto positions = deformedHexahedron(elements[e], u);
      const auto before = previous ? deformedHexahedron(elements[e], *previous) : positions;
      auto &item = results[e];
      for (int point = 0; point < 8; ++point) {
        const auto endpoint = certifyGaussPoint(elements[e], point, positions, nullptr);
        const double minimum = Hex::down(gaussDeterminantFloor + Hex::down(endpoint.lower / elements[e].gaussRest[point].hi));
        item.endpointPassed = item.endpointPassed && endpoint.certified;
        item.endpointBelowFloor += endpoint.nonpositiveSample;
        const auto raw = gaussJacobian(positions, point, elements[e].gaussOrientation[point]);
        item.endpointNonpositive += !finiteInterval(raw) || raw.hi <= 0;
        if (minimum < item.minimum) {
          item.minimum = minimum;
          item.worstPoint = point;
        }
        if (previous) {
          const auto segment = certifyGaussPoint(elements[e], point, positions, &before);
          item.segmentPassed = item.segmentPassed && segment.certified;
          item.segmentBelowFloor += segment.nonpositiveSample;
          item.sweptMinimum = std::min(item.sweptMinimum,
            Hex::down(gaussDeterminantFloor + Hex::down(segment.lower / elements[e].gaussRest[point].hi)));
        }
      }
    });
    double minimum = INFINITY, sweptMinimum = INFINITY;
    int failed = 0, sweptFailed = 0, nonpositive = 0, belowFloor = 0, sweptBelowFloor = 0, worst = -1, point = -1;
    for (std::size_t e = 0; e < results.size(); ++e) {
      const auto &item = results[e];
      failed += !item.endpointPassed;
      sweptFailed += !item.segmentPassed;
      nonpositive += item.endpointNonpositive != 0;
      belowFloor += item.endpointBelowFloor != 0;
      sweptBelowFloor += item.segmentBelowFloor != 0;
      sweptMinimum = std::min(sweptMinimum, item.sweptMinimum);
      if (item.minimum < minimum) {
        minimum = item.minimum;
        worst = static_cast<int>(e);
        point = item.worstPoint;
      }
    }
    result = { { "passed", failed == 0 }, { "minimum", minimum }, { "worst_hexahedron", worst },
      { "worst_gauss_point", point }, { "criterion", "gauss" }, { "determinant_floor", gaussDeterminantFloor },
      { "gauss_points_per_element", 8 }, { "minimum_is_certified_lower_bound", failed == 0 },
      { "nonpositive_or_nonfinite_count", nonpositive }, { "at_or_below_floor_or_nonfinite_count", belowFloor },
      { "uncertified_count", failed },
      { "method", "outward interval Jacobians at the original eight Gauss points; cubic Bernstein subdivision certifies det(F)>1e-8 throughout each original linear time segment" } };
    if (previous) {
      result["linear_segment_passed"] = sweptFailed == 0;
      result["linear_segment_minimum"] = sweptMinimum;
      result["linear_segment_minimum_is_certified_lower_bound"] = sweptFailed == 0;
      result["linear_segment_at_or_below_floor_or_nonfinite_count"] = sweptBelowFloor;
      result["linear_segment_uncertified_count"] = sweptFailed;
    }
  }
  result["whole_cell_diagnostics"] = std::move(whole);
  return result;
}

namespace TriangleArea = pgo::RunSim::TriangleAreaValidation;
struct ShellTriangle
{
  std::array<int, 3> vertices;
  TriangleArea::Positions rest;
  int edgeScalePower = 0;
  double restAreaSquaredUpper = 0;
};

std::vector<ShellTriangle> prepareShellTriangles(const pgo::SolidDeformationModel::SimulationMesh &mesh)
{
  if (mesh.getElementType() != pgo::SolidDeformationModel::SimulationMeshType::SHELL || mesh.getNumElementVertices() != 6)
    throw std::invalid_argument("Expected six-vertex shell patches.");
  if (mesh.getNumElements() == 0)
    throw std::invalid_argument("Shell trajectory validation requires at least one triangle.");
  std::vector<ShellTriangle> triangles(mesh.getNumElements());
  for (int e = 0; e < mesh.getNumElements(); ++e) {
    auto &triangle = triangles[e];
    for (int v = 0; v < 3; ++v) {
      triangle.vertices[v] = mesh.getVertexIndex(e, v);
      mesh.getVertex(e, v, triangle.rest[v].data());
    }
    triangle.edgeScalePower = TriangleArea::edgeScalePower(triangle.rest);
    const auto certificate = TriangleArea::certify(triangle.rest, nullptr, triangle.edgeScalePower);
    if (!certificate.certified)
      throw std::invalid_argument("Could not certify a nonzero rest triangle area.");
    triangle.restAreaSquaredUpper = certificate.upper;
  }
  return triangles;
}

TriangleArea::Positions deformedTriangle(const ShellTriangle &triangle, ES::ConstRefVecXd u)
{
  auto positions = triangle.rest;
  // Shell collision embedding is the identity. These are the same represented
  // world endpoints checked by collision queries, including rest+u rounding.
  for (int v = 0; v < 3; ++v)
    for (int axis = 0; axis < 3; ++axis)
      positions[v][axis] += u[3 * triangle.vertices[v] + axis];
  return positions;
}

Json checkTriangleAreas(const std::vector<ShellTriangle> &triangles, ES::ConstRefVecXd u, const ES::VXd *previous)
{
  struct Result
  {
    TriangleArea::Certificate endpoint, segment;
    double minimum = 0, segmentMinimum = 0;
  };
  std::vector<Result> results(triangles.size());
  tbb::parallel_for(std::size_t(0), triangles.size(), [&](std::size_t e) {
    const auto positions = deformedTriangle(triangles[e], u);
    auto &result = results[e];
    result.endpoint = TriangleArea::certify(positions, nullptr, triangles[e].edgeScalePower);
    result.minimum = result.endpoint.certified ? std::max(0.0, TriangleArea::down(result.endpoint.lower / triangles[e].restAreaSquaredUpper)) : 0;
    if (previous) {
      const auto before = deformedTriangle(triangles[e], *previous);
      result.segment = TriangleArea::certify(before, &positions, triangles[e].edgeScalePower);
      result.segmentMinimum = result.segment.certified ? std::max(0.0, TriangleArea::down(result.segment.lower / triangles[e].restAreaSquaredUpper)) : 0;
    }
  });
  double minimum = std::numeric_limits<double>::infinity(), sweptMinimum = minimum;
  int failed = 0, sweptFailed = 0, worst = -1, subdivisions = 0;
  for (std::size_t e = 0; e < results.size(); ++e) {
    const auto &result = results[e];
    failed += !result.endpoint.certified;
    if (result.minimum < minimum) {
      minimum = result.minimum;
      worst = static_cast<int>(e);
    }
    if (previous) {
      sweptFailed += !result.segment.certified;
      sweptMinimum = std::min(sweptMinimum, result.segmentMinimum);
      subdivisions += result.segment.subdivisions;
    }
  }
  Json result = { { "passed", failed == 0 }, { "minimum", minimum }, { "worst_triangle", worst },
    { "metric", "squared area / squared rest area" }, { "minimum_is_certified_lower_bound", failed == 0 },
    { "uncertified_count", failed },
    { "method", "outward-rounded squared cross-product bounds; quartic Bernstein subdivision along the complete original linear segment; no orientation-sign constraint" } };
  if (previous) {
    result["linear_segment_minimum"] = sweptMinimum;
    result["linear_segment_minimum_is_certified_lower_bound"] = sweptFailed == 0;
    result["linear_segment_uncertified_count"] = sweptFailed;
    result["linear_segment_passed"] = sweptFailed == 0;
    result["linear_segment_bernstein_subdivisions"] = subdivisions;
  }
  return result;
}

Json checkSurface(const IPC::SurfaceIPCTopology &topology, ES::ConstRefVecXd positions, double dhat)
{
  std::vector<std::array<int, 2>> crossings;
  IPC::SurfaceIPCSelfBroadPhase broadPhase;
  broadPhase.findIntersectingEdgeTrianglePairs(topology, positions, crossings);
  int selfCrossings = 0, externalCrossings = 0;
  for (const auto &crossing : crossings) {
    if (topology.edgeContainsDeformableVertex(crossing[0]) && topology.triangleContainsDeformableVertex(crossing[1]))
      ++selfCrossings;
    else
      ++externalCrossings;
  }
  std::vector<IPC::PTPair> ptPairs;
  std::vector<IPC::EEPair> eePairs;
  broadPhase.buildPairs(topology, positions, dhat, ptPairs, eePairs);
  int invalidDistances = 0;
  double minimumDistance2 = std::numeric_limits<double>::infinity();
  auto point = [&](int i) -> ES::V3d {
    return positions.segment<3>(3 * i);
  };
  auto consider = [&](double distance2) {
    if (!std::isfinite(distance2) || distance2 <= 0)
      ++invalidDistances;
    minimumDistance2 = std::min(minimumDistance2, distance2);
  };
  for (const auto &pair : ptPairs)
    consider(IPC::distance::computePTSqDist(point(pair.p), point(pair.t0), point(pair.t1), point(pair.t2)));
  for (const auto &pair : eePairs)
    consider(IPC::distance::computeEESqDist(point(pair.ea0), point(pair.ea1), point(pair.eb0), point(pair.eb1)));
  Json result = { { "passed", crossings.empty() && invalidDistances == 0 },
    { "self_edge_triangle_intersections", selfCrossings },
    { "external_edge_triangle_intersections", externalCrossings },
    { "nonpositive_or_nonfinite_distances", invalidDistances },
    { "active_point_triangle_pairs", ptPairs.size() }, { "active_edge_edge_pairs", eePairs.size() },
    { "minimum_active_pair_distance", nullptr } };
  if (std::isfinite(minimumDistance2) && minimumDistance2 >= 0)
    result["minimum_active_pair_distance"] = std::sqrt(minimumDistance2);
  if (!crossings.empty())
    result["first_intersection_edge_triangle"] = crossings.front();
  return result;
}

struct SweepResult
{
  int calls = 0;
  int deepestSubdivision = 0;
  double minimumRawAlpha = 1;
  double unresolvedBegin = 0, unresolvedEnd = 1;
};

bool certifySegment(const IPC::SurfaceIPCTopology &topology, ES::ConstRefVecXd originalFrom, ES::ConstRefVecXd originalTo,
  double dhat, int depth, int maxDepth, SweepResult &result, double begin = 0, double end = 1)
{
  if ((originalTo - originalFrom).isZero(0))
    return true;
  const ES::VXd from = originalFrom + begin * (originalTo - originalFrom);
  const ES::VXd to = originalFrom + end * (originalTo - originalFrom);
  const ES::VXd dx = to - from;
  ++result.calls;
  result.deepestSubdivision = std::max(result.deepestSubdivision, depth);
  const double alpha = IPC::SurfaceIPCMaxStep().compute(topology, from, dx, dhat, 1.0);
  if (!std::isfinite(alpha) || alpha <= 0) {
    result.unresolvedBegin = begin;
    result.unresolvedEnd = end;
    return false;
  }
  result.minimumRawAlpha = std::min(result.minimumRawAlpha, alpha);
  if (alpha >= 1)
    return true;
  if (depth >= maxDepth) {
    result.unresolvedBegin = begin;
    result.unresolvedEnd = end;
    return false;
  }
  const double timeMiddle = begin + 0.5 * (end - begin);
  return certifySegment(topology, originalFrom, originalTo, dhat, depth + 1, maxDepth, result, begin, timeMiddle) &&
    certifySegment(topology, originalFrom, originalTo, dhat, depth + 1, maxDepth, result, timeMiddle, end);
}

Json checkSweep(const IPC::SurfaceIPCTopology &topology, ES::ConstRefVecXd from,
  ES::ConstRefVecXd to, double dhat, int maxDepth, bool endpointsValid)
{
  if (!endpointsValid)
    return { { "certified", false }, { "reason", "invalid_or_missing_endpoint" } };
  SweepResult result;
  const ES::VXd dx = to - from;
  const double fullBound = dx.isZero(0) ? 1.0 : IPC::SurfaceIPCMaxStep().compute(topology, from, dx, dhat, 1.0);
  const bool certified = std::isfinite(fullBound) && fullBound == 1.0;
  result.calls = dx.isZero(0) ? 0 : 1;
  result.minimumRawAlpha = std::isfinite(fullBound) ? fullBound : 0.0;
  // Only this original-chord certificate determines pass/fail. Bisection is
  // diagnostic localization, and can never promote an uncertified original
  // chord to a pass by reinterpreting rounded subinterval endpoints.
  if (!certified)
    certifySegment(topology, from, to, dhat, 0, maxDepth, result);
  Json report = { { "certified", certified }, { "ccd_calls", result.calls },
    { "deepest_subdivision", result.deepestSubdivision }, { "minimum_raw_alpha", result.minimumRawAlpha },
    { "reason", certified ? "complete" : "cpu_ccd_did_not_certify_original_segment" },
    { "subdivision_scope", "diagnostic localization only; pass requires original full-chord certificate" } };
  if (!certified) {
    report["first_unresolved_time_interval"] = { result.unresolvedBegin, result.unresolvedEnd };
    const double time = result.unresolvedBegin + 0.5 * (result.unresolvedEnd - result.unresolvedBegin);
    const ES::VXd positions = from + time * (to - from);
    report["unresolved_midpoint_time"] = time;
    report["unresolved_midpoint_collision"] = checkSurface(topology, positions, dhat);
  }
  return report;
}

struct PhysicalSubstep
{
  int frame = 0, index = 0, count = 0, iterations = 0, pcgIterations = 0, status = -1;
  std::size_t activePairs = 0;
  double fraction = 0, time = 0, gradient = 0, tolerance = 0;
  fs::path stateFile;
  std::string stop;
};

struct SubstepManifest
{
  bool present = false;
  std::map<int, std::vector<PhysicalSubstep>> frames;
  std::vector<std::string> errors;
};

SubstepManifest readSubstepManifest(const fs::path &filename, int numFrames)
{
  SubstepManifest result;
  result.present = fs::exists(filename);
  if (!result.present)
    return result;
  std::ifstream stream(filename);
  std::string line;
  const std::string header = "frame,substep,substeps,fraction,physical_time,state_file,newton_iterations,pcg_iterations,gradient_inf,tolerance,status,stop,active_pairs";
  auto trimCR = [](std::string &value) {
    if (!value.empty() && value.back() == '\r')
      value.pop_back();
  };
  if (!std::getline(stream, line)) {
    result.errors.push_back("Missing accepted-substep CSV header.");
    return result;
  }
  trimCR(line);
  if (line != header) {
    result.errors.push_back("Unexpected accepted-substep CSV schema.");
    return result;
  }
  auto integer = [](const std::string &value) {
    std::size_t consumed = 0;
    const int result = std::stoi(value, &consumed);
    if (consumed != value.size())
      throw std::invalid_argument("Invalid integer field.");
    return result;
  };
  auto scalar = [](const std::string &value) {
    std::size_t consumed = 0;
    const double result = std::stod(value, &consumed);
    if (consumed != value.size() || !std::isfinite(result))
      throw std::invalid_argument("Invalid or nonfinite scalar field.");
    return result;
  };
  int lineNumber = 1, previousFrame = -1, previousIndex = -1;
  while (std::getline(stream, line)) {
    ++lineNumber;
    trimCR(line);
    try {
      std::istringstream row(line);
      std::vector<std::string> fields;
      std::string field;
      while (std::getline(row, field, ','))
        fields.push_back(field);
      if (fields.size() != 13)
        throw std::invalid_argument("Expected 13 CSV fields.");
      PhysicalSubstep step;
      step.frame = integer(fields[0]);
      step.index = integer(fields[1]);
      step.count = integer(fields[2]);
      step.fraction = scalar(fields[3]);
      step.time = scalar(fields[4]);
      step.stateFile = fields[5];
      step.iterations = integer(fields[6]);
      step.pcgIterations = integer(fields[7]);
      step.gradient = scalar(fields[8]);
      step.tolerance = scalar(fields[9]);
      step.status = integer(fields[10]);
      step.stop = fields[11];
      const int pairs = integer(fields[12]);
      if (step.frame < 0 || step.frame >= numFrames || step.index < 0 || step.count < 1 || step.index >= step.count ||
        step.iterations < 0 || step.pcgIterations < 0 || pairs < 0 || step.stop.empty())
        throw std::invalid_argument("Invalid accepted-substep indices or counts.");
      if (step.frame < previousFrame || (step.frame == previousFrame && step.index <= previousIndex))
        throw std::invalid_argument("Accepted-substep rows are duplicated or out of order.");
      previousFrame = step.frame;
      previousIndex = step.index;
      step.activePairs = static_cast<std::size_t>(pairs);
      result.frames[step.frame].push_back(std::move(step));
    }
    catch (const std::exception &error) {
      result.errors.push_back(fmt::format("Line {}: {}", lineNumber, error.what()));
    }
  }
  if (stream.bad())
    result.errors.push_back("Could not read the accepted-substep CSV.");
  return result;
}

bool sameRecordedScalar(double actual, double expected)
{
  return actual == expected || (std::isfinite(actual) && std::isfinite(expected) && std::abs(actual - expected) <= 8 * std::numeric_limits<double>::epsilon() * std::abs(expected));
}

int validate(const fs::path &configPath, const fs::path &reportOverride, const fs::path &csvOverride,
  bool swept, int maxDepth, bool verbose, int firstFrame, int lastFrame, const std::string &cubicOverride)
{
  const auto start = std::chrono::steady_clock::now();
  pgo::ConfigFileJSON config;
  if (!config.open(configPath.string().c_str()))
    throw std::runtime_error("Cannot open trajectory configuration.");
  const bool cubic = config.exist("cubic-mesh");
  const bool shell = !cubic && !config.exist("tet-mesh");
  const bool intervalGeometry = cubic || shell;
  const std::string geometryKey = shell ? "triangle_areas" : cubic ? "cubic_determinants" :
                                                                     "tet_determinants";
  const std::string segmentMinimumKey = shell ? "minimum_linear_segment_area_squared_ratio" : "minimum_linear_segment_determinant";
  const std::string simType = config.getString("sim-type", 1);
  if (simType != "static" && simType != "dynamic")
    throw std::invalid_argument("Expected static or dynamic simulation type.");
  const bool dynamic = simType == "dynamic";
  const int numFrames = dynamic ? config.getInt("num-timestep", 1) : 1;
  if (numFrames < 1 || firstFrame < 0 || firstFrame >= numFrames)
    throw std::invalid_argument("Trajectory validation requires at least one saved frame and a first frame inside its range.");
  if (lastFrame < -1 || lastFrame >= numFrames || (lastFrame >= 0 && lastFrame < firstFrame))
    throw std::invalid_argument("The last frame must follow the first frame inside the configured range, or be -1 for its end.");
  const int endFrame = lastFrame < 0 ? numFrames : lastFrame + 1;
  const int expectedFrames = endFrame - firstFrame;
  const fs::path output = config.getResolvedPath("output", 1);
  auto validCriterion = [](const std::string &value) {
    return value == "gauss" || value == "whole-cell";
  };
  const std::string configuredCriterion = config.exist("cuda-cubic-validity") ? config.getString("cuda-cubic-validity", 1) : "";
  if ((config.exist("cuda-cubic-validity") && !validCriterion(configuredCriterion)) ||
    (!cubicOverride.empty() && !validCriterion(cubicOverride)))
    throw std::invalid_argument("Cubic validity must be gauss or whole-cell.");
  const fs::path validityMetadataPath = output / "material-validity.json";
  Json validityMetadata;
  std::string recordedCriterion;
  if (fs::exists(validityMetadataPath)) {
    std::ifstream stream(validityMetadataPath);
    stream >> validityMetadata;
    recordedCriterion = validityMetadata.at("cubic_validity").get<std::string>();
    const std::string expectedType = cubic ? "cubic" : shell ? "shell" :
                                                               "tet";
    if (validityMetadata.at("schema_version") != 1 || !validCriterion(recordedCriterion) ||
      validityMetadata.at("element_type") != expectedType ||
      !validityMetadata.at("material_max_step_enabled").is_boolean() ||
      validityMetadata.at("determinant_floor").get<double>() != gaussDeterminantFloor)
      throw std::invalid_argument("Invalid or incompatible material-validity metadata.");
    if ((!configuredCriterion.empty() && configuredCriterion != recordedCriterion && cubicOverride.empty()) ||
      (config.exist("enable-material-max-step") &&
        config.getValue<bool>("enable-material-max-step", 1) != validityMetadata.at("material_max_step_enabled").get<bool>()))
      throw std::invalid_argument("Configuration conflicts with recorded material-validity metadata.");
  }
  const std::string cubicCriterion = !cubicOverride.empty() ? cubicOverride : !configuredCriterion.empty() ? configuredCriterion :
    !recordedCriterion.empty()                                                                             ? recordedCriterion :
                                                                                                             "whole-cell";
  const std::string criterionSource = !cubicOverride.empty() ? "CLI" : !configuredCriterion.empty() ? "config" :
    !recordedCriterion.empty()                                                                      ? "material-validity.json" :
                                                                                                      "legacy whole-cell default";
  const fs::path reportPath = reportOverride.empty() ? output / "collision-validation.json" : reportOverride;
  const fs::path csvPath = csvOverride.empty() ? reportPath.parent_path() / "collision-frames.csv" : csvOverride;
  const fs::path substepManifestPath = output / "newton-substeps.csv";
  const auto substepManifest = readSubstepManifest(substepManifestPath, numFrames);
  const fs::path staticPathManifest = output / "static-newton-path.json";
  const bool staticPathPresent = !dynamic && fs::exists(staticPathManifest);
  const bool recordedPathPresent = substepManifest.present || staticPathPresent;
  const bool checkSweeps = swept || recordedPathPresent;
  const double timestep = dynamic ? config.getDouble("timestep", 1) : 0;
  const double solverTolerance = config.getDouble("solver-eps", 1);
  if ((dynamic && (!std::isfinite(timestep) || timestep <= 0)) || !std::isfinite(solverTolerance) || solverTolerance <= 0)
    throw std::invalid_argument("Trajectory timestep and solver tolerance must be valid.");
  pgo::Mesh::initPredicates();
  auto context = shell ? pgo::RunIPCSim::buildShellIpcSimulation(config) : pgo::RunIPCSim::buildVolumeIpcSimulation(config);
  const auto &surface = context.collisionSurface;
  const double dhat = context.collisionHandler->getIPCParameters().dhat;
  IPC::SurfaceIPCTopology topology;
  topology.setMesh(surface.vertices, surface.triangles, surface.vertexIsDeformable);
  const auto tets = cubic || shell ? std::vector<Tet>() : prepareTets(*context.simulationMeshOwner);
  const auto hexes = cubic ? prepareHexahedra(*context.simulationMeshOwner, cubicCriterion) : std::vector<Hexahedron>();
  const auto triangles = shell ? prepareShellTriangles(*context.simulationMeshOwner) : std::vector<ShellTriangle>();
  auto elementGeometry = [&](ES::ConstRefVecXd u, const ES::VXd *previous) {
    return shell ? checkTriangleAreas(triangles, u, previous) : cubic ? checkHexDeterminants(hexes, u, previous, cubicCriterion) :
                                                                        checkDeterminants(tets, u, previous);
  };
  ES::VXd rest(surface.vertices.rows() * 3);
  for (int i = 0; i < surface.vertices.rows(); ++i)
    rest.segment<3>(3 * i) = surface.vertices.row(i).transpose();
  ES::VXd previous = context.initialDisplacement;
  if (firstFrame > 0) {
    ES::MXd preceding;
    if (!readState(output / fmt::format("deform{:04d}.u", firstFrame - 1), static_cast<int>(previous.size()), preceding))
      throw std::invalid_argument("The state preceding --first-frame is missing, malformed or nonfinite.");
    previous = preceding.col(0);
    const double ratio = dynamic ? static_cast<double>(firstFrame - 1) / std::max(1, numFrames - 1) : 0.0;
    for (const auto &object : context.kinematicExternalObjects) {
      const auto positions = object.positionsAt(ratio);
      for (int i = 0; i < positions.rows(); ++i)
        rest.segment<3>(3 * (object.firstSurfaceVertex + i)) = positions.row(i).transpose();
    }
  }
  ES::VXd previousSurface = rest + surface.displacementMap * previous;
  Json initial = { { "collision", checkSurface(topology, previousSurface, dhat) },
    { geometryKey, elementGeometry(previous, nullptr) } };
  bool previousSurfaceValid = initial["collision"]["passed"].get<bool>();
  bool previousAvailable = true;
  bool allEndpointGeometryPassed = initial[geometryKey]["passed"].get<bool>();
  bool allPassed = previousSurfaceValid && allEndpointGeometryPassed && substepManifest.errors.empty();
  Json report = { { "config", fs::absolute(configPath).string() }, { "state_directory", output.string() },
    { "validator", shell ? "independent CPU exact-predicate edge-triangle intersection, IPC primitive distance and triangle squared-area certificates" : "independent CPU exact-predicate edge-triangle intersection, IPC primitive distance and volumetric Jacobian determinant checks" },
    { "collision_filter", "topologically disjoint primitives involving a deformable surface; prescribed-prescribed pairs excluded" },
    { "element_type", shell ? "shell" : cubic ? "cubic" :
                                                "tetrahedral" },
    { "frames_expected", expectedFrames },
    { "first_frame", firstFrame }, { "last_frame", endFrame - 1 }, { "configured_frames", numFrames },
    { "full_frame_range_requested", firstFrame == 0 && endFrame == numFrames }, { "swept_checks_requested", swept },
    { "swept_checks_requested_scope", "CLI --swept flag only; a recorded path manifest independently requires swept checks" },
    { "swept_checks_performed", checkSweeps },
    { "trajectory_path", staticPathPresent ? "recorded_static_newton_iterates" : substepManifest.present ? "recorded_physical_substeps" :
                                                                                                           "legacy_macro_chords" },
    { "sweep_path", staticPathPresent ? "recorded accepted Newton iterates under the full static objective; optimization progress is not physical time; initial-to-final chord is diagnostic only" : substepManifest.present ? "prescribed geometry moves first with deformable state fixed, then every recorded physical substep at fixed prescribed geometry; coarse macroframe interpolation is diagnostic only" :
                                                                                                                                                                                                                               "prescribed geometry moves first with deformable state fixed, then the saved macroframe chord at fixed prescribed geometry" },
    { "initial", initial }, { "frames", Json::array() } };
  report["element_validity_criterion"] = cubic ? cubicCriterion : shell ? "positive-triangle-area" :
                                                                          "positive-tet-determinant";
  if (cubic) {
    report["cubic_validity_source"] = criterionSource;
    report["material_validity_metadata"] = validityMetadata;
    report["config_metadata_criterion_conflict_overridden"] = !configuredCriterion.empty() && !recordedCriterion.empty() && configuredCriterion != recordedCriterion;
    report["material_validity_metadata_file"] = validityMetadataPath.string();
    report["element_determinant_scope"] = cubicCriterion == "gauss" ? "the original eight Gauss points in each hexahedron" : "the entire trilinear hexahedral parameter domain";
    report["determinant_floor"] = cubicCriterion == "gauss" ? gaussDeterminantFloor : 0;
  }
  int checked = 0, failed = 0;
  double globalMinimumGeometry = initial[geometryKey]["minimum"].get<double>();
  std::size_t selfCrossings = 0, externalCrossings = 0;
  std::size_t sweepMidpointSelfCrossings = 0, sweepMidpointExternalCrossings = 0;
  int uncertifiedSweeps = 0, actualSweepsChecked = 0;
  int coarseChecked = 0, coarseFailed = 0, coarseUncertified = 0;
  std::size_t coarseMidpointSelf = 0, coarseMidpointExternal = 0;
  double coarseMinimum = globalMinimumGeometry;
  bool coarseBoundsCertified = allEndpointGeometryPassed;
  int physicalSubstepsChecked = 0, physicalSubstepsFailed = 0, physicalFramesFailed = 0;
  int physicalUncertifiedSweeps = 0;
  std::size_t physicalSelfCrossings = 0, physicalExternalCrossings = 0;
  double physicalMinimum = globalMinimumGeometry, physicalSweptMinimum = globalMinimumGeometry;
  bool physicalEndpointBoundsCertified = allEndpointGeometryPassed, physicalSegmentBoundsCertified = true;
  Json staticPathSummary;
  double staticPathMinimum = globalMinimumGeometry;
  bool staticPathBoundsCertified = allEndpointGeometryPassed;
  auto checkStaticPath = [&](const ES::MXd &finalState, ES::ConstRefVecXd placedSurface) {
    Json result = { { "passed", false }, { "manifest", staticPathManifest.string() },
      { "required_for_validation", true }, { "physical_time", false },
      { "scope", "accepted Newton optimization iterates; only terminal states of each solve must be converged" },
      { "states", Json::array() }, { "errors", Json::array() } };
    int statesChecked = 0, segmentsChecked = 0, failedStates = 0, pathUncertified = 0;
    std::size_t pathSelf = 0, pathExternal = 0;
    bool valid = true;
    auto require = [](bool condition, const std::string &message) {
      if (!condition)
        throw std::invalid_argument(message);
    };
    auto integer = [&](const Json &value) {
      require(value.is_number_integer(), "Static-path indices must be integers.");
      return value.get<int>();
    };
    try {
      Json manifest;
      std::ifstream stream(staticPathManifest);
      require(bool(stream), "Cannot open the static Newton path manifest.");
      stream >> manifest;
      require(integer(manifest.at("schema_version")) == 1 && manifest.at("path_type") == "static_newton_iterates" &&
          manifest.at("parameter") == "accepted_newton_iteration" && manifest.at("physical_time") == false &&
          manifest.at("complete") == true && integer(manifest.at("num_dofs")) == previous.size(),
        "Incomplete or unsupported static Newton path schema.");
      const auto &states = manifest.at("states");
      const auto &solves = manifest.at("solves");
      const int expectedSolves = config.getValue<double>("contact-friction-coeff", 0, 0) > 0 ?
        config.getValue<int>("ipc-friction-iterations", 0, 1) :
        1;
      const int maximumIterations = config.getInt("solver-max-iter", 1);
      require(states.is_array() && !states.empty() && solves.is_array() &&
          expectedSolves >= 1 && solves.size() == static_cast<std::size_t>(expectedSolves),
        "Static Newton path does not cover the configured solves.");
      result["recorded_states"] = states.size();
      result["recorded_solves"] = solves.size();
      ES::MXd precedingState = ES::MXd::Zero(previous.size(), 3);
      precedingState.col(0) = context.initialDisplacement;
      ES::VXd precedingSurface = placedSurface;
      bool precedingValid = previousSurfaceValid;
      std::size_t nextState = 0;
      for (int solveIndex = 0; solveIndex < expectedSolves; ++solveIndex) {
        const auto &solve = solves.at(solveIndex);
        const int first = integer(solve.at("first_state")), last = integer(solve.at("last_state"));
        const double tolerance = solve.at("tolerance").get<double>(), terminalGradient = solve.at("gradient_inf").get<double>();
        require(integer(solve.at("solve")) == solveIndex && first >= 0 && static_cast<std::size_t>(first) == nextState &&
            last >= first && static_cast<std::size_t>(last) < states.size() && last - first <= maximumIterations &&
            integer(solve.at("newton_iterations")) == last - first &&
            integer(solve.at("status")) == 0 && solve.at("stop") == "gradient_tolerance" &&
            tolerance == solverTolerance && std::isfinite(terminalGradient) && terminalGradient >= 0 && terminalGradient < tolerance,
          "Static Newton solve indices, convergence metadata or iteration limit are invalid.");
        for (int index = first; index <= last; ++index) {
          const auto &knot = states.at(index);
          const int iteration = integer(knot.at("iteration"));
          const double gradient = knot.at("gradient_inf").get<double>(), energy = knot.at("energy").get<double>();
          require(integer(knot.at("solve")) == solveIndex && iteration == index - first &&
              std::isfinite(gradient) && gradient >= 0 && std::isfinite(energy),
            "Static Newton iterates are out of order or have invalid energy/gradient metadata.");
          if (index == last)
            require(gradient < solverTolerance && sameRecordedScalar(gradient, terminalGradient),
              "Terminal static Newton state does not match the converged solve.");
          const fs::path filename(knot.at("state_file").get<std::string>());
          require(!filename.empty() && !filename.is_absolute() && filename == filename.lexically_normal(),
            "Static Newton state must use a normalized relative file path.");
          for (const auto &component : filename)
            require(component != "..", "Static Newton state cannot escape the output directory.");
          ES::MXd state;
          require(readState(output / filename, static_cast<int>(previous.size()), state) &&
              state.rightCols(2).isZero(0),
            "Static Newton state is missing, malformed, nonfinite or contains physical velocity/acceleration.");
          if (index == first)
            require((state.array() == precedingState.array()).all(),
              "Static Newton initial state or lagged-solve chain does not match its predecessor.");
          const ES::VXd current = state.col(0), preceding = precedingState.col(0);
          const ES::VXd positions = rest + surface.displacementMap * current;
          Json item = { { "solve", solveIndex }, { "iteration", iteration }, { "state_file", filename.generic_string() },
            { "energy", energy }, { "gradient_inf", gradient },
            { "collision", checkSurface(topology, positions, dhat) },
            { geometryKey, elementGeometry(current, index > 0 ? &preceding : nullptr) } };
          const bool endpointValid = item["collision"]["passed"].get<bool>();
          const bool geometryValid = item[geometryKey]["passed"].get<bool>();
          bool segmentValid = true, sweepValid = true;
          if (index > 0) {
            segmentValid = intervalGeometry ? item[geometryKey]["linear_segment_passed"].get<bool>() :
                                              item[geometryKey]["linear_segment_nonpositive_or_nonfinite_count"].get<int>() == 0;
            item["motion_sweep"] = geometryValid && segmentValid ?
              checkSweep(topology, precedingSurface, positions, dhat, maxDepth, precedingValid && endpointValid) :
              Json{ { "certified", false }, { "reason", "invalid_or_uncertified_element_segment" } };
            item["motion_sweep"]["required_for_validation"] = true;
            item["motion_sweep"]["scope"] = "accepted Newton optimization transition; not physical time";
            sweepValid = item["motion_sweep"]["certified"].get<bool>();
            ++segmentsChecked;
            ++actualSweepsChecked;
            uncertifiedSweeps += !sweepValid;
            pathUncertified += !sweepValid;
            staticPathMinimum = std::min(staticPathMinimum, item[geometryKey]["linear_segment_minimum"].get<double>());
            if (item["motion_sweep"].contains("unresolved_midpoint_collision")) {
              const auto &collision = item["motion_sweep"]["unresolved_midpoint_collision"];
              sweepMidpointSelfCrossings += collision["self_edge_triangle_intersections"].get<std::size_t>();
              sweepMidpointExternalCrossings += collision["external_edge_triangle_intersections"].get<std::size_t>();
            }
          }
          const bool statePassed = endpointValid && geometryValid && segmentValid && sweepValid;
          item["passed"] = statePassed;
          valid = valid && statePassed;
          failedStates += !statePassed;
          ++statesChecked;
          pathSelf += item["collision"]["self_edge_triangle_intersections"].get<std::size_t>();
          pathExternal += item["collision"]["external_edge_triangle_intersections"].get<std::size_t>();
          staticPathMinimum = std::min(staticPathMinimum, item[geometryKey]["minimum"].get<double>());
          staticPathBoundsCertified = staticPathBoundsCertified && geometryValid && segmentValid;
          precedingState = state;
          precedingSurface = positions;
          precedingValid = endpointValid;
          result["states"].push_back(std::move(item));
        }
        nextState = static_cast<std::size_t>(last) + 1;
      }
      require(nextState == states.size(), "Static Newton path contains states outside the declared solves.");
      require((precedingState.array() == finalState.array()).all(),
        "Final static Newton state does not exactly match the exported equilibrium.");
      if (substepManifest.present) {
        const auto found = substepManifest.frames.find(0);
        require(substepManifest.errors.empty() && substepManifest.frames.size() == 1 &&
            found != substepManifest.frames.end() && found->second.size() == 1,
          "Static terminal metadata must contain exactly one final state.");
        const auto &step = found->second.front();
        require(step.index == 0 && step.count == 1 && step.fraction == 1 && step.time == 0 &&
            step.stateFile == "deform0000.u" && step.status == 0 && step.stop == "gradient_tolerance" &&
            step.tolerance == solverTolerance && step.gradient >= 0 && step.gradient < solverTolerance &&
            sameRecordedScalar(step.gradient, solves.back().at("gradient_inf").get<double>()),
          "Static terminal metadata conflicts with the recorded Newton path.");
      }
    }
    catch (const std::exception &error) {
      valid = false;
      staticPathBoundsCertified = false;
      result["errors"].push_back(error.what());
    }
    result["passed"] = valid;
    result["states_checked"] = statesChecked;
    result["segments_checked"] = segmentsChecked;
    result["failed_states"] = failedStates;
    result["uncertified_sweep_segments"] = pathUncertified;
    result["self_edge_triangle_intersections"] = pathSelf;
    result["external_edge_triangle_intersections"] = pathExternal;
    result[segmentMinimumKey] = staticPathMinimum;
    result["minimum_is_certified_lower_bound"] = intervalGeometry && valid && staticPathBoundsCertified;
    return result;
  };
  for (int frame = firstFrame; frame < endFrame; ++frame) {
    const fs::path statePath = output / fmt::format("deform{:04d}.u", frame);
    Json entry = { { "frame", frame }, { "state_file", statePath.string() }, { "passed", false } };
    ES::MXd state;
    if (!readState(statePath, static_cast<int>(previous.size()), state)) {
      entry["error"] = "missing, malformed or nonfinite state; expected displacement, velocity and acceleration columns";
      previousAvailable = false;
    }
    else {
      ++checked;
      const ES::VXd u = state.col(0);
      const double objectRatio = dynamic ? static_cast<double>(frame) / std::max(1, numFrames - 1) : 0.0;
      for (const auto &object : context.kinematicExternalObjects) {
        const ES::MXd positions = object.positionsAt(objectRatio);
        for (int i = 0; i < positions.rows(); ++i)
          rest.segment<3>(3 * (object.firstSurfaceVertex + i)) = positions.row(i).transpose();
      }
      const ES::VXd beforeSolve = rest + surface.displacementMap * previous;
      const ES::VXd currentSurface = rest + surface.displacementMap * u;
      entry["collision"] = checkSurface(topology, currentSurface, dhat);
      entry[geometryKey] = elementGeometry(u, checkSweeps && previousAvailable ? &previous : nullptr);
      globalMinimumGeometry = std::min(globalMinimumGeometry, entry[geometryKey]["minimum"].get<double>());
      allEndpointGeometryPassed = allEndpointGeometryPassed && entry[geometryKey]["passed"].get<bool>();
      selfCrossings += entry["collision"]["self_edge_triangle_intersections"].get<std::size_t>();
      externalCrossings += entry["collision"]["external_edge_triangle_intersections"].get<std::size_t>();
      bool passed = entry["collision"]["passed"].get<bool>() && entry[geometryKey]["passed"].get<bool>();
      bool placementValid = false;
      if (previousAvailable) {
        entry["after_prescribed_placement"] = checkSurface(topology, beforeSolve, dhat);
        placementValid = entry["after_prescribed_placement"]["passed"].get<bool>();
        passed = passed && placementValid;
      }
      if (checkSweeps) {
        entry["prescribed_motion_sweep"] = checkSweep(topology, previousSurface, beforeSolve, dhat, maxDepth,
          previousAvailable && previousSurfaceValid && placementValid);
        entry["deformable_motion_sweep"] = checkSweep(topology, beforeSolve, currentSurface, dhat, maxDepth,
          previousAvailable && placementValid && entry["collision"]["passed"].get<bool>());
        const bool coarseGeometryPassed = previousAvailable &&
          (intervalGeometry ? entry[geometryKey]["linear_segment_passed"].get<bool>() :
                              entry[geometryKey]["linear_segment_nonpositive_or_nonfinite_count"].get<int>() == 0);
        const bool coarseCertified = entry["deformable_motion_sweep"]["certified"].get<bool>();
        const bool coarsePassed = previousAvailable && coarseCertified && coarseGeometryPassed;
        entry["deformable_motion_sweep"]["required_for_validation"] = !recordedPathPresent;
        entry["deformable_motion_sweep"]["scope"] = "direct coarse macroframe shortcut";
        entry["coarse_chord_diagnostics"] = { { "passed", coarsePassed },
          { "required_for_validation", !recordedPathPresent },
          { "collision_certified", coarseCertified }, { shell ? "triangle_areas_passed" : "determinants_passed", coarseGeometryPassed } };
        ++coarseChecked;
        coarseFailed += !coarsePassed;
        coarseUncertified += !coarseCertified;
        coarseBoundsCertified = coarseBoundsCertified && coarseGeometryPassed;
        if (previousAvailable) {
          const double minimum = entry[geometryKey]["linear_segment_minimum"].get<double>();
          coarseMinimum = std::min(coarseMinimum, minimum);
          entry["coarse_chord_diagnostics"][segmentMinimumKey] = minimum;
        }
        passed = passed && previousAvailable && entry["prescribed_motion_sweep"]["certified"].get<bool>();
        if (!recordedPathPresent)
          passed = passed && coarsePassed;
        for (const char *key : { "prescribed_motion_sweep", "deformable_motion_sweep" }) {
          const auto &sweepReport = entry[key];
          const bool isCoarse = std::string(key) == "deformable_motion_sweep";
          const bool authoritative = !isCoarse || !recordedPathPresent;
          if (authoritative) {
            ++actualSweepsChecked;
            uncertifiedSweeps += !sweepReport["certified"].get<bool>();
          }
          if (sweepReport.contains("unresolved_midpoint_collision")) {
            const auto &collision = sweepReport["unresolved_midpoint_collision"];
            const auto self = collision["self_edge_triangle_intersections"].get<std::size_t>();
            const auto external = collision["external_edge_triangle_intersections"].get<std::size_t>();
            if (authoritative) {
              sweepMidpointSelfCrossings += self;
              sweepMidpointExternalCrossings += external;
            }
            if (isCoarse) {
              coarseMidpointSelf += self;
              coarseMidpointExternal += external;
            }
          }
        }
      }
      if (substepManifest.present && !staticPathPresent) {
        Json physical = { { "passed", false }, { "steps", Json::array() } };
        const auto found = substepManifest.frames.find(frame);
        bool physicalPassed = substepManifest.errors.empty() && previousAvailable && placementValid;
        if (found == substepManifest.frames.end() || found->second.empty()) {
          physical["error"] = "Accepted-substep manifest does not cover this frame.";
          physicalPassed = false;
        }
        else {
          const auto &steps = found->second;
          physical["expected_substeps"] = steps.front().count;
          physical["recorded_substeps"] = steps.size();
          physicalPassed = physicalPassed && static_cast<std::size_t>(steps.front().count) == steps.size();
          ES::VXd preceding = previous, precedingSurface = beforeSolve;
          bool precedingValid = previousAvailable && placementValid;
          for (std::size_t i = 0; i < steps.size(); ++i) {
            const auto &step = steps[i];
            const double expectedFraction = double(i + 1) / steps.size();
            const double expectedTime = dynamic ? (frame + expectedFraction) * timestep : 0;
            const bool finalStep = i + 1 == steps.size();
            const fs::path expectedPath = finalStep ? fs::path(fmt::format("deform{:04d}.u", frame)) :
                                                      fs::path("substeps") / fmt::format("frame{:04d}-step{:04d}.u", frame, i);
            bool stepPassed = step.index == static_cast<int>(i) && step.count == static_cast<int>(steps.size()) &&
              sameRecordedScalar(step.fraction, expectedFraction) && sameRecordedScalar(step.time, expectedTime) &&
              step.stateFile == expectedPath && step.status == 0 && step.gradient >= 0 &&
              step.tolerance == solverTolerance && step.gradient < step.tolerance;
            Json item = { { "substep", step.index }, { "state_file", step.stateFile.generic_string() },
              { "fraction", step.fraction }, { "physical_time", step.time }, { "gradient_inf", step.gradient },
              { "tolerance", step.tolerance }, { "status", step.status }, { "manifest_row_valid", stepPassed } };
            ES::MXd substate;
            if (step.stateFile != expectedPath || !readState(output / expectedPath, static_cast<int>(previous.size()), substate)) {
              item["error"] = "Unexpected state path or missing, malformed or nonfinite physical state.";
              stepPassed = false;
              precedingValid = false;
            }
            else {
              ++physicalSubstepsChecked;
              const ES::VXd current = substate.col(0);
              const ES::VXd positions = rest + surface.displacementMap * current;
              const bool reuseMacro = steps.size() == 1 && checkSweeps && previousAvailable;
              item["collision"] = finalStep ? entry["collision"] : checkSurface(topology, positions, dhat);
              item[geometryKey] = reuseMacro ? entry[geometryKey] : elementGeometry(current, &preceding);
              const bool endpointValid = item["collision"]["passed"].get<bool>();
              const bool geometryPassed = item[geometryKey]["passed"].get<bool>();
              const bool segmentPassed = intervalGeometry ? item[geometryKey]["linear_segment_passed"].get<bool>() :
                                                            item[geometryKey]["linear_segment_nonpositive_or_nonfinite_count"].get<int>() == 0;
              // Invalid element geometry already fails the trajectory.
              // Endpoint collisions remain checked; the contact sweep is not
              // computed through a segment with no geometry certificate.
              if (!geometryPassed || !segmentPassed)
                item["motion_sweep"] = { { "certified", false }, { "ccd_calls", 0 }, { "skipped", true },
                  { "reason", shell ? "degenerate_or_uncertified_triangle_segment" : "nonpositive_or_uncertified_volumetric_segment" } };
              else
                item["motion_sweep"] = reuseMacro ? entry["deformable_motion_sweep"] :
                                                    checkSweep(topology, precedingSurface, positions, dhat, maxDepth, precedingValid && endpointValid);
              item["motion_sweep"]["required_for_validation"] = true;
              item["motion_sweep"]["scope"] = "recorded physical substep";
              const bool sameMacroState = !finalStep || (substate.array() == state.array()).all();
              item["matches_macro_state"] = sameMacroState;
              stepPassed = stepPassed && endpointValid && geometryPassed && segmentPassed &&
                item["motion_sweep"]["certified"].get<bool>() && precedingValid && sameMacroState;
              physicalSelfCrossings += item["collision"]["self_edge_triangle_intersections"].get<std::size_t>();
              physicalExternalCrossings += item["collision"]["external_edge_triangle_intersections"].get<std::size_t>();
              physicalUncertifiedSweeps += !item["motion_sweep"]["certified"].get<bool>();
              ++actualSweepsChecked;
              uncertifiedSweeps += !item["motion_sweep"]["certified"].get<bool>();
              if (item["motion_sweep"].contains("unresolved_midpoint_collision")) {
                const auto &collision = item["motion_sweep"]["unresolved_midpoint_collision"];
                sweepMidpointSelfCrossings += collision["self_edge_triangle_intersections"].get<std::size_t>();
                sweepMidpointExternalCrossings += collision["external_edge_triangle_intersections"].get<std::size_t>();
              }
              physicalMinimum = std::min(physicalMinimum, item[geometryKey]["minimum"].get<double>());
              physicalSweptMinimum = std::min(physicalSweptMinimum, item[geometryKey]["linear_segment_minimum"].get<double>());
              physicalEndpointBoundsCertified = physicalEndpointBoundsCertified && geometryPassed;
              physicalSegmentBoundsCertified = physicalSegmentBoundsCertified && segmentPassed;
              preceding = current;
              precedingSurface = positions;
              precedingValid = endpointValid;
            }
            item["passed"] = stepPassed;
            physicalPassed = physicalPassed && stepPassed;
            physicalSubstepsFailed += !stepPassed;
            physical["steps"].push_back(std::move(item));
          }
        }
        physical["passed"] = physicalPassed;
        physicalFramesFailed += !physicalPassed;
        passed = passed && physicalPassed;
        entry["physical_substeps"] = std::move(physical);
      }
      if (staticPathPresent) {
        Json path = checkStaticPath(state, beforeSolve);
        passed = passed && path["passed"].get<bool>();
        staticPathSummary = path;
        staticPathSummary.erase("states");
        staticPathSummary["manifest_present"] = true;
        entry["static_newton_path"] = std::move(path);
        if (substepManifest.present)
          entry["physical_substeps"] = { { "passed", entry["static_newton_path"]["errors"].empty() },
            { "scope", "terminal equilibrium metadata only; motion is the recorded static Newton path" },
            { "authoritative_for_motion", false } };
      }
      entry["passed"] = passed;
      previous = u;
      previousSurface = currentSurface;
      previousSurfaceValid = entry["collision"]["passed"].get<bool>();
      previousAvailable = true;
    }
    const bool passed = entry["passed"].get<bool>();
    allPassed = allPassed && passed;
    failed += !passed;
    if (verbose || !passed || (frame + 1) % 100 == 0 || frame + 1 == endFrame)
      std::cout << "CPU validation frame " << frame << ": " << (passed ? "passed" : "FAILED") << std::endl;
    report["frames"].push_back(std::move(entry));
  }
  report["passed"] = allPassed;
  report["frames_checked"] = checked;
  report["frames_failed"] = failed;
  report["self_edge_triangle_intersections"] = selfCrossings;
  report["external_edge_triangle_intersections"] = externalCrossings;
  report["intersection_count_scope"] = "saved frame endpoints; separate sweep midpoint counts are diagnostic samples within uncertified intervals";
  report["actual_sweep_segments_checked"] = actualSweepsChecked;
  report["uncertified_sweep_segments"] = uncertifiedSweeps;
  report["sweep_count_scope"] = "authoritative prescribed motion and recorded physical substeps, or legacy macroframe chords when no manifest exists";
  report["coarse_chord_diagnostics"] = { { "required_for_validation", !recordedPathPresent },
    { "checked", coarseChecked }, { "failed", coarseFailed },
    { "passed", checkSweeps ? Json(coarseChecked == expectedFrames && coarseFailed == 0) : Json(nullptr) },
    { "uncertified_sweep_segments", coarseUncertified },
    { "sweep_midpoint_self_edge_triangle_intersections", coarseMidpointSelf },
    { "sweep_midpoint_external_edge_triangle_intersections", coarseMidpointExternal },
    { segmentMinimumKey, checkSweeps ? Json(coarseMinimum) : Json(nullptr) },
    { "minimum_is_certified_lower_bound", intervalGeometry && checkSweeps && coarseChecked == expectedFrames && coarseBoundsCertified },
    { "scope", "direct interpolation between saved macroframes; diagnostic only when a physical-substep manifest is present" } };
  report["sweep_midpoint_self_edge_triangle_intersections"] = sweepMidpointSelfCrossings;
  report["sweep_midpoint_external_edge_triangle_intersections"] = sweepMidpointExternalCrossings;
  const bool physicalCoverageComplete = physicalFramesFailed == 0 && checked == expectedFrames && substepManifest.errors.empty();
  const double trajectoryMinimum = staticPathPresent ? std::min(globalMinimumGeometry, staticPathMinimum) :
    substepManifest.present                          ? std::min(globalMinimumGeometry, physicalSweptMinimum) :
    checkSweeps                                      ? std::min(globalMinimumGeometry, coarseMinimum) :
                                                       globalMinimumGeometry;
  report[shell ? "minimum_triangle_area_squared_ratio" : "minimum_element_determinant"] = trajectoryMinimum;
  if (shell)
    report["element_geometry_metric"] = "squared area / squared rest area; normals may rotate without an orientation-sign constraint";
  report["minimum_scope"] = staticPathPresent ? "initial state, final equilibrium and recorded accepted Newton-iterate interpolation; optimization progress is not physical time; excludes diagnostic coarse shortcut" :
    substepManifest.present                   ? "initial state, saved endpoints and recorded physical-substep interpolation; excludes diagnostic coarse shortcuts" :
    checkSweeps                               ? "initial state, saved endpoints and legacy macroframe interpolation" :
                                                "initial state and saved macroframe endpoints";
  report["minimum_is_certified_lower_bound"] = intervalGeometry && allEndpointGeometryPassed && checked == expectedFrames &&
    (staticPathPresent          ? !staticPathSummary.is_null() && staticPathSummary["passed"].get<bool>() && staticPathBoundsCertified :
        substepManifest.present ? physicalCoverageComplete && physicalEndpointBoundsCertified && physicalSegmentBoundsCertified :
                                  !checkSweeps || coarseBoundsCertified);
  if (!cubic && !shell)
    report["minimum_tet_determinant"] = trajectoryMinimum;
  report["physical_substeps"] = { { "manifest_present", substepManifest.present },
    { "manifest", substepManifestPath.string() }, { "manifest_errors", substepManifest.errors },
    { "checked", physicalSubstepsChecked }, { "failed", physicalSubstepsFailed },
    { "frames_failed", physicalFramesFailed }, { "swept_checks_required_when_present", true },
    { "passed", substepManifest.present ? Json(physicalCoverageComplete) : Json(nullptr) },
    { shell ? "minimum_endpoint_area_squared_ratio" : "minimum_endpoint_determinant", substepManifest.present ? Json(physicalMinimum) : Json(nullptr) },
    { segmentMinimumKey, substepManifest.present ? Json(physicalSweptMinimum) : Json(nullptr) },
    { "endpoint_minimum_is_certified_lower_bound", substepManifest.present && intervalGeometry && physicalCoverageComplete && physicalEndpointBoundsCertified },
    { "segment_minimum_is_certified_lower_bound", substepManifest.present && intervalGeometry && physicalCoverageComplete && physicalSegmentBoundsCertified },
    { "self_edge_triangle_intersections", physicalSelfCrossings }, { "external_edge_triangle_intersections", physicalExternalCrossings },
    { "uncertified_sweep_segments", physicalUncertifiedSweeps },
    { "motion", "authoritative recorded physical substeps at their timestamps, with prescribed objects frozen at the macroframe placement; coarse shortcut checks are diagnostic only" } };
  if (staticPathPresent) {
    if (staticPathSummary.is_null())
      staticPathSummary = { { "manifest_present", true }, { "required_for_validation", true }, { "passed", false },
        { "errors", Json::array({ "Final equilibrium is unavailable; the static Newton path cannot be completed." }) } };
    report["static_newton_path"] = staticPathSummary;
    report["sweep_count_scope"] = "authoritative prescribed motion and recorded accepted Newton transitions; optimization progress is not physical time";
    report["coarse_chord_diagnostics"]["scope"] = "initial-to-equilibrium straight shortcut; diagnostic only because the accepted Newton path is recorded";
    report["physical_substeps"] = { { "manifest_present", substepManifest.present }, { "manifest", substepManifestPath.string() },
      { "manifest_errors", substepManifest.errors }, { "authoritative_for_motion", false },
      { "passed", substepManifest.present ? Json(staticPathSummary["errors"].empty()) : Json(nullptr) },
      { "scope", "terminal equilibrium metadata only; motion is the recorded static Newton path" } };
  }
  if (cubic) {
    int endpoints = 0, endpointFailures = 0, segments = 0, segmentFailures = 0, invalidRest = 0;
    int nonpositiveEndpoints = 0, nonpositiveSegments = 0;
    int invertedEndpointEvaluations = 0, invertedSegmentEvaluations = 0;
    std::set<int> invertedElements, unresolvedElements, invertedFrames;
    double endpointMinimum = INFINITY, segmentMinimum = INFINITY;
    bool endpointBounds = true, segmentBounds = true;
    auto collectWholeCell = [&](const Json &geometry, bool includeSegment, int frame) {
      const auto &whole = geometry.at("whole_cell_diagnostics");
      auto collectIds = [&](const char *invertedKey, const char *unresolvedKey) {
        const auto ids = whole.at(invertedKey).get<std::vector<int>>();
        invertedElements.insert(ids.begin(), ids.end());
        if (!ids.empty() && frame >= firstFrame)
          invertedFrames.insert(frame);
        const auto uncertain = whole.at(unresolvedKey).get<std::vector<int>>();
        unresolvedElements.insert(uncertain.begin(), uncertain.end());
        return static_cast<int>(ids.size());
      };
      invertedEndpointEvaluations += collectIds("confirmed_inverted_element_ids", "uncertified_without_negative_witness_element_ids");
      ++endpoints;
      endpointFailures += !whole.at("passed").get<bool>();
      nonpositiveEndpoints += whole.at("nonpositive_or_nonfinite_count").get<int>();
      invalidRest = std::max(invalidRest, whole.at("uncertified_rest_count").get<int>());
      endpointBounds = endpointBounds && whole.at("minimum_is_certified_lower_bound").get<bool>();
      if (whole.at("minimum").is_number())
        endpointMinimum = std::min(endpointMinimum, whole.at("minimum").get<double>());
      if (includeSegment && whole.contains("linear_segment_passed")) {
        invertedSegmentEvaluations += collectIds("linear_segment_confirmed_inverted_element_ids", "linear_segment_uncertified_without_negative_witness_element_ids");
        ++segments;
        segmentFailures += !whole.at("linear_segment_passed").get<bool>();
        nonpositiveSegments += whole.at("linear_segment_nonpositive_or_nonfinite_count").get<int>();
        segmentBounds = segmentBounds && whole.at("linear_segment_minimum_is_certified_lower_bound").get<bool>();
        if (whole.at("linear_segment_minimum").is_number())
          segmentMinimum = std::min(segmentMinimum, whole.at("linear_segment_minimum").get<double>());
      }
    };
    collectWholeCell(report["initial"][geometryKey], false, firstFrame - 1);
    for (const auto &entry : report["frames"]) {
      const int frame = entry["frame"].get<int>();
      if (entry.contains(geometryKey))
        collectWholeCell(entry[geometryKey], !recordedPathPresent && checkSweeps, frame);
      if (!staticPathPresent && entry.contains("physical_substeps") && entry["physical_substeps"].contains("steps"))
        for (const auto &step : entry["physical_substeps"]["steps"])
          if (step.contains(geometryKey))
            collectWholeCell(step[geometryKey], true, frame);
      if (staticPathPresent && entry.contains("static_newton_path") && entry["static_newton_path"].contains("states"))
        for (const auto &state : entry["static_newton_path"]["states"])
          if (state.contains(geometryKey))
            collectWholeCell(state[geometryKey], true, frame);
    }
    for (int element : invertedElements)
      unresolvedElements.erase(element);
    report["whole_cell_diagnostics"] = { { "criterion", "whole-cell" },
      { "required_for_validation", cubicCriterion == "whole-cell" },
      { "passed", endpointFailures == 0 && segmentFailures == 0 && checked == expectedFrames },
      { "endpoint_checks", endpoints }, { "failed_endpoint_checks", endpointFailures },
      { "segment_checks", segments }, { "failed_segment_checks", segmentFailures },
      { "nonpositive_endpoint_element_evaluations", nonpositiveEndpoints },
      { "nonpositive_segment_element_evaluations", nonpositiveSegments }, { "uncertified_rest_elements", invalidRest },
      { "confirmed_inverted_element_ids", invertedElements }, { "confirmed_inverted_element_count", invertedElements.size() },
      { "confirmed_inverted_endpoint_evaluations", invertedEndpointEvaluations },
      { "confirmed_inverted_segment_evaluations", invertedSegmentEvaluations },
      { "frames_with_confirmed_inversion", invertedFrames }, { "frames_with_confirmed_inversion_count", invertedFrames.size() },
      { "initial_state_confirmed_inverted_element_ids", report["initial"][geometryKey]["whole_cell_diagnostics"]["confirmed_inverted_element_ids"] },
      { "uncertified_without_negative_witness_element_ids", unresolvedElements },
      { "uncertified_without_negative_witness_element_count", unresolvedElements.size() },
      { "inversion_count_scope", "unique element/frame IDs from strict negative corner or certificate witnesses on available actual-path states/segments; evaluation totals may repeat endpoints; global unresolved IDs exclude every element with a negative witness anywhere in this checked path; witness count is not exhaustive for unproved interior inversions" },
      { "minimum_endpoint_determinant", std::isfinite(endpointMinimum) && !invalidRest ? Json(endpointMinimum) : Json(nullptr) },
      { "minimum_linear_segment_determinant", std::isfinite(segmentMinimum) && !invalidRest ? Json(segmentMinimum) : Json(nullptr) },
      { "minimum_is_certified_lower_bound", endpointBounds && segmentBounds && checked == expectedFrames },
      { "scope", "available initial/saved/recorded-path endpoint evaluations and authoritative interpolation segments; repeated endpoint evaluations are counted; diagnostic coarse shortcuts excluded" } };
  }
  report["elapsed_seconds"] = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
  if (!reportPath.parent_path().empty())
    fs::create_directories(reportPath.parent_path());
  std::ofstream stream(reportPath);
  if (!stream || !(stream << report.dump(2) << '\n'))
    throw std::runtime_error("Failed to write trajectory validation report: " + reportPath.string());
  if (!csvPath.parent_path().empty())
    fs::create_directories(csvPath.parent_path());
  std::ofstream csv(csvPath);
  csv << std::setprecision(17);
  csv << "frame,passed,self_edge_triangle_intersections,external_edge_triangle_intersections,invalid_primitive_distances,"
      << (shell ? "min_area_squared_ratio,min_linear_segment_area_squared_ratio" : "min_det_F,min_linear_segment_det_F")
      << ",prescribed_sweep_certified,deformable_sweep_certified,minimum_is_certified_lower_bound,deformable_sweep_required_for_validation,physical_substeps_passed"
      << ",element_validity_criterion";
  if (cubic)
    csv << ",whole_cell_endpoint_passed,whole_cell_linear_segment_passed,whole_cell_min_det_F,whole_cell_confirmed_inverted_elements,whole_cell_linear_segment_confirmed_inverted_elements";
  csv << '\n';
  for (const auto &entry : report["frames"]) {
    csv << entry["frame"].get<int>() << ',' << entry["passed"].get<bool>();
    if (entry.contains("collision")) {
      const auto &collision = entry["collision"];
      csv << ',' << collision["self_edge_triangle_intersections"].get<int>()
          << ',' << collision["external_edge_triangle_intersections"].get<int>()
          << ',' << collision["nonpositive_or_nonfinite_distances"].get<int>()
          << ',' << entry[geometryKey]["minimum"];
      csv << ',';
      if (entry[geometryKey].contains("linear_segment_minimum"))
        csv << entry[geometryKey]["linear_segment_minimum"];
      csv << ',';
      if (entry.contains("prescribed_motion_sweep"))
        csv << entry["prescribed_motion_sweep"]["certified"].get<bool>();
      csv << ',';
      if (entry.contains("deformable_motion_sweep"))
        csv << entry["deformable_motion_sweep"]["certified"].get<bool>();
    }
    else {
      csv << ",,,,,,,";
    }
    csv << ',';
    if (entry.contains(geometryKey))
      csv << (intervalGeometry && entry[geometryKey]["minimum_is_certified_lower_bound"].get<bool>());
    csv << ',' << !recordedPathPresent << ',';
    if (entry.contains("physical_substeps"))
      csv << entry["physical_substeps"]["passed"].get<bool>();
    csv << ',' << report["element_validity_criterion"].get<std::string>();
    if (cubic) {
      csv << ',';
      if (entry.contains(geometryKey)) {
        const auto &whole = entry[geometryKey]["whole_cell_diagnostics"];
        csv << whole["passed"].get<bool>() << ',';
        if (whole.contains("linear_segment_passed"))
          csv << whole["linear_segment_passed"].get<bool>();
        csv << ',' << whole["minimum"] << ',' << whole["confirmed_inverted_element_count"] << ',';
        if (whole.contains("linear_segment_confirmed_inverted_element_count"))
          csv << whole["linear_segment_confirmed_inverted_element_count"];
      }
      else
        csv << ",,,,";
    }
    csv << '\n';
  }
  if (!csv)
    throw std::runtime_error("Failed to write trajectory validation CSV: " + csvPath.string());
  std::cout << "CPU trajectory validation " << (allPassed ? "PASSED" : "FAILED") << ": " << checked << '/' << expectedFrames
            << " frames; self intersections=" << selfCrossings << "; external intersections=" << externalCrossings
            << "; sweep midpoint intersections=" << sweepMidpointSelfCrossings + sweepMidpointExternalCrossings
            << (shell ? "; min area squared/rest area squared=" : "; min det(F)=") << trajectoryMinimum << "; report=" << reportPath << std::endl;
  return allPassed ? 0 : 2;
}
}  // namespace

int main(int argc, char *argv[])
{
  pgo::Logging::init();
  argparse::ArgumentParser program("validateIPCTrajectory");
  program.add_argument("config").help("IPC simulation config whose saved .u frames are checked").required();
  program.add_argument("--report").default_value(std::string()).help("JSON report path; default: output/collision-validation.json");
  program.add_argument("--csv").default_value(std::string()).help("Per-frame CSV path; default: collision-frames.csv next to the JSON report");
  program.add_argument("--swept").default_value(false).implicit_value(true).help("Also check piecewise linear frame motion with CPU CCD and volume/triangle geometry certificates");
  program.add_argument("--max-sweep-depth").default_value(16).scan<'i', int>().help("Maximum CPU CCD subdivision depth");
  program.add_argument("--first-frame").default_value(0).scan<'i', int>().help("Start at this frame, using its preceding saved state; default checks the entire trajectory");
  program.add_argument("--last-frame").default_value(-1).scan<'i', int>().help("Stop after this frame without changing prescribed-motion timing; -1 checks through the configured final frame");
  program.add_argument("--verbose").default_value(false).implicit_value(true).help("Print a result for every frame");
  program.add_argument("--cubic-validity").default_value(std::string()).help("gauss or whole-cell; overrides explicit config/recorded metadata; legacy default is whole-cell");
  try {
    program.parse_args(argc, argv);
    const int depth = program.get<int>("--max-sweep-depth");
    if (depth < 0 || depth > 24)
      throw std::invalid_argument("Maximum sweep depth must be between 0 and 24.");
    return validate(program.get<std::string>("config"), program.get<std::string>("--report"), program.get<std::string>("--csv"),
      program.get<bool>("--swept"), depth, program.get<bool>("--verbose"), program.get<int>("--first-frame"), program.get<int>("--last-frame"), program.get<std::string>("--cubic-validity"));
  }
  catch (const std::exception &error) {
    std::cerr << "IPC trajectory validation error: " << error.what() << std::endl;
    return 1;
  }
}
