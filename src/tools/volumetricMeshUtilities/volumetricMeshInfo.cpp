#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
struct Vec3
{
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;
};

Vec3 operator-(const Vec3 &a, const Vec3 &b)
{
  return Vec3{ a.x - b.x, a.y - b.y, a.z - b.z };
}

Vec3 operator*(double s, const Vec3 &v)
{
  return Vec3{ s * v.x, s * v.y, s * v.z };
}

Vec3 &operator+=(Vec3 &a, const Vec3 &b)
{
  a.x += b.x;
  a.y += b.y;
  a.z += b.z;
  return a;
}

Vec3 &operator/=(Vec3 &a, double s)
{
  a.x /= s;
  a.y /= s;
  a.z /= s;
  return a;
}

double dot(const Vec3 &a, const Vec3 &b)
{
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

Vec3 cross(const Vec3 &a, const Vec3 &b)
{
  return Vec3{
    a.y * b.z - a.z * b.y,
    a.z * b.x - a.x * b.z,
    a.x * b.y - a.y * b.x
  };
}

std::string trim(const std::string &s)
{
  const auto begin = std::find_if_not(s.begin(), s.end(), [](unsigned char c) { return std::isspace(c); });
  const auto end = std::find_if_not(s.rbegin(), s.rend(), [](unsigned char c) { return std::isspace(c); }).base();
  if (begin >= end)
    return "";
  return std::string(begin, end);
}

std::string uppercase(std::string s)
{
  std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return std::toupper(c); });
  return s;
}

bool readDataLine(std::istream &in, std::string &line)
{
  while (std::getline(in, line)) {
    line = trim(line);
    if (line.empty() || line[0] == '#')
      continue;
    return true;
  }

  return false;
}

struct VolumetricMeshData
{
  std::vector<Vec3> vertices;
  std::vector<std::vector<int>> elements;
  int numElementVertices = 0;
};

void readVerticesSection(std::istream &in, VolumetricMeshData &mesh)
{
  std::string line;
  if (readDataLine(in, line) == false)
    throw std::runtime_error("Unexpected end of file after *VERTICES.");

  int numVertices = 0;
  int dimension = 0;
  std::istringstream header(line);
  if (!(header >> numVertices >> dimension) || numVertices < 0 || dimension != 3)
    throw std::runtime_error("Invalid *VERTICES header.");

  mesh.vertices.assign(numVertices, Vec3{});

  for (int i = 0; i < numVertices; ++i) {
    if (readDataLine(in, line) == false)
      throw std::runtime_error("Unexpected end of file while reading vertices.");

    int vertexID = 0;
    Vec3 vertex;
    std::istringstream vertexLine(line);
    if (!(vertexLine >> vertexID >> vertex.x >> vertex.y >> vertex.z))
      throw std::runtime_error("Invalid vertex line: " + line);
    if (vertexID < 1 || vertexID > numVertices)
      throw std::runtime_error("Vertex index out of range: " + std::to_string(vertexID));

    mesh.vertices[vertexID - 1] = vertex;
  }
}

void readElementsSection(std::istream &in, VolumetricMeshData &mesh)
{
  std::string line;
  if (readDataLine(in, line) == false)
    throw std::runtime_error("Unexpected end of file after *ELEMENTS.");

  const std::string elementType = uppercase(line);
  if (elementType != "TET" && elementType != "CUBIC")
    throw std::runtime_error("Unsupported *ELEMENTS type: " + line);

  if (readDataLine(in, line) == false)
    throw std::runtime_error("Unexpected end of file in *ELEMENTS header.");

  int numElements = 0;
  int numElementVertices = 0;
  std::istringstream header(line);
  if (!(header >> numElements >> numElementVertices) || numElements < 0)
    throw std::runtime_error("Invalid *ELEMENTS header.");

  const int expectedElementVertices = (elementType == "TET") ? 4 : 8;
  if (numElementVertices != expectedElementVertices)
    throw std::runtime_error("Unexpected number of vertices per element in *ELEMENTS.");

  mesh.numElementVertices = numElementVertices;
  mesh.elements.assign(numElements, std::vector<int>(numElementVertices, 0));

  for (int i = 0; i < numElements; ++i) {
    if (readDataLine(in, line) == false)
      throw std::runtime_error("Unexpected end of file while reading elements.");

    int elementID = 0;
    std::istringstream elementLine(line);
    if (!(elementLine >> elementID))
      throw std::runtime_error("Invalid element line: " + line);
    if (elementID < 1 || elementID > numElements)
      throw std::runtime_error("Element index out of range: " + std::to_string(elementID));

    for (int j = 0; j < numElementVertices; ++j) {
      int vertexID = 0;
      if (!(elementLine >> vertexID))
        throw std::runtime_error("Invalid element line: " + line);
      if (vertexID < 1 || vertexID > static_cast<int>(mesh.vertices.size()))
        throw std::runtime_error("Element vertex index out of range: " + std::to_string(vertexID));

      mesh.elements[elementID - 1][j] = vertexID - 1;
    }
  }
}

VolumetricMeshData loadVegFile(const std::string &filename)
{
  std::ifstream in(filename);
  if (!in)
    throw std::runtime_error("Failed to open file: " + filename);

  VolumetricMeshData mesh;
  bool readVertices = false;
  bool readElements = false;

  std::string line;
  while (readDataLine(in, line)) {
    if (line[0] != '*')
      continue;

    const std::string section = uppercase(line);
    if (section == "*VERTICES") {
      readVerticesSection(in, mesh);
      readVertices = true;
    }
    else if (section == "*ELEMENTS") {
      if (readVertices == false)
        throw std::runtime_error("*ELEMENTS appeared before *VERTICES.");

      readElementsSection(in, mesh);
      readElements = true;
    }
  }

  if (readVertices == false)
    throw std::runtime_error("Missing *VERTICES section.");
  if (readElements == false)
    throw std::runtime_error("Missing *ELEMENTS section.");

  return mesh;
}

double tetVolume(const Vec3 &a, const Vec3 &b, const Vec3 &c, const Vec3 &d)
{
  return std::abs(dot(d - a, cross(b - a, c - a))) / 6.0;
}

double cubicVolume(const std::array<Vec3, 8> &v)
{
  return std::abs(dot(v[1] - v[0], cross(v[3] - v[0], v[4] - v[0])));
}

Vec3 elementCenter(const std::vector<Vec3> &vertices, const std::vector<int> &element)
{
  Vec3 center;
  for (int vertexID : element)
    center += vertices[vertexID];

  center /= static_cast<double>(element.size());
  return center;
}

double elementVolume(const VolumetricMeshData &mesh, const std::vector<int> &element)
{
  if (mesh.numElementVertices == 4) {
    return tetVolume(
      mesh.vertices[element[0]],
      mesh.vertices[element[1]],
      mesh.vertices[element[2]],
      mesh.vertices[element[3]]);
  }

  std::array<Vec3, 8> vertices;
  for (int i = 0; i < 8; ++i)
    vertices[i] = mesh.vertices[element[i]];

  return cubicVolume(vertices);
}

Vec3 computeCenterOfMass(const VolumetricMeshData &mesh, double &totalVolume)
{
  Vec3 centerOfMass;
  totalVolume = 0.0;

  for (const std::vector<int> &element : mesh.elements) {
    const double volume = elementVolume(mesh, element);
    totalVolume += volume;
    centerOfMass += volume * elementCenter(mesh.vertices, element);
  }

  if (totalVolume == 0.0)
    throw std::runtime_error("Cannot compute center of mass because total volume is zero.");

  centerOfMass /= totalVolume;
  return centerOfMass;
}

void printUsage(const char *programName)
{
  std::cerr << "Usage: " << programName << " <input.veg>\n";
}
}  // namespace

int main(int argc, char **argv)
{
  if (argc != 2) {
    printUsage(argv[0]);
    return 1;
  }

  try {
    const VolumetricMeshData mesh = loadVegFile(argv[1]);

    double totalVolume = 0.0;
    const Vec3 centerOfMass = computeCenterOfMass(mesh, totalVolume);

    std::cout << std::setprecision(17);
    std::cout << "#vtx: " << mesh.vertices.size() << '\n';
    std::cout << "#elements: " << mesh.elements.size() << '\n';
    std::cout << "#element vertices: " << mesh.numElementVertices << '\n';
    std::cout << "total volume: " << totalVolume << '\n';
    std::cout << "center of mass: "
              << centerOfMass.x << ' '
              << centerOfMass.y << ' '
              << centerOfMass.z << '\n';
  }
  catch (const std::exception &e) {
    std::cerr << "volumetricMeshInfo error: " << e.what() << '\n';
    return 1;
  }

  return 0;
}
