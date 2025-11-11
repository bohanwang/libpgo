/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "volumetricMesh" library , Copyright (C) 2007 CMU, 2009 MIT, 2018 USC *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code author: Jernej Barbic                                            *
 * http://www.jernejbarbic.com/vega                                      *
 *                                                                       *
 * Research: Jernej Barbic, Hongyi Xu, Yijing Li,                        *
 *           Danyong Zhao, Bohan Wang,                                   *
 *           Fun Shing Sin, Daniel Schroeder,                            *
 *           Doug L. James, Jovan Popovic                                *
 *                                                                       *
 * Funding: National Science Foundation, Link Foundation,                *
 *          Singapore-MIT GAMBIT Game Lab,                               *
 *          Zumberge Research and Innovation Fund at USC,                *
 *          Sloan Foundation, Okawa Foundation,                          *
 *          USC Annenberg Foundation                                     *
 *                                                                       *
 * This library is free software; you can redistribute it and/or         *
 * modify it under the terms of the BSD-style license that is            *
 * included with this library in the file LICENSE.txt                    *
 *                                                                       *
 * This library is distributed in the hope that it will be useful,       *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the file     *
 * LICENSE.TXT for more details.                                         *
 *                                                                       *
 *************************************************************************/

#include "tetMesh.h"
#include "volumetricMeshParser.h"
#include "geometryQuery.h"
#include "tetMeshGeo.h"
#include "predicates.h"

namespace pgo
{
namespace VolumetricMeshes
{

TetMesh::TetMesh(const char *filename, fileFormatType fileFormat, int verbose):
  VolumetricMesh(filename, fileFormat, 4, &temp, verbose)
{
  if (temp != elementType()) {
    printf("Error: mesh is not a tet mesh.\n");
    throw 11;
  }
}

TetMesh::TetMesh(void *binaryStream, int memoryLoad):
  VolumetricMesh(binaryStream, 4, &temp, memoryLoad)
{
  if (temp != elementType()) {
    printf("Error: mesh is not a tet mesh.\n");
    throw 11;
  }
}

TetMesh::TetMesh(const char *filename, int specialFileType, int verbose):
  VolumetricMesh(4)
{
  if (specialFileType != 0) {
    printf("Unknown special file type %d requested.\n", specialFileType);
    throw 1;
  }

  char lineBuffer[1024];
  VolumetricMeshParser parser;

  // first, read the vertices
  sprintf(lineBuffer, "%s.node", filename);
  if (parser.open(lineBuffer) != 0)
    throw 2;

  parser.getNextLine(lineBuffer, 1);
  int dim;
  sscanf(lineBuffer, "%d %d", &numVertices, &dim);
  if (dim != 3)
    throw 3;

  vertices.resize(numVertices);
  for (int i = 0; i < numVertices; i++) {
    parser.getNextLine(lineBuffer, 1);
    int index;
    double x, y, z;
    sscanf(lineBuffer, "%d %lf %lf %lf", &index, &x, &y, &z);
    if (index != (i + 1))
      throw 3;
    vertices[i] = Vec3d(x, y, z);
  }

  parser.close();

  // next, read the elements
  sprintf(lineBuffer, "%s.ele", filename);
  if (parser.open(lineBuffer) != 0)
    throw 4;

  parser.getNextLine(lineBuffer, 1);
  sscanf(lineBuffer, "%d %d", &numElements, &dim);
  if (dim != 4) {
    printf("Error: not a tet mesh file (%d vertices per tet encountered).\n", dim);
    throw 5;
  }

  elements.resize(4 * numElements);
  elementMaterial.resize(numElements);

  for (int i = 0; i < numElements; i++) {
    parser.getNextLine(lineBuffer, 1);
    int index;
    int v[4];
    sscanf(lineBuffer, "%d %d %d %d %d", &index, &v[0], &v[1], &v[2], &v[3]);
    if (index != (i + 1))
      throw 6;
    for (int j = 0; j < 4; j++)  // vertices are 1-indexed in .ele files
      elements[i * 4 + j] = v[j] - 1;
  }

  parser.close();

  setSingleMaterial(E_default, nu_default, density_default);
}

TetMesh::TetMesh(const Vec3d &p0, const Vec3d &p1, const Vec3d &p2, const Vec3d &p3):
  VolumetricMesh(4)
{
  numVertices = 4;
  vertices = { p0, p1, p2, p3 };
  //  vertices[0] = p0;
  //  vertices[1] = p1;
  //  vertices[2] = p2;
  //  vertices[3] = p3;

  numElements = 1;
  elements.resize(4 * numElements);
  elementMaterial.resize(numElements);
  for (int i = 0; i < 4; i++)
    elements[i] = i;

  setSingleMaterial(E_default, nu_default, density_default);
}

TetMesh::TetMesh(int numVertices_, const double *vertices_, int numElements_, const int *elements_, double E, double nu, double density):
  VolumetricMesh(numVertices_, vertices_, numElements_, 4, elements_, E, nu, density)
{
}

TetMesh::TetMesh(int numVertices_, const double *vertices_, int numElements_, const int *elements_,
  int numMaterials_, const Material *const *materials_, int numSets_, const Set *sets_, int numRegions_, const Region *regions_):
  VolumetricMesh(numVertices_, vertices_, numElements_, 4, elements_,
    numMaterials_, materials_, numSets_, sets_, numRegions_, regions_)
{
}

TetMesh::TetMesh(const std::vector<Vec3d> &vertices, const std::vector<Vec4i> &elements, double E, double nu, double density):
  TetMesh(vertices.size(), (double *)vertices.data(), elements.size(), (int *)elements.data(), E, nu, density)
{
}

TetMesh::TetMesh(const TetMesh &source):
  VolumetricMesh(source)
{
}

VolumetricMesh *TetMesh::clone()
{
  TetMesh *mesh = new TetMesh(*this);
  return mesh;
}

TetMesh::TetMesh(const TetMesh &tetMesh, int numElements_, int *elements_, std::map<int, int> *vertexMap_):
  VolumetricMesh(tetMesh, numElements_, elements_, vertexMap_)
{
}

TetMesh::~TetMesh()
{
}

int TetMesh::saveToAscii(const char *filename) const
{
  return VolumetricMesh::saveToAscii(filename, elementType());
}

int TetMesh::saveToBinary(const char *filename, unsigned int *bytesWritten) const
{
  return VolumetricMesh::saveToBinary(filename, bytesWritten, elementType());
}

int TetMesh::saveToBinary(FILE *binaryOutputStream, unsigned int *bytesWritten, bool countBytesOnly) const
{
  return VolumetricMesh::saveToBinary(binaryOutputStream, bytesWritten, elementType(), countBytesOnly);
}

void TetMesh::exportMeshGeometry(std::vector<Vec3d> &vertices, std::vector<Vec4i> &tets) const
{
  exportMeshGeometry(vertices);
  tets.resize(getNumElements());
  for (int i = 0; i < getNumElements(); i++) {
    const int *ii = getVertexIndices(i);
    tets[i] = Vec4i(ii[0], ii[1], ii[2], ii[3]);
  }
}

void TetMesh::exportMeshGeometry(Mesh::TetMeshGeo &geo) const
{
  exportMeshGeometry(geo.positions(), geo.tets());
}

void TetMesh::computeElementMassMatrix(int el, double *massMatrix) const
{
  /*
    Consistent mass matrix of a tetrahedron =

                   [ 2  1  1  1  ]
                   [ 1  2  1  1  ]
       mass / 20 * [ 1  1  2  1  ]
                   [ 1  1  1  2  ]

    Note: mass = density * volume. Other than via the mass, the
    consistent mass matrix does not depend on the shape of the tetrahedron.
    (This can be seen after a long algebraic derivation; see:
     Singiresu S. Rao: The finite element method in engineering, 2004)
  */

  const double mtx[16] = { 2, 1, 1, 1,
    1, 2, 1, 1,
    1, 1, 2, 1,
    1, 1, 1, 2 };

  double density = getElementDensity(el);
  double factor = density * getElementVolume(el) / 20;

  for (int i = 0; i < 16; i++)
    massMatrix[i] = factor * mtx[i];

  // lumped mass
  /*
    double mass = element(el)->density() * getElementVolume(el);
    massMatrix[0] = massMatrix[5] = massMatrix[10] = massMatrix[15] = mass / 4.0;
  */
}

double TetMesh::getSignedTetVolume(const Vec3d &a, const Vec3d &b, const Vec3d &c, const Vec3d &d)
{
  return 1.0 / 6 * getTetDeterminant(a, b, c, d);
}

double TetMesh::getTetVolume(const Vec3d &a, const Vec3d &b, const Vec3d &c, const Vec3d &d)
{
  // volume = 1/6 * | (d-a) . ((b-a) x (c-a)) |
  return (1.0 / 6 * std::abs(getTetDeterminant(a, b, c, d)));
}

double TetMesh::getElementVolume(int el) const
{
  const Vec3d &a = getVertex(el, 0);
  const Vec3d &b = getVertex(el, 1);
  const Vec3d &c = getVertex(el, 2);
  const Vec3d &d = getVertex(el, 3);
  return getTetVolume(a, b, c, d);
}

void TetMesh::getElementInertiaTensor(int el, Mat3d &inertiaTensor) const
{
  Vec3d a = getVertex(el, 0);
  Vec3d b = getVertex(el, 1);
  Vec3d c = getVertex(el, 2);
  Vec3d d = getVertex(el, 3);

  Vec3d center = getElementCenter(el);
  a -= center;
  b -= center;
  c -= center;
  d -= center;

  double absdetJ = fabs(getTetDeterminant(a, b, c, d));

  double x1 = a[0], x2 = b[0], x3 = c[0], x4 = d[0];
  double y1 = a[1], y2 = b[1], y3 = c[1], y4 = d[1];
  double z1 = a[2], z2 = b[2], z3 = c[2], z4 = d[2];

  double A = absdetJ * (y1 * y1 + y1 * y2 + y2 * y2 + y1 * y3 + y2 * y3 + y3 * y3 + y1 * y4 + y2 * y4 + y3 * y4 + y4 * y4 + z1 * z1 + z1 * z2 + z2 * z2 + z1 * z3 + z2 * z3 + z3 * z3 + z1 * z4 + z2 * z4 + z3 * z4 + z4 * z4) / 60.0;

  double B = absdetJ * (x1 * x1 + x1 * x2 + x2 * x2 + x1 * x3 + x2 * x3 + x3 * x3 + x1 * x4 + x2 * x4 + x3 * x4 + x4 * x4 + z1 * z1 + z1 * z2 + z2 * z2 + z1 * z3 + z2 * z3 + z3 * z3 + z1 * z4 + z2 * z4 + z3 * z4 + z4 * z4) / 60.0;

  double C = absdetJ * (x1 * x1 + x1 * x2 + x2 * x2 + x1 * x3 + x2 * x3 + x3 * x3 + x1 * x4 + x2 * x4 + x3 * x4 + x4 * x4 + y1 * y1 + y1 * y2 + y2 * y2 + y1 * y3 + y2 * y3 + y3 * y3 + y1 * y4 + y2 * y4 + y3 * y4 + y4 * y4) / 60.0;

  double Ap = absdetJ * (2 * y1 * z1 + y2 * z1 + y3 * z1 + y4 * z1 + y1 * z2 + 2 * y2 * z2 + y3 * z2 + y4 * z2 + y1 * z3 + y2 * z3 + 2 * y3 * z3 + y4 * z3 + y1 * z4 + y2 * z4 + y3 * z4 + 2 * y4 * z4) / 120.0;

  double Bp = absdetJ * (2 * x1 * z1 + x2 * z1 + x3 * z1 + x4 * z1 + x1 * z2 + 2 * x2 * z2 + x3 * z2 + x4 * z2 + x1 * z3 + x2 * z3 + 2 * x3 * z3 + x4 * z3 + x1 * z4 + x2 * z4 + x3 * z4 + 2 * x4 * z4) / 120.0;

  double Cp = absdetJ * (2 * x1 * y1 + x2 * y1 + x3 * y1 + x4 * y1 + x1 * y2 + 2 * x2 * y2 + x3 * y2 + x4 * y2 + x1 * y3 + x2 * y3 + 2 * x3 * y3 + x4 * y3 + x1 * y4 + x2 * y4 + x3 * y4 + 2 * x4 * y4) / 120.0;

  inertiaTensor = asMat3d(A, -Bp, -Cp, -Bp, B, -Ap, -Cp, -Ap, C);
}

bool TetMesh::containsVertex(int el, Vec3d pos) const  // true if given element contain given position, false otherwise
{
  double weights[4];
  return Mesh::pointInTet(pos.data(),
    getVertex(el, 0).data(),
    getVertex(el, 1).data(),
    getVertex(el, 2).data(),
    getVertex(el, 3).data());

  // computeBarycentricWeights(el, pos, weights);
  // all weights must be non-negative
  // return ((weights[0] >= 0) && (weights[1] >= 0) && (weights[2] >= 0) && (weights[3] >= 0));
}

void TetMesh::computeBarycentricWeights(const Vec3d tetVtxPos[4], const Vec3d &pos, double weights[4])
{
  return computeBarycentricWeights(tetVtxPos[0], tetVtxPos[1], tetVtxPos[2], tetVtxPos[3], pos, weights);
}

void TetMesh::computeBarycentricWeights(const Vec3d &tetVtxPos0, const Vec3d &tetVtxPos1, const Vec3d &tetVtxPos2, const Vec3d &tetVtxPos3,
  const Vec3d &pos, double weights[4])
{
  Mesh::getTetBarycentricWeights(pos, tetVtxPos0, tetVtxPos1, tetVtxPos2, tetVtxPos3, weights);
}

void TetMesh::computeBarycentricWeights(int el, const Vec3d &pos, double *weights) const
{
  Vec3d vtx[4];
  for (int i = 0; i < 4; i++)
    vtx[i] = getVertex(el, i);

  computeBarycentricWeights(vtx, pos, weights);
}

double TetMesh::getTetDeterminant(const Vec3d &a, const Vec3d &b, const Vec3d &c, const Vec3d &d)
{
  return Mesh::getTetDeterminant(a, b, c, d);
}

void TetMesh::interpolateGradient(int element, const double *U, int numFields, Vec3d pos, double *grad) const
{
  computeGradient(element, U, numFields, grad);
}

void TetMesh::computeGradient(int element, const double *U, int numFields, double *grad) const
{
  // grad is 9 x numFields
  // grad is constant inside a tet
  Vec3d vtx[4];
  for (int i = 0; i < 4; i++)
    vtx[i] = getVertex(element, i);

  // form M =
  // [b - a]
  // [c - a]
  // [d - a]

  Mat3d M = asMat3d(vtx[1] - vtx[0], vtx[2] - vtx[0], vtx[3] - vtx[0]);
  Mat3d MInv = M.fullPivLu().inverse();
  // printf("M=\n");
  // M.print();

  for (int field = 0; field < numFields; field++) {
    // form rhs =
    // [U1 - U0]
    // [U2 - U0]
    // [U3 - U0]
    const double *u[4];
    for (int i = 0; i < 4; i++)
      u[i] = &U[3 * numVertices * field + 3 * getVertexIndex(element, i)];

    Vec3d rows[3];
    for (int i = 0; i < 3; i++)
      rows[i] = asVec3d(u[i + 1]) - asVec3d(u[0]);

    Mat3d rhs = asMat3d(rows[0], rows[1], rows[2]);
    // printf("rhs=\n");
    // rhs.print();

    Mat3d gradM = (MInv * rhs).transpose();
    (EigenSupport::Mp<Mat3d>)(&grad[9 * field]) = gradM;

    /*
        // test gradient
        if (field == 0)
        {
          printf("----\n");
          printf("0: pos: %.15f %.15f %.15f | uExact: %.15f %.15f %.15f\n", vtx[0][0], vtx[0][1], vtx[0][2], u[0][0], u[0][1], u[0][2]);
          for(int vertex=0; vertex<3; vertex++)
          {
            Vec3d u1 = asVec3d(u[vertex+1]);
            printf("%d: ", vertex+1);
            printf("pos: %.15f %.15f %.15f | uExact: %.15f %.15f %.15f | ", vtx[vertex+1][0], vtx[vertex+1][1], vtx[vertex+1][2], u1[0], u1[1], u1[2]);
            Vec3d u1approx = asVec3d(u[0]) + gradM * (vtx[vertex+1] - vtx[0]);
            printf("uApprox: %.15f %.15f %.15f\n", u1approx[0], u1approx[1], u1approx[2]);
          }
          printf("----\n");
        }
    */
  }
}

int TetMesh::getNumElementEdges() const
{
  return 6;
}

void TetMesh::getElementEdges(int el, int *edgeBuffer) const
{
  int v[4];
  for (int i = 0; i < 4; i++)
    v[i] = getVertexIndex(el, i);

  int edgeMask[6][2] = {
    { 0, 1 }, { 1, 2 }, { 2, 0 },
    { 0, 3 }, { 1, 3 }, { 2, 3 }
  };

  for (int edge = 0; edge < 6; edge++) {
    edgeBuffer[2 * edge + 0] = v[edgeMask[edge][0]];
    edgeBuffer[2 * edge + 1] = v[edgeMask[edge][1]];
  }
}

void TetMesh::orient()
{
  for (int el = 0; el < numElements; el++) {
    // a, b, c, d
    // dot(d - a, cross(b - a, c - a))
    double det = getTetDeterminant(getVertex(el, 0), getVertex(el, 1), getVertex(el, 2), getVertex(el, 3));

    if (det < 0) {
      // reverse tet
      int *elementVertices = &elements[el * 4];
      // swap 2 and 3
      std::swap(elementVertices[2], elementVertices[3]);
    }
  }
}

int TetMesh::getClosestElement(const Vec3d &pos) const
{
  // linear scan
  double closestDist = DBL_MAX;
  int closestElement = 0;
  for (int element = 0; element < numElements; element++) {
    double dist = Mesh::getSquaredDistanceToTet(pos,
      getVertex(element, 0),
      getVertex(element, 1),
      getVertex(element, 2),
      getVertex(element, 3));

    if (dist < closestDist) {
      closestDist = dist;
      closestElement = element;
    }
  }

  return closestElement;
}
}  // namespace VolumetricMeshes
}  // namespace pgo
