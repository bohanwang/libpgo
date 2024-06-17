/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "volumetricMesh" library , Copyright (C) 2007 CMU, 2009 MIT, 2018 USC *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code authors: Jernej Barbic, Yijing Li                                *
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

#include "generateSurfaceMesh.h"
#include "cubicMesh.h"
#include "triKey.h"
#include "tetMesh.h"
#include "rectKey.h"
#include "pgoLogging.h"

#include <unordered_map>
#include <cfloat>

// the main routine
void pgo::VolumetricMeshes::GenerateSurfaceMesh::computeMesh(const VolumetricMesh *volumetricMesh, std::vector<EigenSupport::V3d> &vertices, std::vector<std::vector<int>> &faces,
  bool triangulate, bool allElementFaces)
{
  int numElementVertices = volumetricMesh->getNumElementVertices();
  int faceDegree = 0;

  if (numElementVertices == 4) {
    faceDegree = 3;
    triangulate = false;
  }

  if (numElementVertices == 8)
    faceDegree = 4;

  if (faceDegree == 0) {
    printf("Error: unsupported volumetricMesh type encountered.\n");
    return;
  }

  // create an empty surface volumetricMesh
  vertices.clear();
  faces.clear();

  // add all vertices
  for (int i = 0; i < volumetricMesh->getNumVertices(); i++)
    vertices.emplace_back(volumetricMesh->getVertex(i));

  // build unique list of all surface faces

  if (volumetricMesh->getElementType() == VolumetricMesh::TET)  // tet volumetricMesh
  {
    const TetMesh *tetMesh = dynamic_cast<const TetMesh *>(volumetricMesh);
    PGO_ALOG(tetMesh != nullptr);

    std::unordered_map<Mesh::OTriKey, int> surfaceFaces;
    for (int i = 0; i < volumetricMesh->getNumElements(); i++) {
      // compute determinant to establish orientation
      double det = tetMesh->getTetDeterminant(i);

      auto processFace = [&](int q0, int q1, int q2) {
        Mesh::OTriKey key(volumetricMesh->getVertexIndex(i, q0), volumetricMesh->getVertexIndex(i, q1), volumetricMesh->getVertexIndex(i, q2));
        if (allElementFaces)  // get all faces
        {
          faces.emplace_back(std::vector<int>{ key[0], key[1], key[2] });
          return;
        }

        auto it = surfaceFaces.find(key);
        if (it != surfaceFaces.end())
          it->second++;
        else {
          auto revKey = key.getReversedTriKey();
          it = surfaceFaces.find(revKey);
          if (it != surfaceFaces.end())
            it->second--;
          else
            surfaceFaces.emplace(key, 1);
        }
      };

      if (det >= 0) {
        processFace(1, 2, 3);
        processFace(2, 0, 3);
        processFace(3, 0, 1);
        processFace(1, 0, 2);
      }
      else {
        processFace(3, 2, 1);
        processFace(3, 0, 2);
        processFace(1, 0, 3);
        processFace(2, 0, 1);
      }
    }

    if (allElementFaces == false)  // we build surface volumetricMesh
    {
      for (const auto &p : surfaceFaces) {
        if (p.second == 0)
          continue;  // inner face
        auto key = p.first;
        int numFaces = abs(p.second);
        if (p.second < 0)
          key.reverse();
        for (int i = 0; i < numFaces; i++) {
          faces.emplace_back(std::vector<int>{ key[0], key[1], key[2] });
        }
      }
    }
  }
  else if (volumetricMesh->getElementType() == VolumetricMesh::CUBIC)  // cubic volumetricMesh
  {
    std::unordered_map<Mesh::ORectKey, int> surfaceFaces;
    for (int i = 0; i < volumetricMesh->getNumElements(); i++) {
      auto processFace = [&](int q0, int q1, int q2, int q3) {
        Mesh::ORectKey key(volumetricMesh->getVertexIndex(i, q0), volumetricMesh->getVertexIndex(i, q1), volumetricMesh->getVertexIndex(i, q2), volumetricMesh->getVertexIndex(i, q3));
        if (allElementFaces) {
          faces.emplace_back(std::vector<int>{ key[0], key[1], key[2], key[3] });
          return;
        }
        auto it = surfaceFaces.find(key);
        if (it != surfaceFaces.end())
          it->second++;
        else {
          auto revKey = key.getReversedRectKey();
          it = surfaceFaces.find(revKey);
          if (it != surfaceFaces.end())
            it->second--;
          else
            surfaceFaces.emplace(key, 1);
        }
      };

      processFace(0, 3, 2, 1);
      processFace(4, 5, 6, 7);
      processFace(0, 1, 5, 4);
      processFace(3, 7, 6, 2);
      processFace(1, 2, 6, 5);
      processFace(0, 4, 7, 3);
    }
    if (allElementFaces == false)  // we build surface volumetricMesh
    {
      for (const auto &p : surfaceFaces) {
        if (p.second == 0)
          continue;  // inner face
        auto key = p.first;
        int numFaces = abs(p.second);
        if (p.second < 0)
          key.reverse();
        for (int i = 0; i < numFaces; i++) {
          if (triangulate) {
            faces.emplace_back(std::vector<int>{ key[0], key[1], key[2] });
            faces.emplace_back(std::vector<int>{ key[2], key[3], key[0] });
          }
          else
            faces.emplace_back(std::vector<int>{ key[0], key[1], key[2], key[3] });
        }
      }
    }
  }
  else {
    std::cerr << "Error: unknown VolumetricMesh element type in GenerateSurfaceMesh" << std::endl;
    return;
  }
}