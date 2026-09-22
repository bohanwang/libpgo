/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "interpolationCoordinates" library , Copyright (C) 2018 USC           *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code authors: Yijing Li, Jernej Barbic                                *
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

#pragma once

#include "interpolationCoordinatesBase.h"
#include "EigenDef.h"
#include "volumetricMesh.h"

#include <vector>
#include <string>

namespace pgo
{
namespace InterpolationCoordinates
{

// Compute and store barycentric coordinates.
// Optionally, compute the indices of the elements that the embedded vertices belong to.

class BarycentricCoordinates : public InterpolationCoordinatesBase
{
public:
  BarycentricCoordinates() {}

  // Load from file. No element indices are available. The file stores only
  // (vertex index, weight) pairs, so pass the cage's vertex count when it has
  // vertices beyond the largest referenced index; otherwise it is inferred as
  // the largest referenced index + 1.
  BarycentricCoordinates(const std::string &filename, int numCageVertices = -1);

  // Compute barycentric coordinates and element indices.
  BarycentricCoordinates(int numLocations, const double *locations, const VolumetricMeshes::VolumetricMesh *volumetricMesh);

  // Initialize with only element vertex indices and weights. No element indices are available.
  // numCageVertices < 0 infers the cage size as the largest referenced index + 1.
  BarycentricCoordinates(int numLocations, int numElementVertices, const int *elementVertexIndices, const double *weights,
    const int *elementIndices = nullptr, int numCageVertices = -1);

  BarycentricCoordinates(const std::vector<EigenSupport::V4i> &tetVtxIndices, const std::vector<EigenSupport::V4d> &tetWeights,
    const int *elementIndices = nullptr, int numCageVertices = -1);

  virtual ~BarycentricCoordinates() {}

  virtual void deform(const double *verticesDisp, double *locationDisp) const override;

  int getNumLocations() const { return numLocations; }
  int getNumElementVertices() const { return numElementVertices; }

  const int *getEmbeddingVertexIndices(int embeddedVtx) const { return &indices[numElementVertices * embeddedVtx]; }
  const double *getEmbeddingWeights(int embeddedVtx) const { return &weights[numElementVertices * embeddedVtx]; }

  // Get internal data.
  std::vector<int> &getEmbeddingVertexIndices() { return indices; }
  std::vector<double> &getEmbeddingWeights() { return weights; }
  std::vector<int> &getElements() { return elements; }

  const std::vector<int> &getEmbeddingVertexIndices() const { return indices; }
  const std::vector<double> &getEmbeddingWeights() const { return weights; }
  const std::vector<int> &getElements() const { return elements; }

  // Return embedding element index if it's available; -1 otherwise.
  int getEmbeddingElement(int embeddedVtx) const { return elements[embeddedVtx]; }

  // Save weights to file.
  int saveInterpolationWeights(const std::string &filename) const;

  // Number of columns (vertices) of the cage the weights refer to.
  int getNumCageVertices() const { return numCageVertices; }
  void setNumCageVertices(int numCageVertices_);

  // Throws std::invalid_argument when the cage size is unknown or an index exceeds it.
  EigenSupport::SpMatD generateInterpolationMatrix() const;

protected:
  void initializeInterpolationWeights(int numLocations_, const double *locations, const VolumetricMeshes::VolumetricMesh *volumetricMesh);
  void resolveNumCageVertices(int requested);
  int numLocations = 0;
  int numElementVertices = 0;
  int numCageVertices = 0;

  std::vector<int> indices;
  std::vector<double> weights;
  std::vector<int> elements;
};
}  // namespace InterpolationCoordinates

}  // namespace pgo
