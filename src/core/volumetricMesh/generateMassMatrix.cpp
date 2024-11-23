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

#include "generateMassMatrix.h"

#include <tbb/parallel_for.h>
#include <tbb/spin_mutex.h>

using namespace pgo;
using namespace pgo::VolumetricMeshes;

namespace ES = pgo::EigenSupport;

void GenerateMassMatrix::computeMassMatrix(const VolumetricMesh *volumetricMesh, EigenSupport::SpMatD &massMatrix,
  bool inflate3Dim, const double *elementWeight)
{
  int n = volumetricMesh->getNumVertices();
  int numElementVertices = volumetricMesh->getNumElementVertices();

  if (inflate3Dim)
    massMatrix.resize(3 * n, 3 * n);
  else
    massMatrix.resize(n, n);

  std::vector<ES::TripletD> entries;

  if (inflate3Dim) {
    entries.resize(numElementVertices * numElementVertices * 3 * volumetricMesh->getNumElements());
  }
  else {
    entries.resize(numElementVertices * numElementVertices * volumetricMesh->getNumElements());
  }

  tbb::parallel_for(0, volumetricMesh->getNumElements(), [&](int el) {
    thread_local ES::MXd elementMass;

    if (elementMass.rows() == 0) {
      elementMass.resize(numElementVertices, numElementVertices);
    }

    volumetricMesh->computeElementMassMatrix(el, elementMass.data());
    for (int i = 0; i < numElementVertices; i++) {
      int vtxi = volumetricMesh->getVertexIndex(el, i);
      for (int j = 0; j < numElementVertices; j++) {
        int vtxj = volumetricMesh->getVertexIndex(el, j);
        double w = 1.0;
        if (elementWeight)
          w = elementWeight[el];

        double entry = elementMass(i, j) * w;  // since element mass matrix is symmetric
        if (inflate3Dim == false) {
          entries[el * numElementVertices * numElementVertices + i * numElementVertices + j] = ES::TripletD(vtxi, vtxj, entry);
        }
        else {
          for (int d = 0; d < 3; d++) {
            entries[(el * numElementVertices * numElementVertices + i * numElementVertices + j) * 3 + d] = ES::TripletD(vtxi * 3 + d, vtxj * 3 + d, entry);
          }
        }
      }
    }
  });

  massMatrix.setFromTriplets(entries.begin(), entries.end());
}

void GenerateMassMatrix::computeVertexMasses(const VolumetricMesh *volumetricMesh, double *masses, bool inflate3Dim)
{
  int n = volumetricMesh->getNumVertices();
  int numElementVertices = volumetricMesh->getNumElementVertices();
  memset(masses, 0, sizeof(double) * n * (inflate3Dim ? 3 : 1));

  std::vector<tbb::spin_mutex> vtxLocks(n);
  tbb::parallel_for(0, volumetricMesh->getNumElements(), [&](int el) {
    thread_local ES::MXd elementMass;

    if (elementMass.rows() == 0) {
      elementMass.resize(numElementVertices, numElementVertices);
    }

    volumetricMesh->computeElementMassMatrix(el, elementMass.data());
    for (int i = 0; i < numElementVertices; i++) {
      int vtxi = volumetricMesh->getVertexIndex(el, i);
      double vtxMass = 0.0;
      for (int j = 0; j < numElementVertices; j++) {
        vtxMass += elementMass(i, j);  // since element mass matrix is symmetric
      }

      double *massBufferPtr = (inflate3Dim ? &masses[3 * vtxi] : &masses[vtxi]);

      vtxLocks[vtxi].lock();
      *massBufferPtr += vtxMass;
      vtxLocks[vtxi].unlock();
    }
  });

  if (inflate3Dim)
    for (int i = 0; i < n; i++)
      masses[3 * i + 1] = masses[3 * i + 2] = masses[3 * i];
}

void GenerateMassMatrix::computeVertexMassesByAveragingNeighboringElements(const VolumetricMesh *volumetricMesh, double *masses, bool inflate3Dim)
{
  int n = volumetricMesh->getNumVertices();
  int numElementVertices = volumetricMesh->getNumElementVertices();
  double invNumEleVtx = 1.0 / numElementVertices;
  memset(masses, 0, sizeof(double) * n * (inflate3Dim ? 3 : 1));
  std::vector<tbb::spin_mutex> vtxLocks(n);
  tbb::parallel_for(0, volumetricMesh->getNumElements(), [&](int el) {
    double vtxMass = volumetricMesh->getElementVolume(el) * volumetricMesh->getElementDensity(el) * invNumEleVtx;
    for (int i = 0; i < numElementVertices; i++) {
      int vtxi = volumetricMesh->getVertexIndex(el, i);

      double *massBufferPtr = (inflate3Dim ? &masses[3 * vtxi] : &masses[vtxi]);
      vtxLocks[vtxi].lock();
      *massBufferPtr += vtxMass;
      vtxLocks[vtxi].unlock();
    }
  });

  if (inflate3Dim)
    for (int i = 0; i < n; i++)
      masses[3 * i + 1] = masses[3 * i + 2] = masses[3 * i];
}
