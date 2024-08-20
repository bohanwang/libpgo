/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "mesh" library , Copyright (C) 2018 USC                               *
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

#include "labelOuterTets.h"

#include "pgoLogging.h"

#include <tbb/parallel_for.h>

std::vector<bool> pgo::Mesh::labelOuterTets(const TetMeshRef &tetMesh, const TetNeighbor &tetNeighbor,
  std::function<bool(int tetID)> isTetBoundary, std::function<bool(int tetID)> istetOuter)
{
  enum class TetLabel : int
  {
    TL_IN,
    TL_OUT,
    UNKNOWN
  };

  auto tetBoundaries = tetNeighbor.findTetBoundaries(tetMesh.numTets(), tetMesh.tets());
  std::vector<int> airTetIDs;

  std::vector<TetLabel> tetLabel(tetMesh.numTets(), TetLabel::UNKNOWN);

  //  airProfiler.startTimer("initialAirTetLabel");
  tbb::parallel_for(0, tetMesh.numTets(), [&](int tetID) {
    if (isTetBoundary(tetID)) {
      tetLabel[tetID] = TetLabel::TL_IN;
    }
  });

  for (auto p : tetBoundaries) {
    int tetID = p.first;
    PGO_ALOG(tetID >= 0 && tetID < tetMesh.numTets());
    if (tetLabel[tetID] == TetLabel::UNKNOWN) {
      airTetIDs.push_back(tetID);
      tetLabel[tetID] = TetLabel::TL_OUT;
    }
  }
  //  airProfiler.stopLastTimer();

  //  cout << "Try finding air tets..." << endl;
  size_t candidateBegin = 0, candidateEnd = airTetIDs.size();

  auto floodFillAir = [&]() {
    while (candidateBegin != candidateEnd) {
      for (size_t i = candidateBegin; i < candidateEnd; i++) {
        int tetID = airTetIDs[i];
        for (int nbr : tetNeighbor.getTetNeighbors(tetID)) {
          if (nbr < 0)
            continue;
          PGO_ALOG(nbr < tetMesh.numTets());
          if (tetLabel[nbr] != TetLabel::UNKNOWN)
            continue;
          tetLabel[nbr] = TetLabel::TL_OUT;
          airTetIDs.push_back(nbr);
        }
      }
      candidateBegin = candidateEnd;
      candidateEnd = airTetIDs.size();
    }
  };

  auto floodFillInterior = [&](int seedTetID) {
    std::vector<int> buffer = { seedTetID }, nextBuffer;
    tetLabel[seedTetID] = TetLabel::TL_IN;
    while (buffer.size() > 0) {
      for (int tetID : buffer) {
        for (int nbr : tetNeighbor.getTetNeighbors(tetID)) {
          if (nbr < 0)
            continue;
          PGO_ALOG(nbr < tetMesh.numTets());
          if (tetLabel[nbr] != TetLabel::UNKNOWN)
            continue;
          tetLabel[nbr] = TetLabel::TL_IN;
          nextBuffer.push_back(nbr);
        }
      }
      buffer.swap(nextBuffer);
      nextBuffer.clear();
    }
  };
  //  airProfiler.startTimer("initialFloodFill");
  floodFillAir();

  // now process tets that are in the interior holes
  for (int tetID = 0; tetID < tetMesh.numTets(); tetID++) {
    if (tetLabel[tetID] != TetLabel::UNKNOWN)
      continue;
    if (istetOuter(tetID))  // outside
    {
      candidateBegin = airTetIDs.size();
      airTetIDs.push_back(tetID);
      tetLabel[tetID] = TetLabel::TL_OUT;
      candidateEnd = airTetIDs.size();
      floodFillAir();
    }
    else  // inside
    {
      floodFillInterior(tetID);
    }
  }

  std::vector<bool> ret(tetMesh.numTets());
  for (int tetID = 0; tetID < tetMesh.numTets(); tetID++)
    ret[tetID] = (tetLabel[tetID] == TetLabel::TL_OUT);

  return ret;
}
