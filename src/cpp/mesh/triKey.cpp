/*
  This code is based on code from the Geometric Tools library,
  which is licensed under a boost license.
  Such usage is permitted by the boost license; for details,
  please see the boost license below.
*/

// Geometric Tools, LLC
// Copyright (c) 1998-2014
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt

/*************************************************************************
 *                                                                       *
 * We release our improvements to the wildMagic code under our standard  *
 * Vega FEM license, as follows:                                         *
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "improvements to the wildMagic library" , Copyright (C) 2018 USC      *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code author: Yijing Li                                                *
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

#include "triKey.h"

const int pgo::Mesh::UTriKey::triEdgeIndex[3][2] = { { 1, 2 }, { 0, 2 }, { 0, 1 } };

const int pgo::Mesh::OTriKey::triEdgeIndex[3][2] = { { 1, 2 }, { 2, 0 }, { 0, 1 } };

void pgo::Mesh::OTriKey::permute(int v0, int v1, int v2, int &r0, int &r1, int &r2)
{
  if (v0 < v1)
  {
    if (v0 < v2)
    {
      r0 = v0; r1 = v1; r2 = v2;
    }
    else
    { // v2 <= v0
      r0 = v2; r1 = v0; r2 = v1;
    }
  }
  else
  { // v1 <= v0
    if (v1 < v2)
    {
      r0 = v1; r1 = v2; r2 = v0;
    }
    else
    { // v2 <= v1
      r0 = v2; r1 = v0; r2 = v1;
    }
  }
}

