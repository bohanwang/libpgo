/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "macros" library , Copyright (C) 2018 USC                             *
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

#include <cfloat>

namespace pgo
{

// a simple struct to find the index of an entry with min/max value in a search algorithm

struct MaxValueIndex
{
  MaxValueIndex() : value(-DBL_MAX), index(-1) {}
  // record the new index if the newer value is larger
  void update(double newValue, int newIndex) { if (newValue > value) { value = newValue; index = newIndex; } }
  void reset() { value = -DBL_MAX; index = -1; }
  double value;
  int index;
};

struct MinValueIndex
{
  MinValueIndex() : value(DBL_MAX), index(-1) {}
  // record the new index if the newer value is smaller
  void update(double newValue, int newIndex) { if (newValue < value) { value = newValue; index = newIndex; } }
  void reset() { value = DBL_MAX; index = -1; }
  double value;
  int index;
};

template<class T>
struct MaxValueKey
{
  MaxValueKey() : value(-DBL_MAX) {} 
  MaxValueKey(const T & initialKey) : value(-DBL_MAX), key(initialKey) {} 
  void update(double newValue, const T & newKey) { if (newValue > value) { value = newValue; key = newKey; } }
  void reset(const T & initialKey) { value = -DBL_MAX; key = initialKey; }
  double value;
  T key;
};

template<class T>
struct MinValueKey
{
  MinValueKey() : value(DBL_MAX) {} 
  MinValueKey(const T & initialKey) : value(DBL_MAX), key(initialKey) {} 
  void update(double newValue, const T & newKey) { if (newValue < value) { value = newValue; key = newKey; } }
  void reset(const T & initialKey) { value = DBL_MAX; key = initialKey; }
  double value;
  T key;
};

}
