/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "basicAlgorithms" library , Copyright (C) 2018 USC                    *
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

#include "averagingBuffer.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <cmath>
#include <algorithm>

using namespace std;

namespace pgo::BasicAlgorithms
{

AveragingBuffer::AveragingBuffer(int size_):
  size(size_), index(0)
{
  assert(size_ > 0);
  buffer.resize(size);
  average = 0.0;
  numValuesAdded = 0;
}

AveragingBuffer::~AveragingBuffer()
{
}

void AveragingBuffer::addValue(double value)
{
  if (numValuesAdded < size) {
    average = (average * numValuesAdded + value) / (numValuesAdded + 1);
    numValuesAdded++;
  }
  else {
    average -= buffer[index] / size;
    average += value / size;
  }

  buffer[index] = value;
  index = (index + 1) % size;
}

void AveragingBuffer::setBufferSize(int newSize)
{
  assert(newSize > 0);
  if (newSize == size)
    return;

  // create a new buffer, copy enough values from old buffer to the new buffer
  vector<double> newBuffer(newSize);
  int numToCopy = min(newSize, numValuesAdded);  // max number of values that can be copied
  for (int i = 0; i < numToCopy; i++)
    index = (index + (size - 1)) % size;         // move index backward #numToCopy times

  average = 0.0;                                 // recompute average value
  for (int i = 0; i < numToCopy; i++)            // move index forward, copy data from old buffer to new buffer
  {
    newBuffer[i] = buffer[index];
    average += buffer[index];
    index = (index + 1) % size;
  }

  // replace old vars with new vars
  buffer = newBuffer;
  index = numToCopy % newSize;
  size = newSize;
  numValuesAdded = numToCopy;
  average /= numValuesAdded;
}

}  // namespace pgo::BasicAlgorithms
