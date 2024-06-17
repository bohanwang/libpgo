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

#pragma once

#include <vector>
#include <string>

/*
  A parser for the volumetric mesh text file format.
  Note: the end user never needs to use this class directly.
  See volumetricMesh.h .
*/

namespace pgo
{
namespace VolumetricMeshes
{

class VolumetricMeshParser
{
public:
  VolumetricMeshParser(const char *includeToken = nullptr);  // pass nullptr for normal usage
  ~VolumetricMeshParser();

  int open(const char *filename);

  // return the next line, s must be externally allocated string
  // if last line, return will be nullptr
  char *getNextLine(char *s, int numRetainedSpaces = 0, int removeWhitespace = 1);

  void rewindToStart();
  void close();

  static void upperCase(char *s);
  static void removeWhitespace(char *s, int numRetainedSpaces = 0);                    // any whitespace equal in length or longer to "numRetainedSpaces" is shrunk to "numRetainedSpaces" and retained
  static void beautifyLine(char *s, int numRetainedSpaces, int removeWhitespace = 1);  // strip whitespace + removes trailing "\n"

protected:
  FILE *fin = nullptr;
  std::vector<FILE *> fileStack;
  // int fileStackDepth;

  std::string directoryName;

  char includeToken[96];   // normally "*INCLUDE "
  int includeTokenLength;  // normally 9
};

}  // namespace VolumetricMeshes
}  // namespace pgo