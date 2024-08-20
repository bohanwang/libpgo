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

#pragma once

#include <fstream>
#include <functional>
#include <string>
#include <vector>
#include <set>
#include <cassert>

namespace pgo
{
namespace BasicAlgorithms
{

// read each line from a file or std::istream
// empty line will be skipped, so input parameter line to processLine is always not empty
// the input function processLine should return 0 when success
// return 0 when success, return the error code of processLine otherwise
// Note: when using filename as input parameter, the function returns 1 when it fails to open the file
int readEachLine(const std::string &filename, std::function<int(std::string &line)> processLine, const char *comment = "#");
int readEachLine(std::istream &in, std::function<int(std::string &line)> processLine, const char *comment = "#");

// use in >> s to read each string s, skipping white spaces
int readEachString(std::istream &in, std::function<int(std::string &s)> processString);

// get the directory part of a path
// examples: a/b/c -> a/b     a -> .
std::string getPathDirectoryName(const std::string &path);

// return each line from fin
// each line will be stripped
// empty line will be skipped
class LineParser
{
public:
  LineParser() {}
  LineParser(std::istream &in) { setStream(in); }

  // defulat comment string: "#"
  void setCommentString(const char *comment) { commentString = comment; }

  void setStream(std::istream &newStream) { this->in = &newStream; }

  bool eof() const { return in == nullptr || in->eof(); }

  std::string &getNextLine();

protected:
  std::istream *in = nullptr;
  std::string line;

  std::string commentString = "#";
};

// =========================================================================
// low-level helper functions to read/write data into binary files
// =========================================================================

inline void writeBinary(std::ostream &o, const int &v)
{
  o.write((char *)&v, sizeof(int));
}

inline void writeBinary(std::ostream &o, const double &v)
{
  o.write((char *)&v, sizeof(double));
}

inline void writeBinary(std::ostream &o, const float &v)
{
  o.write((char *)&v, sizeof(float));
}

template<class T>
void writeBinary(std::ostream &o, int num, const T *array)
{
  writeBinary(o, num);
  for (int i = 0; i < num; i++)
    writeBinary(o, array[i]);
}

template<class T, class A>
void writeBinary(std::ostream &o, const std::vector<T, A> &s)
{
  int size = s.size();
  writeBinary(o, size, s.data());
}

template<class T, class A>
void writeBinary(std::ostream &o, const std::set<T, A> &container)
{
  int size = container.size();
  writeBinary(o, size);
  for (const auto &v : container)
    writeBinary(o, v);
}

inline void readBinary(std::istream &is, int &v)
{
  is.read((char *)&v, sizeof(int));
}

inline void readBinary(std::istream &is, double &v)
{
  is.read((char *)&v, sizeof(double));
}

inline void readBinary(std::istream &is, float &v)
{
  is.read((char *)&v, sizeof(float));
}

template<class T, class A>
void readBinary(std::istream &is, std::set<T, A> &s)
{
  int size = 0;
  readBinary(is, size);
  assert(size >= 0);
  for (int i = 0; i < size; i++) {
    T v = 0.0;
    is.read((char *)&v, sizeof(T));
    s.insert(v);
  }
}

template<class T, class A>
void readBinary(std::istream &is, std::vector<T, A> &s)
{
  int size = 0;
  readBinary(is, size);
  assert(size >= 0);
  for (int i = 0; i < size; i++) {
    T v = 0.0;
    is.read((char *)&v, sizeof(T));
    s.push_back(v);
  }
}

}  // namespace BasicAlgorithms
}  // namespace pgo