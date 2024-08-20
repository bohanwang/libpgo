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

#include "containerHelper.h"

using namespace std;

namespace pgo::BasicAlgorithms
{

template<class T>
bool saveToAscii(const std::vector<T> &v, std::ostream &out)
{
  out << v.size() << " ";
  if (out.fail())
    return false;
  for (size_t i = 0; i < v.size(); i++) {
    out << v[i] << " ";
    if (out.fail())
      return false;
  }
  return true;
}

template<class T>
bool loadFromAscii(std::vector<T> &v, std::istream &in)
{
  size_t num = 0;
  in >> num;
  if (in.fail())
    return false;
  v.reserve(v.size() + num);
  for (size_t i = 0; i < num; i++) {
    T k = T();
    in >> k;
    if (in.fail())
      return false;
    v.push_back(k);
  }
  return true;
}

template bool saveToAscii<int>(const std::vector<int> &v, std::ostream &out);
template bool saveToAscii<double>(const std::vector<double> &v, std::ostream &out);

template bool loadFromAscii<int>(std::vector<int> &v, std::istream &in);
template bool loadFromAscii<double>(std::vector<double> &v, std::istream &in);

}  // namespace pgo::BasicAlgorithms
