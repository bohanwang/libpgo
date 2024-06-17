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

#include <string>
#include <cstdarg>

namespace pgo
{
namespace BasicAlgorithms
{

// check with str ends with substr
bool endWith(const std::string &str, const std::string &substr);
// version with case-insensitive comparison
bool iendWith(const std::string &str, const std::string &substr);
bool iendWith(const std::string &str, const char *substr);

// strip a string, removing white-space characters at the beginning and end of the string
// if the string contains only white-space characters, return an empty string
void stripSelf(std::string &s);
std::string strip(const std::string &s);

// strip a c-string s with minimal modification to the memory
// remove the end white-space and return the beginning of the stripped s
char *stripLight(char *s);

inline void skipSpace(char *&s)
{
  while (std::isspace(*s))
    s++;
}

// converts a string to upper case
void upperCase(char *s);
void upperCase(std::string &s);

// convert the args in printf(fmt, ...) to a string
std::string format2string(const char *fmt, ...);
// the second input is a va_list which represents a group of arguments for fmt
std::string format2string(const char *fmt, va_list args);

}  // namespace BasicAlgorithms
}  // namespace pgo
