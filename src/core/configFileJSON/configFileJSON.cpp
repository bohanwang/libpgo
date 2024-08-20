/*
author: Bohan Wang
copyright to MIT, USC
*/

#include "configFileJSON.h"

#include "pgoLogging.h"

#include <fmt/format.h>

#include <regex>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace pgo;

ConfigFileJSON::ConfigFileJSON()
{
}

ConfigFileJSON::ConfigFileJSON(const char *filename)
{
  if (ConfigFileJSON::open(filename) == false) {
    std::string msg = fmt::format("Cannot find file {}.", filename);
    SPDLOG_LOGGER_ERROR(Logging::lgr(), "{}", msg);

    throw std::domain_error(msg);
  }
}

std::ptrdiff_t ConfigFileJSON::findMatchedToken(const std::string &str, std::ptrdiff_t start, char tkLeft, char tkRight) const
{
  if (str.size() <= (size_t)start)
    return -1;

  if (str[start] != tkLeft)
    return -1;

  std::ptrdiff_t endIdx = (std::ptrdiff_t)str.size();
  if (tkLeft == tkRight) {
    start++;

    // we traverse the string
    while (start < endIdx) {
      // if the char matches the right
      if (str[start] == tkRight && start - 1 >= 0 && str[start - 1] != '\\') {
        return start;
      }

      start++;
    }

    return -1;
  }
  else {
    std::vector<char> tokenStack;

    // we traverse the string
    while (start < endIdx) {
      // if the char matches the right
      if (str[start] == tkRight) {
        // if there are more than one items in the stack, we pop the top of the stack
        if (tokenStack.size() > 1 && tokenStack.back() == tkLeft) {
          tokenStack.pop_back();
        }
        // if there are only one left, then we found the matched token
        else if (tokenStack.size() == 1 && tokenStack.back() == tkLeft) {
          return start;
        }
        // any case else are invalid
        else {
          return -1;
        }
      }
      // if the char matches the left, we push it into the stack
      else if (str[start] == tkLeft) {
        tokenStack.push_back(str[start]);
      }

      start++;
    }

    // any case else are invalid
    return -1;
  }
}

std::string ConfigFileJSON::removeComments(const std::string &line) const
{
  std::ptrdiff_t idx = 0;
  std::ptrdiff_t endIdx = (std::ptrdiff_t)line.size();

  while (idx < endIdx) {
    // if it is a string
    if (line[idx] == '\"') {
      std::ptrdiff_t next = findMatchedToken(line, idx, '\"', '\"');

      if (next <= 0) {
        throw std::domain_error("unmatched \'\"\'");
      }

      idx = next + 1;
    }
    // if it is '/'
    else if (line[idx] == '/') {
      // we found the comment
      if (idx + 1 < endIdx && line[idx + 1] == '/') {
        if (idx > 0)
          return line.substr(0, idx);
        else
          return "";
      }
      // it is not a comment
      else
        idx++;
    }
    else
      idx++;
  }

  return std::string(line);
}

bool ConfigFileJSON::open(const char *filename)
{
  // static std::regex removeComments("\\s*\\/\\/.*");
  std::ifstream infile(filename);
  if (!infile) {
    SPDLOG_LOGGER_ERROR(Logging::lgr(), "Cannot open file {}", filename);
    return false;
  }

  try {
    std::string line;
    std::stringstream ss;

    while (getline(infile, line)) {
      line = removeComments(line);
      ss << line << '\n';
    }

    infile.close();

    ss >> j;

    // LGI << ss.str();
  }
  catch (std::exception &e) {
    SPDLOG_LOGGER_ERROR(Logging::lgr(), "Exception: {}", e.what());
    return false;
  }

  return true;
}

std::string ConfigFileJSON::getPathFromKey(const char *key, const char *token, const char *dir)
{
  std::regex re(fmt::format("\\{{{}\\}}", token));

  std::string path = getValue<std::string>(key, 1);
  return std::regex_replace(path, re, dir);
}

std::string ConfigFileJSON::getPathFromKey(const char *key, const char *workingDir)
{
  static std::regex re("\\{work\\}");

  std::string path = getValue<std::string>(key, 1);
  return std::regex_replace(path, re, workingDir);
}

std::string ConfigFileJSON::getPathFromInput(const char *pathString, const char *token, const char *dir)
{
  std::regex re(fmt::format("\\{{{}\\}}", token));
  return std::regex_replace(pathString, re, dir);
}

std::vector<std::string> ConfigFileJSON::getVectorPath(const char *key, int forceExistance) const
{
  std::vector<std::string> paths = getVectorString(key, forceExistance);

  return paths;
}

void ConfigFileJSON::printError(const char *msg) const
{
  SPDLOG_LOGGER_ERROR(Logging::lgr(), msg);
}

void ConfigFileJSON::printWarning(const char *msg) const
{
  SPDLOG_LOGGER_WARN(Logging::lgr(), msg);
}
