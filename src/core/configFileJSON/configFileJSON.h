/*
author: Bohan Wang
copyright to MIT, USC
*/

#pragma once

#include <nlohmann/json.hpp>

#include <string>
#include <vector>
#include <stdexcept>
#include <sstream>

namespace pgo
{
class ConfigFileJSON
{
public:
  ConfigFileJSON();
  ConfigFileJSON(const char *filename);

  bool open(const char *filename);

  template<typename T>
  T getValue(const char *key, int forceExistance = 0, const T &default_value = T()) const;

  int getInt(const char *key, int forceExistance = 0) const { return getValue<int>(key, forceExistance); }
  float getFloat(const char *key, int forceExistance = 0) const { return getValue<float>(key, forceExistance); }
  double getDouble(const char *key, int forceExistance = 0) const { return getValue<double>(key, forceExistance); }
  std::string getString(const char *key, int forceExistance = 0) const { return getValue<std::string>(key, forceExistance); }
  std::vector<int> getVectorInt(const char *key, int forceExistance = 0) const { return getValue<std::vector<int>>(key, forceExistance); }
  std::vector<double> getVectorDouble(const char *key, int forceExistance = 0) const { return getValue<std::vector<double>>(key, forceExistance); }
  std::vector<std::string> getVectorString(const char *key, int forceExistance = 0) const { return getValue<std::vector<std::string>>(key, forceExistance); }
  std::vector<std::string> getVectorPath(const char *key, int forceExistance = 0) const;

  template<typename T>
  nlohmann::json::reference operator[](T &&key);
  template<typename T>
  nlohmann::json::const_reference operator[](T &&key) const;

  nlohmann::json::reference handle() { return j; }
  nlohmann::json::const_reference handle() const { return j; }

  bool exist(const char *key) const { return j.find(key) != j.end(); }

  std::string getPathFromKey(const char *key, const char *token, const char *dir);
  std::string getPathFromKey(const char *key, const char *workingDir);
  static std::string getPathFromInput(const char *pathString, const char *token, const char *dir);

protected:
  void printError(const char *msg) const;
  void printWarning(const char *msg) const;
  std::string removeComments(const std::string &line) const;
  std::ptrdiff_t findMatchedToken(const std::string &str, std::ptrdiff_t start, char tkLeft, char tkRight) const;

  nlohmann::json j;
};

template<typename T>
inline T ConfigFileJSON::getValue(const char *key, int forceExistance, const T &default_value) const
{
  if (j.find(key) != j.end()) {
    try {
      return j[key].get<T>();
    }
    catch (std::exception &ex) {
      printError(ex.what());
      throw std::invalid_argument("failed to parse the json.");
    }
  }
  else {
    if (forceExistance) {
      std::ostringstream oss;
      oss << key << " does not exist. Throw exception.";
      printError(oss.str().c_str());
      throw std::invalid_argument("Key does not exist");
    }
    else {
      std::ostringstream oss;
      oss << key << " does not exist. Use default value.";
      printWarning(oss.str().c_str());
      return default_value;
    }
  }
}

template<typename T>
inline nlohmann::json::reference ConfigFileJSON::operator[](T &&key)
{
  return j[key];
}

template<typename T>
inline nlohmann::json::const_reference ConfigFileJSON::operator[](T &&key) const
{
  return j[key];
}
}  // namespace pgo
