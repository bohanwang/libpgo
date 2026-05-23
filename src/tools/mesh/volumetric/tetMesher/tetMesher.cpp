#include "tetMesherBackend.h"

#include "configFileJSON.h"
#include "pgoLogging.h"

#include <argparse/argparse.hpp>

#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

namespace
{

using json = nlohmann::json;

std::string requireString(const json &j, const char *key)
{
  if (j.contains(key) == false)
    throw std::runtime_error(std::string("Missing required JSON field: ") + key);

  if (j.at(key).is_string() == false)
    throw std::runtime_error(std::string("JSON field must be a string: ") + key);

  return j.at(key).get<std::string>();
}

int readVersion(const json &j)
{
  if (j.contains("version") == false)
    return 1;

  if (j.at("version").is_number_integer() == false)
    throw std::runtime_error("JSON field must be an integer: version");

  return j.at("version").get<int>();
}

bool readBool(const json &j, const char *key, bool defaultValue)
{
  if (j.contains(key) == false)
    return defaultValue;

  if (j.at(key).is_boolean() == false)
    throw std::runtime_error(std::string("JSON field must be a boolean: ") + key);

  return j.at(key).get<bool>();
}

double readDouble(const json &j, const char *key, double defaultValue)
{
  if (j.contains(key) == false)
    return defaultValue;

  if (j.at(key).is_number() == false)
    throw std::runtime_error(std::string("JSON field must be a number: ") + key);

  return j.at(key).get<double>();
}

int readInt(const json &j, const char *key, int defaultValue)
{
  if (j.contains(key) == false)
    return defaultValue;

  if (j.at(key).is_number_integer() == false)
    throw std::runtime_error(std::string("JSON field must be an integer: ") + key);

  return j.at(key).get<int>();
}

const json &readObject(const json &j, const char *key, bool required)
{
  static const json empty = json::object();

  if (j.contains(key) == false) {
    if (required)
      throw std::runtime_error(std::string("Missing required JSON object: ") + key);
    return empty;
  }

  if (j.at(key).is_object() == false)
    throw std::runtime_error(std::string("JSON field must be an object: ") + key);

  return j.at(key);
}

tet_mesher::MaterialOptions readMaterialOptions(const json &j)
{
  tet_mesher::MaterialOptions material;
  if (j.contains("material") == false)
    return material;

  const json &materialJson = readObject(j, "material", false);
  material.enabled = true;
  material.density = readDouble(materialJson, "density", material.density);
  material.youngModulus = readDouble(materialJson, "young_modulus", material.youngModulus);
  material.poissonRatio = readDouble(materialJson, "poisson_ratio", material.poissonRatio);

  if (std::isfinite(material.density) == false || material.density <= 0.0)
    throw std::runtime_error("material.density must be positive");

  if (std::isfinite(material.youngModulus) == false || material.youngModulus <= 0.0)
    throw std::runtime_error("material.young_modulus must be positive");

  if (std::isfinite(material.poissonRatio) == false ||
    material.poissonRatio <= -1.0 || material.poissonRatio >= 0.5) {
    throw std::runtime_error("material.poisson_ratio must be in (-1, 0.5)");
  }

  return material;
}

tet_mesher::CommonOptions readCommonOptions(const pgo::ConfigFileJSON &config)
{
  const json &j = config.handle();

  tet_mesher::CommonOptions options;
  options.inputMesh = config.resolvePath(requireString(j, "input_mesh"));
  options.outputMesh = config.resolvePath(requireString(j, "output_mesh"));
  if (j.contains("output_surface"))
    options.outputSurface = config.resolvePath(requireString(j, "output_surface"));
  options.material = readMaterialOptions(j);
  options.printStats = readBool(j, "print_stats", false);
  options.quiet = readBool(j, "quiet", false);
  return options;
}

int runTetgen(const pgo::ConfigFileJSON &config)
{
  const json &tetgen = readObject(config.handle(), "tetgen", true);

  tet_mesher::TetgenOptions options;
  options.common = readCommonOptions(config);
  options.command = requireString(tetgen, "command");

  std::unique_ptr<pgo::VolumetricMeshes::TetMesh> tetMesh = tet_mesher::generateTetgenMesh(options);
  if (options.common.material.enabled)
    tetMesh->setSingleMaterial(options.common.material.youngModulus,
      options.common.material.poissonRatio, options.common.material.density);
  tet_mesher::saveTetMeshOutputs(*tetMesh, options.common);
  return 0;
}

int runTetwild(const pgo::ConfigFileJSON &config)
{
  const json &tetwild = readObject(config.handle(), "tetwild", false);

  tet_mesher::TetwildOptions options;
  options.common = readCommonOptions(config);
  options.lr = readDouble(tetwild, "lr", options.lr);
  options.epsr = readDouble(tetwild, "epsr", options.epsr);
  options.stopEnergy = readDouble(tetwild, "stop_energy", options.stopEnergy);
  options.maxThreads = readInt(tetwild, "max_threads", options.maxThreads);

  if (tetwild.contains("la")) {
    options.la = readDouble(tetwild, "la", options.la);
    options.hasLa = true;
  }

  if (options.hasLa && tetwild.contains("lr"))
    throw std::runtime_error("tetwild.la and tetwild.lr are mutually exclusive");

  std::unique_ptr<pgo::VolumetricMeshes::TetMesh> tetMesh = tet_mesher::generateTetwildMesh(options);
  if (options.common.material.enabled)
    tetMesh->setSingleMaterial(options.common.material.youngModulus,
      options.common.material.poissonRatio, options.common.material.density);
  tet_mesher::saveTetMeshOutputs(*tetMesh, options.common);
  return 0;
}

}  // namespace

int main(int argc, char *argv[])
{
  argparse::ArgumentParser program("tetMesher");
  program.add_description("Tetrahedralize a surface mesh into a .veg simulation mesh from a JSON config");
  program.add_argument("--config")
    .help("Tet meshing JSON config filename")
    .required()
    .metavar("PATH");

  try {
    program.parse_args(argc, argv);
  }
  catch (const std::exception &err) {
    std::cerr << err.what() << std::endl;
    std::cerr << program;
    return 1;
  }

  pgo::Logging::init();

  try {
    pgo::ConfigFileJSON config;
    const std::string configPath = program.get<std::string>("--config");
    if (config.open(configPath.c_str()) == false)
      throw std::runtime_error("Failed to open tet mesher config: " + configPath);

    const int version = readVersion(config.handle());
    if (version != 1)
      throw std::runtime_error("Unsupported tet mesher config version: " + std::to_string(version));

    const std::string backend = requireString(config.handle(), "backend");
    if (backend == "tetgen")
      return runTetgen(config);

    if (backend == "tetwild")
      return runTetwild(config);

    throw std::runtime_error("Unsupported tet mesher backend: " + backend);
  }
  catch (const std::exception &err) {
    std::cerr << err.what() << std::endl;
    return 1;
  }
}
