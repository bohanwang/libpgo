#include "surfaceQuality.h"

#include <argparse/argparse.hpp>
#include <nlohmann/json.hpp>

#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

namespace
{

void writeJsonReport(const std::string &path, const nlohmann::json &report)
{
  const std::filesystem::path outputPath(path);
  if (outputPath.has_parent_path())
    std::filesystem::create_directories(outputPath.parent_path());

  std::ofstream out(outputPath);
  if (out.is_open() == false)
    throw std::runtime_error("Failed to open JSON report for writing: " + path);

  out << report.dump(2) << '\n';
}

}  // namespace

int main(int argc, char *argv[])
{
  argparse::ArgumentParser program("meshQualityCheck");

  argparse::ArgumentParser surfaceCommand("surface");
  surfaceCommand.add_description("Check surface triangle mesh quality");
  surfaceCommand.add_argument("--input")
    .help("Input surface mesh filename")
    .required()
    .metavar("PATH");
  surfaceCommand.add_argument("--json")
    .help("Output JSON report filename")
    .required()
    .metavar("PATH");
  surfaceCommand.add_argument("--check-level")
    .help("Check level: full, raw, or degenerate-only")
    .default_value(std::string("full"))
    .metavar("LEVEL");
  surfaceCommand.add_argument("--self-intersection-triangle-limit")
    .help("Maximum triangle count for exact self-intersection checking")
    .default_value(200000)
    .metavar("INT")
    .scan<'i', int>();
  surfaceCommand.add_argument("--expected-components")
    .help("Expected number of edge-connected triangle components when topology is checked")
    .default_value(-1)
    .metavar("INT")
    .scan<'i', int>();
  surfaceCommand.add_argument("--invalid-triangles-policy")
    .help("How invalid or degenerate triangles affect pass/fail: fail or warn. Omit to preserve the check-level default.")
    .default_value(std::string(""))
    .metavar("POLICY");

  program.add_subparser(surfaceCommand);

  try {
    program.parse_args(argc, argv);
  }
  catch (const std::exception &err) {
    std::cerr << err.what() << std::endl;
    std::cerr << program;
    return 1;
  }

  if (program.is_subcommand_used(surfaceCommand) == false) {
    std::cerr << program;
    return 1;
  }

  try {
    const std::string inputMesh = surfaceCommand.get<std::string>("--input");
    const std::string jsonReport = surfaceCommand.get<std::string>("--json");
    const std::string checkLevel = surfaceCommand.get<std::string>("--check-level");
    const int selfIntersectionTriangleLimit = surfaceCommand.get<int>("--self-intersection-triangle-limit");
    const int expectedComponents = surfaceCommand.get<int>("--expected-components");
    const std::string invalidTrianglesPolicy = surfaceCommand.get<std::string>("--invalid-triangles-policy");

    mesh_quality_check::SurfaceQualityOptions options;
    if (selfIntersectionTriangleLimit < 0)
      throw std::runtime_error("--self-intersection-triangle-limit must be non-negative");
    options.selfIntersectionTriangleLimit = selfIntersectionTriangleLimit;
    if (expectedComponents < -1)
      throw std::runtime_error("--expected-components must be non-negative");
    if (expectedComponents >= 0)
      options.expectedComponents = expectedComponents;
    if (checkLevel == "degenerate-only") {
      options.checkTopology = false;
      options.checkWinding = false;
      options.checkSelfIntersection = false;
    }
    else if (checkLevel == "raw") {
      options.checkTopology = true;
      options.checkWinding = false;
      options.checkSelfIntersection = false;
      options.failOnInvalidTriangles = false;
    }
    else if (checkLevel != "full") {
      throw std::runtime_error("--check-level must be one of: full, raw, degenerate-only");
    }
    if (invalidTrianglesPolicy == "fail") {
      options.failOnInvalidTriangles = true;
    }
    else if (invalidTrianglesPolicy == "warn") {
      options.failOnInvalidTriangles = false;
    }
    else if (invalidTrianglesPolicy.empty() == false) {
      throw std::runtime_error("--invalid-triangles-policy must be one of: fail, warn");
    }

    const mesh_quality_check::SurfaceQualityReport report = mesh_quality_check::checkSurfaceMesh(inputMesh, options);
    writeJsonReport(jsonReport, mesh_quality_check::surfaceQualityReportToJson(report));

    if (report.passed == false) {
      for (const std::string &error : report.errors)
        std::cerr << error << std::endl;
      return 1;
    }

    return 0;
  }
  catch (const std::exception &err) {
    std::cerr << err.what() << std::endl;
    return 1;
  }
}
