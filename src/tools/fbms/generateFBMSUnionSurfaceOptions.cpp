#include "generateFBMSUnionSurfaceOptions.h"

#include <argparse/argparse.hpp>

#include <cmath>
#include <limits>
#include <stdexcept>

namespace pgo::Tools::FBMSUnionSurface
{
namespace
{

void addCommonOptions(argparse::ArgumentParser &program)
{
  program.add_argument("--fbms")
    .help("Input FBMS triangle surface OBJ")
    .required()
    .metavar("PATH");
  program.add_argument("--sphere")
    .help("Input bounding sphere OBJ")
    .required()
    .metavar("PATH");
  program.add_argument("--fbms-thickness")
    .help("FBMS unsigned-distance thickening width. Positive: use as-is. "
          "Negative V: treat |V| as a target volume budget and binary-search "
          "the largest thickness whose union mesh has volume <= |V|.")
    .required()
    .metavar("FLOAT")
    .scan<'g', double>();
  program.add_argument("--sphere-thickness")
    .help("Bounding sphere shell thickening width")
    .required()
    .metavar("FLOAT")
    .scan<'g', double>();
  program.add_argument("--resolution")
    .help("Uniform SDF grid resolution per axis")
    .required()
    .metavar("INT")
    .scan<'i', int>();
  program.add_argument("--padding-ratio")
    .help("BBox diagonal ratio added to each side after thickness expansion")
    .required()
    .metavar("FLOAT")
    .scan<'g', double>();
  program.add_argument("--output-surface")
    .help("Output raw union surface OBJ")
    .required()
    .metavar("PATH");
  program.add_argument("--enable-truncating")
    .help("Truncate the FBMS shell by the bounding sphere mid-surface; "
          "only FBMS structure inside the sphere radius is kept.")
    .default_value(false)
    .implicit_value(true);
  program.add_argument("--project-fbms-boundary-to-sphere")
    .help("Before thickening, project open FBMS boundary vertices to the "
          "bounding sphere radius so the sphere shell consumes the thickened "
          "boundary side wall.")
    .default_value(false)
    .implicit_value(true);
  program.add_argument("--debug-field-mode")
    .help("Diagnostic field to polygonize: union, fbms, sphere, or union-minus-sphere")
    .default_value(std::string("union"))
    .metavar("MODE");
  program.add_argument("--surface-mode")
    .help("Surface field to polygonize. Overrides --debug-field-mode when set: union, fbms, sphere, or union-minus-sphere")
    .default_value(std::string(""))
    .metavar("MODE");
  program.add_argument("--filter-small-components")
    .help("Remove tiny edge-connected triangle components after surface extraction")
    .default_value(false)
    .implicit_value(true);
  program.add_argument("--min-component-triangles")
    .help("Minimum triangle count kept by --filter-small-components")
    .default_value(100)
    .metavar("INT")
    .scan<'i', int>();
  program.add_argument("--keep-largest-components")
    .help("With --filter-small-components, optionally keep only the largest N components after thresholding; -1 keeps all components above threshold")
    .default_value(-1)
    .metavar("INT")
    .scan<'i', int>();
}

void addMarchingCubesOptions(argparse::ArgumentParser &program)
{
  program.add_argument("--iso-offset-mode")
    .help("Iso-surface extraction mode: zero, fixed, or auto. zero keeps the historical 0-level extraction. fixed uses --iso-offset. auto chooses a small nearby iso level away from sampled grid values.")
    .default_value(std::string("zero"))
    .metavar("MODE");
  program.add_argument("--iso-offset")
    .help("Iso level used when --iso-offset-mode fixed is selected")
    .default_value(0.0)
    .metavar("FLOAT")
    .scan<'g', double>();
}

void addOpenVDBOptions(argparse::ArgumentParser &program)
{
  program.add_argument("--vdb-voxel-size")
    .help("OpenVDB voxel size. If 0, derive it from --resolution and the expanded bbox.")
    .default_value(0.0)
    .metavar("FLOAT")
    .scan<'g', double>();
  program.add_argument("--vdb-half-width")
    .help("OpenVDB level-set half width in voxels")
    .default_value(3.0)
    .metavar("FLOAT")
    .scan<'g', double>();
  program.add_argument("--vdb-adaptivity")
    .help("OpenVDB volumeToMesh adaptivity in [0, 1]")
    .default_value(0.0)
    .metavar("FLOAT")
    .scan<'g', double>();
  program.add_argument("--vdb-smooth-steps")
    .help("Number of OpenVDB mean-curvature smoothing iterations before meshing")
    .default_value(0)
    .metavar("INT")
    .scan<'i', int>();
}

void readCommonOptions(const argparse::ArgumentParser &program, Options &options)
{
  options.fbmsPath = program.get<std::string>("--fbms");
  options.spherePath = program.get<std::string>("--sphere");
  options.fbmsThickness = program.get<double>("--fbms-thickness");
  options.sphereThickness = program.get<double>("--sphere-thickness");
  options.resolution = program.get<int>("--resolution");
  options.paddingRatio = program.get<double>("--padding-ratio");
  options.outputSurfacePath = program.get<std::string>("--output-surface");
  options.enableTruncating = program.get<bool>("--enable-truncating");
  options.projectFBMSBoundaryToSphere = program.get<bool>("--project-fbms-boundary-to-sphere");
  options.debugFieldMode = program.get<std::string>("--debug-field-mode");
  const std::string surfaceMode = program.get<std::string>("--surface-mode");
  if (!surfaceMode.empty())
    options.debugFieldMode = surfaceMode;
  options.filterSmallComponents = program.get<bool>("--filter-small-components");
  options.minComponentTriangles = program.get<int>("--min-component-triangles");
  options.keepLargestComponents = program.get<int>("--keep-largest-components");
}

void readMarchingCubesOptions(const argparse::ArgumentParser &program, Options &options)
{
  options.isoOffsetMode = program.get<std::string>("--iso-offset-mode");
  options.isoOffset = program.get<double>("--iso-offset");
}

void readOpenVDBOptions(const argparse::ArgumentParser &program, Options &options)
{
  options.vdbVoxelSize = program.get<double>("--vdb-voxel-size");
  options.vdbHalfWidth = program.get<double>("--vdb-half-width");
  options.vdbAdaptivity = program.get<double>("--vdb-adaptivity");
  options.vdbSmoothSteps = program.get<int>("--vdb-smooth-steps");
}

bool isBackendSubcommand(const char *arg)
{
  const std::string value(arg == nullptr ? "" : arg);
  return value == "marching-cubes" || value == "openvdb";
}

Options parseSubcommandOptions(int argc, char *argv[])
{
  argparse::ArgumentParser program("generateFBMSUnionSurface");

  argparse::ArgumentParser marchingCubesCommand("marching-cubes");
  marchingCubesCommand.add_description("Extract the union surface with the historical dense grid marching-cubes backend");
  addCommonOptions(marchingCubesCommand);
  addMarchingCubesOptions(marchingCubesCommand);

  argparse::ArgumentParser openVDBCommand("openvdb");
  openVDBCommand.add_description("Extract the union surface with OpenVDB sparse level-set meshing");
  addCommonOptions(openVDBCommand);
  addOpenVDBOptions(openVDBCommand);

  program.add_subparser(marchingCubesCommand);
  program.add_subparser(openVDBCommand);
  program.parse_args(argc, argv);

  Options options;
  if (program.is_subcommand_used(marchingCubesCommand)) {
    readCommonOptions(marchingCubesCommand, options);
    readMarchingCubesOptions(marchingCubesCommand, options);
    options.extractionBackend = "marching-cubes";
  }
  else if (program.is_subcommand_used(openVDBCommand)) {
    readCommonOptions(openVDBCommand, options);
    readOpenVDBOptions(openVDBCommand, options);
    options.extractionBackend = "openvdb";
  }
  else {
    throw std::runtime_error("Expected subcommand: marching-cubes or openvdb");
  }

  return options;
}

Options parseLegacyOptions(int argc, char *argv[])
{
  argparse::ArgumentParser program("generateFBMSUnionSurface");
  addCommonOptions(program);
  program.add_argument("--extraction-backend")
    .help("Surface extraction backend: marching-cubes or openvdb")
    .default_value(std::string("marching-cubes"))
    .metavar("BACKEND");
  addMarchingCubesOptions(program);
  addOpenVDBOptions(program);

  program.parse_args(argc, argv);

  Options options;
  readCommonOptions(program, options);
  options.extractionBackend = program.get<std::string>("--extraction-backend");
  readMarchingCubesOptions(program, options);
  readOpenVDBOptions(program, options);
  return options;
}

}  // namespace

Options parseOptions(int argc, char *argv[])
{
  if (argc > 1 && isBackendSubcommand(argv[1]))
    return parseSubcommandOptions(argc, argv);

  return parseLegacyOptions(argc, argv);
}

void validateOptions(const Options &options)
{
  if (options.fbmsThickness == 0.0 || !std::isfinite(options.fbmsThickness))
    throw std::runtime_error("--fbms-thickness must be a finite non-zero value "
                             "(positive = explicit thickness; negative = volume budget)");
  if (options.sphereThickness <= 0.0 || !std::isfinite(options.sphereThickness))
    throw std::runtime_error("--sphere-thickness must be positive");
  if (options.resolution < 2)
    throw std::runtime_error("--resolution must be at least 2");
  if (options.paddingRatio < 0.0 || !std::isfinite(options.paddingRatio))
    throw std::runtime_error("--padding-ratio must be finite and non-negative");
  if (options.debugFieldMode != "union" && options.debugFieldMode != "fbms" &&
    options.debugFieldMode != "sphere" && options.debugFieldMode != "union-minus-sphere")
    throw std::runtime_error("--debug-field-mode must be one of: union, fbms, sphere, union-minus-sphere");
  if (options.fbmsThickness < 0.0 && options.debugFieldMode != "union")
    throw std::runtime_error("Volume-budget mode only supports --debug-field-mode union");
  if (options.minComponentTriangles < 1)
    throw std::runtime_error("--min-component-triangles must be at least 1");
  if (options.keepLargestComponents != -1 && options.keepLargestComponents < 1)
    throw std::runtime_error("--keep-largest-components must be -1 or a positive integer");
  if (options.extractionBackend != "marching-cubes" && options.extractionBackend != "openvdb")
    throw std::runtime_error("--extraction-backend must be one of: marching-cubes, openvdb");
  if (options.isoOffsetMode != "zero" && options.isoOffsetMode != "fixed" && options.isoOffsetMode != "auto")
    throw std::runtime_error("--iso-offset-mode must be one of: zero, fixed, auto");
  if (!std::isfinite(options.isoOffset))
    throw std::runtime_error("--iso-offset must be finite");
  if (options.isoOffsetMode != "fixed" && options.isoOffset != 0.0)
    throw std::runtime_error("--iso-offset is only used with --iso-offset-mode fixed");
  if (options.extractionBackend == "openvdb" && options.isoOffsetMode != "zero")
    throw std::runtime_error("--iso-offset-mode is only supported by the marching-cubes backend");
  if (options.vdbVoxelSize < 0.0 || !std::isfinite(options.vdbVoxelSize))
    throw std::runtime_error("--vdb-voxel-size must be finite and non-negative");
  if (options.vdbHalfWidth <= 0.0 || !std::isfinite(options.vdbHalfWidth))
    throw std::runtime_error("--vdb-half-width must be positive");
  if (options.vdbAdaptivity < 0.0 || options.vdbAdaptivity > 1.0 || !std::isfinite(options.vdbAdaptivity))
    throw std::runtime_error("--vdb-adaptivity must be finite and in [0, 1]");
  if (options.vdbSmoothSteps < 0)
    throw std::runtime_error("--vdb-smooth-steps must be non-negative");

  const long long gridCount = 1LL * options.resolution * options.resolution * options.resolution;
  if (gridCount > std::numeric_limits<int>::max())
    throw std::runtime_error("--resolution is too large for the current marching cubes interface");
}

}  // namespace pgo::Tools::FBMSUnionSurface
