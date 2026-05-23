#include "app/output.h"

#include "deformationModelAssembler.h"
#include "simulationMesh.h"

#include <fmt/format.h>
#include <nlohmann/json.hpp>

#include <fstream>
#include <stdexcept>
#include <utility>
#include <vector>

namespace pgo::RunIPCSim
{
namespace ES = pgo::EigenSupport;

namespace
{
const double *dataOrNull(const ES::VXd &values)
{
  return values.size() > 0 ? values.data() : nullptr;
}

void createOutputSubdirectories(const OutputDirectories &outputDirs)
{
  std::error_code ec;
  std::filesystem::create_directories(outputDirs.root, ec);
  if (ec)
    throw std::runtime_error("Failed to create output folder `" + outputDirs.root.string() + "`: " + ec.message());

  for (const auto &dir : { outputDirs.states, outputDirs.surface, outputDirs.stress }) {
    std::filesystem::create_directories(dir, ec);
    if (ec)
      throw std::runtime_error("Failed to create output subfolder `" + dir.string() + "`: " + ec.message());
  }
}

void clearOutputDirectory(const std::filesystem::path &outputFolder)
{
  std::error_code ec;
  std::filesystem::remove_all(outputFolder, ec);
  if (ec)
    throw std::runtime_error("Failed to clear output folder `" + outputFolder.string() + "`: " + ec.message());

  std::filesystem::create_directories(outputFolder, ec);
  if (ec)
    throw std::runtime_error("Failed to create output folder `" + outputFolder.string() + "`: " + ec.message());
}

const char *vonMisesStressLocation(const pgo::SolidDeformationModel::SimulationMesh &mesh)
{
  using pgo::SolidDeformationModel::SimulationMeshType;
  switch (mesh.getElementType()) {
    case SimulationMeshType::TET:
      return "tet_element";
    case SimulationMeshType::CUBIC:
      return "cubic_element";
    default:
      return "element";
  }
}
}  // namespace

std::filesystem::path framePath(const std::filesystem::path &dir, const char *prefix, int frame, const char *extension)
{
  return dir / fmt::format("{}{:04d}{}", prefix, frame, extension);
}

RunIPCSimOutput::RunIPCSimOutput(std::filesystem::path outputFolder):
  outputDirs_{
    std::move(outputFolder),
    {},
    {},
    {},
  }
{
  outputDirs_.states = outputDirs_.root / "states";
  outputDirs_.surface = outputDirs_.root / "surface";
  outputDirs_.stress = outputDirs_.root / "stress";
}

std::filesystem::path RunIPCSimOutput::logPath() const
{
  return outputDirs_.root / "runIPCSim.log";
}

std::filesystem::path RunIPCSimOutput::statePath(int frame) const
{
  return framePath(outputDirs_.states, "deform", frame, ".u");
}

std::filesystem::path RunIPCSimOutput::surfacePath(int outputFrame) const
{
  return framePath(outputDirs_.surface, "ret", outputFrame, ".obj");
}

std::filesystem::path RunIPCSimOutput::stressPath(int frame) const
{
  return framePath(outputDirs_.stress, "von_mises", frame, ".json");
}

void RunIPCSimOutput::prepare(bool restartFromU) const
{
  if (!restartFromU)
    clearOutputDirectory(outputDirs_.root);
  createOutputSubdirectories(outputDirs_);
}

int RunIPCSimOutput::loadLatestRestartState(int numSimSteps, int n3,
  ES::VXd &u, ES::VXd &uvel, ES::VXd &uacc) const
{
  for (int framei = numSimSteps - 1; framei >= 0; --framei) {
    const std::filesystem::path deformFilename = statePath(framei);
    if (!std::filesystem::exists(deformFilename))
      continue;

    ES::MXd uMat(n3, 3);
    if (ES::readMatrix(deformFilename.string().c_str(), uMat) == 0) {
      u.noalias() = uMat.col(0);
      uvel.noalias() = uMat.col(1);
      uacc.noalias() = uMat.col(2);
      return framei;
    }
  }

  return -1;
}

void RunIPCSimOutput::writeState(int frame, const ES::VXd &u, const ES::VXd &uvel, const ES::VXd &uacc) const
{
  ES::MXd uMat(u.size(), 3);
  uMat.col(0) = u;
  uMat.col(1) = uvel;
  uMat.col(2) = uacc;
  ES::writeMatrix(statePath(frame).string().c_str(), uMat);
}

void RunIPCSimOutput::writeSurface(int outputFrame, const pgo::Mesh::TriMeshGeo &mesh) const
{
  mesh.save(surfacePath(outputFrame).string());
}

void RunIPCSimOutput::writeStateAndSurfaceFrame(
  int frame,
  int outputFrame,
  const IpcSimulationContext &context,
  const ES::VXd &u,
  const ES::VXd &uvel,
  const ES::VXd &uacc,
  double scale,
  bool writeStateFile,
  bool writeSurfaceFile) const
{
  if (writeStateFile)
    writeState(frame, u, uvel, uacc);

  if (!writeSurfaceFile)
    return;

  pgo::Mesh::TriMeshGeo mesh = context.surfaceMesh;
  ES::VXd usurf(context.surfaceRestPositions.size());
  ES::mv(context.surfaceFromSimulationDispMap, u, usurf);
  const ES::VXd psurf = context.surfaceRestPositions + usurf;
  for (int vi = 0; vi < mesh.numVertices(); ++vi)
    mesh.pos(vi) = psurf.segment<3>(vi * 3) / scale;
  writeSurface(outputFrame, mesh);
}

void RunIPCSimOutput::writeVonMisesStressJson(int frame, double timestep,
  const IpcSimulationContext &context, const ES::VXd &displacement) const
{
  if (!context.deformationModelAssemblerOwner || !context.simulationMeshOwner)
    throw std::runtime_error("runIPCSim cannot output von Mises stresses without a simulation mesh and assembler.");

  const int elementCount = context.simulationMeshOwner->getNumElements();
  std::vector<double> elementStresses(elementCount, 0.0);
  const ES::VXd absolutePositions = context.simulationRestPosition + displacement;
  context.deformationModelAssemblerOwner->computeVonMisesStresses(
    absolutePositions.data(),
    dataOrNull(context.plasticParams),
    dataOrNull(context.elasticParams),
    elementStresses.data());

  nlohmann::json stressJson;
  stressJson["frame"] = frame;
  stressJson["time"] = static_cast<double>(frame) * timestep;
  stressJson["stress_type"] = "von_mises";
  stressJson["location"] = vonMisesStressLocation(*context.simulationMeshOwner);
  stressJson["values"] = elementStresses;

  const std::filesystem::path outputPath = stressPath(frame);
  std::ofstream out(outputPath);
  if (!out.is_open())
    throw std::runtime_error("Failed to write von Mises stress JSON: " + outputPath.string());
  out << stressJson.dump(2) << '\n';
}
}  // namespace pgo::RunIPCSim
