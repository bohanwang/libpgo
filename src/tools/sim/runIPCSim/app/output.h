#pragma once

#include "EigenSupport.h"
#include "setup/setup.h"
#include "triMeshGeo.h"

#include <filesystem>

namespace pgo::RunIPCSim
{
struct OutputDirectories
{
  std::filesystem::path root;
  std::filesystem::path states;
  std::filesystem::path surface;
  std::filesystem::path stress;
};

class RunIPCSimOutput
{
public:
  explicit RunIPCSimOutput(std::filesystem::path outputFolder);

  const OutputDirectories &directories() const { return outputDirs_; }
  std::filesystem::path logPath() const;
  std::filesystem::path statePath(int frame) const;
  std::filesystem::path surfacePath(int outputFrame) const;
  std::filesystem::path stressPath(int frame) const;

  void prepare(bool restartFromU) const;
  int loadLatestRestartState(int numSimSteps, int n3,
    EigenSupport::VXd &u, EigenSupport::VXd &uvel, EigenSupport::VXd &uacc) const;
  void writeState(int frame, const EigenSupport::VXd &u,
    const EigenSupport::VXd &uvel, const EigenSupport::VXd &uacc) const;
  void writeSurface(int outputFrame, const pgo::Mesh::TriMeshGeo &mesh) const;
  void writeStateAndSurfaceFrame(
    int frame,
    int outputFrame,
    const IpcSimulationContext &context,
    const EigenSupport::VXd &u,
    const EigenSupport::VXd &uvel,
    const EigenSupport::VXd &uacc,
    double scale,
    bool writeStateFile,
    bool writeSurfaceFile) const;
  void writeVonMisesStressJson(int frame, double timestep,
    const IpcSimulationContext &context, const EigenSupport::VXd &displacement) const;

private:
  OutputDirectories outputDirs_;
};

std::filesystem::path framePath(const std::filesystem::path &dir, const char *prefix, int frame, const char *extension);
}  // namespace pgo::RunIPCSim
