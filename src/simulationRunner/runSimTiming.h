#pragma once

#include "simulationRunner.h"

#include <chrono>
#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace pgo::RunSim
{
// Optional wall timings of a frame invocation. Compute includes runner/logging
// and every attempted solve, but excludes the explicitly measured export spans.
class IPCSimulationTiming
{
public:
  explicit IPCSimulationTiming(const SimulationRunner::IPCSimulationRunOptions &options)
  {
    if (options.maxFrames != -1 && options.maxFrames <= 0)
      throw std::invalid_argument("`maxFrames` must be -1 or a positive integer.");
    if (!options.timingOutput.empty()) {
      if (!options.timingOutput.parent_path().empty())
        std::filesystem::create_directories(options.timingOutput.parent_path());
      stream.open(options.timingOutput);
      if (!stream)
        throw std::runtime_error("Unable to open IPC frame timing CSV.");
      stream << "frame,compute_seconds,output_seconds,total_seconds,terminal_status,all_solver_status_converged,accepted_solvers_converged,committed,newton_solves\n"
             << std::setprecision(17);
    }
  }
  bool enabled() const { return stream.is_open(); }

private:
  friend class IPCSimulationFrameTimer;
  std::ofstream stream;
};

class IPCSimulationFrameTimer
{
public:
  using Clock = std::chrono::steady_clock;
  IPCSimulationFrameTimer(IPCSimulationTiming &timing, int frame, double precedingComputeSeconds = 0):
    timing(timing), frame(frame), enabled(timing.stream.is_open()), start(enabled ? Clock::now() : Clock::time_point{}), precedingComputeSeconds(precedingComputeSeconds) {}

  IPCSimulationFrameTimer(const IPCSimulationFrameTimer &) = delete;
  IPCSimulationFrameTimer &operator=(const IPCSimulationFrameTimer &) = delete;

  ~IPCSimulationFrameTimer()
  {
    // Preserve the timing of a failed invocation without replacing its original
    // exception or fatal solver status with an I/O exception from a destructor.
    if (!finished) {
      try {
        finish(false, false);
      }
      catch (...) {
      }
    }
  }

  void recordSolve(int status)
  {
    terminalStatus = status;
    allConverged = allConverged && status == 0;
    ++solves;
  }

  template<class Callback>
  void output(Callback &&callback)
  {
    if (!enabled) {
      callback();
      return;
    }
    const auto begin = Clock::now();
    try {
      callback();
    }
    catch (...) {
      outputSeconds += seconds(begin, Clock::now());
      throw;
    }
    outputSeconds += seconds(begin, Clock::now());
  }

  void finish(bool committed, bool acceptedConverged)
  {
    if (finished)
      return;
    finished = true;
    if (!enabled)
      return;
    const double total = precedingComputeSeconds + seconds(start, Clock::now());
    timing.stream << frame << ',' << total - outputSeconds << ',' << outputSeconds << ',' << total
                  << ',' << terminalStatus << ',' << (solves > 0 && allConverged)
                  << ',' << acceptedConverged << ',' << committed << ',' << solves << std::endl;
    if (!timing.stream)
      throw std::runtime_error("Failed to write IPC frame timing CSV.");
  }

  bool allSolvesConverged() const { return solves > 0 && allConverged; }

private:
  static double seconds(Clock::time_point begin, Clock::time_point end)
  {
    return std::chrono::duration<double>(end - begin).count();
  }
  IPCSimulationTiming &timing;
  int frame;
  bool enabled;
  Clock::time_point start;
  double precedingComputeSeconds;
  double outputSeconds = 0;
  int terminalStatus = 2, solves = 0;
  bool allConverged = true, finished = false;
};
}  // namespace pgo::RunSim
