#pragma once

#include "averagingBuffer.h"

#include <chrono>
#include <string>
#include <unordered_map>
#include <mutex>

namespace pgo
{
namespace BasicAlgorithms
{

class TimeTable
{
public:
  // Start timing for a named statistic
  void start(const std::string &name)
  {
    std::lock_guard<std::mutex> lock(mutex_);
    startTimes[name] = std::chrono::high_resolution_clock::now();
  }

  // End timing for a named statistic and record the elapsed time
  void end(const std::string &name)
  {
    auto endTime = std::chrono::high_resolution_clock::now();
    std::lock_guard<std::mutex> lock(mutex_);
    if (startTimes.find(name) != startTimes.end()) {
      auto elapsed = (double)std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTimes[name]).count() * 1e-6;
      auto it = statistics.find(name);
      if (it == statistics.end()) {
        auto tp = std::make_tuple((int)statistics.size(), std::make_shared<AveragingBuffer>(1024));
        it = statistics.emplace_hint(it, name, tp);
      }

      std::get<1>(it->second)->addValue(elapsed);
      startTimes.erase(name);
    }
  }

  // Get the average time for a named statistic
  double getAverage(const std::string &name) const
  {
    std::lock_guard<std::mutex> lock(mutex_);
    if (statistics.find(name) != statistics.end()) {
      return std::get<1>(statistics.at(name))->getAverage();
    }

    return 0.0;
  }

  // Print all recorded time statistics and their total sum
  void printStatistics() const
  {
    std::lock_guard<std::mutex> lock(mutex_);
    double totalSum = 0.0;

    std::cout << "Time Statistics:" << std::endl;

    std::vector<std::tuple<std::string, int, std::shared_ptr<AveragingBuffer>>> temp;
    for (const auto &entry : statistics) {
      temp.emplace_back(std::make_tuple(entry.first, std::get<0>(entry.second), std::get<1>(entry.second)));
    }

    std::sort(temp.begin(), temp.end(), [](const auto &v1, const auto &v2) {
      return std::get<1>(v1) < std::get<1>(v2);
    });

    for (const auto &entry : temp) {
      double average = std::get<2>(entry)->getAverage();
      std::cout << "  " << std::get<0>(entry) << ": " << average << " seconds" << std::endl;
      totalSum += average;
    }

    std::cout << "Total Sum: " << totalSum << " seconds" << std::endl;
  }

private:
  mutable std::mutex mutex_;                                                                      // Mutex for thread safety
  std::unordered_map<std::string, std::chrono::high_resolution_clock::time_point> startTimes;     // Start times
  std::unordered_map<std::string, std::tuple<int, std::shared_ptr<AveragingBuffer>>> statistics;  // Named statistics
};
}  // namespace BasicAlgorithms
}  // namespace pgo