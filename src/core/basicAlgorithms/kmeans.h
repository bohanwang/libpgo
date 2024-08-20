#pragma once

#include <random>

namespace pgo
{
namespace BasicAlgorithms
{

template<int dim, typename EleType>
class KMeans
{
protected:
  template<typename PointType>
  static PointType add(const PointType &p1, const PointType &p2)
  {
    PointType pret;
    for (int i = 0; i < dim; i++)
      pret[i] = p1[i] + p2[i];

    return pret;
  }

  template<typename PointType>
  static PointType zero()
  {
    PointType zeroPt;
    for (int i = 0; i < dim; i++)
      zeroPt[i] = static_cast<EleType>(0.0);

    return zeroPt;
  }

  template<typename PointType>
  static void div(PointType &p, EleType val)
  {
    for (int i = 0; i < dim; i++)
      p[i] /= val;
  }

public:
  template<typename PointType, typename Index, typename Index1, typename DistFunc>
  static void compute(Index numPoints, const PointType *points, Index1 numCenters, DistFunc distFunc, int *pointCenterIDs)
  {
    if (numCenters < 1)
      return;

    if (numCenters == 1) {
      for (Index i = 0; i < numPoints; i++)
        pointCenterIDs[i] = 0;

      return;
    }

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> distrib(0, static_cast<int>(numCenters) - 1);
    PointType *centers = new PointType[numCenters];
    EleType *counters = new EleType[numCenters];

    for (Index i = 0; i < numPoints; i++)
      pointCenterIDs[i] = distrib(gen);

    int *pointCenterIDs1 = new int[numPoints];
    for (Index i = 0; i < numPoints; i++) {
      pointCenterIDs1[i] = -1;
    }

    while (1) {
      // check whether it terminates
      bool stop = true;
      for (Index i = 0; i < numPoints; i++) {
        if (pointCenterIDs1[i] != pointCenterIDs[i]) {
          stop = false;
          break;
        }
      }
      if (stop)
        break;

      for (Index i = 0; i < numPoints; i++) {
        pointCenterIDs1[i] = pointCenterIDs[i];
      }

      // compute centers
      for (Index1 i = 0; i < numCenters; i++) {
        centers[i] = zero<PointType>();
        counters[i] = static_cast<EleType>(0.0);
      }

      for (Index i = 0; i < numPoints; i++) {
        centers[pointCenterIDs[i]] = add<PointType>(centers[pointCenterIDs[i]], points[i]);
        counters[pointCenterIDs[i]] += static_cast<EleType>(1.0);
      }

      for (Index1 i = 0; i < numCenters; i++) {
        div<PointType>(centers[i], counters[i]);
      }

      // assign new centers
      for (Index pi = 0; pi < numPoints; pi++) {
        EleType minDist = 1e100;
        int centerID = 0;
        for (Index1 ci = 0; ci < numCenters; ci++) {
          EleType dist = distFunc(points[pi], centers[ci]);
          if (dist < minDist) {
            centerID = static_cast<int>(ci);
            minDist = dist;
          }
        }
        pointCenterIDs[pi] = centerID;
      }
    }

    delete[] pointCenterIDs1;
    delete[] counters;
    delete[] centers;
  }
};
}  // namespace BasicAlgorithms
}  // namespace pgo