#include "ipc/cuda/cudaIPCContinuousCollision.h"
#include "cudaDeviceBuffer.h"

namespace
{
__global__ void intervalCertificate(const double *positions, const double *motion,
  const double *error, int kind, double firstTime, int *result)
{
  using namespace pgo::Contact::CIPC::CudaIPCContinuousCollision;
  Vec3 p[4], v[4], e[4];
  for (int i = 0; i < 4; ++i)
    for (int k = 0; k < 3; ++k) {
      p[i][k] = positions[3 * i + k];
      v[i][k] = motion[3 * i + k];
      e[i][k] = error[3 * i + k];
    }
  *result = certifyRemainingInterval(p, v, e, kind, firstTime);
}
}  // namespace

bool cudaIPCIntervalCertificateForTest(const double *positions, const double *motion,
  const double *error, int kind, double firstTime)
{
  using pgo::CudaUtilities::DeviceBuffer;
  DeviceBuffer<double> p, v, e;
  p.upload(positions, 12);
  v.upload(motion, 12);
  e.upload(error, 12);
  DeviceBuffer<int> result(1);
  intervalCertificate<<<1, 1>>>(p.data(), v.data(), e.data(), kind, firstTime, result.data());
  pgo::CudaUtilities::checkCuda(cudaGetLastError(), "launch interval CCD test");
  int certified = 0;
  result.download(&certified, 1);
  return certified != 0;
}
