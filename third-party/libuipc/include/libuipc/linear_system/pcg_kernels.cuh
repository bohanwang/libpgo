// SPDX-License-Identifier: Apache-2.0
// Adapted from libuipc src/backends/cuda/linear_system/linear_pcg.cu,
// commit 9c748a759bb9952a7aeedd59076a70591dc61d6c.
// Local modifications: double/raw device pointers, PGO namespace, explicit
// bounds, and removal of upstream SimSystem, cuda_tool views and registration.
// See PROVENANCE.txt and the libuipc LICENSE/NOTICE at the import root.
#pragma once

#include <cuda_runtime.h>

namespace pgo::NonlinearOptimization::detail
{

__global__ void pcgUpdateXR(double alpha, double *x, const double *p,
  double *r, const double *Ap, int n)
{
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= n)
    return;
  x[i] += alpha * p[i];
  r[i] -= alpha * Ap[i];
}

__global__ void pcgUpdateP(double *p, const double *z, double beta, int n)
{
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= n)
    return;
  p[i] = z[i] + beta * p[i];
}

}  // namespace pgo::NonlinearOptimization::detail
