/*
copyright to Bohan Wang
*/

#include "ipcBarrier.h"

#include <cmath>

namespace pgo {
namespace Contact {
namespace IPC {
// =========================================================================
//  Barrier function
// =========================================================================
namespace barrier {

// Barrier on squared distance (matching Codim-IPC elastic formulation):
//   b(s, shat) = -(s/shat - 1)^2 * ln(s/shat)    for 0 < s < shat
// where s = d^2 (squared distance), shat = dhat^2
double b(double s, double shat)
{
  if (s >= shat || s <= 0.0)
    return 0.0;
  double r = s / shat;
  double rm1 = r - 1.0;
  return -rm1 * rm1 * std::log(r);
}

// db/ds = -(1/shat) * (r - 1) * (2 ln(r) + (r-1)/r)
double dbds(double s, double shat)
{
  if (s >= shat || s <= 0.0)
    return 0.0;
  double r = s / shat;
  double rm1 = r - 1.0;
  return -(rm1 * (2.0 * std::log(r) + rm1 / r)) / shat;
}

// d2b/ds2 = -(1/shat^2) * (2 ln(r) + 4(r-1)/r - (r-1)^2/r^2)
double d2bds2(double s, double shat)
{
  if (s >= shat || s <= 0.0)
    return 0.0;
  double r = s / shat;
  double rm1 = r - 1.0;
  double invr = 1.0 / r;
  double shat2 = shat * shat;
  return -(2.0 * std::log(r) + 4.0 * rm1 * invr - rm1 * rm1 * invr * invr) / shat2;
}

}  // namespace barrier

}  // namespace IPC
}  // namespace Contact
}  // namespace pgo
