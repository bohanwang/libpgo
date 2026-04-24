/*
copyright to Bohan Wang
*/

#pragma once

#include "../topology/surfaceIPCTopology.h"

namespace pgo
{
namespace Contact
{
namespace CIPC
{

class SurfaceIPCMaxStep
{
public:
  double compute(
    const SurfaceIPCTopology &topology,
    EigenSupport::ConstRefVecXd x,
    EigenSupport::ConstRefVecXd dx,
    double dhat,
    double slackness) const;
};

}  // namespace CIPC
}  // namespace Contact
}  // namespace pgo
