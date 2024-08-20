#pragma once

namespace pgo
{
namespace Simulation
{
enum class TimeIntegratorSolverOption
{
  SO_IPOPT,
  SO_NEWTON,
  SO_KNITRO
};

}
}