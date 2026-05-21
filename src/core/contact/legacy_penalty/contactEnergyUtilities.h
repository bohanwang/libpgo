#pragma once

#include "EigenDef.h"

#include <functional>

namespace pgo
{
namespace Contact
{

using PosFunction = std::function<void(const EigenSupport::V3d &, EigenSupport::V3d &, int dofStart)>;

}  // namespace Contact
}  // namespace pgo