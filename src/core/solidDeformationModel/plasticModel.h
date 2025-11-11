#pragma once

namespace pgo
{
namespace SolidDeformationModel
{
class PlasticModel
{
public:
  PlasticModel() {}
  virtual ~PlasticModel() {}

  virtual int getNumParameters() const { return 0; }
};
}  // namespace SolidDeformationModel
}  // namespace pgo