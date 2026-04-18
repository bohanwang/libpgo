#pragma once

namespace pgo
{
namespace SolidDeformationModel
{
class ElasticModel
{
public:
  ElasticModel() {}

  virtual ~ElasticModel() {}

  virtual int getNumParameters() const { return 0; }

  virtual void enableSPD(int enable) {}
};
}  // namespace SolidDeformationModel
}  // namespace pgo