/*
copyright to Bohan Wang
*/

#pragma once

#include "EigenDef.h"

namespace pgo
{
namespace Contact
{
namespace IPC
{

using namespace pgo::EigenSupport;

// Point-Triangle: which feature of the triangle is closest to the point?
enum class PTDistType
{
  PP_PT0,
  PP_PT1,
  PP_PT2,
  PE_PT0T1,
  PE_PT1T2,
  PE_PT2T0,
  PT
};

// Edge-Edge: which sub-features are closest?
enum class EEDistType
{
  PP_Ea0Eb0,
  PP_Ea0Eb1,
  PP_Ea1Eb0,
  PP_Ea1Eb1,
  PE_Ea0_Eb,
  PE_Ea1_Eb,
  PE_Eb0_Ea,
  PE_Eb1_Ea,
  EE
};

namespace distance
{

PTDistType classifyPT(const V3d &p, const V3d &t0,
  const V3d &t1, const V3d &t2);

EEDistType classifyEE(const V3d &ea0, const V3d &ea1,
  const V3d &eb0, const V3d &eb1);

double ppSqDist(const V3d &a, const V3d &b);
V6d ppSqDistGrad(const V3d &a, const V3d &b);
M6d ppSqDistHess(const V3d &a, const V3d &b);

double peSqDist(const V3d &p, const V3d &e0, const V3d &e1);
V9d peSqDistGrad(const V3d &p, const V3d &e0, const V3d &e1);
M9d peSqDistHess(const V3d &p, const V3d &e0, const V3d &e1);

double ptSqDist(const V3d &p, const V3d &t0,
  const V3d &t1, const V3d &t2);
V12d ptSqDistGrad(const V3d &p, const V3d &t0,
  const V3d &t1, const V3d &t2);
M12d ptSqDistHess(const V3d &p, const V3d &t0,
  const V3d &t1, const V3d &t2);

double eeSqDist(const V3d &ea0, const V3d &ea1,
  const V3d &eb0, const V3d &eb1);
V12d eeSqDistGrad(const V3d &ea0, const V3d &ea1,
  const V3d &eb0, const V3d &eb1);
M12d eeSqDistHess(const V3d &ea0, const V3d &ea1,
  const V3d &eb0, const V3d &eb1);

double eeMollifier(const V3d &ea0, const V3d &ea1,
  const V3d &eb0, const V3d &eb1,
  double eps_x);
V12d eeMollifierGrad(const V3d &ea0, const V3d &ea1,
  const V3d &eb0, const V3d &eb1,
  double eps_x);
M12d eeMollifierHess(const V3d &ea0, const V3d &ea1,
  const V3d &eb0, const V3d &eb1,
  double eps_x);

double computePTSqDist(const V3d &p, const V3d &t0,
  const V3d &t1, const V3d &t2);
V12d computePTSqDistGrad(const V3d &p, const V3d &t0,
  const V3d &t1, const V3d &t2);
M12d computePTSqDistHess(const V3d &p, const V3d &t0,
  const V3d &t1, const V3d &t2);

double computeEESqDist(const V3d &ea0, const V3d &ea1,
  const V3d &eb0, const V3d &eb1);
V12d computeEESqDistGrad(const V3d &ea0, const V3d &ea1,
  const V3d &eb0, const V3d &eb1);
M12d computeEESqDistHess(const V3d &ea0, const V3d &ea1,
  const V3d &eb0, const V3d &eb1);

struct EEDistAll
{
  double d2;
  V12d grad;
  M12d hess;
};

EEDistAll computeEESqDistAll(const V3d &ea0, const V3d &ea1,
  const V3d &eb0, const V3d &eb1);

struct PTDistAll
{
  double d2;
  V12d grad;
  M12d hess;
};

PTDistAll computePTSqDistAll(const V3d &p, const V3d &t0,
  const V3d &t1, const V3d &t2);

}  // namespace distance

}  // namespace IPC
}  // namespace Contact
}  // namespace pgo
