/*
author: Bohan Wang
copyright to USC,MIT
*/

#pragma once

namespace pgo
{
namespace SolidDeformationModel
{
class PlasticModel
{
public:
  PlasticModel(int np):
    numParams(np) {}
  virtual ~PlasticModel() {}

  virtual int getNumParameters() const { return 0; }
  virtual void computeA(const double *param, double A[9]) const;
  virtual void computeAInv(const double *param, double AInv[9]) const;
  virtual double compute_detA(const double *param) const { return 1.0; }
  // here ddetA_da (and d2detA_da2)'s dimension is matched with param.
  // For example, if param is \in R^3, then ddetA_da is \in R^3 as well.
  // And d2detA_da2 is \in R^3 x R^3
  virtual void compute_ddetA_da(const double *param, double *ddetaA_da, int leadingDim) const;
  virtual void compute_d2detA_da2(const double *param, double *d2detA_da2, int leadingDim) const;
  virtual void compute_dAInv_da(const double *param, int pi, double ret[9]) const;
  virtual void compute_d2AInv_da2(const double *param, int pi, int pj, double ret[9]) const;

  virtual void defaultFp(double *Fp) const;
  virtual void projectParam(double *param, double zeroThreshold) const;
  virtual void toParam(const double *Fp, double *param) const;

  virtual void computeR(const double *param, double R[9]) const;
  virtual void compute_dparamfull_dparamsub(const double *param, const double *basis, int numHandles, double *dpf_dps) const;
  virtual bool has2OrderDeriv() const { return false; }
  virtual void compute_d2paramfull_dparamsub2(const double *param, const double *basis, int numHandles, int pi, double *d2a_dz2) const;

  virtual void compute_paramfull(double *param) const;

protected:
  int numParams;
};

inline void PlasticModel::computeA(const double *, double A[9]) const
{
  A[0] = 1.0, A[3] = 0.0, A[6] = 0.0;
  A[1] = 0.0, A[4] = 1.0, A[7] = 0.0;
  A[2] = 0.0, A[5] = 0.0, A[8] = 1.0;
}

inline void PlasticModel::computeAInv(const double *, double AInv[9]) const
{
  AInv[0] = 1.0, AInv[3] = 0.0, AInv[6] = 0.0;
  AInv[1] = 0.0, AInv[4] = 1.0, AInv[7] = 0.0;
  AInv[2] = 0.0, AInv[5] = 0.0, AInv[8] = 1.0;
}

inline void PlasticModel::compute_ddetA_da(const double *, double *, int) const
{
}

inline void PlasticModel::compute_d2detA_da2(const double *, double *, int) const
{
}

inline void PlasticModel::compute_dAInv_da(const double *, int, double[9]) const
{
}

inline void PlasticModel::compute_d2AInv_da2(const double *, int, int, double[9]) const
{
}

inline void PlasticModel::defaultFp(double *) const
{
}

inline void PlasticModel::projectParam(double *, double) const
{
}

inline void PlasticModel::toParam(const double *, double *) const
{
}

inline void PlasticModel::computeR(const double *, double R[9]) const
{
  R[0] = 1.0, R[1] = 0.0, R[2] = 0.0;
  R[3] = 0.0, R[4] = 1.0, R[5] = 0.0;
  R[6] = 0.0, R[7] = 0.0, R[8] = 1.0;
}

inline void PlasticModel::compute_dparamfull_dparamsub(const double *, const double *, int, double *) const
{
}

inline void PlasticModel::compute_d2paramfull_dparamsub2(const double *, const double *, int, int, double *) const
{
}

inline void PlasticModel::compute_paramfull(double *) const
{
}
}  // namespace SolidDeformationModel
}  // namespace pgo