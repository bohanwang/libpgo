#pragma once

#include "contactEnergyUtilities.h"

#include "potentialEnergy.h"

namespace pgo {
namespace Contact {
struct PointPenetrationBarrierEnergyBuffer;

class PointPenetrationBarrierEnergy : public NonlinearOptimization::PotentialEnergy {
public:
    static double normalizePositiveZero(double dhat, double positiveZero);

    PointPenetrationBarrierEnergy(int numPoints, int nAll, const double* constraintCoeffs,
                                  const double* constraintTargetPositions, const double* constraintNormals,
                                  const std::vector<std::vector<int>>&    barycentricIdx,
                                  const std::vector<std::vector<double>>& barycentricWeights, double dhat,
                                  double positiveZero);

    PointPenetrationBarrierEnergyBuffer* allocateBuffer() const;
    void                                 freeBuffer(PointPenetrationBarrierEnergyBuffer* buf) const;
    void                                 setBuffer(PointPenetrationBarrierEnergyBuffer* b) { buf = b; }

    void setComputePosFunction(PosFunction f) { posFunc = f; }

    double func(EigenSupport::ConstRefVecXd u) const override;
    void   gradient(EigenSupport::ConstRefVecXd u, EigenSupport::RefVecXd grad) const override;
    void   hessian(EigenSupport::ConstRefVecXd u, EigenSupport::SpMatD& hess) const override;

    void setCoeff(double v) { coeffAll = v; }

    void createHessian(EigenSupport::SpMatD& h) const override { h = hessianTemplate; }
    void getDOFs(std::vector<int>& dofs) const override { dofs = allDOFs; }
    int  getNumDOFs() const override { return static_cast<int>(hessianTemplate.rows()); }

protected:
    int                                     nAll;
    int                                     numPoints;
    const double*                           constraintCoeffs;
    const double*                           constraintTargetPositions;
    const double*                           constraintNormals;
    const std::vector<std::vector<int>>&    barycentricIdx;
    const std::vector<std::vector<double>>& barycentricWeights;

    double      coeffAll = 1.0;
    double      dhat;
    double      positiveZero;
    PosFunction posFunc;

    std::vector<int>       allDOFs;
    EigenSupport::SpMatD   hessianTemplate;
    EigenSupport::EntryMap hessianEntryMap;

    PointPenetrationBarrierEnergyBuffer* buf = nullptr;
};

}  // namespace Contact
}  // namespace pgo
