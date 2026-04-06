#pragma once

#include "contactEnergyUtilities.h"

#include "potentialEnergy.h"

#include <array>
#include <vector>

namespace pgo {
namespace Contact {
struct PointTrianglePairBarrierEnergyBuffer;

class PointTrianglePairBarrierEnergy : public NonlinearOptimization::PotentialEnergy {
public:
    static double normalizePositiveZero(double dhat, double positiveZero);

    PointTrianglePairBarrierEnergy(int numPairs, int numObjects, const std::array<int, 4>* objectIDs,
                                   const std::array<int, 4>* pointTrianglePairs, const int* objectDOFOffsets,
                                   const EigenSupport::SpMatD* sampleEmbeddingWeights,
                                   const std::vector<double>* sampleWeights, const double* constraintNormals,
                                   const double* constraintBarycentricWeights, double dhat, double positiveZero);

    PointTrianglePairBarrierEnergyBuffer* allocateBuffer() const;
    void                                  freeBuffer(PointTrianglePairBarrierEnergyBuffer* buf) const;
    void                                  setBuffer(PointTrianglePairBarrierEnergyBuffer* b) { buf = b; }

    void setToPosFunction(PosFunction f) { toPosFunc = f; }
    void setCoeff(double c) { coeffAll = c; }

    double func(EigenSupport::ConstRefVecXd x) const override;
    void   gradient(EigenSupport::ConstRefVecXd x, EigenSupport::RefVecXd grad) const override;
    void   hessian(EigenSupport::ConstRefVecXd x, EigenSupport::SpMatD& hess) const override;

    void createHessian(EigenSupport::SpMatD& hess) const override { hess = hessTemplate; }
    void getDOFs(std::vector<int>& dofs) const override { dofs = allDOFs; }
    int  getNumDOFs() const override { return static_cast<int>(hessTemplate.rows()); }

protected:
    EigenSupport::V3d computePosition(EigenSupport::ConstRefVecXd x, int objID, int sampleID) const;
    double            getPairWeight(int pairID) const;

    int                               numPairs;
    int                               numObjects;
    const std::array<int, 4>* const   objectIDs;
    const std::array<int, 4>* const   pointTrianglePairs;
    const int* const                  objectDOFOffsets;
    const EigenSupport::SpMatD* const sampleEmbeddingWeights;
    const std::vector<double>* const  sampleWeights;
    const double* const               constraintNormals;
    const double* const               constraintBarycentricWeights;

    double      coeffAll     = 1.0;
    double      dhat         = 0.0;
    double      positiveZero = 0.0;
    PosFunction toPosFunc    = [](const EigenSupport::V3d& x, EigenSupport::V3d& p, int) { p = x; };

    std::vector<int>       allDOFs;
    EigenSupport::SpMatD   hessTemplate;
    EigenSupport::EntryMap entryMap;

    PointTrianglePairBarrierEnergyBuffer* buf = nullptr;
};

}  // namespace Contact
}  // namespace pgo
