#pragma once

#include "potentialEnergy.h"

#if defined(PGO_HAS_MKL)
#include "EigenMKLPardisoSupport.h"
#endif

#include <cfloat>
#include <memory>

namespace pgo {
namespace NonlinearOptimization {

class LineSearchHandle;

class NewtonRaphsonSolver {
public:
    enum SolverSubiterationType { SST_SUBITERATION_LINE_SEARCH, SST_SUBITERATION_STATIC_DAMPING, SST_SUBITERATION_ONE };

    enum SolveStatus {
        SR_CONVERGED = 0,
        SR_STALLED_SMALL_STEP = 1,
        SR_MAX_ITER_REACHED = 2,
        SR_FEASIBLE_STEP_ZERO = 3,
        SR_STOP_AFTER_INCREASE = 4,

        SR_LINE_SEARCH_FAILED = -1,
        SR_NUMERICAL_FAILURE = -2,
    };

    enum LineSearchMethod {
        LSM_GOLDEN,
        LSM_BRENTS,
        LSM_BACKTRACK,
        LSM_SIMPLE,
    };

    struct SolverParam {
        double                 alpha             = 0.5;
        SolverSubiterationType sst               = SST_SUBITERATION_LINE_SEARCH;
        LineSearchMethod       lsm               = LSM_SIMPLE;
        int                    stopAfterIncrease = 1;
        int                    addDamping        = 0;
    };

    NewtonRaphsonSolver(const double* x, SolverParam sp, PotentialEnergy_const_p energy_,
                        const std::vector<int>& fixedDOFs, const double* fixedValues_ = nullptr);

    void setFixedDOFs(const std::vector<int>& fixedDOFs, const double* fixedValues);
    int  solve(double* x, int numIter, double epsilon, int verbose);

    using AlphaTestFunc = std::function<double(const EigenSupport::VXd&, const EigenSupport::VXd&)>;
    void setAlphaTestFunc(AlphaTestFunc func) { alphaTestFunc = func; }
    void clearAlphaTestFunc() { alphaTestFunc = nullptr; }

    using StepFunc = std::function<void(const EigenSupport::VXd&, int)>;
    void setStepFunc(StepFunc func) { stepFunc = func; }

    const EigenSupport::VXd& getx() const { return x; }
    static bool              isHardFailure(int status) { return status < 0; }
    static bool              isSoftTermination(int status) { return status > 0; }
    static const char*       statusToString(int status);

protected:
    void filterVector(EigenSupport::VXd& v);

    PotentialEnergy_const_p           energy;
    SolverParam                       solverParam;
    std::shared_ptr<LineSearchHandle> lineSearchHandle;

    EigenSupport::VXd    x, grad, deltax, deltaxSmall, lineSearchx;
    EigenSupport::SpMatD sysFull, A11, A12;
    EigenSupport::SpMatI A11Mapping, A12Mapping;

#if defined(PGO_HAS_MKL)
    std::shared_ptr<EigenSupport::EigenMKLPardisoSupport> solver;
#else
    std::shared_ptr<EigenSupport::SymSolver> solver;
#endif

    std::vector<int>  allDOFs, fixedDOFs;
    std::vector<int>  rhss2b, rhsb2s;
    EigenSupport::VXd rhs;
    EigenSupport::VXd fixedValues;
    int               n3;

    EigenSupport::VXd historyx;
    double            historyGradNormMin;

    AlphaTestFunc alphaTestFunc;
    StepFunc      stepFunc;
};
}  // namespace NonlinearOptimization
}  // namespace pgo
