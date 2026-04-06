#include "NewtonRaphsonSolver.h"

#include "EigenSupport.h"
#include "lineSearch.h"
#include "pgoLogging.h"

#include <iostream>
#include <numeric>
#include <chrono>
#include <cmath>

using namespace pgo;
using namespace pgo::NonlinearOptimization;
namespace ES = pgo::EigenSupport;
using hclock = std::chrono::high_resolution_clock;

namespace pgo::NonlinearOptimization {
class LineSearchHandle {
public:
    std::shared_ptr<LineSearch> nativeLineSearch;
};
}  // namespace pgo::NonlinearOptimization

inline double dura(const hclock::time_point& t1, const hclock::time_point& t2) {
    return std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1e6;
}

namespace {

constexpr double kSmallAcceptedStepSizeTol = 1e-15;
constexpr double kTinyDirectionTol         = 1e-12;

bool shouldTreatLineSearchFailureAsStall(double alpha, double deltaxMax, double gradNorm, double epsilon, double eng,
                                         double eng1) {
    const double trialStepSize = std::abs(alpha) * deltaxMax;
    const double energyTol     = std::max(1e-12, std::abs(eng) * 1e-12);
    const double gradTol       = std::max(1e-12, epsilon * 10.0);

    return trialStepSize < kSmallAcceptedStepSizeTol || deltaxMax < kTinyDirectionTol || gradNorm <= gradTol ||
           std::abs(eng1 - eng) <= energyTol;
}

}  // namespace

NewtonRaphsonSolver::NewtonRaphsonSolver(const double* x_, SolverParam sp, PotentialEnergy_const_p energy_,
                                         const std::vector<int>& fixedDOFs_, const double* fixedValues_)
    : energy(energy_), solverParam(sp) {
    n3 = (int)energy->getNumDOFs();
    allDOFs.resize(energy->getNumDOFs());
    std::iota(allDOFs.begin(), allDOFs.end(), 0);

    grad.resize(energy->getNumDOFs());
    x = Eigen::Map<const ES::VXd>(x_, energy->getNumDOFs());
    deltax.resize(energy->getNumDOFs());
    lineSearchx.resize(energy->getNumDOFs());
    historyx.resize(energy->getNumDOFs());
    historyGradNormMin = 1e100;

    setFixedDOFs(fixedDOFs_, fixedValues_);

    if (solverParam.sst == SST_SUBITERATION_LINE_SEARCH) {
        lineSearchHandle = std::make_shared<LineSearchHandle>();

        LineSearch::EvaluateFunction evalFunc = [this](const double* x, double* f, double* grad) -> int {
            if (f)
                *f = energy->func(Eigen::Map<const ES::VXd>(x, n3));

            if (grad) {
                memset(grad, 0, sizeof(double) * n3);
                energy->gradient(Eigen::Map<const ES::VXd>(x, n3), Eigen::Map<ES::VXd>(grad, n3));
            }

            return 0;
        };

        lineSearchHandle->nativeLineSearch = std::make_shared<LineSearch>(n3, evalFunc);
    }
}

void NewtonRaphsonSolver::setFixedDOFs(const std::vector<int>& fixedDOFs_, const double* fixedValues_) {
    if (fixedDOFs_.size() != 0 && fixedDOFs.size() == fixedDOFs_.size() &&
        std::memcmp(fixedDOFs.data(), fixedDOFs_.data(), sizeof(int) * fixedDOFs.size()) == 0) {
    } else {
        fixedDOFs = fixedDOFs_;

        // dofs
        rhsb2s.clear();
        rhss2b.clear();
        ES::removeRows(n3, fixedDOFs, rhsb2s, rhss2b);

        rhs.resize(n3 - (int)fixedDOFs.size());
        deltaxSmall.resize(n3 - (int)fixedDOFs.size());

        // sparse matrix
        energy->createHessian(sysFull);
        energy->hessian(x, sysFull);

        ES::removeRowsCols(sysFull, fixedDOFs, A11);
        ES::removeRowsCols(sysFull, A11, fixedDOFs, A11Mapping);

        // std::vector<int> flexibleDOFs;
        // std::set_difference(allDOFs.begin(), allDOFs.end(), fixedDOFs.begin(), fixedDOFs.end(),
        // std::back_inserter(flexibleDOFs)); ES::SelectRowsCols(sysFull, flexibleDOFs, fixedDOFs, A12);
        // ES::Big2Small(sysFull, A12, flexibleDOFs, fixedDOFs, A12Mapping, 0);

#if defined(PGO_HAS_MKL)
        solver = std::make_shared<ES::EigenMKLPardisoSupport>(
            A11, ES::EigenMKLPardisoSupport::MatrixType::REAL_SYM_INDEFINITE,
            ES::EigenMKLPardisoSupport::ReorderingType::NESTED_DISSECTION, 0, 0, 2, 0, 0, 0);
        solver->analyze(A11);
#else
        solver = std::make_shared<EigenSupport::SymSolver>();
        solver->analyzePattern(A11);
#endif
    }

    if (fixedDOFs_.empty()) {
        fixedValues.resize(0);
    } else {
        PGO_ALOG(fixedValues_ != nullptr);
        fixedValues = ES::Mp<const ES::VXd>(fixedValues_, fixedDOFs_.size());
    }
}

const char* NewtonRaphsonSolver::statusToString(int status) {
    switch (status) {
    case SR_CONVERGED:
        return "converged";
    case SR_STALLED_SMALL_STEP:
        return "stalled_small_step";
    case SR_MAX_ITER_REACHED:
        return "max_iter_reached";
    case SR_FEASIBLE_STEP_ZERO:
        return "feasible_step_zero";
    case SR_STOP_AFTER_INCREASE:
        return "stop_after_increase";
    case SR_LINE_SEARCH_FAILED:
        return "line_search_failed";
    case SR_NUMERICAL_FAILURE:
        return "numerical_failure";
    default:
        return "unknown_status";
    }
}

int NewtonRaphsonSolver::solve(double* x_, int numIter, double epsilon, int verbose) {
    hclock::time_point t1 = hclock::now();

    x.noalias() = Eigen::Map<ES::VXd>(x_, energy->getNumDOFs());

    if (solverParam.sst == SST_SUBITERATION_ONE || solverParam.sst == SST_SUBITERATION_STATIC_DAMPING) {
        historyGradNormMin = 1e100;
    }

    for (size_t i = 0; i < fixedDOFs.size(); i++) {
        x[fixedDOFs[i]] = fixedValues[i];
    }

    int printGap = 10;
    if (verbose == 2) {
        printGap = 10;
    } else if (verbose == 3) {
        printGap = 1;
    }

    int    iter        = 0;
    double lambdaScale = 1.0;
    double lambda0     = 1.0;
    int    status      = SR_MAX_ITER_REACHED;
    // compute lambda initial
    memset(grad.data(), 0, sizeof(double) * grad.size());
    energy->gradient(x, grad);
    filterVector(grad);
    lambda0 = grad.cwiseAbs().maxCoeff();

    double gradNormLast = lambda0;
    for (; iter < numIter; iter++) {
        if (verbose >= 2 && iter % printGap == 0)
            std::cout << "    Iter=" << iter << std::endl;

        double eng = energy->func(x);

        // we solve f(x_i) + K(x_i) deltax = 0
        memset(grad.data(), 0, sizeof(double) * grad.size());
        energy->gradient(x, grad);
        filterVector(grad);

        double gradNorm = grad.cwiseAbs().maxCoeff();
        if (verbose >= 2 && iter % printGap == 0)
            std::cout << "        E= " << eng << "; ||grad||_max=" << grad.cwiseAbs().maxCoeff()
                      << "; ||x||=" << x.norm() << "; ||grad||=" << gradNorm << std::endl;

        if (!std::isfinite(eng) || !std::isfinite(gradNorm)) {
            if (verbose >= 1) {
                std::cout << "    Iter=" << iter << "; encountered non-finite energy or gradient." << std::endl;
            }

            status = SR_NUMERICAL_FAILURE;
            break;
        }

        if (solverParam.sst == SST_SUBITERATION_ONE || solverParam.sst == SST_SUBITERATION_STATIC_DAMPING) {
            if (gradNorm < historyGradNormMin) {
                historyx.noalias() = x;
                historyGradNormMin = gradNorm;
            } else {
                if (solverParam.stopAfterIncrease) {
                    if (verbose >= 1) {
                        std::cout << "    Iter=" << iter << "; stopping after gradient increase." << std::endl;
                    }

                    status = SR_STOP_AFTER_INCREASE;
                    break;
                }
            }
        }

        if (gradNorm < epsilon) {
            if (verbose >= 1) {
                std::cout << "    Iter=" << iter << "; ||grad||=" << gradNorm << " < eps. Done." << std::endl;
            }

            status = SR_CONVERGED;
            break;
        }

        memset(sysFull.valuePtr(), 0, sizeof(double) * sysFull.nonZeros());
        energy->hessian(x, sysFull);

        // grad too small, we don't need damping
        if (gradNorm < 1e-4) {
            lambdaScale = 0.0;
        } else {
            // if lambda too small, we set it to zero
            if (lambdaScale < 1e-8) {
                lambdaScale = 0;
            } else {
                // if gradient increasing, we increase lambda
                if (gradNorm > gradNormLast) {
                } else {
                    // otherwise we decrease lambda
                    lambdaScale *= 0.9;
                }

                // clamp
                lambdaScale = std::min(lambdaScale, 1.0);
            }
        }
        gradNormLast = gradNorm;

        // std::cout << "        Damping lambda=" << lambda << std::endl;

        // remove column rows
        ES::transferBigToSmall(sysFull, A11, A11Mapping, 1);

        if (solverParam.addDamping) {
            for (int i = 0; i < A11.rows(); i++) {
                A11.coeffRef(i, i) += lambdaScale * lambda0;
            }
        }

        // ES::transferBigToSmall(sysFull, A12, A12Mapping, 1);
        ES::transferBigToSmall(grad, rhs, rhsb2s, 1);

        // std::cout << "@@:\n"
        //   << x.norm() << '\n'
        //   << sysFull.coeff(0, 0) << '\n'
        //   << sysFull.coeff(1, 0) << '\n'
        //   << sysFull.coeff(1, 1) << '\n'
        //   << sysFull.coeff(2, 0) << '\n'
        //   << sysFull.coeff(2, 1) << '\n'
        //   << sysFull.coeff(2, 2) << std::endl;

        // rhs = -grad - A12 * fixedvalue
        // ES::mv(A12, fixedValues, rhs, -1.0, -1.0, 0);
        rhs *= -1.0;
#if defined(PGO_HAS_MKL)
        solver->factorize(A11);
        solver->solve(A11, deltaxSmall.data(), rhs.data(), 1);
#else
        solver->factorize(A11);
        deltaxSmall.noalias() = solver->solve(rhs);
#endif

        if (verbose >= 3 && iter % printGap == 0)
            std::cout << (A11 * deltaxSmall - rhs).norm() << ' ' << rhs.norm() << std::endl;

        memset(deltax.data(), 0, sizeof(double) * n3);
        ES::transferSmallToBig(deltaxSmall, deltax, rhss2b);

        if (!deltax.allFinite()) {
            if (verbose >= 1) {
                std::cout << "    Iter=" << iter << "; Newton direction contains non-finite entries." << std::endl;
            }

            status = SR_NUMERICAL_FAILURE;
            break;
        }

        // if (iter == 0)
        //   error0 = deltax.norm();
        // else {
        //   if (deltax.norm() < epsilon) {
        //     if (verbose >= 1)
        //       std::cout << "    dir too small." << std::endl;

        //     break;
        //   }
        // }

        const double deltaxMax = deltax.cwiseAbs().maxCoeff();
        if (!std::isfinite(deltaxMax)) {
            if (verbose >= 1) {
                std::cout << "    Iter=" << iter << "; Newton direction norm is non-finite." << std::endl;
            }

            status = SR_NUMERICAL_FAILURE;
            break;
        }

        if (verbose >= 2 && iter % printGap == 0)
            std::cout << "        ||deltax||_max=" << deltaxMax << std::endl;

        // x += alpha delta x ?
        if (solverParam.sst == SST_SUBITERATION_LINE_SEARCH) {
            double alpha    = 1;
            double stepSize = 0;

            alpha = alphaTestFunc ? alphaTestFunc(x, deltax) : 1.0;

            if (alpha <= 0.0) {
                if (verbose >= 1) {
                    std::cout << "    Iter=" << iter << "; feasible step upper bound is zero." << std::endl;
                }

                status = SR_FEASIBLE_STEP_ZERO;
                break;
            }

            lineSearchx.noalias() = x + deltax * alpha;
            double eng1           = energy->func(lineSearchx);
            int    maxIter        = 50;
            if (eng1 < eng) {
                maxIter = 3;
            }

            if (solverParam.lsm == LSM_GOLDEN) {
                lineSearchHandle->nativeLineSearch->setMaxIterations(maxIter);
                LineSearch::Result ret = lineSearchHandle->nativeLineSearch->golden(x.data(), deltax.data(), eng);
                alpha                  = ret.alpha;
                eng1                   = ret.f;

                if (verbose >= 2 && iter % printGap == 0) {
                    std::cout << "        f=" << ret.f << ";alpha=" << alpha << std::endl;
                }
            } else if (solverParam.lsm == LSM_BRENTS) {
                lineSearchHandle->nativeLineSearch->setMaxIterations(maxIter);
                LineSearch::Result ret = lineSearchHandle->nativeLineSearch->BrentsMethod(x.data(), deltax.data(), eng);
                alpha                  = ret.alpha;
                eng1                   = ret.f;

                if (verbose >= 2 && iter % printGap == 0) {
                    std::cout << "        f=" << ret.f << ";alpha=" << alpha << std::endl;
                }
            } else if (solverParam.lsm == LSM_BACKTRACK) {
                lineSearchHandle->nativeLineSearch->setMaxIterations(maxIter);
                LineSearch::Result ret = lineSearchHandle->nativeLineSearch->backtracking(x.data(), deltax.data(), eng,
                                                                                          grad.data(), 0.0001, 0.5);
                alpha                  = ret.alpha;
                eng1                   = ret.f;

                if (verbose >= 2 && iter % printGap == 0) {
                    std::cout << "        f=" << ret.f << ";alpha=" << alpha << std::endl;
                }
            } else if (solverParam.lsm == LSM_SIMPLE) {
                eng1 = eng;
                for (int i = 0; i < 100; i++) {
                    lineSearchx.noalias() = x + deltax * alpha;
                    eng1                  = energy->func(lineSearchx);

                    if (eng1 < eng) {
                        break;
                    }

                    alpha *= 0.75;
                }

                if (verbose >= 2 && iter % printGap == 0) {
                    std::cout << "        f=" << eng1 << ";alpha=" << alpha << std::endl;
                }
            }

            if (!std::isfinite(alpha) || !std::isfinite(eng1)) {
                if (verbose >= 1) {
                    std::cout << "    Iter=" << iter << "; line search produced non-finite values." << std::endl;
                }

                status = SR_NUMERICAL_FAILURE;
                break;
            }

            if (eng1 > eng) {
                if (shouldTreatLineSearchFailureAsStall(alpha, deltaxMax, gradNorm, epsilon, eng, eng1)) {
                    if (verbose >= 1) {
                        std::cout << "    Iter=" << iter
                                  << "; line search stalled near a stationary point." << std::endl;
                    }

                    status = SR_STALLED_SMALL_STEP;
                } else {
                    if (verbose >= 1) {
                        std::cout << "    Iter=" << iter << "; line search failed." << std::endl;
                    }

                    status = SR_LINE_SEARCH_FAILED;
                }

                break;
            }

            x += deltax * alpha;

            stepSize = std::abs(alpha * deltaxMax);
            if (stepSize < kSmallAcceptedStepSizeTol) {
                if (verbose >= 1) {
                    std::cout << "    Iter=" << iter << "; dx = " << stepSize << "; dx too small." << std::endl;
                }

                status = SR_STALLED_SMALL_STEP;
                break;
            }
        } else if (solverParam.sst == SST_SUBITERATION_ONE) {
            x += deltax;
            historyx.noalias() = x;

            memset(grad.data(), 0, sizeof(double) * grad.size());
            energy->gradient(x, grad);
            filterVector(grad);

            historyGradNormMin = grad.norm();

            if (verbose >= 2 && iter % printGap == 0)
                std::cout << "    f=" << energy->func(x) << std::endl;
        } else if (solverParam.sst == SST_SUBITERATION_STATIC_DAMPING) {
            if (verbose >= 2 && iter % printGap == 0)
                std::cout << "        E= " << eng << "; ||grad||_max=" << grad.cwiseAbs().maxCoeff()
                          << "; ||grad||=" << grad.norm() << std::endl;

            x += deltax * solverParam.alpha;
        }

        if (stepFunc) {
            stepFunc(x, iter);
        }
    }

    if (iter == numIter && status == SR_MAX_ITER_REACHED && verbose >= 1) {
        std::cout << "    Reached max iterations (" << numIter << ")." << std::endl;
    }

    if (solverParam.sst == SST_SUBITERATION_ONE || solverParam.sst == SST_SUBITERATION_STATIC_DAMPING) {
        if (verbose >= 1)
            std::cout << "        Final ||grad||=" << historyGradNormMin << std::endl;

        Eigen::Map<ES::VXd>(x_, energy->getNumDOFs()) = historyx;
    } else
        Eigen::Map<ES::VXd>(x_, energy->getNumDOFs()) = x;

    hclock::time_point t2 = hclock::now();

    double timeCost = dura(t1, t2);

    if (auto logger = Logging::lgr()) {
        SPDLOG_LOGGER_INFO(logger, "Newton solve time: {}", timeCost);
    }
    return status;
}

void NewtonRaphsonSolver::filterVector(ES::VXd& v) {
    for (int dof : fixedDOFs)
        v[dof] = 0;
}

//
// OPK_DROP(lineSearchHandle);
