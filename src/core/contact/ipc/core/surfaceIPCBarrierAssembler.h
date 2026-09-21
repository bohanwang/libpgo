/*
copyright to Bohan Wang
*/

#pragma once

#include "EigenDef.h"
#include "ipc/core/surfaceIPCPairs.h"

#include <cstddef>
#include <vector>

namespace pgo
{
namespace Contact
{
namespace CIPC
{

// Destination-major map from surface vertices to the per-pair gradient blocks
// that touch them. Pair index i < numPT refers to ptPairs[i] and numPT + j to
// eePairs[j]; block 4 * pair + slot is the 3-vector that pair writes for its
// slot-th vertex. Summing a vertex's blocks in the stored order makes the
// gradient assembly bit-reproducible while the vertices are still gathered in
// parallel.
struct SurfaceIPCGradientGatherMap
{
  std::vector<std::ptrdiff_t> start;   // numVerts + 1 offsets into blocks
  std::vector<std::ptrdiff_t> blocks;  // block indices grouped by vertex

  void build(const std::vector<PTPair> &ptPairs, const std::vector<EEPair> &eePairs, int numVerts);
  int numVertices() const { return start.empty() ? 0 : static_cast<int>(start.size()) - 1; }
  std::size_t numBlocks() const { return blocks.size(); }
};

class SurfaceIPCBarrierAssembler
{
public:
  // deterministic = true sums the per-pair energies in pair order instead of
  // a parallel reduction, so the result depends only on the pair order.
  double computeEnergy(
    EigenSupport::ConstRefVecXd pos,
    const std::vector<PTPair> &ptPairs,
    const std::vector<EEPair> &eePairs,
    int numVerts,
    double dhat,
    double kappa,
    double eps_ee,
    bool deterministic = false) const;

  // A non-null gather map selects the deterministic path: per-pair gradient
  // blocks are stored and then gathered per vertex in the map's fixed order.
  // The map must have been built for exactly these pairs. Without a map the
  // blocks are accumulated with atomic adds in scheduling order.
  void computeGradient(
    EigenSupport::ConstRefVecXd pos,
    const std::vector<PTPair> &ptPairs,
    const std::vector<EEPair> &eePairs,
    int numVerts,
    double dhat,
    double kappa,
    double eps_ee,
    EigenSupport::RefVecXd grad,
    const SurfaceIPCGradientGatherMap *gather = nullptr) const;

  // Each pair owns a fixed slice of the triplet list, which setFromTriplets
  // sums sequentially, so the Hessian is bit-reproducible whenever the pair
  // order is.
  void computeHessian(
    EigenSupport::ConstRefVecXd pos,
    const std::vector<PTPair> &ptPairs,
    const std::vector<EEPair> &eePairs,
    int numVerts,
    double dhat,
    double kappa,
    double eps_ee,
    bool projectHessianToPSD,
    EigenSupport::SpMatD &hess) const;

  void computeAll(
    EigenSupport::ConstRefVecXd x,
    const std::vector<PTPair> &ptPairs,
    const std::vector<EEPair> &eePairs,
    int numVerts,
    double dhat,
    double kappa,
    double eps_ee,
    bool projectHessianToPSD,
    double &energy,
    EigenSupport::VXd &grad,
    EigenSupport::SpMatD &hess,
    const SurfaceIPCGradientGatherMap *gather = nullptr) const;
};

}  // namespace CIPC
}  // namespace Contact
}  // namespace pgo
