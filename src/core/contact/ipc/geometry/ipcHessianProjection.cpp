/*
copyright to Bohan Wang
*/

#include "ipcHessianProjection.h"

namespace pgo {
namespace Contact {
namespace IPC {
using namespace pgo::EigenSupport;

// =========================================================================
//  PSD projection of 12x12 symmetric matrix
// =========================================================================
M12d projectToPSD(const M12d &H)
{
  Eigen::SelfAdjointEigenSolver<M12d> es(H);
  const auto &evals = es.eigenvalues();

  // Early exit: eigenvalues are sorted ascending, so if the smallest is >= 0
  // the matrix is already PSD.
  if (evals(0) >= 0.0)
    return H;

  // Clamp negative eigenvalues to zero
  Eigen::DiagonalMatrix<double, 12> D(evals);
  for (int i = 0; i < 12; ++i) {
    if (D.diagonal()(i) < 0.0)
      D.diagonal()(i) = 0.0;
    else
      break;  // remaining eigenvalues are >= 0
  }

  return es.eigenvectors() * D * es.eigenvectors().transpose();
}

}  // namespace IPC
}  // namespace Contact
}  // namespace pgo
