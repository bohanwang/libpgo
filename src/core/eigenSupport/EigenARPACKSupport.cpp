#include "EigenARPACKSupport.h"
#include "EigenMKLPardisoSupport.h"

#include <iostream>

using namespace pgo;
using namespace pgo::EigenSupport;

extern "C"
{
  void dsaupd_(int *IDO, char *BMAT, int *N, char *WHICH, int *NEV,
    double *TOL, double *RESID, int *NCV, double *V, int *LDV, int *IPARAM,
    int *IPNTR, double *WORKD, double *WORKL, int *LWORKL, int *INFO);

  // call DSAUPD
  // c     ( IDO, BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM,
  // c       IPNTR, WORKD, WORKL, LWORKL, INFO )

  void dseupd_(
    int *RVEC, char *HOWMNY, int *SELECT, double *D, double *Z, int *LDZ, double *SIGMA,
    char *BMAT, int *N, char *WHICH, int *NEV,
    double *TOL, double *RESID, int *NCV, double *V, int *LDV, int *IPARAM,
    int *IPNTR, double *WORKD, double *WORKL, int *LWORKL, int *INFO);

  // c  call DSEUPD
  // c     ( RVEC, HOWMNY, SELECT, D, Z, LDZ, SIGMA, BMAT, N, WHICH, NEV, TOL,
  // c       RESID, NCV, V, LDV, IPARAM, IPNTR, WORKD, WORKL, LWORKL, INFO )
}

int EigenSupport::solveGenEigShInv(const SpMatD &A, const SpMatD &B, int numEigenvalues, VXd &eigenvalues, MXd &eigenvectors, double sigma, int verbose)
{
  if (A.rows() != A.cols() || B.rows() != B.cols() || A.rows() != B.rows()) {
    std::cerr << "A or B are either not square matrix or their sizes are mismatched." << std::endl;
    return 1;
  }

  // solve Ax = lambda B x with ARPACK, shift-invert mode (mode number 3)
  // need multiplication with (A-sigma*B)^{-1}, (A-sigma*B)^{-1} B, and with B

  const SpMatD *AMinusSigmaB = nullptr;
  if (sigma != 0) {
    SpMatD *ptr = new SpMatD;
    *ptr = A - sigma * B;

    AMinusSigmaB = ptr;
  }
  else {
    AMinusSigmaB = &A;
  }

  // create (A-sigma*B)^{-1} solver
  // EigenPardisoSupport invASolver(*AMinusSigmaB, EigenPardisoSupport::MatrixType::REAL_SYM_INDEFINITE,
  //   EigenPardisoSupport::ReorderingType::NESTED_DISSECTION, 0, 0, 0, 0, 0, 0);

  EigenMKLPardisoSupport invASolver(*AMinusSigmaB, EigenMKLPardisoSupport::MatrixType::REAL_SYM_INDEFINITE,
    EigenMKLPardisoSupport::ReorderingType::NESTED_DISSECTION, 0, 0, 0, 0, 0, 0);

  // Eigen::PardisoLDLT<SpMatD> invASolver;
  // invASolver.analyze(*AMinusSigmaB);
  // invASolver.analyzePattern(*AMinusSigmaB);
  // invASolver.factorize(*AMinusSigmaB);
  invASolver.analyze(*AMinusSigmaB);
  invASolver.factorize(*AMinusSigmaB);

  // create (A-sigma*B)^{-1}*B solver
  const int maxIter = 10000;

  // call ARPACK
  int IDO = 0;
  char BMAT = 'G';
  int N = int(A.rows());
  char WHICH[3] = "LM";
  int NEV = numEigenvalues;
  double TOL = 0.0;

  VXd RESID(N);

  int NCV = std::min(3 * NEV, N);
  VXd V(N * NCV);

  int LDV = N;
  int IPARAM[11] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
  IPARAM[0] = 1;  // no user-provided shifts
  IPARAM[2] = maxIter;
  IPARAM[3] = 1;
  IPARAM[6] = 3;  // the mode (shift-inverted generalized eigenproblem)

  int IPNTR[11];
  VXd WORKD(3 * N);

  int LWORKL = 2 * (NCV * NCV + 8 * NCV);  // at least NCV**2 + 8*NCV;
  VXd WORKL(LWORKL);

  int INFO = 0;
  VXd buffer(N);

  if (verbose >= 1) {
    std::cout << "Dimension of the system            : " << A.rows() << std::endl;
    std::cout << "Number of 'requested' eigenvalues  : " << NEV << std::endl;
    std::cout << "Entering the ARPACK eigenvalue routine..." << std::endl;
  }

  do {
    dsaupd_(&IDO, &BMAT, &N, WHICH, &NEV, &TOL, RESID.data(), &NCV, V.data(), &LDV, IPARAM, IPNTR, WORKD.data(), WORKL.data(), &LWORKL, &INFO);
    if (INFO != 0) {
      std::cerr << "Error: DSAUPD returned non-zero exit code " << INFO << " .\n";
      if (sigma != 0)
        delete AMinusSigmaB;

      return 1;
    }

    // printf("IDO=%d\n", (int)IDO); fflush(NULL);
    // c           ===> OP = (inv[K - sigma*M])*M  and  B = M.

    switch (IDO) {
    case -1:
      // c          IDO = -1: compute  Y = OP * X  where
      // c                    IPNTR(1) is the pointer into WORKD for X,
      // c                    IPNTR(2) is the pointer into WORKD for Y.
      // c                    This is for the initialization phase to force the
      // c                    starting vector into the range of OP.
      // printf("IDO = -1\n");
      {
        double *X = WORKD.data() + IPNTR[0] - 1;
        double *Y = WORKD.data() + IPNTR[1] - 1;
        mv(B, Mp<const VXd>(X, B.cols()), buffer);
        invASolver.solve(Y, buffer.data(), 1);
        // Eigen::Map<VXd>(Y, N) = invASolver.solve(buffer);
      }
      break;
    case 1:
      // c          IDO =  1: compute  Y = (K - sigma*M)^-1 * Z
      // c                    IPNTR(2) is the pointer into WORKD for Y,
      // c                    IPNTR(3) is the pointer into WORKD for Z.
      // c           (see dsdrv4.f example)
      // printf("IDO = 1\n");
      {
        double *Y = WORKD.data() + IPNTR[1] - 1;
        double *Z = WORKD.data() + IPNTR[2] - 1;
        invASolver.solve(Y, Z, 1);
        // Eigen::Map<VXd>(Y, N) = invASolver.solve(Eigen::Map<VXd>(Z, N));
      }
      break;
    case 2:
      // c          IDO =  2: compute  Y = B * X  where
      // c                    IPNTR(1) is the pointer into WORKD for X,
      // c                    IPNTR(2) is the pointer into WORKD for Y.
      // printf("IDO = 2\n");
      {
        double *X = WORKD.data() + IPNTR[0] - 1;
        double *Y = WORKD.data() + IPNTR[1] - 1;
        mv(B, Mp<const VXd>(X, B.cols()), Mp<VXd>(Y, B.cols()));
      }
      break;
    case 3:
      std::cerr << "Error: case IDO=3 should have never happened.\n";
      if (sigma != 0)
        delete AMinusSigmaB;

      return 1;
    case 99:
      break;

    default:
      std::cerr << "Error: unknown case.\n";
      if (sigma != 0)
        delete AMinusSigmaB;
      return 1;
    };
  } while (IDO != 99);

  // obtain the eigenvalues and eigenvectors
  int RVEC = 1;
  char HOWMNY = 'A';  // all eigenvectors
  VXi SELECT(NCV);
  double *D = eigenvalues.data();
  double *Z = eigenvectors.data();
  int LDZ = N;
  double SIGMA = sigma;
  dseupd_(&RVEC, &HOWMNY, SELECT.data(), D, Z, &LDZ, &SIGMA,
    &BMAT, &N, WHICH, &NEV, &TOL, RESID.data(), &NCV, V.data(), &LDV, IPARAM, IPNTR, WORKD.data(), WORKL.data(), &LWORKL, &INFO);

  if (sigma != 0) {
    delete AMinusSigmaB;
  }

  int nconv = IPARAM[4];

  if (verbose >= 1) {
    std::cout << "ARPACK solver is done." << std::endl;
    // cout << "Dimension of the system            : " << np              << endl;
    // cout << "Number of 'requested' eigenvalues  : " << NEV << endl;
    std::cout << "Number of 'converged' eigenvalues  : " << nconv << std::endl;
    std::cout << "Number of Arnoldi vectors generated: " << NCV << std::endl;
    std::cout << "Number of iterations taken         : " << IPARAM[2] << std::endl;
    std::cout << std::endl;

    // Printing eigenvalues.

    std::cout << "Eigenvalues:" << std::endl;
    // for (int i=numEigenvalues-nconv; i<numEigenvalues; i++)
    for (int i = 0; i < nconv; i++) {
      std::cout << "  lambda[" << (i + 1) << "]: " << eigenvalues[i] << std::endl;
    }
    std::cout << std::endl;

    // Printing the residual norm || A*x - lambda*B*x ||
    // for the nconv accurately computed eigenvectors.

    // double *Ax = new double[np];
    // double *Bx = new double[np];
    // double *ResNorm = new double[nconv];

    // double infinityNormK = K->GetInfinityNorm();
    // double infinityNormM = M->GetInfinityNorm();
    // double infNorm = infinityNormK;
    // if (infinityNormM > infNorm)
    //   infNorm = infinityNormM;

    // for (int i = 0; i < nconv; i++) {
    //   K->MultiplyVector(&eigenvectors[np * i], Ax);
    //   M->MultiplyVector(&eigenvectors[np * i], Bx);
    //   for (int j = 0; j < np; j++)
    //     Ax[j] = Ax[j] - eigenvalues[i] * Bx[j];
    //   // cblas_daxpy(np, -eigenvalues[i], Bx, 1, Ax, 1);
    //   ResNorm[i] = 0;
    //   for (int j = 0; j < np; j++)
    //     ResNorm[i] += Ax[j] * Ax[j];
    //   ResNorm[i] = sqrt(ResNorm[i]) / fabs(eigenvalues[i]) / infNorm;
    //   // ResNorm[i] = cblas_dnrm2(np, Ax, 1) / fabs(eigenvalues[i]);
    // }

    /*
    for (int i=0; i<nconv; i++)
    {
      cout << "||A*x(" << (i+1) << ") - lambda(" << (i+1);
      cout << ")*B*x(" << (i+1) << ")|| / |lambda| / max(||A||,||B||) = " << ResNorm[i] << endl;
    }
    cout << endl;
*/

    // printf("Cleaning up ARPACK workspace.\n");fflush(NULL);
  }

  return 0;
}
