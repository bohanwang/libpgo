#include <math.h>

typedef struct
{
  float real;
  float imag;
} complex_float;

extern void cgemm_(const char *transa, const char *transb,
  const int *m, const int *n, const int *k,
  const complex_float *alpha, const complex_float *a, const int *lda,
  const complex_float *b, const int *ldb,
  const complex_float *beta, complex_float *c, const int *ldc);

extern void dgemm_(const char *transa, const char *transb,
  const int *m, const int *n, const int *k,
  const double *alpha, const double *a, const int *lda,
  const double *b, const int *ldb,
  const double *beta, double *c, const int *ldc);

int main(void)
{
  const char noTranspose = 'N';
  const int one = 1;

  const double dAlpha = 1.0;
  const double dBeta = 0.0;
  const double dA = 2.0;
  const double dB = 3.0;
  double dC = 0.0;
  dgemm_(&noTranspose, &noTranspose, &one, &one, &one,
    &dAlpha, &dA, &one, &dB, &one, &dBeta, &dC, &one);
  if(fabs(dC - 6.0) > 1e-12)
    return 1;

  const complex_float cAlpha = { 1.0f, 0.0f };
  const complex_float cBeta = { 0.0f, 0.0f };
  const complex_float cA = { 2.0f, 1.0f };
  const complex_float cB = { 3.0f, -1.0f };
  complex_float cC = { 0.0f, 0.0f };
  cgemm_(&noTranspose, &noTranspose, &one, &one, &one,
    &cAlpha, &cA, &one, &cB, &one, &cBeta, &cC, &one);
  if(fabsf(cC.real - 7.0f) > 1e-5f || fabsf(cC.imag - 1.0f) > 1e-5f)
    return 2;

  return 0;
}
