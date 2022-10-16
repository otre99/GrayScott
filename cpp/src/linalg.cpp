#include "linalg.h"
#include <cstring>
#include <algorithm>
#include <cstring>

void Tridiag(const double *a, const double *b, const double *c, const int n,
             const double *r, double *x) {
  int j;
  double bet;
  double gam[n];
  if (b[0] == 0.0) throw("Singular matrix in Tridiag2()");
  x[0] = r[0] / (bet = b[0]);
  for (j = 1; j < n; j++) {
    gam[j] = c[j - 1] / bet;
    bet = b[j] - a[j] * gam[j];
    if (bet == 0.0) throw("Singular matrix in Tridiag2()");
    x[j] = (r[j] - a[j] * x[j - 1]) / bet;
  }
  for (j = (n - 2); j >= 0; j--) {
    x[j] -= gam[j + 1] * x[j + 1];
  }
}

void Tridiag2(const double *a, const double *b, const double *c, const int n,
              const double *r1, const double *r2, double *x1, double *u2) {
  int j;
  double bet;
  double gam[n];
  if (b[0] == 0.0) throw("Singular matrix in Tridiag2()");
  x1[0] = r1[0] / (bet = b[0]);
  u2[0] = r2[0] / bet;
  for (j = 1; j < n; j++) {
    gam[j] = c[j - 1] / bet;
    bet = b[j] - a[j] * gam[j];
    if (bet == 0.0) throw("Singular matrix in Tridiag2()");
    x1[j] = (r1[j] - a[j] * x1[j - 1]) / bet;
    u2[j] = (r2[j] - a[j] * u2[j - 1]) / bet;
  }
  for (j = (n - 2); j >= 0; j--) {
    x1[j] -= gam[j + 1] * x1[j + 1];
    u2[j] -= gam[j + 1] * u2[j + 1];
  }
}

void CyclicPack(const double *a, double *b, const double *c, int n,
                const double *r, double *x) {
  int i;
  double fact, gamma;
  const double alpha = c[n - 1];
  const double beta = a[0];
  if (n <= 2) throw("System dimension is too small in cyclic");
  double *z = new double[n];

  gamma = -b[0];
  b[0] -= gamma;
  b[n - 1] -= alpha * beta / gamma;
  z[0] = gamma;
  z[n - 1] = alpha;
  //std::fill(z + 1, z + n - 1, 0.0);
  std::memset(z+1, 0, (n-2)*sizeof(double));
  Tridiag2(a, b, c, n, r, z, x, z);
  fact =
      (x[0] + beta * x[n - 1] / gamma) / (1.0 + z[0] + beta * z[n - 1] / gamma);
  for (i = 0; i < n; i++) x[i] -= fact * z[i];
  delete[] z;
}

void CyclicPackV1(const double *a, const double *b, const double *c, int n,
                  const double *r, double *x) {
  double *bb = new double[n];
  //std::copy(b, b + n, bb);
  std::memcpy(bb,b,n*sizeof(double));
  CyclicPack(a, bb, c, n, r, x);
  delete[] bb;
}
