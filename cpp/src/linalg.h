#ifndef LINALG_H
#define LINALG_H

void Tridiag(const double *a, const double *b, const double *c, const int n,
             const double *r, double *x);
void CyclicPackV1(const double *a, const double *b, const double *c, int n,
                  const double *r, double *x);

#endif  // LINALG_H
