#include "gs_solver.h"
#include <algorithm>
#include <execution>
#include <iostream>

void EulerInOneRow(const double *u, const double *v, double *uOut, double *vOut,
                   const GrayScottParams &gs) {
  double uv2;
  size_t i = 0, j, k;

  uv2 = u[i] * v[i] * v[i];
  uOut[i] = u[i] +
            gs.mu * (u[i + gs.ncols - 1] + u[i + 1] + u[i - gs.ncols] +
                     u[i + gs.ncols] - 4 * u[i]) -
            uv2 + gs.f * (1 - u[i]);
  vOut[i] = v[i] +
            gs.mv * (v[i + gs.ncols - 1] + v[i + 1] + v[i - gs.ncols] +
                     v[i + gs.ncols] - 4 * v[i]) +
            uv2 - (gs.f + gs.k) * v[i];
  for (j = 1; j < gs.ncols - 1; ++j) {
    k = i + j;
    uv2 = u[k] * v[k] * v[k];
    uOut[k] = u[k] +
              gs.mu * (u[k - 1] + u[k + 1] + u[k - gs.ncols] + u[k + gs.ncols] -
                       4 * u[k]) -
              uv2 + gs.f * (1 - u[k]);
    vOut[k] = v[k] +
              gs.mv * (v[k - 1] + v[k + 1] + v[k - gs.ncols] + v[k + gs.ncols] -
                       4 * v[k]) +
              uv2 - (gs.f + gs.k) * v[k];
  }
  k = i + j;
  uv2 = u[k] * v[k] * v[k];
  uOut[k] =
      u[k] +
      gs.mu * (u[k - 1] + u[i] + u[k - gs.ncols] + u[k + gs.ncols] - 4 * u[k]) -
      uv2 + gs.f * (1 - u[k]);
  vOut[k] =
      v[k] +
      gs.mv * (v[k - 1] + v[i] + v[k - gs.ncols] + v[k + gs.ncols] - 4 * v[k]) +
      uv2 - (gs.f + gs.k) * v[k];
}

void TimeStepEuler(const double *u, const double *v, double *uOut, double *vOut,
                   const GrayScottParams &gs) {

  const auto exec_policy = std::execution::par_unseq;
  std::for_each(exec_policy, gs.index.begin() + 1, gs.index.end() - 1,
                [&](int i) {
                  const size_t o = i * gs.ncols;
                  EulerInOneRow(u + o, v + o, uOut + o, vOut + o, gs);
                });

  double uv2;
  const size_t ll = (gs.nrows - 1) * (gs.ncols);
  size_t it = 1, ib = 1 + ll;
  for (; it < gs.ncols - 1; ++it, ++ib) {
    uv2 = u[it] * v[it] * v[it];
    uOut[it] =
        u[it] +
        gs.mu * (u[it - 1] + u[it + 1] + u[ib] + u[gs.ncols + it] - 4 * u[it]) -
        uv2 + gs.f * (1 - u[it]);
    vOut[it] =
        v[it] +
        gs.mv * (v[it - 1] + v[it + 1] + v[ib] + v[gs.ncols + it] - 4 * v[it]) +
        uv2 - (gs.f + gs.k) * v[it];

    uv2 = u[ib] * v[ib] * v[ib];
    uOut[ib] =
        u[ib] +
        gs.mu * (u[ib - 1] + u[ib + 1] + u[it] + u[ib - gs.ncols] - 4 * u[ib]) -
        uv2 + gs.f * (1 - u[ib]);
    vOut[ib] =
        v[ib] +
        gs.mv * (v[ib - 1] + v[ib + 1] + v[it] + v[ib - gs.ncols] - 4 * v[ib]) +
        uv2 - (gs.f + gs.k) * v[ib];
  }

  // tl
  size_t i = 0;
  uv2 = u[i] * v[i] * v[i];
  uOut[i] = u[i] +
            gs.mu * (u[ll] + u[gs.ncols] + u[1] + u[gs.ncols - 1] - 4 * u[i]) -
            uv2 + gs.f * (1 - u[i]);
  vOut[i] = v[i] +
            gs.mv * (v[ll] + v[gs.ncols] + v[1] + v[gs.ncols - 1] - 4 * v[i]) +
            uv2 - (gs.f + gs.k) * v[i];
  // tr
  i = gs.ncols - 1;
  uv2 = u[i] * v[i] * v[i];
  uOut[i] = u[i] +
            gs.mu * (u[i - 1] + u[0] + u[ll + gs.ncols - 1] + u[gs.ncols + i] -
                     4 * u[i]) -
            uv2 + gs.f * (1 - u[i]);
  vOut[i] = v[i] +
            gs.mv * (v[i - 1] + v[0] + v[ll + gs.ncols - 1] + v[gs.ncols + i] -
                     4 * v[i]) +
            uv2 - (gs.f + gs.k) * v[i];
  // bl
  i = ll;
  uv2 = u[i] * v[i] * v[i];
  uOut[i] = u[i] +
            gs.mu * (u[i + gs.ncols - 1] + u[i + 1] + u[i - gs.ncols] + u[0] -
                     4 * u[i]) -
            uv2 + gs.f * (1 - u[i]);
  vOut[i] = v[i] +
            gs.mv * (v[i + gs.ncols - 1] + v[i + 1] + v[i - gs.ncols] + v[0] -
                     4 * v[i]) +
            uv2 - (gs.f + gs.k) * v[i];
  // br
  i = ll + gs.ncols - 1;
  uv2 = u[i] * v[i] * v[i];
  uOut[i] = u[i] +
            gs.mu * (u[ll] + u[i - 1] + u[gs.ncols - 1] + u[i - gs.ncols] -
                     4 * u[i]) -
            uv2 + gs.f * (1 - u[i]);
  vOut[i] = v[i] +
            gs.mv * (v[ll] + v[i - 1] + v[gs.ncols - 1] + v[i - gs.ncols] -
                     4 * v[i]) +
            uv2 - (gs.f + gs.k) * v[i];
}
