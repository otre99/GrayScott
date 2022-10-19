#include "solver.h"

#include <algorithm>
#include <execution>
#include <iostream>

#include "grayscott_utils.h"
#include "linalg.h"

void GrayScottParams::Init(const double &pk, const double &pf,
                           const double &pmu, const double &pmv,
                           const double &pdt, const double &pn) {
  k = pk;
  f = pf;
  mu = pmu;
  mv = pmv;
  dt = pdt;
  n = pn;

  uf = 0.5 * dt * mu;
  vf = 0.5 * dt * mv;

  au.resize(this->n, -uf);
  bu.resize(this->n, 1.0 + 2.0 * uf);
  cu.resize(this->n, -uf);

  av.resize(this->n, -vf);
  bv.resize(this->n, 1.0 + 2.0 * vf);
  cv.resize(this->n, -vf);

  index.resize(n);
  std::iota(index.begin(), index.end(), 0);
};

void TimeStepSymRK2(std::vector<double> &u, std::vector<double> &v,
                    const GrayScottParams &gs) {
  const auto exec_policy = std::execution::par_unseq;
  // xdir
  std::for_each(exec_policy, gs.index.begin(), gs.index.end(), [&](int i) {
    CNXDir(u.data() + i * gs.n, v.data() + i * gs.n, gs);
  });
  // ydir
  std::for_each(exec_policy, gs.index.begin(), gs.index.end(),
                [&](int i) { CNYDir(u.data() + i, v.data() + i, gs); });
  // reaction
  std::for_each(exec_policy, gs.index.begin(), gs.index.end(), [&](int i) {
    SolveReaction(u.data() + i * gs.n, v.data() + i * gs.n, gs);
  });
  // ydir
  std::for_each(exec_policy, gs.index.begin(), gs.index.end(),
                [&](int i) { CNYDir(u.data() + i, v.data() + i, gs); });
  // xdir
  std::for_each(exec_policy, gs.index.begin(), gs.index.end(), [&](int i) {
    CNXDir(u.data() + i * gs.n, v.data() + i * gs.n, gs);
  });
}

void TimeStepSymRK2Seq(std::vector<double> &u, std::vector<double> &v,
                       const GrayScottParams &gs) {
  std::for_each(gs.index.begin(), gs.index.end(), [&](int i) {
    CNXDir(u.data() + i * gs.n, v.data() + i * gs.n, gs);
  });
  // ydir
  std::for_each(gs.index.begin(), gs.index.end(),
                [&](int i) { CNYDir(u.data() + i, v.data() + i, gs); });
  // reaction
  std::for_each(gs.index.begin(), gs.index.end(), [&](int i) {
    SolveReaction(u.data() + i * gs.n, v.data() + i * gs.n, gs);
  });
  // ydir
  std::for_each(gs.index.begin(), gs.index.end(),
                [&](int i) { CNYDir(u.data() + i, v.data() + i, gs); });
  // xdir
  std::for_each(gs.index.begin(), gs.index.end(), [&](int i) {
    CNXDir(u.data() + i * gs.n, v.data() + i * gs.n, gs);
  });
}

void CNXDir(double *pu, double *pv, const GrayScottParams &gs) {
  const size_t N = gs.n;
  double *ru = new double[N];
  double *rv = new double[N];
  size_t j = 0;
  ru[j] = (pu[j + 1] + pu[N - 1]) * gs.uf + (1 - 2 * gs.uf) * pu[j];
  rv[j] = (pv[j + 1] + pv[N - 1]) * gs.vf + (1 - 2 * gs.vf) * pv[j];
  for (j = 1; j < N - 1; ++j) {
    ru[j] = (pu[j + 1] + pu[j - 1]) * gs.uf + (1 - 2 * gs.uf) * pu[j];
    rv[j] = (pv[j + 1] + pv[j - 1]) * gs.vf + (1 - 2 * gs.vf) * pv[j];
  }
  ru[j] = (pu[0] + pu[j - 1]) * gs.uf + (1 - 2 * gs.uf) * pu[j];
  rv[j] = (pv[0] + pv[j - 1]) * gs.vf + (1 - 2 * gs.vf) * pv[j];

  //  Tridiag(gs.au.data(), gs.bu.data(), gs.cu.data(), N, ru, pu);
  //  Tridiag(gs.av.data(), gs.bv.data(), gs.cv.data(), N, rv, pv);
  CyclicPackV1(gs.au.data(), gs.bu.data(), gs.cu.data(), N, ru, pu);
  CyclicPackV1(gs.av.data(), gs.bv.data(), gs.cv.data(), N, rv, pv);
  delete[] ru;
  delete[] rv;
}

void CNYDir(double *pu, double *pv, const GrayScottParams &gs) {
  const size_t N = gs.n;
  const size_t NL = (N - 1) * N;
  double *ru = new double[N];
  double *rv = new double[N];

  size_t j = 0;
  size_t jj = j * N;
  ru[j] = (pu[jj + N] + pu[NL]) * gs.uf + (1 - 2 * gs.uf) * pu[jj];
  rv[j] = (pv[jj + N] + pv[NL]) * gs.vf + (1 - 2 * gs.vf) * pv[jj];
  for (j = 1; j < N - 1; ++j) {
    jj = j * N;
    ru[j] = (pu[jj + N] + pu[jj - N]) * gs.uf + (1 - 2 * gs.uf) * pu[jj];
    rv[j] = (pv[jj + N] + pv[jj - N]) * gs.vf + (1 - 2 * gs.vf) * pv[jj];
  }
  jj = j * N;
  ru[j] = (pu[jj - N] + pu[0]) * gs.uf + (1 - 2 * gs.uf) * pu[jj];
  rv[j] = (pv[jj - N] + pv[0]) * gs.vf + (1 - 2 * gs.vf) * pv[jj];

  //  Tridiag(gs.au.data(), gs.bu.data(), gs.cu.data(), N, ru, ru);
  //  Tridiag(gs.av.data(), gs.bv.data(), gs.cv.data(), N, rv, rv);
  CyclicPackV1(gs.au.data(), gs.bu.data(), gs.cu.data(), N, ru, ru);
  CyclicPackV1(gs.av.data(), gs.bv.data(), gs.cv.data(), N, rv, rv);
  for (j = 0; j < N; ++j) {
    jj = j * N;
    pu[jj] = ru[j];
    pv[jj] = rv[j];
  }
  delete[] ru;
  delete[] rv;
}

void SolveReaction(double *u, double *v, const GrayScottParams &gs) {
  double uv2, uk1, vk1, uk2, vk2;
  for (size_t k = 0; k < gs.n; ++k) {
    uv2 = u[k] * v[k] * v[k];
    uk1 = u[k] + 0.5 * gs.dt * (-uv2 + gs.f * (1 - u[k]));
    vk1 = v[k] + 0.5 * gs.dt * (uv2 - (gs.f + gs.k) * v[k]);

    uv2 = uk1 * vk1 * vk1;
    uk2 = -uv2 + gs.f * (1 - uk1);
    vk2 = uv2 - (gs.f + gs.k) * vk1;
    u[k] += gs.dt * uk2;
    v[k] += gs.dt * vk2;
    // ----------------------------------------
    uv2 = u[k] * v[k] * v[k];
    uk1 = u[k] + 0.5 * gs.dt * (-uv2 + gs.f * (1 - u[k]));
    vk1 = v[k] + 0.5 * gs.dt * (uv2 - (gs.f + gs.k) * v[k]);

    uv2 = uk1 * vk1 * vk1;
    uk2 = -uv2 + gs.f * (1 - uk1);
    vk2 = uv2 - (gs.f + gs.k) * vk1;
    u[k] += gs.dt * uk2;
    v[k] += gs.dt * vk2;
  }
}

/////////////////////////////////////////////////////////////////
void EulerInOneRow(const double *u, const double *v, double *uOut, double *vOut,
                   const GrayScottParams &gs) {
  double uv2;
  size_t i = 0, j, k;

  uv2 = u[i] * v[i] * v[i];
  uOut[i] = u[i] +
            gs.mu * (u[i + gs.n - 1] + u[i + 1] + u[i - gs.n] + u[i + gs.n] -
                     4 * u[i]) -
            uv2 + gs.f * (1 - u[i]);
  vOut[i] = v[i] +
            gs.mv * (v[i + gs.n - 1] + v[i + 1] + v[i - gs.n] + v[i + gs.n] -
                     4 * v[i]) +
            uv2 - (gs.f + gs.k) * v[i];
  for (j = 1; j < gs.n - 1; ++j) {
    k = i + j;
    uv2 = u[k] * v[k] * v[k];
    uOut[k] =
        u[k] +
        gs.mu * (u[k - 1] + u[k + 1] + u[k - gs.n] + u[k + gs.n] - 4 * u[k]) -
        uv2 + gs.f * (1 - u[k]);
    vOut[k] =
        v[k] +
        gs.mv * (v[k - 1] + v[k + 1] + v[k - gs.n] + v[k + gs.n] - 4 * v[k]) +
        uv2 - (gs.f + gs.k) * v[k];
  }
  k = i + j;
  uv2 = u[k] * v[k] * v[k];
  uOut[k] = u[k] +
            gs.mu * (u[k - 1] + u[i] + u[k - gs.n] + u[k + gs.n] - 4 * u[k]) -
            uv2 + gs.f * (1 - u[k]);
  vOut[k] = v[k] +
            gs.mv * (v[k - 1] + v[i] + v[k - gs.n] + v[k + gs.n] - 4 * v[k]) +
            uv2 - (gs.f + gs.k) * v[k];
}

void TimeStepEuler(std::vector<double> &u, std::vector<double> &v,
                   std::vector<double> &uT, std::vector<double> &vT,
                   const GrayScottParams &gs) {
  const auto exec_policy = std::execution::par_unseq;
  std::for_each(exec_policy, gs.index.begin() + 1, gs.index.end() - 1,
                [&](int i) {
                  const size_t o = i * gs.n;
                  EulerInOneRow(u.data() + o, v.data() + o, uT.data() + o,
                                vT.data() + o, gs);
                });

  double uv2;
  const size_t ll = gs.n * (gs.n - 1);
  size_t it = 1, ib = 1 + ll;
  for (; it < gs.n - 1; ++it, ++ib) {
    uv2 = u[it] * v[it] * v[it];
    uT[it] =
        u[it] +
        gs.mu * (u[it - 1] + u[it + 1] + u[ib] + u[gs.n + it] - 4 * u[it]) -
        uv2 + gs.f * (1 - u[it]);
    vT[it] =
        v[it] +
        gs.mv * (v[it - 1] + v[it + 1] + v[ib] + v[gs.n + it] - 4 * v[it]) +
        uv2 - (gs.f + gs.k) * v[it];

    uv2 = u[ib] * v[ib] * v[ib];
    uT[ib] =
        u[ib] +
        gs.mu * (u[ib - 1] + u[ib + 1] + u[it] + u[ib - gs.n] - 4 * u[ib]) -
        uv2 + gs.f * (1 - u[ib]);
    vT[ib] =
        v[ib] +
        gs.mv * (v[ib - 1] + v[ib + 1] + v[it] + v[ib - gs.n] - 4 * v[ib]) +
        uv2 - (gs.f + gs.k) * v[ib];
  }

  // tl
  size_t i = 0;
  uv2 = u[i] * v[i] * v[i];
  uT[i] = u[i] + gs.mu * (u[ll] + u[gs.n] + u[1] + u[gs.n - 1] - 4 * u[i]) -
          uv2 + gs.f * (1 - u[i]);
  vT[i] = v[i] + gs.mv * (v[ll] + v[gs.n] + v[1] + v[gs.n - 1] - 4 * v[i]) +
          uv2 - (gs.f + gs.k) * v[i];
  // tr
  i = gs.n - 1;
  uv2 = u[i] * v[i] * v[i];
  uT[i] =
      u[i] +
      gs.mu * (u[i - 1] + u[0] + u[ll + gs.n - 1] + u[gs.n + i] - 4 * u[i]) -
      uv2 + gs.f * (1 - u[i]);
  vT[i] =
      v[i] +
      gs.mv * (v[i - 1] + v[0] + v[ll + gs.n - 1] + v[gs.n + i] - 4 * v[i]) +
      uv2 - (gs.f + gs.k) * v[i];
  // bl
  i = ll;
  uv2 = u[i] * v[i] * v[i];
  uT[i] = u[i] +
          gs.mu * (u[i + gs.n - 1] + u[i + 1] + u[i - gs.n] + u[0] - 4 * u[i]) -
          uv2 + gs.f * (1 - u[i]);
  vT[i] = v[i] +
          gs.mv * (v[i + gs.n - 1] + v[i + 1] + v[i - gs.n] + v[0] - 4 * v[i]) +
          uv2 - (gs.f + gs.k) * v[i];
  // br
  i = ll + gs.n - 1;
  uv2 = u[i] * v[i] * v[i];
  uT[i] = u[i] +
          gs.mu * (u[ll] + u[i - 1] + u[gs.n - 1] + u[i - gs.n] - 4 * u[i]) -
          uv2 + gs.f * (1 - u[i]);
  vT[i] = v[i] +
          gs.mv * (v[ll] + v[i - 1] + v[gs.n - 1] + v[i - gs.n] - 4 * v[i]) +
          uv2 - (gs.f + gs.k) * v[i];

  std::swap(u, uT);
  std::swap(v, vT);
}

void TimeStepEulerSeq(std::vector<double> &u, std::vector<double> &v,
                      std::vector<double> &uTemp, std::vector<double> &vTemp,
                      const GrayScottParams &gs) {
  const size_t ll = gs.n * (gs.n - 1);
  double uv2;
  size_t i, j, k;
  for (i = gs.n; i < ll; i += gs.n) {
    uv2 = u[i] * v[i] * v[i];
    uTemp[i] = u[i] +
               gs.mu * (u[i + gs.n - 1] + u[i + 1] + u[i - gs.n] + u[i + gs.n] -
                        4 * u[i]) -
               uv2 + gs.f * (1 - u[i]);
    vTemp[i] = v[i] +
               gs.mv * (v[i + gs.n - 1] + v[i + 1] + v[i - gs.n] + v[i + gs.n] -
                        4 * v[i]) +
               uv2 - (gs.f + gs.k) * v[i];
    for (j = 1; j < gs.n - 1; ++j) {
      k = i + j;
      uv2 = u[k] * v[k] * v[k];
      uTemp[k] =
          u[k] +
          gs.mu * (u[k - 1] + u[k + 1] + u[k - gs.n] + u[k + gs.n] - 4 * u[k]) -
          uv2 + gs.f * (1 - u[k]);
      vTemp[k] =
          v[k] +
          gs.mv * (v[k - 1] + v[k + 1] + v[k - gs.n] + v[k + gs.n] - 4 * v[k]) +
          uv2 - (gs.f + gs.k) * v[k];
    }
    k = i + j;
    uv2 = u[k] * v[k] * v[k];
    uTemp[k] =
        u[k] +
        gs.mu * (u[k - 1] + u[i] + u[k - gs.n] + u[k + gs.n] - 4 * u[k]) - uv2 +
        gs.f * (1 - u[k]);
    vTemp[k] =
        v[k] +
        gs.mv * (v[k - 1] + v[i] + v[k - gs.n] + v[k + gs.n] - 4 * v[k]) + uv2 -
        (gs.f + gs.k) * v[k];
  }

  size_t it = 1, ib = 1 + ll;
  for (; it < gs.n - 1; ++it, ++ib) {
    uv2 = u[it] * v[it] * v[it];
    uTemp[it] =
        u[it] +
        gs.mu * (u[it - 1] + u[it + 1] + u[ib] + u[gs.n + it] - 4 * u[it]) -
        uv2 + gs.f * (1 - u[it]);
    vTemp[it] =
        v[it] +
        gs.mv * (v[it - 1] + v[it + 1] + v[ib] + v[gs.n + it] - 4 * v[it]) +
        uv2 - (gs.f + gs.k) * v[it];

    uv2 = u[ib] * v[ib] * v[ib];
    uTemp[ib] =
        u[ib] +
        gs.mu * (u[ib - 1] + u[ib + 1] + u[it] + u[ib - gs.n] - 4 * u[ib]) -
        uv2 + gs.f * (1 - u[ib]);
    vTemp[ib] =
        v[ib] +
        gs.mv * (v[ib - 1] + v[ib + 1] + v[it] + v[ib - gs.n] - 4 * v[ib]) +
        uv2 - (gs.f + gs.k) * v[ib];
  }

  // tl
  i = 0;
  uv2 = u[i] * v[i] * v[i];
  uTemp[i] = u[i] + gs.mu * (u[ll] + u[gs.n] + u[1] + u[gs.n - 1] - 4 * u[i]) -
             uv2 + gs.f * (1 - u[i]);
  vTemp[i] = v[i] + gs.mv * (v[ll] + v[gs.n] + v[1] + v[gs.n - 1] - 4 * v[i]) +
             uv2 - (gs.f + gs.k) * v[i];
  // tr
  i = gs.n - 1;
  uv2 = u[i] * v[i] * v[i];
  uTemp[i] =
      u[i] +
      gs.mu * (u[i - 1] + u[0] + u[ll + gs.n - 1] + u[gs.n + i] - 4 * u[i]) -
      uv2 + gs.f * (1 - u[i]);
  vTemp[i] =
      v[i] +
      gs.mv * (v[i - 1] + v[0] + v[ll + gs.n - 1] + v[gs.n + i] - 4 * v[i]) +
      uv2 - (gs.f + gs.k) * v[i];
  // bl
  i = ll;
  uv2 = u[i] * v[i] * v[i];
  uTemp[i] =
      u[i] +
      gs.mu * (u[i + gs.n - 1] + u[i + 1] + u[i - gs.n] + u[0] - 4 * u[i]) -
      uv2 + gs.f * (1 - u[i]);
  vTemp[i] =
      v[i] +
      gs.mv * (v[i + gs.n - 1] + v[i + 1] + v[i - gs.n] + v[0] - 4 * v[i]) +
      uv2 - (gs.f + gs.k) * v[i];
  // br
  i = ll + gs.n - 1;
  uv2 = u[i] * v[i] * v[i];
  uTemp[i] = u[i] +
             gs.mu * (u[ll] + u[i - 1] + u[gs.n - 1] + u[i - gs.n] - 4 * u[i]) -
             uv2 + gs.f * (1 - u[i]);
  vTemp[i] = v[i] +
             gs.mv * (v[ll] + v[i - 1] + v[gs.n - 1] + v[i - gs.n] - 4 * v[i]) +
             uv2 - (gs.f + gs.k) * v[i];

  u.swap(uTemp);
  v.swap(vTemp);
}
