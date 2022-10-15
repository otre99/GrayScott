#include "solver.h"

#include <algorithm>
#include <execution>
#include <iostream>

#include "linalg.h"
void PrintVector(double *v, int n, bool add_endl) {
  std::cout << v[0];
  for (int i = 1; i < n; ++i) {
    std::cout << ' ' << v[i];
  }
  if (add_endl) std::cout << std::endl;
}

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

void TimeStep(std::vector<double> &u, std::vector<double> &v,
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
