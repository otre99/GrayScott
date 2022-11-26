#ifndef GS_SOLVER_H
#define GS_SOLVER_H
#include <vector>

struct GrayScottParams {
  double k, f, mu, mv;
  std::size_t nrows, ncols;
  std::vector<int> index;
};

void TimeStepEuler(const double *u, const double *v, double *uOut, double *vOut,
                   const GrayScottParams &gs);
#endif // GS_SOLVER_H
