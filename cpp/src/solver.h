#ifndef SOLVER_H
#define SOLVER_H
#include <cstddef>
#include <vector>

void PrintVector(double *v, int n, bool add_endl = true);

struct GrayScottParams {
  void Init(const double &pk, const double &pf, const double &pmu,
            const double &pmv, const double &pdt, const double &pn);
  size_t n;
  double k, f, mu, mv;
  double uf, vf;
  double dt;
  std::vector<double> au, bu, cu;
  std::vector<double> av, bv, cv;
  std::vector<int> index;
};

void TimeStepSymRK2(std::vector<double> &u, std::vector<double> &v,
                    const GrayScottParams &);
void TimeStepSymRK2Seq(std::vector<double> &u, std::vector<double> &v,
                       const GrayScottParams &);
void CNXDir(double *pu, double *pv, const GrayScottParams &gs);
void CNYDir(double *pu, double *pv, const GrayScottParams &gs);
void SolveReaction(double *u, double *v, const GrayScottParams &gs);

void TimeStepEuler(std::vector<double> &u, std::vector<double> &v,
                   std::vector<double> &uTemp, std::vector<double> &vTemp,
                   const GrayScottParams &gs);
void TimeStepEulerSeq(std::vector<double> &u, std::vector<double> &v,
                      std::vector<double> &uTemp, std::vector<double> &vTemp,
                      const GrayScottParams &gs);
#endif  // SOLVER_H
