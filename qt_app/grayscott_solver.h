#ifndef GRAYSCOTT_SOLVER_H
#define GRAYSCOTT_SOLVER_H
#include <QThread>
#include <vector>

#include "solver.h"

class GrayScottSolver : public QObject {
  Q_OBJECT
  QThread *thread_;

 public:
  GrayScottSolver();
  ~GrayScottSolver() override;
  void initParams(int dim, double du, double dv, double f, double k);
  void initParams(int dim, double du, double dv, double f, double k, double dt, bool reset, int method);
  void initialize(int pattern, bool with_noise);
  void initialize(const QImage &img);
  void updateDpData();
  size_t getSize() { return m_gsParams.n; }
  const std::vector<double> &getData() const { return m_dataDp; }
  double currTime() const {return m_currTime;}
  bool isRunning() { return !m_loopBreak; }
  void selectUOrV(bool select_u) { m_uIsSelected = select_u; }
  void enableParallel(bool enable) { m_enableParallel = enable; }

 public slots:
  void stop() { m_loopBreak = true; }
  void solveSymRK2(int n);
  void solveEuler(int n);

 private:
  std::vector<double> m_u;
  std::vector<double> m_v;
  std::vector<double> m_uTemp;
  std::vector<double> m_vTemp;
  std::vector<double> m_dataDp;

  bool m_loopBreak;
  bool m_uIsSelected;
  bool m_enableParallel{true};
  GrayScottParams m_gsParams;
  double m_currTime{0.0};
 signals:
  void dataReady();
};

#endif  // GRAYSCOTT_SOLVER_H
