#include "grayscott_solver.h"

#include <QApplication>
#include <QDebug>
#include <QImage>
#include <algorithm>
#include <future>
#include <random>

#include "grayscott_utils.h"
#include "solver.h"

GrayScottSolver::GrayScottSolver() : m_loopBreak(true), m_uIsSelected(true) {
  thread_ = new QThread();
  this->moveToThread(thread_);
  thread_->start();
}

GrayScottSolver::~GrayScottSolver() {
  m_loopBreak = true;
  thread_->quit();
  thread_->wait();
  delete thread_;
}

void GrayScottSolver::updateDpData() {
  if (m_uIsSelected)
    m_dataDp = m_u;
  else
    m_dataDp = m_v;
}

void GrayScottSolver::initParams(int dim, double du, double dv, double f,
                                 double k, double dt, bool reset, int method) {
  m_gsParams.Init(k, f, du, dv, dt, dim);
  if (reset) {
    m_u.resize(dim * dim);
    m_v.resize(dim * dim);
  }
  if (method == 0) {
    m_uTemp.resize(m_u.size());
    m_vTemp.resize(m_v.size());
  } else {
    m_uTemp.clear();
    m_vTemp.clear();
  }
}

void GrayScottSolver::initialize(int pattern, bool with_noise) {
  switch (pattern) {
    case 0:
      InitializeP0(m_u, m_v, with_noise);
      break;
    case 1:
      InitializeP1(m_u, m_v, with_noise);
      break;
    case 2:
      InitializeP2(m_u, m_v, with_noise);
      break;
    default:
      break;
  }
  m_currTime = 0.0;
  updateDpData();
}

void GrayScottSolver::solveSymRK2(int n) {
  m_loopBreak = false;
  while (true) {
    if (m_loopBreak) break;
    int i;
    for (i = 0; i < n; ++i) {
      (m_enableParallel) ? TimeStepSymRK2(m_u, m_v, m_gsParams)
                         : TimeStepSymRK2Seq(m_u, m_v, m_gsParams);
      if (m_loopBreak) {
        break;
      }
    }
    m_currTime += 2 * i * m_gsParams.dt;
    updateDpData();
    emit dataReady();
  }
}

void GrayScottSolver::solveEuler(int n) {
  m_loopBreak = false;
  while (true) {
    if (m_loopBreak) break;
    int i;
    for (i = 0; i < n; ++i) {
      (m_enableParallel)
          ? TimeStepEuler(m_u, m_v, m_uTemp, m_vTemp, m_gsParams)
          : TimeStepEulerSeq(m_u, m_v, m_uTemp, m_vTemp, m_gsParams);
      if (m_loopBreak) {
        break;
      }
    }
    m_currTime += i;
    updateDpData();
    emit dataReady();
  }
}
