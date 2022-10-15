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
                                 double k, double dt) {
  m_gsParams.Init(k, f, du, dv, dt, dim);
  m_u.resize(dim * dim);
  m_v.resize(dim * dim);
}

void GrayScottSolver::initialize(int pattern) {
  switch (pattern) {
    case 0:
      InitializeP0(m_u, m_v);
      break;
    case 1:
      InitializeP1(m_u, m_v);
      break;
    default:
      break;
  }
  updateDpData();
}

void GrayScottSolver::solve(int n) {
  m_loopBreak = false;

  //    auto
  while (true) {
    if (m_loopBreak) break;
    for (int i = 0; i < n; ++i) {
      (m_enableParallel) ? TimeStep(m_u, m_v, m_gsParams)
                         : TimeStepSeq(m_u, m_v, m_gsParams);
      if (m_loopBreak) break;
    }
    updateDpData();
    emit dataReady();
  }
}
