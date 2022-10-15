#include "mainwindow.h"

#include <QDebug>
#include <QMap>
#include <QPainter>
#include <cmath>

#include "./ui_mainwindow.h"
#include "imageviewer.h"

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent), ui(new Ui::MainWindow) {
  ui->setupUi(this);
  SetUp();
}

MainWindow::~MainWindow() { delete ui; }

void MainWindow::on_push_button_run__clicked() {
  if (!m_gray_scott_solver.isRunning()) {
    emit StartSolver(ui->spin_box_iters_->value());
    ui->push_button_run_->setText("Stop");
    ui->group_box_params_->setEnabled(false);
  } else {
    m_gray_scott_solver.stop();
    ui->push_button_run_->setText("Run");
    ui->group_box_params_->setEnabled(true);
  }
}

void MainWindow::on_push_button_init__clicked() {
  UpdateParams();
  m_gray_scott_solver.initialize(ui->comboBoxPattern->currentIndex());
  UpdatePlot();
}

void MainWindow::on_combo_box_uov__activated(int index) {
  m_gray_scott_solver.selectUOrV(!index);
  if (!m_gray_scott_solver.isRunning()) {
    m_gray_scott_solver.updateDpData();
    UpdatePlot();
  }
}

void MainWindow::SetUp() {
  m_name2gp["gpGrayscale"] = ColorMapper::gpGrayscale;
  m_name2gp["gpHot"] = ColorMapper::gpHot;
  m_name2gp["gpCold"] = ColorMapper::gpCold;
  m_name2gp["gpNight"] = ColorMapper::gpNight;
  m_name2gp["gpCandy"] = ColorMapper::gpCandy;
  m_name2gp["gpGeography"] = ColorMapper::gpGeography;
  m_name2gp["gpIon"] = ColorMapper::gpIon;
  m_name2gp["gpThermal"] = ColorMapper::gpThermal;
  m_name2gp["gpPolar"] = ColorMapper::gpPolar;
  m_name2gp["gpSpectrum"] = ColorMapper::gpSpectrum;
  m_name2gp["gpJet"] = ColorMapper::gpJet;
  m_name2gp["gpHues"] = ColorMapper::gpHues;

  connect(&m_gray_scott_solver, &GrayScottSolver::dataReady, this,
          &MainWindow::UpdatePlot);
  connect(this, &MainWindow::StartSolver, &m_gray_scott_solver,
          &GrayScottSolver::solve);

  setCentralWidget(m_viewer = new ImageViewer);
  UpdateParams();
}

void MainWindow::UpdateParams() {
  double du = ui->le_du->text().toDouble();
  double dv = ui->le_dv->text().toDouble();
  double f = ui->le_f->text().toDouble();
  double k = ui->le_k->text().toDouble();
  double dt = ui->line_dt->text().toDouble();
  int n = ui->le_res->text().toDouble();
  m_gray_scott_solver.initParams(n, du, dv, f, k, dt);
}

void MainWindow::UpdatePlot() {
  m_viewer->updateData(m_gray_scott_solver.getData());
}

void MainWindow::on_combo_box_cmaps__currentTextChanged(const QString &arg1) {
  const auto p = m_name2gp[arg1];
  m_viewer->setColormap(p);
  if (!m_gray_scott_solver.isRunning()) {
    m_viewer->colorize();
    m_viewer->forceToRepaint();
  }
}
void MainWindow::on_checkBoxParallel_clicked(bool checked) {
  m_gray_scott_solver.enableParallel(checked);
}
