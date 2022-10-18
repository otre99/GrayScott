#include "mainwindow.h"

#include <QDebug>
#include <QFileDialog>
#include <QMap>
#include <QPainter>
#include <cmath>

#include "./ui_mainwindow.h"
#include "grayscott_utils.h"
#include "imageviewer.h"

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent), ui(new Ui::MainWindow) {
  ui->setupUi(this);
  SetUp();
}

MainWindow::~MainWindow() { delete ui; }

void MainWindow::on_push_button_run__clicked() {
  if (!m_grayscottSolver.isRunning()) {
    UpdateParams(false);

    if (ui->comboBoxMethod->currentIndex() == 0) {
      emit StartSolverEuler(ui->spin_box_iters_->value());
    } else {
      emit StartSolverSymRK2(ui->spin_box_iters_->value());
    }
    ui->push_button_run_->setText("Stop");
    ui->group_box_params_->setEnabled(false);
  } else {
    m_grayscottSolver.stop();
    ui->push_button_run_->setText("Run");
    ui->group_box_params_->setEnabled(true);
  }
}

void MainWindow::on_push_button_init__clicked() {
  UpdateParams(true);
  m_grayscottSolver.initialize(ui->comboBoxPattern->currentIndex(),
                               ui->checkBoxNoise->isChecked());
  UpdatePlot();
}

void MainWindow::on_combo_box_uov__activated(int index) {
  m_grayscottSolver.selectUOrV(!index);
  if (!m_grayscottSolver.isRunning()) {
    m_grayscottSolver.updateDpData();
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

  //
  m_gsPatterns["type-alpha-01"] = {0.010, 0.047, 0};
  m_gsPatterns["type-alpha-02"] = {0.014, 0.053, 0};
  m_gsPatterns["type-beta-01"] = {0.014, 0.039, 0};
  m_gsPatterns["type-beta-02"] = {0.026, 0.051, 0};
  m_gsPatterns["type-gamma-01"] = {0.022, 0.051, 0};
  m_gsPatterns["type-gamma-02"] = {0.026, 0.055, 0};
  m_gsPatterns["type-delta-01"] = {0.030, 0.055, 0};
  m_gsPatterns["type-delta-02"] = {0.042, 0.059, 0};
  m_gsPatterns["type-epsilon-01"] = {0.018, 0.055, 0};
  m_gsPatterns["type-epsilon-02"] = {0.022, 0.059, 0};
  m_gsPatterns["type-zeta-01"] = {0.022, 0.061, 0};
  m_gsPatterns["type-zeta-02"] = {0.026, 0.059, 0};
  m_gsPatterns["type-eta-01"] = {0.034, 0.063, 0};
  m_gsPatterns["type-theta-01"] = {0.030, 0.057, 0};
  m_gsPatterns["type-theta-02"] = {0.038, 0.061, 0};
  m_gsPatterns["type-iota-01"] = {0.046, 0.0594, 0};
  m_gsPatterns["type-kappa-01"] = {0.050, 0.063, 0};
  m_gsPatterns["type-kappa-02"] = {0.058, 0.063, 0};
  m_gsPatterns["type-lambda-01"] = {0.026, 0.061, 0};
  m_gsPatterns["type-lambda-02"] = {0.034, 0.065, 0};
  m_gsPatterns["type-mu-01"] = {0.046, 0.065, 0};
  m_gsPatterns["type-mu-02"] = {0.058, 0.065, 0};
  m_gsPatterns["type-xi-01"] = {0.014, 0.047, 0};
  m_gsPatterns["type-pi-01"] = {0.062, 0.061, 0};

  connect(&m_grayscottSolver, &GrayScottSolver::dataReady, this,
          &MainWindow::UpdatePlot);
  connect(this, &MainWindow::StartSolverSymRK2, &m_grayscottSolver,
          &GrayScottSolver::solveSymRK2);
  connect(this, &MainWindow::StartSolverEuler, &m_grayscottSolver,
          &GrayScottSolver::solveEuler);

  setCentralWidget(m_viewer = new ImageViewer);
  UpdateParams(true);
}

void MainWindow::UpdateParams(bool reset) {
  double du = ui->le_du->text().toDouble();
  double dv = ui->le_dv->text().toDouble();
  double f = ui->le_f->text().toDouble();
  double k = ui->le_k->text().toDouble();
  double dt = ui->line_dt->text().toDouble();
  int n = ui->le_res->text().toDouble();

  m_grayscottSolver.initParams(n, du, dv, f, k, dt, reset,
                               ui->comboBoxMethod->currentIndex());
}

void MainWindow::UpdatePlot() {
  m_viewer->updateData(m_grayscottSolver.getData());
  ui->labelCurrTime->setNum(m_grayscottSolver.currTime());
}

void MainWindow::on_combo_box_cmaps__currentTextChanged(const QString &arg1) {
  const auto p = m_name2gp[arg1];
  m_viewer->setColormap(p);
  if (!m_grayscottSolver.isRunning()) {
    m_viewer->colorize();
    m_viewer->forceToRepaint();
  }
}
void MainWindow::on_checkBoxParallel_clicked(bool checked) {
  m_grayscottSolver.enableParallel(checked);
}

void MainWindow::on_comboBox_currentTextChanged(const QString &arg1) {
  if (arg1 == "custom") {
    enableGSParamEdit(true);
    return;
  }
  const auto gs = m_gsPatterns[arg1];
  ui->comboBoxPattern->setCurrentIndex(gs.pattern);
  ui->le_f->setText(QString::number(gs.f));
  ui->le_k->setText(QString::number(gs.k));
  enableGSParamEdit(false);
  on_push_button_init__clicked();
}

void MainWindow::enableGSParamEdit(bool enable) {
  ui->comboBoxPattern->setEnabled(enable);
  ui->le_du->setEnabled(enable);
  ui->le_dv->setEnabled(enable);
  ui->le_f->setEnabled(enable);
  ui->le_k->setEnabled(enable);
}

void MainWindow::on_comboBoxMethod_currentIndexChanged(int index) {
  if (index == 0) {
    ui->line_dt->setText("1.0");
    ui->line_dt->setEnabled(false);
  } else {
    ui->line_dt->setEnabled(true);
  }
}

void MainWindow::on_pushButtonSaveData_clicked() {
  QString fname = QFileDialog::getSaveFileName(this, "Save data");
  if (!fname.isEmpty()) {
    SaveToFile(fname.toStdString(), m_viewer->getData());
  }
}

void MainWindow::on_pushButtonSaveImage_clicked() {
  QString fname = QFileDialog::getSaveFileName(this, "Save data");
  if (!fname.isEmpty()) {
    if (m_viewer->imagePtr()) {
      m_viewer->imagePtr()->save(fname);
    }
  }
}
