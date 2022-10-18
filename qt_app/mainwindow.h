#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "grayscott_solver.h"
#include "colormapping.h"

QT_BEGIN_NAMESPACE
namespace Ui {
class MainWindow;
}
QT_END_NAMESPACE


class ImageViewer;

struct GSPattern {
    double f, k;
    int pattern;
};

class MainWindow : public QMainWindow {
  Q_OBJECT

  enum { kDp2D, kDp3D } m_active_dp = kDp2D;

 public:
  MainWindow(QWidget *parent = nullptr);
  ~MainWindow();

 private slots:
  void on_push_button_run__clicked();
  void on_push_button_init__clicked();
  void on_combo_box_uov__activated(int index);
  void on_combo_box_cmaps__currentTextChanged(const QString &arg1);
  void on_checkBoxParallel_clicked(bool checked);
  void on_comboBox_currentTextChanged(const QString &arg1);
  void on_comboBoxMethod_currentIndexChanged(int index);
  void on_pushButtonSaveData_clicked();
  void on_pushButtonSaveImage_clicked();

  private:
  Ui::MainWindow *ui;
  GrayScottSolver m_grayscottSolver;
  ImageViewer *m_viewer;
  void SetUp();
  void UpdateParams(bool reset);
  QMap<QString, ColorMapper::GradientPreset> m_name2gp;
  QMap<QString, GSPattern> m_gsPatterns;
  void enableGSParamEdit(bool enable);

 signals:
  void StartSolverSymRK2(int);
  void StartSolverEuler(int);
  void StopSolver();

 public slots:
  void UpdatePlot();
};
#endif  // MAINWINDOW_H
