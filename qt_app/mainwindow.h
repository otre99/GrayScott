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

private:
  Ui::MainWindow *ui;
  GrayScottSolver m_gray_scott_solver;
  ImageViewer *m_viewer;
  void SetUp();
  void UpdateParams();
  QMap<QString, ColorMapper::GradientPreset> m_name2gp;

 signals:
  void StartSolver(int);
  void StopSolver();

 public slots:
  void UpdatePlot();
};
#endif  // MAINWINDOW_H
