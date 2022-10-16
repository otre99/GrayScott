#include "mainwindow.h"
#include <QApplication>
#include "linalg.h"
#include <iostream>

//void PrintVector(double *v, int n, bool add_endl) {
//  std::cout << v[0];
//  for (int i = 1; i < n; ++i) {
//    std::cout << ' ' << v[i];
//  }
//  if (add_endl) std::cout << std::endl;
//}

int main(int argc, char *argv[]) {
  QApplication a(argc, argv);
  MainWindow w;


//  double la[4] = {1,1,1,3};
//  double lb[4] = {6,7,8,9};
//  double lc[4] = {2,2,1,2};

//  double bb[4] = {3,4,5,1};
//  double xx[4] = {3,4,5,1};

//  CyclicPackV1(la, lb, lc, 4, bb, xx);
//  PrintVector(xx,4);

//  CyclicPackV1(la, lb, lc, 4, bb, bb);
//  PrintVector(bb,4);

  w.showMaximized();
  return a.exec();
}
