#ifndef IMAGEVIEWER_H
#define IMAGEVIEWER_H

#include <QAbstractScrollArea>
#include <QPixmap>

#include "colormapping.h"
class QImage;

class ImageViewer : public QAbstractScrollArea {
  Q_OBJECT
 public:
  explicit ImageViewer(QWidget *parent = nullptr);
  void updateData(const std::vector<double> &data);
  void setColormap(ColorMapper::GradientPreset p);
  void forceToRepaint();
  void setUseLog(bool use_log) { m_useLog = use_log; }
  void colorize();
  void attachImagePtr(QImage *ptr);
  void init();
  const std::vector<double> &getData() const {return m_rawData;}
  QImage *imagePtr();
  double scale() const { return m_scaleFactor; }
  QRect canvasRect() const;
  QRect imageRect() const;

 protected:
  void paintEvent(QPaintEvent *) override;
  void resizeEvent(QResizeEvent *ev) override;
  void mouseMoveEvent(QMouseEvent *) override;
  void mousePressEvent(QMouseEvent *) override;
  void wheelEvent(QWheelEvent *event) override;
 signals:
  void pixelTrack(int x, int y, double f);

 public slots:
  void setXmov(int x);
  void setYmov(int y);
  void setScf(double);
  void fitWidth();
  void generateCache();
  void adjustAll();

 private:
  QImage *m_imagePtr;
  QImage m_imageData;
  int m_xMov, m_yMov;
  int m_cW, m_cH;
  int m_screenW, m_screenH;
  double m_scaleFactor;
  QPoint m_lastPt;
  void selectScf();
  QTimer *cacheTimer_;
  void queueGenerateCache();
  QPixmap m_cachedPixmap;
  ColorMapper m_colorMapper;
  std::vector<double> m_rawData;
  bool m_useLog{false};
};

#endif  // IMAGEVIEWER_H
