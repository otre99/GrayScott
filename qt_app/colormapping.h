#ifndef COLORMAPPING_H
#define COLORMAPPING_H
#include <QColor>
#include <QMap>
#include <QRgb>
#include <QVector>

class ColorMapper {
 public:
  enum ColorInterpolation { ciRGB, ciHSV };
  enum GradientPreset {
    gpGrayscale  ///< Continuous lightness from black to white (suited for
                 ///< non-biased data representation)
    ,
    gpHot  ///< Continuous lightness from black over firey colors to white
           ///< (suited for non-biased data representation)
    ,
    gpCold  ///< Continuous lightness from black over icey colors to white
            ///< (suited for non-biased data representation)
    ,
    gpNight  ///< Continuous lightness from black over weak blueish colors to
             ///< white (suited for non-biased data representation)
    ,
    gpCandy  ///< Blue over pink to white
    ,
    gpGeography  ///< Colors suitable to represent different elevations on
                 ///< geographical maps
    ,
    gpIon  ///< Half hue spectrum from black over purple to blue and finally
           ///< green (creates banding illusion but allows more precise
           ///< magnitude estimates)
    ,
    gpThermal  ///< Colors suitable for thermal imaging, ranging from dark blue
               ///< over purple to orange, yellow and white
    ,
    gpPolar  ///< Colors suitable to emphasize polarity around the center, with
             ///< blue for negative, black in the middle and red for positive
             ///< values
    ,
    gpSpectrum  ///< An approximation of the visible light spectrum (creates
                ///< banding illusion but allows more precise magnitude
                ///< estimates)
    ,
    gpJet  ///< Hue variation similar to a spectrum, often used in numerical
           ///< visualization (creates banding illusion but allows more precise
           ///< magnitude estimates)
    ,
    gpHues  ///< Full hue cycle, with highest and lowest color red (suitable for
            ///< periodic data, such as angles and phases, see \ref setPeriodic)
  };

  ColorMapper();
  ColorMapper(GradientPreset preset);
  bool operator==(const ColorMapper &other) const;
  bool operator!=(const ColorMapper &other) const { return !(*this == other); }

  // getters:
  int levelCount() const { return mLevelCount; }
  QMap<double, QColor> colorStops() const { return mColorStops; }
  ColorInterpolation colorInterpolation() const { return mColorInterpolation; }
  bool periodic() const { return mPeriodic; }

  // setters:
  void setLevelCount(int n);
  void setColorStops(const QMap<double, QColor> &colorStops);
  void setColorStopAt(double position, const QColor &color);
  void setColorInterpolation(ColorInterpolation interpolation);
  void setPeriodic(bool enabled);

  // non-property methods:
  void colorize(const double *data, const double &lower, const double &upper,
                QRgb *scanLine, int n, int dataIndexFactor = 1,
                bool logarithmic = false);
  void colorize(const double *data, const unsigned char *alpha,
                const double &lower, const double &upper, QRgb *scanLine, int n,
                int dataIndexFactor = 1, bool logarithmic = false);
  QRgb color(double position, const double &lower, const double &upper,
             bool logarithmic = false);
  void loadPreset(GradientPreset preset);
  void clearColorStops();
  ColorMapper inverted() const;

 protected:
  // property members:
  int mLevelCount;
  QMap<double, QColor> mColorStops;
  ColorInterpolation mColorInterpolation;
  bool mPeriodic;

  // non-property members:
  QVector<QRgb>
      mColorBuffer;  // have colors premultiplied with alpha (for usage with
                     // QImage::Format_ARGB32_Premultiplied)
  bool mColorBufferInvalidated;

  // non-virtual methods:
  bool stopsUseAlpha() const;
  void updateColorBuffer();
};

#endif  // COLORMAPPING_H
