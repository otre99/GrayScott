#include "grayscott_utils.h"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <random>

using namespace std;

void PrintMatrix(const double *data, int n, int m, bool new_line) {
  size_t idx = 0;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < m; ++j) {
      std::printf("%21.17f     ", data[idx++]);
    }
    if (i < n - 1) std::printf("\n");
  }
  if (new_line) std::printf("\n");
}

MatrixWrapper::MatrixWrapper(std::vector<double> &data) : data_(&data) {
  n_ = static_cast<int>(std::sqrt(data_->size()));
}
double &MatrixWrapper::operator()(int i, int j) const {
  i %= n_;
  j %= n_;
  if (i < 0) i += n_;
  if (j < 0) j += n_;
  int pos = i * n_ + j;
  return data_->operator[](pos);
}

double *MatrixWrapper::GetRow(int row) { return data_->data() + row * n_; }

// void RandomInit(vector<double> &u, vector<double> &v,
//                uniform_real_distribution<double> &dist, mt19937 &mt) {
//  generate(u.begin(), u.end(), [&]() { return 1 /*+ dist(mt)*/; });
//  generate(v.begin(), v.end(), [&]() { return 0.0/*+dist(mt)*/; });
//}

constexpr double delta1 = 0.0;
constexpr double delta2 = 0.05;
const static double PI = std::acos(-1);
void AddNoise(vector<double> &u, vector<double> &v) {
  random_device rd;
  mt19937 mt(rd());
  uniform_real_distribution<double> dist(delta1, delta2);
  for_each(u.begin(), u.end(), [&](double &v) { v += dist(mt); });
  for_each(v.begin(), v.end(), [&](double &v) { v += dist(mt); });
}

void AddRect(MatrixWrapper &&uu, MatrixWrapper &&vv, int i, int j, int r) {
  for (int ii = i - r; ii < i + r; ++ii) {
    for (int jj = j - r; jj < j + r; ++jj) {
      uu(ii, jj) = 0.5;
      vv(ii, jj) = 0.25;
    }
  }
}

void Initialize(std::vector<double> &u, std::vector<double> &v, int pattern,
                bool with_noise) {
  switch (pattern) {
    case 0:
      InitializeP0(u, v, with_noise);
      break;
    case 1:
      InitializeP1(u, v, with_noise);
      break;
    case 2:
      InitializeP2(u, v, with_noise);
      break;
    default:
      throw("Wrong initial pattern option. Valid options: [0, 1, 2]");
  }
}

void InitializeP0(vector<double> &u, vector<double> &v, bool witth_noise) {
  std::fill(u.begin(), u.end(), 1.0);
  std::fill(v.begin(), v.end(), 0.0);
  size_t dim = sqrt(u.size());
  const size_t radius = dim / (dim / 16);
  AddRect(MatrixWrapper(u), MatrixWrapper(v), dim / 2, dim / 2, radius);
  if (witth_noise) {
    AddNoise(u, v);
  }
}

void InitializeP1(vector<double> &u, vector<double> &v, bool witth_noise) {
  std::fill(u.begin(), u.end(), 1.0);
  std::fill(v.begin(), v.end(), 0.0);
  size_t dim = sqrt(u.size());
  const size_t radius = dim / 4;
  for (size_t i = 0; i < dim; i += radius) {
    for (size_t j = 0; j < dim; j += radius) {
      AddRect(MatrixWrapper(u), MatrixWrapper(v), i + radius / 2,
              j + radius / 2, 8);
    }
  }
  if (witth_noise) {
    AddNoise(u, v);
  }
}

void InitializeP2(vector<double> &u, vector<double> &v, bool witth_noise) {
  size_t dim = sqrt(u.size());
  auto f = [](const double &x, const double &y) {
    bool okx = (x >= 1.0) && (x <= 1.5);
    bool oky = (y >= 1.0) && (y <= 1.5);
    if (okx && oky) {
      double sx = sin(4.0 * PI * x);
      double sy = sin(4.0 * PI * y);
      return 0.25 * sx * sx * sy * sy;
    }
    return 0.0;
  };

  const double pixelToVal = 2.5 / dim;
  MatrixWrapper mU(u), mV(v);
  for (size_t i = 0; i < dim; ++i) {
    const double y = pixelToVal * i;
    for (size_t j = 0; j < dim; ++j) {
      double Val = mU(i, j) = f(j * pixelToVal, y);
      mV(i, j) = 1.0 - 2 * Val;
    }
  }
  if (witth_noise) {
    AddNoise(u, v);
  }
}

void SaveToFile(const string &fname, const vector<double> &values) {
  ofstream ofile(fname, ios::binary | ios::out);
  ofile.write(reinterpret_cast<const char *>(values.data()),
              values.size() * sizeof(double));
}

bool LoadFromFile(const string &fname, vector<double> &values) {
  ifstream ifile(fname, ios::binary | ios::in);
  if (ifile) {
    const long NBYTES = std::filesystem::file_size(fname);
    values.resize(NBYTES / sizeof(double));
    return NBYTES ==
           ifile.read(reinterpret_cast<char *>(values.data()), NBYTES).gcount();
  }
  return false;
}
