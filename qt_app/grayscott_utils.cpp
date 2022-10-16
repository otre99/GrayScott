#include "grayscott_utils.h"
#include <random>
#include <algorithm>
#include <fstream>
#include <filesystem>

using namespace std;

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

double *MatrixWrapper::GetRow(int row){
    return data_->data()+row*n_;
}



void RandomInit(vector<double> &u, vector<double> &v,
                uniform_real_distribution<double> &dist, mt19937 &mt) {
  generate(u.begin(), u.end(), [&]() { return 1 + dist(mt); });
  generate(v.begin(), v.end(), [&]() { return dist(mt); });
}

void AddRect(MatrixWrapper &&uu, MatrixWrapper &&vv, int i, int j, int r,
             uniform_real_distribution<double> &dist, mt19937 &mt) {
  for (int ii = i - r; ii < i + r; ++ii) {
    for (int jj = j - r; jj < j + r; ++jj) {
      uu(ii, jj) = 0.5 + dist(mt);
      vv(ii, jj) = 0.25 + dist(mt);
    }
  }
}

constexpr double delta1 = 0.0;
constexpr double delta2 = 0.05;
void InitializeP0(vector<double> &u, vector<double> &v) {
  random_device rd;
  mt19937 mt(rd());
  uniform_real_distribution<double> dist(delta1, delta2);
  RandomInit(u, v, dist, mt);

  size_t dim = sqrt(u.size());
  const size_t radius = dim / (dim / 16);
  AddRect(MatrixWrapper(u), MatrixWrapper(v), dim / 2, dim / 2, radius, dist,
          mt);
}

void InitializeP1(vector<double> &u, vector<double> &v) {
  random_device rd;
  mt19937 mt(rd());
  uniform_real_distribution<double> dist(delta1, delta2);
  RandomInit(u, v, dist, mt);

  size_t dim = sqrt(u.size());
  const size_t radius = dim / 4;

  for (size_t i = 0; i < dim; i += radius) {
    for (size_t j = 0; j < dim; j += radius) {
      AddRect(MatrixWrapper(u), MatrixWrapper(v), i + radius / 2,
              j + radius / 2, 8, dist, mt);
    }
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
