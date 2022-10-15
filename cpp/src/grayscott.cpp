#include "grayscott.h"

#include <algorithm>
#include <random>
#include <vector>

using namespace std;
class MatrixWrapper {
 public:
  explicit MatrixWrapper(vector<double> &data) : data_(&data) {
    n_ = static_cast<int>(std::sqrt(data_->size()));
  }
  double &operator()(int i, int j) const {
    i %= n_;
    j %= n_;
    if (i < 0) i += n_;
    if (j < 0) j += n_;
    int pos = i * n_ + j;
    return data_->operator[](pos);
  }

 private:
  vector<double> *data_;
  int n_;
};

void RandomInit(vector<double> &u, vector<double> &v,
                uniform_real_distribution<double> &dist, mt19937 &mt) {
  generate(u.begin(), u.end(), [&]() { return 1 + dist(mt); });
  generate(v.begin(), v.end(), [&]() { return dist(mt); });
}

auto AddRect(MatrixWrapper &&uu, MatrixWrapper &&vv, int i, int j, int r,
             uniform_real_distribution<double> &dist, mt19937 &mt) {
  for (int ii = i - r / 2; ii < i + r / 2; ++ii) {
    for (int jj = j - r / 2; jj < j + r / 2; ++jj) {
      uu(ii, jj) = 0.5 + dist(mt);
      vv(ii, jj) = 0.25 + dist(mt);
    }
  }
}

constexpr double delta = 0.05;
void InitializeP1(vector<double> &u, vector<double> &v) {
  random_device rd;
  mt19937 mt(rd());
  uniform_real_distribution<double> dist(-delta, delta);
  RandomInit(u, v, dist, mt);

  size_t dim = sqrt(u.size());
  const size_t radius = dim / 10;
  AddRect(MatrixWrapper(u), MatrixWrapper(v), dim / 2, dim / 2, radius, dist,
          mt);
}

void InitializeP2(vector<double> &u, vector<double> &v) {
  random_device rd;
  mt19937 mt(rd());
  uniform_real_distribution<double> dist(-delta, delta);
  RandomInit(u, v, dist, mt);

  size_t dim = sqrt(u.size());
  const size_t radius = dim / 10;

  for (size_t i = 0; i < dim; i += radius) {
    for (size_t j = 0; j < dim; j += radius) {
      AddRect(MatrixWrapper(u), MatrixWrapper(v), i + radius / 2,
              j + radius / 2, radius / 3, dist, mt);
    }
  }
}
