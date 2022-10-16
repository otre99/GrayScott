#ifndef GRAYSCOTT_UTILS_H
#define GRAYSCOTT_UTILS_H
#include <vector>
#include <cmath>
#include <string>

class MatrixWrapper {
 public:
  explicit MatrixWrapper(std::vector<double> &data);
  double &operator()(int i, int j) const;
  double *GetRow(int row);

 private:
  std::vector<double> *data_;
  int n_;
};

void InitializeP0(std::vector<double> &u, std::vector<double> &v);
void InitializeP1(std::vector<double> &u, std::vector<double> &v);

void SaveToFile(const std::string &fname, const std::vector<double> &values);
bool LoadFromFile(const std::string &fname, std::vector<double> &values);

#endif // GRAYSCOTT_H
