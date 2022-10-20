#ifndef GRAYSCOTT_UTILS_H
#define GRAYSCOTT_UTILS_H
#include <cmath>
#include <string>
#include <vector>

class MatrixWrapper {
 public:
  explicit MatrixWrapper(std::vector<double> &data);
  double &operator()(int i, int j) const;
  double *GetRow(int row);

 private:
  std::vector<double> *data_;
  int n_;
};

void PrintMatrix(const double* data, int n, int m=1, bool new_line=true);
void Initialize(std::vector<double> &u, std::vector<double> &v, int pattern,
                bool with_noise);
void InitializeP0(std::vector<double> &u, std::vector<double> &v,
                  bool witth_noise);
void InitializeP1(std::vector<double> &u, std::vector<double> &v,
                  bool witth_noise);
void InitializeP2(std::vector<double> &u, std::vector<double> &v,
                  bool witth_noise);

void SaveToFile(const std::string &fname, const std::vector<double> &values);
bool LoadFromFile(const std::string &fname, std::vector<double> &values);

#endif  // GRAYSCOTT_H
