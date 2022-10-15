#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>

#include "grayscott.h"
#include "simplecmdparser.hpp"
#include "solver.h"
using namespace std;

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

int main(int argc, char *argv[]) {
  SimpleCmdParser arg_parser;
  arg_parser.SetDefaultValue("--dt", "1.0", "Time step [seconds]");
  arg_parser.SetDefaultValue("--n", "512", "Resolution");
  arg_parser.SetDefaultValue("--nsteps", "10000", "No. steps");
  arg_parser.SetDefaultValue("--f", "0.055", "f-param");
  arg_parser.SetDefaultValue("--k", "0.062", "k-param");
  arg_parser.SetDefaultValue("--sv", "-1",
                             "Save partial solution every `sv` steps");
  arg_parser.SetDefaultValue("--mu", "0.16 ", "U difussion coef");
  arg_parser.SetDefaultValue("--p", "1", "Initial condition pattern");
  arg_parser.SetDefaultValue("--mv", "0.08", "V difussion coef");
  arg_parser.SetDefaultValue("--o", "./ ", "Output folder");

  if (!arg_parser.Parse(argc, argv)) {
    arg_parser.PrintHelp();
    return -1;
  }

  GrayScottParams gsParams;
  gsParams.Init(arg_parser.Get<double>("--k"), arg_parser.Get<double>("--f"),
                arg_parser.Get<double>("--mu"), arg_parser.Get<double>("--mv"),
                arg_parser.Get<double>("--dt"), arg_parser.Get<size_t>("--n"));
  const size_t NSTEPS = arg_parser.Get<size_t>("--nsteps");
  const int SV = arg_parser.Get<int>("--sv");
  const string OUTPUT_FOLDER = arg_parser.Get<string>("--o");
  cout << "Params: k=" << gsParams.k << " f=" << gsParams.f
       << " mu=" << gsParams.mu << " mv=" << gsParams.mv
       << " dt=" << gsParams.dt << " n=" << gsParams.n << " nsteps=" << NSTEPS
       << " sv=" << SV << endl;

  vector<double> u(gsParams.n * gsParams.n), v(gsParams.n * gsParams.n);
  const int pattern = arg_parser.Get<int>("--p");
  switch (pattern) {
    case 1:
      InitializeP1(u, v);
      break;
    case 2:
      InitializeP2(u, v);
      break;
    default:
      cerr << "Wrong initial pattern option. Valid options: [1, 2]" << endl;
      return -1;
  }
  const size_t dp = std::max(NSTEPS / 10, size_t(1));
  char fname[128];
  sprintf(fname, "ITER%08d", 0);
  SaveToFile(OUTPUT_FOLDER + "/U" + fname, u);
  SaveToFile(OUTPUT_FOLDER + "/V" + fname, v);
  for (size_t i = 0; i < NSTEPS; ++i) {
    if (SV > 0 && i % SV == 0) {
      cout << "\nSaving partial solution for step " << i << endl;
      sprintf(fname, "ITER%08zu", i);
      SaveToFile(OUTPUT_FOLDER + "/U" + fname, u);
      SaveToFile(OUTPUT_FOLDER + "/V" + fname, v);
    }
    TimeStep(u, v, gsParams);
    if ((i + 1) % dp == 0) {
      cout << "\rIter " << i + 1;
      cout.flush();
    }
  }
  SaveToFile(OUTPUT_FOLDER + "/UFINAL.dat", u);
  SaveToFile(OUTPUT_FOLDER + "/VFINAL.dat", v);
  cout << endl;
  return 0;
}
