#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>

#include "grayscott_utils.h"
#include "simplecmdparser.hpp"
#include "solver.h"

using namespace std;

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
  arg_parser.SetDefaultValue("--pattern", "0",
                             "Initial condition pattern [0,1,2]");
  arg_parser.SetDefaultValue("--mv", "0.08", "V difussion coef");
  arg_parser.SetDefaultValue("--o", "./ ", "Output folder");
  arg_parser.SetDefaultValue("--method", "0", "Method. 0-Euler 1-SymRK2");
  arg_parser.SetDefaultValue("--noise", "1",
                             "Add noise to the initial condition");

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
  const int METHOD = arg_parser.Get<int>("--method");

  const string OUTPUT_FOLDER = arg_parser.Get<string>("--o");
  cout << "Params: k=" << gsParams.k << " f=" << gsParams.f
       << " mu=" << gsParams.mu << " mv=" << gsParams.mv
       << " dt=" << gsParams.dt << " n=" << gsParams.n << " nsteps=" << NSTEPS
       << " sv=" << SV << " method: " << ((METHOD == 0) ? "EULER" : "SymRK2")
       << endl;

  vector<double> u(gsParams.n * gsParams.n), v(gsParams.n * gsParams.n);
  vector<double> uTemp, vTemp;
  if (METHOD == 0) {
    uTemp.resize(u.size());
    vTemp.resize(v.size());
  }

  const int pattern = arg_parser.Get<int>("--pattern");
  const bool with_noise = arg_parser.Get<bool>("--noise");
  Initialize(u, v, pattern, with_noise);

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
    (METHOD == 1)
        ? TimeStepSymRK2(u, v, gsParams) /*TimeStepSymRK2(u, v, gsParams)*/
        : TimeStepEuler(
              u, v, uTemp, vTemp,
              gsParams); /*TimeStepEuler(u, v, uTemp, vTemp, gsParams)*/
    ;
    if ((i + 1) % dp == 0) {
      cout << "\rIter " << i + 1;
      cout.flush();
    }
  }
  SaveToFile(OUTPUT_FOLDER + "/UFINAL", u);
  SaveToFile(OUTPUT_FOLDER + "/VFINAL", v);
  cout << endl;
  return 0;
}
