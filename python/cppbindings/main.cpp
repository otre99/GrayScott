#include "gs_solver.h"
#include <algorithm>
#include <iostream>
#include <numeric>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

using namespace std;
namespace py = pybind11;

typedef py::array_t<double, py::array::c_style | py::array::forcecast>
    NumpyArrF64;

bool euler_step(NumpyArrF64 u, NumpyArrF64 v, NumpyArrF64 uOut,
                NumpyArrF64 vOut, double f, double k, double mu, double mv, int nsteps) {
  static GrayScottParams gs;

  gs.f = f;
  gs.k = k;
  gs.mu = mu;
  gs.mv = mv;
  gs.nrows = u.shape(0);
  gs.ncols = u.shape(1);

  if (gs.index.size() != gs.nrows) {
    gs.index.resize(gs.nrows);
    std::iota(gs.index.begin(), gs.index.end(), 0);
  }

  auto uPtr = static_cast<double *>(u.request().ptr);
  auto vPtr = static_cast<double *>(v.request().ptr);
  auto uOutPtr = static_cast<double *>(uOut.request().ptr);
  auto vOutPtr = static_cast<double *>(vOut.request().ptr);

  for (int i=0; i<nsteps; ++i){    
    ( i%2 == 0) ? TimeStepEuler(uPtr, vPtr, uOutPtr, vOutPtr, gs) : TimeStepEuler(uOutPtr, vOutPtr, uPtr, vPtr, gs);
  }
  return nsteps%2!=0;
}

PYBIND11_MODULE(libgray_scott_cpp, m) {
  m.doc() = "GrayScott pybind11 module";
  m.def("euler_step", &euler_step, py::call_guard<py::gil_scoped_release>(),
        py::arg("u"), py::arg("v"), py::arg("uOut"), py::arg("vOut"),
        py::arg("f"), py::arg("k"), py::arg("mu"), py::arg("mv"), py::arg("nsteps"));
}
