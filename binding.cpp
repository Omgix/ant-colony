#include <pybind11/pybind11.h>
#include "AntColonyBase.h"

namespace py = pybind11;
using namespace pybind11::literals;

PYBIND11_MODULE(antcolony, m) {
  py::class_<AntColonyBase>(m, "AntColony")
      .def(py::init<const std::string &, double, double, double, double, unsigned>(),
          "filename"_a, "alpha"_a=15.0, "beta"_a = 20.0, "rho"_a = 0.1, "colony_eff"_a = 1.0, "maxiter"_a = 500)
      .def("get_path", &AntColonyBase::get_path)
      .def("get_mintour_each", &AntColonyBase::get_mintour_each)
      .def("get_mintour_global", &AntColonyBase::get_mintour_global)
      .def("calcTSP", &AntColonyBase::calcTSP)
      .def("recalcTSP", &AntColonyBase::recalcTSP);
}
