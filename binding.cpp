#include <pybind11/pybind11.h>
#include "AntColonyBase.h"

namespace py = pybind11;

PYBIND11_MODULE(example, m) {
  py::class_<AntColonyBase>(m, "AntColonyBase")
      .def(py::init<const std::string &>())
      .def("get_path", &AntColonyBase::get_path)
      .def("get_mintour_each", &AntColonyBase::get_mintour_each)
      .def("get_mintour_global", &AntColonyBase::get_mintour_global)
      .def("calcTSP", &AntColonyBase::calcTSP)
      .def("recalcTSP", &AntColonyBase::recalcTSP);
}
