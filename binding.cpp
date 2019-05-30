#include <pybind11/pybind11.h>
#include "AntColonyBase.h"

namespace py = pybind11;

PYBIND11_MODULE(example, m) {
  py::class_<AntColonyBase>(m, "antcolony")
      .def(py::init<const std::string &>())
      .def(py::init<const std::string &, double>())
      .def(py::init<const std::string &, double, double>())
      .def(py::init<const std::string &, double, double, double>())
      .def(py::init<const std::string &, double, double, double, double>())
      .def(py::init<const std::string &, double, double, double, double, unsigned>())
      .def("get_path", &AntColonyBase::get_path)
      .def("get_mintour_each", &AntColonyBase::get_mintour_each)
      .def("get_mintour_global", &AntColonyBase::get_mintour_global)
      .def("calcTSP", &AntColonyBase::calcTSP)
      .def("recalcTSP", &AntColonyBase::recalcTSP);
}
