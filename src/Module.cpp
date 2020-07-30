#include "Exceptions.h"
#include "Integrators.h"
#include "kde.h"
#include "DihedralEntropy.h"
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

namespace py = boost::python;

/****
 * The python modules
 */


BOOST_PYTHON_MODULE(kde) {
  py::class_<Gce>("Kde", py::init<py::list &, int>())
    .def("calculate", &Gce::calculate)
    .def("getResult", &Gce::getDensityEstimation)
    .def("integrate", &Gce::integrate_p)
    .def(py::init<py::list &, py::list &, int>())
    ;
  py::class_<std::vector<double> >("DoubleVec")
    .def(py::vector_indexing_suite<std::vector<double>>() )
    ;
  py::class_<DihedralEntropy>("DihedralEntropy", py::init<py::list &, int>())
    .def(py::init<py::list &>())
    .def(py::init<py::list &, py::str &>())
    .def(py::init<py::list &, int, py::str &>())
    .def("getResult", &DihedralEntropy::getEntropy)
    ;
}