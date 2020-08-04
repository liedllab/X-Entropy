#include "Exceptions.h"
#include "Integrators.h"
#include "PyGCE.h"
#include "DihedralEntropy.h"
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

namespace py = boost::python;

/****
 * The python modules
 */


BOOST_PYTHON_MODULE(kde) {
  py::class_<PyGCE>("Kde", py::init<py::list &, int>())
    .def("calculate", &PyGCE::calculate)
    .def("getResult", &PyGCE::getDensityEstimation)
    .def("integrate", &PyGCE::integrate)
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