#include "Exceptions.h"
#include "Integrators.h"
#include "PyGCE.h"
#include "DihedralEntropy.h"
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/args.hpp>
#include <boost/python/module.hpp>



namespace py = boost::python;

/****
 * The python modules
 */


BOOST_PYTHON_MODULE(kde) {
  py::class_<PyGCE>("GCE", "Docstring for the class", py::init<py::list &, int>(
    py::args("self", "data", "bins"),
    "Constructor of a GCE class object. Uses the data and bins\n"
    "to bin the data according to the bin number. The bins will\n"
    "spread from the minimum to the maximum found in data.\n"
    "\n"
    "Parameters:\n"
    "-----------\n"
    "data: list(float)\n"
    "   The data points to bin and calculate the KDE of.\n"
    "bins: int\n"
    "   The number of bins. This should be to base 2, if not\n"
    "   it will be recalculated as base 2 value, by calculating\n"
    "   the next higher number that is base 2 (i.e., 200 will\n"
    "   become 256).\n"
    "\n"
    "Examples:\n"
    "--------\n"
    ">>> data = np.random.normal(scale=0.001, size=600000)\n"
    ">>> gce = entropy.kde.GCE(list(data), -1)\n"
    ))
    .def("calculate", &PyGCE::calculate, "Test Docstring")
    .def("getResult", &PyGCE::getDensityEstimation)
    .def("integrate", &PyGCE::integrate)
    .def(py::init<py::list &, py::list &, int>(
      py::args("self", "histogram", "axis", "nObservations"), 
      "Constructor of a GCE class object. Here, the parameters for\n"
      "the list histogram are all given. Therefore, no new histogram\n"
      "needs to be calculated by the code.\n"
      "\n"
      "Parameters\n"
      "----------\n"
      "histogram: list(float)\n"
      "   The already binned data points.\n"
      "axis: list(float)\n"
      "   The coordinates of the bins. These will be used for the\n"
      "   integration.\n"
      "nObservation: int\n"
      "   The number of observations in the dataset. Could in theory\n"
      "   be calculated using the binned data, but this is more\n"
      "   convenient.\n"
    ))
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