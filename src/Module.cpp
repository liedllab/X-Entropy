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
    "Parameters\n"
    "----------\n"
    "data : list(float)\n"
    "   The data points to bin and calculate the KDE of.\n"
    "bins : int\n"
    "   The number of bins. This should be to base 2, if not\n"
    "   it will be recalculated as base 2 value, by calculating\n"
    "   the next higher number that is base 2 (i.e., 200 will\n"
    "   become 256).\n"
    "\n"
    "Examples\n"
    "--------\n"
    ">>> data = np.random.normal(scale=0.001, size=600000)\n"
    ">>> gce = entropy.kde.GCE(list(data), -1)\n"
    ))
    .def("calculate", &PyGCE::calculate,
      py::arg("self"),
      "Start the calculation\n"
      "\n"
      "Examples\n"
      "--------\n"
      ">>> data = np.random.normal(scale=0.001, size=600000)\n"
      ">>> gce = entropy.kde.GCE(list(data), -1)\n"
      ">>> gce.calculate()\n")
    .def("getResult", &PyGCE::getDensityEstimation,
      py::arg("self"),
      "Returns the results of the kernel density estimation.\n"
      "\n"
      "Returns\n"
      "-------\n"
      "density : list(float)\n"
      "    The estimated density, using the GCE approach.\n")
    .def("getTStar", &PyGCE::getTStar)
    .def("integrate", &PyGCE::integrate,
      py::args("self", "integrator", "min", "max"),
      "Integrates the estimated density and so calculates the\n"
      "entropy of the underlying data. If the implemented\n"
      "integrators are not enough, one can also do this\n"
      "entirely in Python and does not necessarily need to\n"
      "use the C++ implementation.\n"
      "\n"
      "Parameters\n"
      "----------\n"
      "integrator : string\n"
      "    Either \"Simpson\" or \"Riemann\". This is the\n"
      "    Numerical integrator that will be used.\n"
      "min : float\n"
      "    The start of the integration.\n"
      "max : float\n"
      "    The end of the integration.\n"
      "\n"
      "Returns\n"
      "-------\n"
      "area : float\n"
      "    The area under the underlying histogram. If one\n"
      "    needs the estimated density.\n"
      "\n"
      "Examples\n"
      "--------\n"
      ">>> data = np.random.normal(scale=0.001, size=600000)\n"
      ">>> gce = entropy.kde.GCE(list(data), -1)\n"
      ">>> gce.calculate()\n"
      ">>> gce.integrate(\"Simpson\", 0.0, 360.0)\n"
      "0.2988366\n")
    .def(py::init<py::list &, py::list &, int>(
      py::args("self", "histogram", "axis", "nObservations"), 
      "Constructor of a GCE class object. Here, the parameters for\n"
      "the list histogram are all given. Therefore, no new histogram\n"
      "needs to be calculated by the code.\n"
      "\n"
      "Parameters\n"
      "----------\n"
      "histogram : list(float)\n"
      "   The already binned data points.\n"
      "axis : list(float)\n"
      "   The coordinates of the bins. These will be used for the\n"
      "   integration.\n"
      "nObservation : int\n"
      "   The number of observations in the dataset. Could in theory\n"
      "   be calculated using the binned data, but this is more\n"
      "   convenient.\n"
    ))
    ;
  /*py::class_<std::vector<double> >("DoubleVec")
    .def(py::vector_indexing_suite<std::vector<double>>() )
    ;*/
  py::class_<DihedralEntropy>("DihedralEntropy", py::init<py::list &, int>(
      py::args("angles", "bins"),
      "Constructor for the Dihedral Entropy class, be aware that this\n"
      "class calculates the dihedral entropy on construction, so no calls\n"
      "to calculate or similar are necessary. The underlying data should\n"
      "not be binned.\n"
      "\n"
      "Parameters\n"
      "----------\n"
      "angles : list(float)\n"
      "    The angles for which to calculate the dihedral entropy. The\n"
      "    angles must be given in degrees.\n"
      "bins: int\n"
      "    The number of bins for the calculation of the histogram.\n"
    ))
    .def(
      py::init<py::list &>(
        py::arg("angles"),
        "Constructor using only the angles, uses default bin number\n"
        "for the calculation, default bin number is 2<<13.\n"
        "\n"
        "Parameters\n"
        "----------\n"
        "angles : list(float)\n"
        "    The angles for which to calculate the dihedral entropy. The\n"
        "    angles must be given in degrees.\n"
      )
    )
    .def(py::init<py::list &, py::str &>())
    .def(py::init<py::list &, int, py::str &>())
    .def("getResult", &DihedralEntropy::getEntropy, py::return_value_policy<py::return_by_value>())
    ;
}