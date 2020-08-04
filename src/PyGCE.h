#include "kde.h"

#include<vector>
#include <boost/python.hpp>


namespace py = boost::python;


class PyGCE {
private:
  Gce kde;
public:
  PyGCE(py::list &l, int n)
  : kde{ 
    std::vector<double>(
      py::stl_input_iterator<double>(l),
      py::stl_input_iterator<double>()
    ),
    n }
  {
  }


  PyGCE(py::list &l)
  : kde{ 
    std::vector<double>(
      py::stl_input_iterator<double>(l),
      py::stl_input_iterator<double>()
    ),
    -1 }
  {
  }

  PyGCE(py::list &hist, py::list &grid, int frames) 
  : kde{
    std::vector<double>(
      py::stl_input_iterator<double>(hist),
      py::stl_input_iterator<double>()
    ),
    std::vector<double>(
      py::stl_input_iterator<double>(grid),
      py::stl_input_iterator<double>()
    ),
    frames
  } 
  {
  }

  double integrate(const py::str &numericalIntegrationType, double min, double max)
  {
    std::string integralName = std::string(py::extract<char *>(numericalIntegrationType));
    return kde.integrate(integralName, min, max);
  }

  void calculate()
  {
    kde.calculate();
  }

  std::vector<double> getDensityEstimation()
  {
    return kde.getDensityEstimation();
  }


};