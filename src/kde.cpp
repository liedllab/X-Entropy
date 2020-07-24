/************************************************************************************************************************
 *															                                                                                        *
 *						GCE (1.0)					                                                                                    		*
 *	This is an implementation of the Generalized Cross Entropy Method, based on a matlab script by Y. Botev.	          *
 *	It was reprogrammed from the former GCE, for parallelization and to create a python connectable library.	          *
 *	The GCE Method calculates a density estimation, based on a prior probability density p, a generalized 		          *
 *	Cross Entropy distance D between two probability densities and a finite set of constraints connecting the	          *
 *	probability model with the data. The Cross Entropy Postulate states that, given the values stated before,	          *
 *	the posterior density g can be found uniquely, or more general, given any three of those entities, than 	          *
 *	fourth one can be found uniquely.										                                                                *
 *	Here the Csiszar distance is used as a distance between two probability densities and is minimized to 		          *
 *	calculate the density estimation. And the constraint are all bona fide density functions on the set.		            *
 *															                                                                                        *
 ************************************************************************************************************************
 *															                                                                                        *
 *	Written by Johannes Kraml, based on Code by RG Huber and Z Botev. This program was implemented with the		          *
 *	help of Florian Hofer.												                                                                      *
 *															                                                                                        *
 *	Contact: johannes.kraml@uibk.ac.at										                                                              *
 *															                                                                                        *
 ************************************************************************************************************************/


#include <iostream>
#include <limits>

#include <cfloat>
#include <cmath>
#include <cstdlib>

#include "kde.h"

IIntegration *getFunction(const std::string& type) {
  if (type == "Simpson") {
    return new Simpson();
  } else if (type ==  "Riemann") {
      return new Riemann();
  } else {
    std::cerr << "Error: The function you searched for is not yet implemented." << std::endl;
    return NULL;
  }
}


/****
 * Constructor for the GCE class, to make a python importable library, the
 * explanation for the constructor is given only once.
 * @argument array: The values for which to calculate the density estimation.
 * @argument res: The resolution at which to calculate the density estimation.
 */

Gce::Gce(const std::vector<double> &array, int res) 
  : m_resolution(res), m_nFrames(array.size()), m_tStar(0), m_angles(array) 
{
  if (m_resolution == -1) {
    m_resolution = 2 << 13;
  }
  auto ext{ extrema(m_angles) };
  double range{ ext.second - ext.first };
  double histogram_normalizer{ 1.0 / (double)(m_angles.size()) };
  ext.first -= (range * 0.1);
  range *= 1.2;
  
  
  //Necessary steps for the calculation of the histogram.
  double stepsize{ range / (m_resolution - 1.0) };
  for (double i = ext.first; i < (ext.first + range); i += stepsize) {
    m_xgrid.push_back(i);
  }
  m_histogram.resize(m_xgrid.size());

  /*
   * Here a first density is approximated via a histogram. This is the
   * first step for the Cross Entropy Postulate, i.e., a prior
   * probability density p.
   */
  #pragma omp parallel for
  for (int i = 0; i < (int) array.size(); ++i) {
    for (int j = 0; j < (int) (m_xgrid.size() - 1); ++j) {
      if ( (array.at(i) > m_xgrid.at(j)) && (array.at(i) < m_xgrid.at(j + 1)) ) {
        #pragma omp critical
        m_histogram[j] += histogram_normalizer;
        break;
      }
    }
  }
}

/****
 * Constructor for the GCE class, to make a python importable library, the
 * explanation for the constructor is given only once.
 * @argument array: The array for which to calculate the density estimation.
 * @argument res: The resolution at which to calculate the density estimation.
 * @argument length: The length of the array.
 */

Gce::Gce(double *array, int length, int res = -1) 
  : m_resolution(res), m_nFrames(length), m_tStar(0)
{
//	if (length < m_resolution) {
//		std::cerr << "Error: Resolution too high" << std::endl;
//		m_resolution = 0;
//		return;
//	}
  if (m_resolution == -1) {
    m_resolution = 2 << 13;
  }
  for (int i = 0; i < length; ++i) {
    m_angles.push_back(array[i]);
  }
  auto ext{ extrema(m_angles) };
  double range{ ext.second - ext.first };
  double histogram_normalizer{ 1.0 / static_cast<double>(length) };
  ext.first -= (range / 10);
  range *= 1.2;



  double stepsize = range / (m_resolution - 1);
  for (double i = ext.first; i < range; i += stepsize) {
    m_xgrid.push_back(i);
  }



  m_histogram.resize(m_xgrid.size());
  #pragma omp parallel for
  for (int i = 0; i < length; ++i) {
    for (int j = 0; j < (int) (m_xgrid.size() - 1); ++j) {
      if ( (array[i] > m_xgrid.at(j)) && (array[i] < m_xgrid.at(j + 1)) ) {
        #pragma omp critical
        m_histogram[j] += histogram_normalizer;
        break;
      }
    }
  }
}

/****
 * The copy constructor is not possible with the pyhton requirements of boost,
 * thus, it is not explained at this point.
 */
#ifdef NO_PYTHON
Gce::Gce(Gce &other) : resolution(other.resolution), t_star(other.t_star), n_frames(other.n_frames)
  , xgrid(other.xgrid), angles(other.angles), densityEstimation(other.densityEstimation)

{
  std::vector<double> hist = other.getHistogram();
  m_histogram.resize(m_xgrid.size());
  for (int i = 0; i < (int) hist.size(); ++i) {
    m_histogram[i] = hist.at(i);
  }
}
#endif

/****
 * Calculation of the minimum and maximum value of a given data set.
 * @argument array: The dataset to calculate the minimum and maximum from.
 * @argument maximum: A pointer to the value, where the maximum should
 * 		      be stored
 * @return: The minimum value in the array.
 */
std::pair<double, double> Gce::extrema(const std::vector<double> &array) const noexcept {
  double minimum{ array.at(0) };
  double maximum{ array.at(0) };

  // Just go over the entire thing and check if the values are maximum
  // or minimum
  for (int i = 0 ; i < (int) array.size(); ++i) {
    if (array.at(i) < minimum) {
      minimum = array.at(i);
    } else if ( array.at(i) > (maximum) ) {
      maximum = array.at(i);
    }
  }
  return {minimum, maximum};
}

/****
 * Starts the calculation, could have guessed, right?
 */
void Gce::calculate() {
  std::vector<double> i_arr(m_xgrid.size());
  std::vector<double> dct_data(m_xgrid.size());
  double error = std::numeric_limits<double>::max();

  // Precalculate the squared indices for the matrix
  for (int i = 1; i < (int) m_xgrid.size(); ++i)
    i_arr[i - 1] = i*i;

  
  // Calculate a discrete cosine transformation
  dct(&m_histogram[0], &dct_data[0]);


  std::vector<double> a2(m_xgrid.size());
//	#pragma omp parallel for
  for (int i = 1; i < (int) m_xgrid.size(); ++i) {
    a2[i - 1] = (dct_data[i] * 0.5);
    a2[i - 1] *= a2[i - 1];
  }
  // Minimization of the error iteratively, until the convergence criterion is met.
  while (abs(error) > (DBL_EPSILON)) {
    error = fixedpoint(a2, i_arr);
  }
//	#pragma omp parallel for
  for (int i = 0; i < (int) m_xgrid.size(); ++i) {
    dct_data[i] *= std::exp(-0.5 * i * i * M_PI * M_PI * m_tStar);
  }
  std::vector<double> density(m_xgrid.size());

  // Get the data back into real space with the inverse discrete cosine transform
  idct(&dct_data[0], &density[0]);


  double inv_range{ 1 / (m_xgrid.at(m_xgrid.size() - 1) - m_xgrid.at(0)) };
  for (int i = 0; i < (int) m_xgrid.size(); ++i) {
    m_densityEstimation.push_back(density[i] * (0.5 * inv_range));
  }
}

/****
 * Calculates the direct cosine transforamtion using the fftw3 C library.
 */
void Gce::dct(double *before_dct, double *after_dct){
  fftw_plan p = fftw_plan_r2r_1d((int) m_xgrid.size(), before_dct, after_dct, FFTW_REDFT10, FFTW_ESTIMATE);
  fftw_execute(p);
  fftw_destroy_plan(p);
}

/****
 * Calculates the inverse direct cosine transformation using the
 * fftw3 C library.
 */
void Gce::idct(double *before_idct, double *after_idct) {
  fftw_plan p = fftw_plan_r2r_1d((int) m_xgrid.size(), before_idct, after_idct, FFTW_REDFT01, FFTW_ESTIMATE);
  fftw_execute(p);
  fftw_destroy_plan(p);
}

/****
 * Calculates the new t_star value. And the error.
 * @param data: An array holding the data, or the initial density.
 * @param i_arr: An array holding the index number, squared.
 * @return: Sets the t_star value of the class to a new value and returns the
 * 	        error.
 */
double Gce::fixedpoint(const std::vector<double> &data, const std::vector<double> &i_arr) {
  // Used for parallelization purposes.
  double f_helper [omp_get_max_threads()] = { 0 };
  double f{ 0 };
  double error{ 0 };
  double time{ 0 };
  // This part of the code is used to calculate the inital functional, only that no new t_star will be calculated.
  #pragma omp parallel for
  for (int i = 0; i < (int) (m_xgrid.size() - 1); ++i) {
    f_helper[omp_get_thread_num()] += i_arr[i] * i_arr[i] * i_arr[i] * i_arr[i] * 
        i_arr[i] * data[i] * exp(-1 * i_arr[i] * M_PI_2 * m_tStar);
  }
  for (int i = 0; i < omp_get_max_threads(); ++i) {
    f += f_helper[i];
  }

  f *= 2.0 * M_PI_4 * M_PI_4 * M_PI_2;
  // This function uses the old functional (f) to calculate a new one. The functional at this point is the 
  // Csiszar Cross entropy
  for (int s = 4; s >= 2; --s) {
    f = csiszar(f, s, i_arr, data);
  }

  // The actual formula is a bit different, but to optimize away the root,
  // I squared the first part in the pow function. The root is then calculated
  // via the exponent in the pow function (this might be unnecessary optimization,
  // but I think not).
  time = pow(4.0 * m_nFrames * m_nFrames * M_PI * f * f, -0.2);
  // Calculate the relative error of the new t_star (time) and the old t_star.
  error = (m_tStar - time) / time;
  m_tStar = time;
  return error;
}

/****
 * The correctly written functional function. Calculates the f value.
 * Which is a Csiszar measure of Cross Entropy.
 * @param f_before: The last functional, which is updated here.
 * @param s: The s'th step.
 * @param i_arr: An array holding its squared index.
 * @param data: The prior density.
 * @return: The new value for the functional (the new Csiszar measure).
 */
double Gce::csiszar(double f_before, int s, const std::vector<double> &i_arr, const std::vector<double> &data) const noexcept {
  // The first constraint.
  double K0 = 1;
  // A helper for parallelization reasons
  double f[omp_get_max_threads()] = { 0 };
  // A helper that one does not need to calculate this over and over again.
  double helper = 0;
  // The new f to be returned (hence ret).
  double ret = 0;
  // One could parallelize this, but this is not worth the effort (maxmum 4 steps)
  for (int i = 1; i <= ((2 * s) - 1); i += 2) {
    K0 *= i;
  }
  K0 /= sqrt(2 * M_PI);
  helper = pow(2.0 * K0 / (m_nFrames * f_before), 2.0 / (3.0 + 2.0 * s));
  // Calculates the Csiszar measure via Gaussians over the entire grid.
  #pragma omp parallel for
  for (int i = 0; i < (int) (m_xgrid.size() - 1); ++i) {
    f[omp_get_thread_num()] += pow(i_arr[i], s) * data[i] * exp(-i_arr[i] * M_PI_2 * helper);
  }
  for (int i = 0; i < omp_get_max_threads(); ++i) {
    ret += f[i];
  }
  ret *= 2 * pow(M_PI, 2 * s);
  return ret;
}

/*Fancy*/
double Gce::integrate_c(const std::string &type, double min, double max) {
  double norm = 0;
  IIntegration *inte = getFunction(type);
  std::vector<double> grid;
  std::vector<double> fDens;
  for (int i = 0; i < (int) m_xgrid.size(); ++i) {
    if ((m_xgrid.at(i) >= min) && (m_xgrid.at(i) <= max)) {
      grid.push_back(m_xgrid.at(i));
      fDens.push_back(m_densityEstimation.at(i));
    } else if (m_xgrid.at(i) > max) {
      break;
    }
  }
  norm = 0;
#ifdef SHANNON_ENTROPY
  for (int i = 0; i < (int) fDens.size(); ++i) {
    norm += fDens.at(i);
  }
#else
  norm = (*inte)(fDens, grid.size(), grid.at(grid.size() - 1) - grid.at(0));
#endif
  #pragma omp parallel for
  for (int i = 0; i < (int) fDens.size(); ++i) {
    fDens.at(i) /= norm;
    if (fDens.at(i) > 0) {
#ifdef SHANNON_ENTROPY
      #pragma omp critical
      entropy += fDens.at(i) * log(fDens.at(i));
#else
      fDens.at(i) = fDens.at(i) * log(fDens.at(i));
#endif
    } else if (fDens.at(i) != fDens.at(i)) {
      std::cout << "Nan from somewhere!" << std::endl;
      fDens.at(i) = 0;
    }
  }
  delete inte;
#ifdef SHANNON_ENTROPY
  return entropy;
#else
  return (*inte)(fDens, grid.size(), grid.at(grid.size() - 1) - grid.at(0));
#endif
}


/****
 * Returns the calculated density values.
 */
std::vector<double> Gce::getDensityEstimation() {
  return m_densityEstimation;
}


/****
 * Getters
 */
const std::vector<double>& Gce::getAngles() const {
  return m_angles;
}

const std::vector<double>& Gce::getGrid() const {
  return m_xgrid;
}

const std::vector<double>& Gce::getHistogram() const {
  return m_histogram;
}

double Gce::getTStar() const {
  return m_tStar;
}

int Gce::getResolution() const {
  return m_resolution;
}

int Gce::getGridLength() const {
  return (int) m_xgrid.size();
}

double Simpson::operator() (const std::vector<double> &function, int steps, double range) const noexcept{
  if (steps % 2) {
    steps -= 1;
  }
  double ret = function.at(0);
  double h = range / (3 * steps);
  ret += function.at(function.size() - 1);
  double pcalc[omp_get_max_threads()] = { 0 };
  double calc = 0;
  #pragma omp parallel for
  for (int j = 1; j < (steps / 2); ++j) {
    pcalc[omp_get_thread_num()] += function.at(2 * j - 1);
  }
  for (int i = 0; i < omp_get_max_threads(); ++i) {
    calc += pcalc[i];
    pcalc[i] = 0;
  }
  ret += 4 * calc;
  calc = 0;
  #pragma omp parallel for
  for (int j = 1; j < ((steps / 2) - 1); ++j) {
    pcalc[omp_get_thread_num()] += function.at(2 * j);
  }
  for (int i = 0; i < omp_get_max_threads(); ++i) {
    calc += pcalc[i];
  }
  ret += 2 * calc;
  return ret * h;
}

double Riemann::operator() (const std::vector<double> &function, int steps, double range) const noexcept{
		double dx = range / function.size();
		double ret = 0;
		for (int j = 0; j < (int)function.size(); ++j) {
			ret += function.at(j) * dx;
		}
		return ret;
	}


/********************************************************************************
 * 										*
 * Python ports.								*
 *										*
 ********************************************************************************/

namespace py = boost::python;

DihedralEntropy::DihedralEntropy(py::list &l, int n) : m_entropy(0), m_res(n) {
  m_res =std::pow(2, (int) (log(m_res)/log(2)) + 1);
  for (int i = 0; i < (int) py::len(l); ++i) {
    m_angles.push_back(py::extract<double>(l[i]));
  }
  m_length = m_angles.size();
  integrate();
}

DihedralEntropy::DihedralEntropy(py::list &l) : DihedralEntropy(l, (2 << 12)) {}

DihedralEntropy::DihedralEntropy(py::list &l, int n, py::str &numericalIntegral) : m_entropy(0), m_res(n) {
  m_res =std::pow(2, (int) (log(m_res)/log(2)) + 1);
  for (int i = 0; i < (int) py::len(l); ++i) {
    m_angles.push_back(py::extract<double>(l[i]));
  }
  m_length = m_angles.size();
  std::string numInt{ py::extract<char *>(numericalIntegral) };
  integrate(getFunction(numInt));
}

DihedralEntropy::DihedralEntropy(py::list &l, py::str &numericalIntegral) : 
  DihedralEntropy(l, (2 << 12), numericalIntegral) {}

/****
 * Used to calculate the integral of an dihedral angle distribution. If another
 * Entropy is desired, the GCE (or kde in python) has to be used directly and
 * everything encoded here has to be done by python. Or you are also welcome to
 * just change the Code here and send me the copy ;)
 */
void DihedralEntropy::integrate(IIntegration *t) {
  std::vector<double> mirrored(m_angles.size() * 3);
  double dx = 0;
  double norm = 0;
  // Mirrors the dihedrals to get rid of boundary problems
  #pragma omp parallel for
  for (int i = 0; i < (int) m_angles.size(); ++i) {
    mirrored.at(i) = m_angles.at(i) - 360;
    mirrored.at(i + m_angles.size()) = m_angles.at(i);
    mirrored.at(i + 2 * m_angles.size()) = m_angles.at(i) + 360;
  }
  // Uses the GCE to calculate the Density.
  Gce kde(mirrored, m_res);
  kde.calculate();
  dx = kde.getGrid().at(1) - kde.getGrid().at(0);
  // This is defined here, because the handling of those is easier.
  std::vector<double> density = kde.getDensityEstimation();
  std::vector<double> xgrid = kde.getGrid();
  std::vector<double> finalDens;
  // Integrate only over the angles between -180 and 180 degrees.
  for (int i = 0; i < (int) density.size(); ++i) {
    if (xgrid.at(i) >= -180.0 && xgrid.at(i) <= 180.0) {
      finalDens.push_back(density.at(i));
    }
    if (xgrid.at(i) > 180.0) {
      break;
    }
  }
  // Integration part.
  norm = 0;
  //norm = t->operator()(finalDens, finalDens.size(), 360);
  for (int i = 0; i < (int) finalDens.size(); ++i) {
    norm += finalDens.at(i);
  }
  #pragma omp parallel for
  for (int i = 0; i < (int) finalDens.size(); ++i) {
    finalDens.at(i) /= norm;
    if (finalDens.at(i) > 0) {
      finalDens.at(i) = finalDens.at(i) * log(finalDens.at(i));
    } else if (finalDens.at(i) != finalDens.at(i) ) {
      std::cout << "Nan von irgendwo" << std::endl;
      finalDens.at(i) = 0;
    }
  }
  m_entropy = t->operator()(finalDens, finalDens.size(), 360);
  // Finally done, just multiply with the universal gasconstant (which you will find to be defined in the header,
  // where it belongs!
  m_entropy *= GASCONSTANT;

}

void DihedralEntropy::integrate() {
  integrate(new Riemann());
}

/****
 * Get the entropy.
 * @return: The value for the entropy.
 */
double DihedralEntropy::getEntropy(void) {
  return m_entropy;
}


//py::list Gce::getResult(void) {
//	return std_vector_to_py_list<double>(m_densityEstimation);
//}

/****
 * The same as before, but uses python lists now, which it will convert to an C
 * array. Is als a constructor.
 * @param l: The python list object holding the data.
 * @param n: The resolution of the inital grid, standard (if set to -1) is
 * 		2 ^ 12.
 */
Gce::Gce(py::list &l, int n) 
  : m_resolution(n), m_tStar(0) {
  
  for (int i = 0; i < py::len(l); ++i) {
    m_angles.push_back(py::extract<double>(l[i]));
  }
  if (m_resolution == -1) {
    m_resolution = 2 << 13;
  }
  m_nFrames = m_angles.size();
  auto ext{ extrema(m_angles) };
  double range{ ext.second - ext.first };
  double histogram_normalizer{ 1.0 / (double)(m_angles.size()) };
  ext.first -= (range / 10.0);
  range *= 1.2;

  double stepsize = range / (m_resolution - 1);
  for (double i = ext.first; i < range; i += stepsize) {
    m_xgrid.push_back(i);
  }

  m_histogram.resize(m_xgrid.size());
  #pragma omp parallel for
  for (int i = 0; i < (int) m_angles.size(); ++i) {
    for (int j = 0; j < (int) (m_xgrid.size() - 1); ++j) {
      if ( (m_angles.at(i) > m_xgrid.at(j)) && (m_angles.at(i) < m_xgrid.at(j + 1)) ) {
        #pragma omp critical
        m_histogram[j] += histogram_normalizer;
        break;
      }
    }
  }

}




Gce::Gce(py::list &hist, py::list &grid, int frames) : m_resolution(py::len(hist)), m_tStar(0) {
  m_histogram.resize(py::len(hist));
  for (int i = 0; i < py::len(hist); ++i) {
    m_histogram[i] = py::extract<double>(hist[i]);
  }
  for (int i = 0; i < py::len(grid); ++i) {
    m_xgrid.push_back(py::extract<double>(grid[i]));
  }
  m_nFrames = frames;
}

Gce::Gce(py::list &l) : Gce(l, -1) {}

double Gce::integrate_p(py::str &intName, double min, double max) {
  std::string integralName = std::string(py::extract<char *>(intName));
  return integrate_c(integralName, min, max);
}
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

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
