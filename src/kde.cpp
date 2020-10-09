/************************************************************************************************************************
 *															                                                                                        *
 *						GCE (1.0)					                                                                                    		*
 *	This is an implementation of the Generalized Cross Entropy Method, based on a matlab script by Z. Botev.	          *
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

#include "kde.h"

/****
 * @brief Constructor
 * Constructor for the GCE class, to make a python importable library, the
 * explanation for the constructor is given only once.
 * @param array: The values for which to calculate the density estimation.
 * @param res: The resolution at which to calculate the density estimation.
 */

Gce::Gce(const std::vector<double> &array, int res) 
  : m_resolution{ res }, m_nFrames{ static_cast<int>( array.size() ) }, m_tStar{ 0 }, m_angles{ array } 
{
  if (m_resolution == -1) {
    m_resolution = 2 << 13;
  }
  if (m_nFrames <= 0)
  {
    throw EmptyListError("Could not construct the GCE object.");
  }
  auto ext{ extrema(m_angles) };
  m_range = ext.second - ext.first;
  double histogram_normalizer{ 1.0 / static_cast<double>(m_angles.size()) };
  
  ext.first -= (m_range * 0.1);
  m_range *= 1.2;
  
  
  // Necessary steps for the calculation of the histogram.
  double stepsize{ m_range / static_cast<double>(m_resolution) };
  for (int i = 0; i <= m_resolution; i++) {
    m_xgrid.push_back(ext.first + i * stepsize);
  }
  m_histogram.resize(m_xgrid.size() - 1);

  /*
   * Here a first density is approximated via a histogram. This is the
   * first step for the Cross Entropy Postulate, i.e., a prior
   * probability density p.
   */
  //#pragma omp parallel for
  for (int i = 0; i < (int) array.size(); ++i) {
    for (int j = 0; j < (int) (m_xgrid.size() - 1); ++j) {
      if ( (array.at(i) > m_xgrid.at(j)) && (array.at(i) <= m_xgrid.at(j + 1)) ) {
        //#pragma omp critical
        m_histogram[j] += histogram_normalizer;
        break;
      }
    }
  }
}

/*****
 * @brief Constructor
 * Constructor for the GCE class, uses grid and histogram. Also takes number Frames.
 * @param histogram The histogram of the data points
 * @param xGrid The grid in x Dimension
 * @param nFrames The number of observations
 */
Gce::Gce(const std::vector<double>& histogram, const std::vector<double>& xGrid, int nFrames)
: m_histogram{ histogram }, m_xgrid{ xGrid }, m_resolution{ xGrid.size() }, m_nFrames(nFrames)
{

}

/****
 * @brief Constructor
 * Constructor for the GCE class, to make a python importable library, the
 * explanation for the constructor is given only once.
 * @param array: The array for which to calculate the density estimation.
 * @param res: The resolution at which to calculate the density estimation.
 * @param length: The length of the array.
 * @deprecated
 */

Gce::Gce(double *array, int length, int res = -1) 
  : m_resolution{ res }, m_nFrames{ length }, m_tStar{ 0 }
{
  if (m_resolution == -1) {
    m_resolution = 2 << 13;
  }
  auto ext{ extrema(m_angles) };
  double range{ ext.second - ext.first };
  double histogram_normalizer{ 1.0 / static_cast<double>(m_angles.size()) };
  
  ext.first -= (range * 0.1);
  range *= 1.2;
  
  
  // Necessary steps for the calculation of the histogram.
  double stepsize{ range / static_cast<double>(m_resolution - 1) };
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
  for (int i = 0; i < length; ++i) {
    for (int j = 0; j < static_cast<int>(m_xgrid.size() - 1); ++j) {
      if ( (array[i] > m_xgrid.at(j)) && (array[i] < m_xgrid.at(j + 1)) ) {
        #pragma omp critical
        m_histogram[j] += histogram_normalizer;
        break;
      }
    }
  }
}

/****
 * @brief Minimum maximum calculation
 * Calculation of the minimum and maximum value of a given data set. Both values
 * are always calculated, as this does not cost much more and returned as a pair
 * of values.
 * @param array: The dataset to calculate the minimum and maximum from.
 * @param maximum: A pointer to the value, where the maximum should
 * 		      be stored
 * @return: The minimum and maximum value in the array.
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
  return { minimum, maximum };
}

/****
 * @brief Starts the calculation, could have guessed, right?
 */
void Gce::calculate() {
  std::vector<double> i_arr(m_xgrid.size());
  std::vector<double> dct_data(m_xgrid.size());
  double error = std::numeric_limits<double>::max();

  // Precalculate the squared indices for the matrix
  for (int i = 1; i < static_cast<int>( m_xgrid.size() ); ++i)
    i_arr[i - 1] = i*i;
  
  // Calculate a discrete cosine transformation
  dct(&m_histogram[0], &dct_data[0]);


  std::vector<double> a2(m_xgrid.size());
//	#pragma omp parallel for
  for (int i = 1; i < static_cast<int>( m_xgrid.size() ); ++i) {
    a2[i - 1] = (dct_data[i] * 0.5);
    a2[i - 1] *= a2[i - 1];
  }
  // Minimize the error iteratively, until the convergence criterion is met.
  while (abs(error) > (DBL_EPSILON)) {
    error = fixedpoint(a2, i_arr);
  }
//	#pragma omp parallel for
  for (int i = 0; i < static_cast<int>( m_xgrid.size() ); ++i) {
    dct_data[i] *= std::exp(-0.5 * i * i * M_PI* M_PI * m_tStar);
  }
  std::vector<double> density(m_xgrid.size());

  // Get the data back into real space with the inverse discrete cosine transform
  idct(&dct_data[0], &density[0]);


  double inv_range{ 1.0 / (m_xgrid.at(m_xgrid.size() - 1) - m_xgrid.at(0)) };
  for (int i = 0; i < (int) m_xgrid.size(); ++i) {
    m_densityEstimation.push_back(density[i] * (0.5 * inv_range));
  }
}

/****
 * @brief Discrete Cosine Transform
 * Calculates the direct cosine transforamtion using the fftw3 C library.
 * @param before_dct The values of the function before fourier transformation.
 * @param after_dct The values of the function after fourier transformation
 */
void Gce::dct(double *before_dct, double *after_dct){
  fftw_plan p{ fftw_plan_r2r_1d(static_cast<int>(m_xgrid.size()), before_dct, after_dct, FFTW_REDFT10, FFTW_ESTIMATE) };
  fftw_execute(p);
  fftw_destroy_plan(p);
}

/****
 * @brief Inverse Discrete Cosine Transform
 * Calculates the inverse direct cosine transformation using the
 * fftw3 C library.
 * @param before_idct The values of the function before inverse fourier transformation.
 * @param after_idct The values of the function after inverse fourier transformation
 */
void Gce::idct(double *before_idct, double *after_idct) {
  fftw_plan p{ fftw_plan_r2r_1d( static_cast<int>(m_xgrid.size()), before_idct, after_idct, FFTW_REDFT01, FFTW_ESTIMATE) };
  fftw_execute(p);
  fftw_destroy_plan(p);
}

/****
 * @brief Calculates new t_star.
 * Calculates the new t_star value. And the error.
 * @param data: An array holding the data, or the initial density.
 * @param i_arr: An array holding the index number, squared.
 * @return: Sets the t_star value of the class to a new value and returns the
 * 	        error.
 */
double Gce::fixedpoint(const std::vector<double> &data, const std::vector<double> &i_arr) {
  // Used for parallelization purposes.
  double f{ 0 };
  double error{ 0 };
  double time{ 0 };


  // This part of the code is used to calculate the inital functional, only that no new t_star will be calculated.
  #pragma omp parallel for reduction(+: f)
  for (int i = 0; i < (int) (m_xgrid.size() - 1); ++i) {
    f += calcIntPow(i_arr[i], 5) * data[i] * exp(-1 * i_arr[i] * M_PI * M_PI * m_tStar);
  }

  f *= 2.0 * calcIntPow(M_PI, 10);
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
  // Set the new t* to the appropriate value.
  m_tStar = time;
  
  return error;
}

/****
 * @brief Csiszar distance
 * Calculates the Csiszar distance between two different values. This is
 * function is iteratively optimized.
 * @param f_before: The last functional, which is updated here.
 * @param s: The s'th step.
 * @param i_arr: An array holding its squared index.
 * @param data: The prior density.
 * @return: The new value for the functional (the new Csiszar measure).
 */
double Gce::csiszar(
  double f_before,
  int s,
  const std::vector<double> &i_arr,
  const std::vector<double> &data
) const noexcept {
  // The first constraint.
  double K0{ 1 };
  // A helper for parallelization reasons
  double f{ 0 };
  // A helper that one does not need to calculate this particular value over and over again.
  double helper{ 0 };
  // One could parallelize this, but this is not worth the effort (max 4 steps)
  for (int i = 1; i <= (2*s - 1); i += 2) {
    K0 *= static_cast<double>( i );
  }
  K0 *= M_2_SQRTPI * M_SQRT1_2 * 0.5;
  helper = pow(2.0 * K0 / (m_nFrames * f_before), 1.0 / (1.5 + s));

  // Calculates the Csiszar measure via Gaussians over the entire grid.
  #pragma omp parallel for reduction(+: f)
  for (int i = 0; i < static_cast<int>(m_xgrid.size() - 1); ++i) {
    f += calcIntPow(i_arr[i], s) * data[i] * exp(-M_PI * M_PI * i_arr[i] * helper);
  }

  return f *  2.0 * calcIntPow(M_PI, 2 * s);
}

/*****
 * @brief Calculate the power function for an integer exponent.
 * To avoid the use of the pow function, we implement this calculation
 * scheme for the power function. This should be faster than pow.
 * @param value The base of the power function
 * @param exponent The exponent of the power function
 * @return The result of the multiplications (power function)
 */
double Gce::calcIntPow(double value, int exponent) const noexcept {
  double ret{ 1 };
  for (int i{ 0 }; i < exponent; ++i)
  {
    ret *= value;
  }
  return ret;
}

/*****
 * @brief Integrate the values.
 * The integration is done between min and max of a given histogram.
 * The integration scheme is gathered from the getFunction function, which
 * returns the appropriate integration function as a IIntegrator pointer.
 * @param type The type of the integration scheme, can be either Simpson
 *             or Riemann integral.
 * @param min The minimal value in the integration (on the x coordinate)
 * @param max The maximum value in the integration (on the x coordinate)
 * @return The integrated value within the given bounds.
 */
double Gce::integrate(const std::string &type, double min, double max) {
  double norm{ 0 };
  auto inte{ getFunction(type) };
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
  norm = (*inte)(fDens, grid.at(grid.size() - 1) - grid.at(0));
  #pragma omp parallel for
  for (int i = 0; i < (int) fDens.size(); ++i) {
    fDens.at(i) /= norm;

    if (fDens.at(i) > 0) 
    {
      fDens.at(i) = fDens.at(i) * log(fDens.at(i));
    } 
    // Proper error handling here
    else if (fDens.at(i) != fDens.at(i)) 
    {
      throw IntegrationError("Created a NAN during integration!");
    }

  }

  return (*inte)(fDens, grid.at(grid.size() - 1) - grid.at(0));
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

double Gce::getBandwidth() const {
  return sqrt(m_tStar) * m_range;
}
