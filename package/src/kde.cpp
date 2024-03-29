#include "kde.h"

double calcHistogramNormalizer(const std::vector<double> &weights)
{
    return 1.0 / std::accumulate(
        std::begin(weights),
        std::end(weights),
        0.0
    );
}

double calcHistogramNormalizer(int size)
{
    return 1.0 / static_cast<double>(size);
}

void KDE::calculateHistogram(double histogramNormalizer, const std::vector<double> &weights)
{
    auto ext{extrema(m_angles)};
    m_range = ext.second - ext.first;

    ext.first -= (m_range * 0.1);
    m_range *= 1.2;

    // Necessary steps for the calculation of the histogram.
    double stepsize{m_range / static_cast<double>(m_resolution)};
    for (int i = 0; i <= m_resolution; i++)
    {
        m_xgrid.push_back(ext.first + i * stepsize);
        if (i >= 1)
        {
            m_centers.push_back(m_xgrid.at(i - 1) + 0.5 * stepsize);
        }
    }

    m_histogram.resize(m_xgrid.size() - 1);

    /*
   * Here a first density is approximated via a histogram. This is the
   * first step for the Cross Entropy Postulate, i.e., a prior
   * probability density p.
   */
    //#pragma omp parallel for
    for (int i = 0; i < (int)m_angles.size(); ++i)
    {
        for (int j = 0; j < (int)(m_xgrid.size() - 1); ++j)
        {
            if ((m_angles.at(i) > m_xgrid.at(j)) && (m_angles.at(i) <= m_xgrid.at(j + 1)))
            {
                if (static_cast<int>(weights.size()) != 0)
                {
                    m_histogram[j] += histogramNormalizer * weights[i];
                }
                else
                {
                    m_histogram[j] += histogramNormalizer;
                }
                break;
            }
        }
    }
}



KDE::KDE(const std::vector<double> &array, int res)
    : m_resolution{res}, m_nFrames{static_cast<int>(array.size())}, m_star_t{0}, m_angles{array}
{
    if (m_resolution == -1)
    {
        m_resolution = 2 << 13;
    }
    if (m_nFrames <= 0)
    {
        throw EmptyListError("Could not construct the GCE object.");
    }

    double histogramNormalizer{calcHistogramNormalizer(m_angles.size())};
    calculateHistogram(histogramNormalizer, std::vector<double>());
}

/*****
 * @brief Constructor
 * Constructor for the GCE class, uses grid and histogram. Also takes number Frames.
 * @param data The histogram of the data points
 * @param weigths The grid in x Dimension
 * @param resolution The number of observations
 */
KDE::KDE(const std::vector<double> &data, const std::vector<double> &weights, int resolution)
    : m_resolution{resolution}, m_nFrames{static_cast<int>(data.size())}, m_angles{data}, m_star_t{0.0}
{
    if (weights.size() != m_angles.size())
    {
        throw ValueError("Lists weigths and data are of different size.");
    }
    if (m_resolution == -1)
    {
        m_resolution = 2 << 11;
    }
    if (m_nFrames <= 0)
    {
        throw EmptyListError("Could not construct the GCE object.");
    }

    double histogramNormalizer{calcHistogramNormalizer(weights)};
    calculateHistogram(histogramNormalizer, weights);
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

KDE::KDE(double *array, int length, int res = -1)
    : m_resolution{res}, m_nFrames{length}, m_star_t{0}
{
    if (m_resolution == -1)
    {
        m_resolution = 2 << 13;
    }
    auto ext{extrema(m_angles)};
    double range{ext.second - ext.first};
    double histogram_normalizer{1.0 / static_cast<double>(m_angles.size())};

    ext.first -= (range * 0.1);
    range *= 1.2;

    // Necessary steps for the calculation of the histogram.
    double stepsize{range / static_cast<double>(m_resolution - 1)};
    for (double i = ext.first; i < (ext.first + range); i += stepsize)
    {
        m_xgrid.push_back(i);
    }
    m_histogram.resize(m_xgrid.size());

/*
   * Here a first density is approximated via a histogram. This is the
   * first step for the Cross Entropy Postulate, i.e., a prior
   * probability density p.
   */
#pragma omp parallel for
    for (int i = 0; i < length; ++i)
    {
        for (int j = 0; j < static_cast<int>(m_xgrid.size() - 1); ++j)
        {
            if ((array[i] > m_xgrid.at(j)) && (array[i] < m_xgrid.at(j + 1)))
            {
#pragma omp atomic
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
std::pair<double, double> KDE::extrema(const std::vector<double> &array) const noexcept
{
    double minimum{array.at(0)};
    double maximum{array.at(0)};

    // Just go over the entire thing and check if the values are maximum
    // or minimum
    for (int i = 0; i < (int)array.size(); ++i)
    {
        if (array.at(i) < minimum)
        {
            minimum = array.at(i);
        }
        else if (array.at(i) > (maximum))
        {
            maximum = array.at(i);
        }
    }
    return {minimum, maximum};
}

/****
 * @brief Starts the calculation, could have guessed, right?
 */
void KDE::calculate()
{
    std::vector<double> i_arr(m_histogram.size());
    std::vector<double> dct_data(m_histogram.size());
    double error = std::numeric_limits<double>::max();

    // Precalculate the squared indices for the matrix
    for (int i = 1; i < static_cast<int>(m_xgrid.size()); ++i)
        i_arr[i - 1] = i * i;

    // Calculate a discrete cosine transformation
    dct(&m_histogram[0], &dct_data[0]);

    std::vector<double> a2(dct_data.size());
    //	#pragma omp parallel for
    for (int i = 0; i < static_cast<int>(dct_data.size()); ++i)
    {
        a2[i] = (dct_data[i] * 0.5);
        a2[i] *= a2[i];
    }
    // Minimize the error iteratively, until the convergence criterion is met.
    while (abs(error) > (DBL_EPSILON))
    {
        error = fixedpoint(a2, i_arr);
    }
    //	#pragma omp parallel for
    for (int i = 0; i < static_cast<int>(dct_data.size()); ++i)
    {
        dct_data[i] *= std::exp(-0.5 * M_PI * M_PI * i_arr[i] * m_star_t);
    }
    std::vector<double> density(dct_data.size());

    // Get the data back into real space with the inverse discrete cosine transform
    idct(&dct_data[0], &density[0]);

    double inv_range{0.5 / m_resolution};

    for (int i = 0; i < static_cast<int>(dct_data.size()); ++i)
    {
        m_densityEstimation.push_back(density[i] * inv_range);
    }

    Simpson simps;
    double stepsize_half{(m_xgrid.at(1) - m_xgrid.at(0)) / 2};
    double norm{simps(m_densityEstimation, m_xgrid.at(m_xgrid.size() - 2) - m_xgrid.at(0))};
    for ( auto &val : m_densityEstimation)
    {
        val /= norm;
    }
}

/****
 * @brief Discrete Cosine Transform
 * Calculates the direct cosine transforamtion using the fftw3 C library.
 * @param before_dct The values of the function before fourier transformation.
 * @param after_dct The values of the function after fourier transformation
 */
void KDE::dct(double *before_dct, double *after_dct)
{
    fftw_plan p{fftw_plan_r2r_1d(static_cast<int>(m_histogram.size()), before_dct, after_dct, FFTW_REDFT10, FFTW_ESTIMATE)};
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
void KDE::idct(double *before_idct, double *after_idct)
{
    fftw_plan p{fftw_plan_r2r_1d(static_cast<int>(m_histogram.size()), before_idct, after_idct, FFTW_REDFT01, FFTW_ESTIMATE)};
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
double KDE::fixedpoint(const std::vector<double> &data, const std::vector<double> &i_arr)
{
    // Used for parallelization purposes.
    double f{0};
    double error{0};
    double time{0};

// This part of the code is used to calculate the inital functional, only that no new t_star will be calculated.
#pragma omp parallel for reduction(+ \
                                   : f)
    for (int i = 0; i < (int)(m_xgrid.size() - 1); ++i)
    {
        f += calcIntPow(i_arr[i], 5) * data[i] * exp(-1 * i_arr[i] * M_PI * M_PI * m_star_t);
    }

    f *= 2.0 * calcIntPow(M_PI, 10);
    // This function uses the old functional (f) to calculate a new one. The functional at this point is the
    // Csiszar Cross entropy
    for (int s = 4; s >= 2; --s)
    {
        f = derivateKDE(f, s, i_arr, data);
    }

    // The actual formula is a bit different, but to optimize away the root,
    // I squared the first part in the pow function. The root is then calculated
    // via the exponent in the pow function (this might be unnecessary optimization,
    // but I think not).
    time = pow(4.0 * m_nFrames * m_nFrames * M_PI * f * f, -0.2);

    // Calculate the relative error of the new t_star (time) and the old t_star.
    error = (m_star_t - time) / time;
    // Set the new t* to the appropriate value.
    m_star_t = time;

    return error;
}

/****
 * @brief Calculate ||f^(s)||^2
 * This function calculates the s'th derivative of the KDE, in order to estimate
 * t*.
 * @param f_before: The last value for ||f^(s+1)||.
 * @param s: The s'th step.
 * @param i_arr: An array holding its squared index.
 * @param data: The prior density.
 * @return: The new value for the derivative.
 */
double KDE::derivateKDE(
    double f_before,
    int s,
    const std::vector<double> &i_arr,
    const std::vector<double> &data) const noexcept
{
    // The first constraint.
    double K0{1};
    // A helper for parallelization reasons
    double f{0};
    // A helper that one does not need to calculate this particular value over and over again.
    double helper{0};
    // One could parallelize this, but this is not worth the effort (max 4
    // steps), with s < 5.
    for (int i = 1; i <= (2 * s - 1); i += 2)
    {
        K0 *= static_cast<double>(i);
    }
    K0 *= M_2_SQRTPI * M_SQRT1_2 * 0.5;
    helper = pow(2.0 * K0 / (m_nFrames * f_before), 1.0 / (1.5 + s));

// Calculates ||f^(s)||^2.
#pragma omp parallel for reduction(+ \
                                   : f)
    for (int i = 0; i < static_cast<int>(m_xgrid.size() - 1); ++i)
    {
        f += calcIntPow(i_arr[i], s) * data[i] * exp(-M_PI * M_PI * i_arr[i] * helper);
    }

    return f * 2.0 * calcIntPow(M_PI, 2 * s);
}

/*****
 * @brief Calculate the power function for an integer exponent.
 * To avoid the use of the pow function, we implement this calculation
 * scheme for the power function. This should be faster than pow.
 * @param value The base of the power function
 * @param exponent The exponent of the power function
 * @return The result of the multiplications (power function)
 */
double KDE::calcIntPow(double value, int exponent) const noexcept
{
    double ret{1};
    for (int i{0}; i < exponent; ++i)
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
double KDE::integrate(const std::string &type, double min, double max)
{
    auto inte{getFunction(type)};
    std::vector<double> grid;
    std::vector<double> fDens;

    double stepsize_half{(m_xgrid.at(1) - m_xgrid.at(0)) / 2};

    for (int i = 0; i < (int)m_xgrid.size() - 1; ++i)
    {
        if ((m_xgrid.at(i) >= min) && (m_xgrid.at(i + 1) <= max))
        {
            grid.push_back(m_xgrid.at(i) + stepsize_half);
            fDens.push_back(m_densityEstimation.at(i));
        }
        else if (m_xgrid.at(i + 1) > max)
        {
            break;
        }
    }
    return (*inte)(fDens, grid.at(grid.size() - 1) - grid.at(0));
}

/*****
 * @brief Calculate the entropy for a given histogram.
 * The integration is done between min and max of a given histogram.
 * The integration scheme is gathered from the getFunction function, which
 * returns the appropriate integration function as a IIntegrator pointer.
 * The data is integrated twice, first to normalize the area under the curve
 * to 1, then from this, p log(p) is calculated, which is finally integrated
 * using a numerical integration scheme to yield the continuous entropy of the
 * curve.
 * @param type The type of the integration scheme, can be either Simpson
 *             or Riemann integral.
 * @param min The minimal value in the integration (on the x coordinate)
 * @param max The maximum value in the integration (on the x coordinate)
 * @return The integrated value within the given bounds.
 */
double KDE::entropy(const std::string &type, double min, double max)
{
    auto inte{getFunction(type)};
    std::vector<double> grid;
    std::vector<double> fDens;

    double stepsize_half{(m_xgrid.at(1) - m_xgrid.at(0)) / 2};

    for (int i = 0; i < (int)m_xgrid.size() - 1; ++i)
    {
        if ((m_xgrid.at(i) >= min) && (m_xgrid.at(i + 1) <= max))
        {
            grid.push_back(m_xgrid.at(i) + stepsize_half);
            fDens.push_back(m_densityEstimation.at(i));
        }
        else if (m_xgrid.at(i + 1) > max)
        {
            break;
        }
    }

    OMPExceptionHandler except;
    double norm{integrate(type, min, max)};
#pragma omp parallel for
    for (int i = 0; i < static_cast<int>(fDens.size()); ++i)
    {
        except.Run([&] {
            if (fDens.at(i) > 0)
            {
                fDens.at(i) /= norm;
                fDens.at(i) = fDens.at(i) * log(fDens.at(i));
            }
            // Proper error handling here
            else if (fDens.at(i) != fDens.at(i))
            {
                throw IntegrationError("Created a NAN during integration!");
            }
        });
    }

    except.Rethrow();

    return (*inte)(fDens, max - min);
}

/****
 * Returns the calculated density values.
 */
std::vector<double> KDE::getDensityEstimation()
{
    return m_densityEstimation;
}

/****
 * Getters
 */
const std::vector<double> &KDE::getAngles() const
{
    return m_angles;
}

const std::vector<double> &KDE::getGrid() const
{
    return m_xgrid;
}

const std::vector<double> &KDE::getCenters() const
{
    return m_centers;
}

const std::vector<double> &KDE::getHistogram() const
{
    return m_histogram;
}

double KDE::getTStar() const
{
    return m_star_t;
}

int KDE::getResolution() const
{
    return m_resolution;
}

int KDE::getGridLength() const
{
    return (int)m_xgrid.size();
}

double KDE::getBandwidth() const
{
    return sqrt(m_star_t) * m_range;
}
