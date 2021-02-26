#pragma once

#include <vector>
#include <string>
#include <iostream>
#include <limits>
#include <numeric>

#include <cfloat>
#include <cmath>
#include <cstdlib>

#include <fftw3.h>

#include "Exceptions.h"
#include "Integrators.h"
#include "OMPExceptionHandler.h"

#define GASCONSTANT 8.314462618

class Gce
{
private:
    int m_resolution;
    int m_nFrames;
    double m_range;
    double m_tStar;
    std::vector<double> m_xgrid;
    std::vector<double> m_centers;
    std::vector<double> m_histogram;
    std::vector<double> m_angles;
    std::vector<double> m_densityEstimation;

    /****
     * @brief Calculates new t_star.
     * Calculates the new t_star value. And the error.
     * @param data: An array holding the data, or the initial density.
     * @param i_arr: An array holding the index number, squared.
     * @return: Sets the t_star value of the class to a new value and returns the
     * 	        error.
     */
    double fixedpoint(const std::vector<double> &, const std::vector<double> &);

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
    std::pair<double, double> extrema(const std::vector<double> &) const noexcept;
    
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
    double derivateKDE(double, int, const std::vector<double> &, const std::vector<double> &) const noexcept;
    
    /*****
     * @brief Calculate the power function for an integer exponent.
     * To avoid the use of the pow function, we implement this calculation
     * scheme for the power function. This should be faster than pow.
     * @param value The base of the power function
     * @param exponent The exponent of the power function
     * @return The result of the multiplications (power function)
     */
    double calcIntPow(double, int) const noexcept;
    double calcHistogramNormalizer(int) const;
    double calcHistogramNormalizer(const std::vector<double> &) const;
    void calculateHistogram(double histogramNormalizer, const std::vector<double> &weights);

public:
    /****
     * @brief Constructor
     * Constructor for the GCE class, to make a python importable library, the
     * explanation for the constructor is given only once.
     * @param array: The array for which to calculate the density estimation.
     * @param res: The resolution at which to calculate the density estimation.
     * @param length: The length of the array.
     * @deprecated
     */
    Gce(double *array, int res, int length);

    /****
     * @brief Constructor
     * Constructor for the GCE class, to make a python importable library, the
     * explanation for the constructor is given only once.
     * @param array: The values for which to calculate the density estimation.
     * @param res: The resolution at which to calculate the density estimation.
     */
    Gce(const std::vector<double> &array, int res);
    
    /*****
     * @brief Constructor
     * Constructor for the GCE class, uses grid and histogram. Also takes number Frames.
     * @param data The histogram of the data points
     * @param weigths The grid in x Dimension
     * @param resolution The number of observations
     */
    Gce(const std::vector<double> &data, const std::vector<double> &weigths, int resolution);

    /****
     * @brief Discrete Cosine Transform
     * Calculates the direct cosine transforamtion using the fftw3 C library.
     * @param before_dct The values of the function before fourier transformation.
     * @param after_dct The values of the function after fourier transformation
     */
    void dct(double *before_dct, double *after_dct);

    /****
     * @brief Inverse Discrete Cosine Transform
     * Calculates the inverse direct cosine transformation using the
     * fftw3 C library.
     * @param before_idct The values of the function before inverse fourier transformation.
     * @param after_idct The values of the function after inverse fourier transformation
     */
    void idct(double *before_idct, double *after_idct);

    /****
     * @brief Starts the calculation, could have guessed, right?
     */
    void calculate(void);

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
    double integrate(const std::string &, double, double);

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
    double entropy(const std::string &type, double min, double max);

    // Getters
    int getGridLength() const;
    int getResolution() const;
    double getTStar() const;
    const std::vector<double> &getHistogram() const;
    const std::vector<double> &getGrid() const;
    const std::vector<double> &getCenters() const;
    std::vector<double> getDensityEstimation();
    const std::vector<double> &getAngles() const;
    double getBandwidth() const;
};
