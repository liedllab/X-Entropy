#pragma once


#include <vector>
#include <memory>
#include <algorithm>

#include <omp.h>

#include "Exceptions.h"




class IIntegration {
public:
  virtual ~IIntegration() {}
	virtual double operator() (const std::vector<double> &function, double range) const noexcept = 0;
};

class Simpson : public IIntegration {
public:
	Simpson() {}
  /*****
   * @brief Uses the composite Simpson's Rule.
   * This function implements the Simpson's rule of numerical integration.
   * In the following is the formula, it is split into different lines in
   * the source code, which represent the different summands. The comments
   * in the code (stated line 1 to 4), are to be understood in the sense
   * of these summands, where the first line is the first summand.
   * @f{eqnarray}
   * \int_{a}^{b}f(x)dx \approx \frac{h}{3}\left[
   * f(x_0) +
   * 2 \Sum_{j = 1}^{n/2 - 1} f(x_{2j} + 
   * 4 \Sum_{j = 1}^{n / 2} f(x_{2j - 1} + 
   * f(x_n)
   * \right]
   * @f}
   * 
   * @param function The values of the function for which the area should be
   *                 calcualted.
   * @param range The range for the integration, i.e., a - b.
   * 
   * @returns The integral for the function, assuming range as a - b.
   */
	virtual double operator() (const std::vector<double> &function, double range) const noexcept override;
};

class Riemann : public IIntegration {
public:
	Riemann() {};
  /*****
   * @brief Uses the Riemann type integral.
   * This function implements the Riemann type integral, which uses rectangles
   * to approximate the area under the curve of a dataset.
   * 
   * @param function The values of the function for which the area should be
   *                 calcualted.
   * @param range The range for the integration, i.e., a - b.
   * 
   * @returns The integral for the function, assuming range as a - b.
   */
	virtual double operator() (const std::vector<double> &function, double range) const noexcept override;
};


/****
 * @brief Returns the integration functor.
 * Get the function, so the numerical integration scheme associated with the type string.
 * If the chosen function does not exist, throws an error, hinting the user towards that
 * problem.
 * @param type The type of the numerical integration scheme, i.e., Simpson or Riemann.
 * @return An object that implements the IIntegration interface.
 * @throw An UnknownIntegral exception, if the chosen integral does not exist (or is)
 *        not implemented.
 * 
 * @see Simpson
 * @see Riemann
 */
std::unique_ptr<IIntegration> getFunction(const std::string& type);