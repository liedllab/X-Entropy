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
	virtual double operator() (const std::vector<double> &function, double range) const noexcept override;
};

class Riemann : public IIntegration {
public:
	Riemann() {};
	virtual double operator() (const std::vector<double> &function, double range) const noexcept override;
};


/****
 * @brief Returns the integration functor.
 * Get the function, so the numerical integration scheme associated with the type string.
 * If the chosen function does not exist, throws an error, hinting the user towards that
 * problem.
 * @param type The type of the numerical integration scheme, i.e., Simpson or Riemann.
 * @return An object that implements the IIntegration interface.
 */
std::unique_ptr<IIntegration> getFunction(const std::string& type);