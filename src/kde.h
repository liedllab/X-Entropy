#ifndef KDE_H
#define KDE_H

#include <vector>
#include <string>
#include <iostream>
#include <limits>

#include <cfloat>
#include <cmath>
#include <cstdlib>

#include <fftw3.h>
#include <omp.h>

#include "Exceptions.h"
#include "Integrators.h"

#define GASCONSTANT 8.314462618

class Gce {
private:
	int m_resolution;
	int m_nFrames;
	double m_tStar;
	std::vector<double> m_xgrid;
	std::vector<double> m_histogram;
	std::vector<double> m_angles;
	std::vector<double> m_densityEstimation;
	void dct(double *, double  *);
	void idct(double *, double  *);
	double fixedpoint(const std::vector<double> &, const std::vector<double> &);
	std::pair<double, double> extrema(const std::vector<double> &) const noexcept;
	double csiszar(double, int, const std::vector<double>&, const std::vector<double>&) const noexcept;
  double calcIntPow(double, int) const noexcept;
public:
	Gce(double *, int, int);
	Gce(const std::vector<double>&, int);
  Gce(const std::vector<double>& histogram, const std::vector<double>& xGrid, int res);

	void calculate(void);
	double integrate(const std::string &, double, double);
	int getGridLength() const;
	int getResolution() const;
	double getTStar() const;
	const std::vector<double>& getHistogram() const;
	const std::vector<double>& getGrid() const;
	std::vector<double> getDensityEstimation();
	const std::vector<double>& getAngles() const;
//	boost::python::list getResult(void);
};

#endif
