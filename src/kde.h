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
#include <boost/python.hpp>

#define GASCONSTANT 8.314462618


class IIntegration {
public:
	virtual double operator() (const std::vector<double> &function, int steps, double range) const noexcept = 0;
};

class Simpson : public IIntegration {
public:
	Simpson() {}
	virtual double operator() (const std::vector<double> &function, int steps, double range) const noexcept;
};

class Riemann : public IIntegration {
public:
	Riemann() {};
	virtual double operator() (const std::vector<double> &function, int steps, double range) const noexcept;
};


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
public:
	Gce(double *, int);
	Gce(double *, int, int);
	Gce(const std::vector<double>&, int);

	Gce(boost::python::list &, int);
	Gce(boost::python::list &);
	Gce(boost::python::list &, boost::python::list &, int);

	void calculate(void);
	double integrate_c(const std::string &, double, double);
	double integrate_p(boost::python::str &, double, double);
	int getGridLength() const;
	int getResolution() const;
	double getTStar() const;
	const std::vector<double>& getHistogram() const;
	const std::vector<double>& getGrid() const;
	std::vector<double> getDensityEstimation();
	const std::vector<double>& getAngles() const;
//	boost::python::list getResult(void);
};


class DihedralEntropy {
private:
	double m_entropy = 0;
	double m_res = 0;
	int m_length = 0;
	std::vector<double> m_angles;
	void integrate(IIntegration *t);
	void integrate();
public:
	DihedralEntropy(void) {}
	DihedralEntropy(boost::python::list &, int);
	DihedralEntropy(boost::python::list &, boost::python::str &);
	DihedralEntropy(boost::python::list &, int, boost::python::str &);
	DihedralEntropy(boost::python::list &);
	double getEntropy(void);
};

#endif
