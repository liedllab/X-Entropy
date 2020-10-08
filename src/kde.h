#ifndef KDE_H
#define KDE_H

#include <vector>
#include <string>
#include <fftw3.h>
#include <omp.h>
#include <boost/python.hpp>

#define GASCONSTANT 8.3144598


class Integration {
public:
	Integration(){};
	virtual double operator() (std::vector<double> function, int steps, double range) = 0;
};

class Simpson : public Integration {
public:
	Simpson();
	virtual double operator() (std::vector<double> function, int steps, double range);
};


class Gce {
private:
	int resolution = 0;
	int n_frames = 0;
	double t_star = 0;
	std::vector<double> xgrid;
	double *histogram = NULL;
	std::vector<double> angles;
	std::vector<double> densityEstimation;
	void dct(double *, double  *);
	void idct(double *, double  *);
	double fixedpoint(double *, double *);
	double minMax(std::vector<double>, double *);
	double functional(double, int, double *, double *);
public:
	Gce(void) {this->histogram = NULL;};
	Gce(double *, int);
	Gce(double *, int, int);
	Gce(std::vector<double>, int);
	Gce(boost::python::list &, int);
	Gce(boost::python::list &);
	Gce(boost::python::list &, boost::python::list &, int);
	~Gce(void);
	void calculate(void);
	double integrate_c(std::string, double, double);
	double integrate_p(boost::python::str &, double, double);
	int getGridLength(void);
	int getResolution(void);
	double getTStar(void);
	std::vector<double> getHistogram(void);
	std::vector<double> getGrid(void);
	std::vector<double> getDensityEstimation(void);
	std::vector<double> getAngles(void);
//	boost::python::list getResult(void);
};

class DihedralEntropy {
private:
	double entropy = 0;
	double res = 0;
	int length = 0;
	std::vector<double> angles;
	void integrate(Integration *t);
	void integrate();
public:
	DihedralEntropy(void) {}
	DihedralEntropy(boost::python::list &, int);
	DihedralEntropy(boost::python::list &, boost::python::str &);
	DihedralEntropy(boost::python::list &, int, boost::python::str &);
	DihedralEntropy(boost::python::list &);
	double getEntropy(void);
};

class Riemann : public Integration {
public:
	Riemann() {};
	virtual double operator() (std::vector<double> function, int steps, double range) {
		double dx = range / function.size();
		double ret = 0;
		for (int j = 0; j < (int)function.size(); ++j) {
			ret += function.at(j) * dx;
		}
		return ret;
	}
};

#endif
