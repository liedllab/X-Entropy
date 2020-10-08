/************************************************************************************************************************
 *															*
 *						GCE (1.0)								*
 *	This is an implementation of the Generalized Cross Entropy Method, based on a matlab script by Y. Botev.	*
 *	It was reprogrammed from the former GCE, for parallelization and to create a python connectable library.	*
 *	The GCE Method calculates a density estimation, based on a prior probability density p, a generalized 		*
 *	Cross Entropy distance D between two probability densities and a finite set of constraints connecting the	*
 *	probability model with the data. The Cross Entropy Postulate states that, given the values stated before,	*
 *	the posterior density g can be found uniquely, or more general, given any three of those entities, than 	*
 *	fourth one can be found uniquely.										*
 *	Here the Csiszar distance is used as a distance between two probability densities and is minimized to 		*
 *	calculate the density estimation. And the constraint are all bona fide density functions on the set.		*
 *															*
 ************************************************************************************************************************
 *															*
 *	Written by Johannes Kraml, based on Code by RG Huber and Z Botev. This program was implemented with the		*
 *	help of Florian Hofer.												*
 *															*
 *	Contact: johannes.kraml@uibk.ac.at										*
 *															*
 ************************************************************************************************************************/









#include <iostream>
#include <stdlib.h>
#include <limits>
#include <cfloat>
#include <math.h>

#include "kde.h"

Integration *getFunction(std::string type) {
	if (type == "Simpson") {
		return new Simpson();
	} else if (type ==  "Riemann") {
			return new Riemann();
	} else {
		std::cerr << "Error: The function you searched for is not yet implemented." << std::endl;
		return NULL;
	}
}


Gce::Gce(double *array, int length) : Gce(array, length, -1) {}

/********************************************************************************
 * 										*
 * Constructor for the GCE class, to make a python importable library, the 	*
 * explanation for the constructor is given only once.				*
 * @argument array: The values for which to calculate the density estimation.	*
 * @argument res: The resolution at which to calculate the density estimation.	*
 * 										*
 ********************************************************************************/

Gce::Gce(std::vector<double> array, int res) 
	: resolution(res), n_frames(array.size()), t_star(0), angles(array) 
{
	if (this->resolution == -1) {
		this->resolution = 2 << 13;
	}
	double maximum = 0;
	double minimum = this->minMax(this->angles, &maximum);
	double range = maximum - minimum;
	double histogram_normalizer = 1.0 / (double)(this->angles.size());
	minimum -= (range / 10);
	range *= 1.2;
	/************************************************************************
	 *									*
	 * Necessary steps for the calculation of the histogram.		*
	 *									*
	 ************************************************************************/
	double stepsize = range / (this->resolution - 1.0);
	for (double i = minimum; i < (minimum + range); i += stepsize) {
		this->xgrid.push_back(i);
	}
	this->histogram = new double [this->xgrid.size()] { 0 };
	/************************************************************************
	 *									*
	 * Here a first density is approximated via a histogram. This is the	*
	 * first step for the Cross Entropy Postulate, i.e., a prior 		*
	 * probability density p.						*
	 *									*
	 ************************************************************************/
	#pragma omp parallel for
	for (int i = 0; i < (int) array.size(); ++i) {
		for (int j = 0; j < (int) (this->xgrid.size() - 1); ++j) {
			if ( (array.at(i) > this->xgrid.at(j)) && (array.at(i) < this->xgrid.at(j + 1)) ) {
				#pragma omp critical
				this->histogram[j] += histogram_normalizer;
				break;
			}
		}
	}
}

/********************************************************************************
 *                                                                              *
 * Constructor for the GCE class, to make a python importable library, the      *
 * explanation for the constructor is given only once.                          *
 * @argument array: The array for which to calculate the density estimation.    *
 * @argument res: The resolution at which to calculate the density estimation.  *
 * @argument length: The length of the array.					*
 *                                                                              *
 *********************************************************************************/

Gce::Gce(double *array, int length, int res) 
	: resolution(res), n_frames(length), t_star(0)
{
//	if (length < this->resolution) {
//		std::cerr << "Error: Resolution too high" << std::endl;
//		this->resolution = 0;
//		return;
//	}
	if (this->resolution == -1) {
		this->resolution = 2 << 13;
	}
	for (int i = 0; i < length; ++i) {
		this->angles.push_back(array[i]);
	}
	double maximum = 0;
	double minimum = this->minMax(this->angles, &maximum);
	double range = maximum - minimum;
	double histogram_normalizer = 1 / (double)length;
	minimum -= (range / 10);
	range *= 1.2;
	double stepsize = range / (this->resolution - 1);
	for (double i = minimum; i < range; i += stepsize) {
		this->xgrid.push_back(i);
	}
	this->histogram = new double [this->xgrid.size()] ();
	#pragma omp parallel for
	for (int i = 0; i < length; ++i) {
		for (int j = 0; j < (int) (this->xgrid.size() - 1); ++j) {
			if ( (array[i] > this->xgrid.at(j)) && (array[i] < this->xgrid.at(j + 1)) ) {
				#pragma omp critical
				this->histogram[j] += histogram_normalizer;
				break;
			}
		}
	}
}

/********************************************************************************
 *										*
 * The copy constructor is not possible with the pyhton requirements of boost,	*
 * thus, it is not explained at this point.					*
 *										*
 ********************************************************************************/
#ifdef NO_PYTHON
Gce::Gce(Gce &other) : resolution(other.resolution), t_star(other.t_star), n_frames(other.n_frames)
	, xgrid(other.xgrid), angles(other.angles), densityEstimation(other.densityEstimation)

{
	std::vector<double> hist = other.getHistogram();
	this->histogram = new double [this->xgrid.size()] ();
	for (int i = 0; i < (int) hist.size(); ++i) {
		this->histogram[i] = hist.at(i);
	}
}
#endif


/********************************************************************************
 *										*
 * This is the Destructor, it frees the only allocated memory, the histogram.	*
 * 										*
 ********************************************************************************/
Gce::~Gce(void) {
	delete[] this->histogram;
}

/********************************************************************************
 *										*
 * Calculation of the minimum and maximum value of a given data set.		*
 * @argument array: The dataset to calculate the minimum and maximum from.	*
 * @argument maximum: A pointer to the value, where the maximum should		*
 * 		      be stored							*
 * @return: The minimum value in the array.					*
 * 										*
 ********************************************************************************/
double Gce::minMax(std::vector<double> array, double *maximum) {
	double ret = array.at(0);
	// If maximum is set to NULL, don't find maximum
	if (maximum != NULL) {
		*maximum = ret;
	}
	for (int i = 0 ; i < (int) array.size(); ++i) {
		if (array.at(i) < ret) {
			ret = array.at(i);
		} else if ( (array.at(i) > (*maximum)) && (maximum != NULL) ) {
			*maximum = array.at(i);
		}
	}
	return ret;
}

/********************************************************************************
 *										*
 * Starts the calculation, could have guessed, right?				*
 *										*
 ********************************************************************************/
void Gce::calculate(void) {
	double *I = new double [this->xgrid.size()] ();
	double *dct_data = new double [this->xgrid.size()] ();
	double error = std::numeric_limits<double>::max();
//	#pragma omp parallel for
	for (int i = 1; i < (int) this->xgrid.size(); ++i) {
		I[i - 1] = i*i;
	}
	// Calculate a discrete cosine transformation
	this->dct(this->histogram, dct_data);
	double *a2 = new double[this->xgrid.size()] ();
//	#pragma omp parallel for
	for (int i = 1; i < (int) this->xgrid.size(); ++i) {
		a2[i - 1] = (dct_data[i] / 2);
		a2[i - 1] *= a2[i - 1];
	}
	// Minimization of the error iteratively, until the convergence criterion is met.
	while (abs(error) > (DBL_EPSILON)) {
		error = this->fixedpoint(a2, I);
	}
	delete[] a2;
//	#pragma omp parallel for
	for (int i = 0; i < (int) this->xgrid.size(); ++i) {
		dct_data[i] *= std::exp(-1 * i * i * M_PI * M_PI * this->t_star / 2);
	}
	double * density = new double[this->xgrid.size()] ();
	this->idct(dct_data, density);
	double range = this->xgrid.at(this->xgrid.size() - 1) - this->xgrid.at(0);
	for (int i = 0; i < (int) this->xgrid.size(); ++i) {
		this->densityEstimation.push_back(density[i]/ (2 * range));
	}
	delete[] dct_data; delete[] I; delete[] density;
}

/********************************************************************************
 *										*
 * Calculates the direct cosine transforamtion using the fftw3 C library.	*
 *										*
 ********************************************************************************/
void Gce::dct(double *before_dct, double *after_dct){
	fftw_plan p = fftw_plan_r2r_1d((int) this->xgrid.size(), before_dct, after_dct, FFTW_REDFT10, FFTW_ESTIMATE);
	fftw_execute(p);
	fftw_destroy_plan(p);
}

/********************************************************************************
 *										*
 * Calculates the inverse direct cosine transformation using the 		*
 * fftw3 C library.								*
 *										*
 ********************************************************************************/
void Gce::idct(double *before_idct, double *after_idct) {
	fftw_plan p = fftw_plan_r2r_1d((int) this->xgrid.size(), before_idct, after_idct, FFTW_REDFT01, FFTW_ESTIMATE);
	fftw_execute(p);
	fftw_destroy_plan(p);
}

/********************************************************************************
 *										*
 * Calculates the new t_star value. And the error.				*
 * @argument data: An array holding the data, or the initial density.		*
 * @argument I: An array holding the index number, squared.			*
 * @return: Sets the t_star value of the class to a new value and returns the	*
 * 	    error.								*
 *										*
 ********************************************************************************/
double Gce::fixedpoint(double *data, double *I) {
	// Used for parallelization purposes.
	double f_helper [omp_get_max_threads()] = { 0 };
	double f = 0;
	double error = 0;
	double time = 0;
	// This part of the code is used to calculate the inital functional, only that no new t_star will be calculated.
	#pragma omp parallel for
	for (int i = 0; i < (int) (this->xgrid.size() - 1); ++i) {
		f_helper[omp_get_thread_num()] += pow(I[i], 
				(double)5) * data[i] * exp(-1 * I[i] * M_PI * M_PI * this->t_star);
	}
	for (int i = 0; i < omp_get_max_threads(); ++i) {
		f += f_helper[i];
	}
	f *= 2 * pow(M_PI, 10);
	// This function uses the old functional (f) to calculate a new one. The functional at this point is the 
	// Csiszar Cross entropy
	for (int s = 4; s >= 2; --s) {
		f = this->functional(f, s, I, data);
	}
	time = pow(2 * this->n_frames * sqrt(M_PI) * f, -2 / (double)5.0);
	// Calculate the relative error of the new t_star (time) and the old t_star.
	error = (this->t_star - time) / time;
	this->t_star = time;
	return error;
}

/********************************************************************************
 *										*
 * The correctly written functional function. Calculates the f value.		*
 * Which is a Csiszar measure of Cross Entropy.					*
 * @argument f_before: The last functional, which is updated here.		*
 * @argument s: The s'th step.							*
 * @argument I: An array holding its squared index.				*
 * @argument data: The prior density.						*
 * @return: The new value for the functional (the new Csiszar measure).		*
 *										*
 ********************************************************************************/
double Gce::functional(double f_before, int s, double *I, double *data) {
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
	helper = pow(2 * K0 / (this->n_frames * f_before), (double)2.0/ (double)(3.0 + 2 * s));
	// Calculates the Csiszar measure via Gaussians over the entire grid.
	#pragma omp parallel for
	for (int i = 0; i < (int) (this->xgrid.size() - 1); ++i) {
		f[omp_get_thread_num()] += pow(I[i], s) * data[i] * exp(-I[i] * M_PI * M_PI * helper);
	}
	for (int i = 0; i < omp_get_max_threads(); ++i) {
		ret += f[i];
	}
	ret *= 2 * pow(M_PI, 2 * s);
	return ret;
}

/*Fancy*/
double Gce::integrate_c(std::string type, double min, double max) {
	double norm = 0;
	Integration *inte = getFunction(type);
	std::vector<double> grid;
	std::vector<double> fDens;
	for (int i = 0; i < (int) this->xgrid.size(); ++i) {
		if ((this->xgrid.at(i) >= min) && (this->xgrid.at(i) <= max)) {
			grid.push_back(this->xgrid.at(i));
			fDens.push_back(this->densityEstimation.at(i));
		} else if (this->xgrid.at(i) > max) {
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
#ifdef SHANNON_ENTROPY
	return entropy;
#else
	return (*inte)(fDens, grid.size(), grid.at(grid.size() - 1) - grid.at(0));
#endif
}


/********************************************************************************
 *										*
 * Returns the calculated density values.					*
 *										*
 ********************************************************************************/
std::vector<double> Gce::getDensityEstimation(void) {
	return this->densityEstimation;
}


/********************************************************************************
 *										*
 * Getters									*
 *										*
 ********************************************************************************/
std::vector<double> Gce::getAngles(void) {
	return this->angles;
}

std::vector<double> Gce::getGrid(void) {
	return this->xgrid;
}

std::vector<double> Gce::getHistogram(void) {
	std::vector<double> ret;
	for (int i = 0; i < (int) this->xgrid.size(); ++i) {
		ret.push_back(this->histogram[i]);
	}
	return ret;
}

double Gce::getTStar(void) {
	return this->t_star;
}

int Gce::getResolution(void) {
	return this->resolution;
}

int Gce::getGridLength(void) {
	return (int) this->xgrid.size();
}

Simpson::Simpson(void) {}

double Simpson::operator() (std::vector<double> function, int steps, double range) {
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



/********************************************************************************
 * 										*
 * Python ports.								*
 *										*
 ********************************************************************************/

namespace py = boost::python;

DihedralEntropy::DihedralEntropy(py::list &l, int n) : entropy(0), res(n) {
	this->res =std::pow(2, (int) (log(this->res)/log(2)) + 1);
	for (int i = 0; i < (int) py::len(l); ++i) {
		this->angles.push_back(py::extract<double>(l[i]));
	}
	this->length = this->angles.size();
	this->integrate();
}

DihedralEntropy::DihedralEntropy(py::list &l) : DihedralEntropy(l, (2 << 12)) {}

DihedralEntropy::DihedralEntropy(py::list &l, int n, py::str &numericalIntegral) : entropy(0), res(n) {
	this->res =std::pow(2, (int) (log(this->res)/log(2)) + 1);
	for (int i = 0; i < (int) py::len(l); ++i) {
		this->angles.push_back(py::extract<double>(l[i]));
	}
	this->length = this->angles.size();
	std::string numInt = std::string(py::extract<char *>(numericalIntegral));
	this->integrate(getFunction(numInt));
}

DihedralEntropy::DihedralEntropy(py::list &l, py::str &numericalIntegral) : 
	DihedralEntropy(l, (2 << 12), numericalIntegral) {}

/********************************************************************************
 *										*
 * Used to calculate the integral of an dihedral angle distribution. If another	*
 * Entropy is desired, the GCE (or kde in python) has to be used directly and	*
 * everything encoded here has to be done by python. Or you are als welcome to	*
 * just change the Code here and send me the copy ;)				*
 * 										*
 ********************************************************************************/
void DihedralEntropy::integrate(Integration *t) {
	std::vector<double> mirrored(this->angles.size() * 3);
	double dx = 0;
	double norm = 0;
	// Mirrors the dihedrals to get rid of boundary problems
	#pragma omp parallel for
	for (int i = 0; i < (int) this->angles.size(); ++i) {
		mirrored.at(i) = this->angles.at(i) - 360;
		mirrored.at(i + this->angles.size()) = this->angles.at(i);
		mirrored.at(i + 2 * this->angles.size()) = this->angles.at(i) + 360;
	}
	// Uses the GCE to calculate the Density.
	Gce kde(mirrored, this->res);
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
	this->entropy = t->operator()(finalDens, finalDens.size(), 360);
	// Finally done, just multiply with the universal gasconstant (which you will find to be defined in the header,
	// where it belongs!
	this->entropy *= GASCONSTANT;

}

void DihedralEntropy::integrate() {
	this->integrate(new Riemann());
}

/********************************************************************************
 * 										*
 * Get the entropy.								*
 * @return: The value for the entropy.						*
 *										*
 ********************************************************************************/
double DihedralEntropy::getEntropy(void) {
	return this->entropy;
}

//py::list Gce::getResult(void) {
//	return std_vector_to_py_list<double>(this->densityEstimation);
//}

/********************************************************************************
 * 										*
 * The same as before, but uses python lists now, which it will convert to an C	*
 * array. Is als a constructor.							*
 * @argument l: The python list object holding the data.			*
 * @argument n: The resolution of the inital grid, standard (if set to -1) is 	*
 * 		2 ^ 12.								*
 * 										*
 ********************************************************************************/
Gce::Gce(py::list &l, int n) : resolution(n), t_star(0) {
	for (int i = 0; i < py::len(l); ++i) {
		this->angles.push_back(py::extract<double>(l[i]));
	}
	if (this->resolution == -1) {
		this->resolution = 2 << 13;
	}
	this->n_frames = angles.size();
	double maximum = 0;
	double minimum = this->minMax(this->angles, &maximum);
	double range = maximum - minimum;
	double histogram_normalizer = 1 / (double)(this->angles.size());
	minimum -= (range / 10);
	range *= 1.2;
	double stepsize = range / (this->resolution - 1);
	for (double i = minimum; i < range; i += stepsize) {
		this->xgrid.push_back(i);
	}
	this->histogram = new double[xgrid.size()] ();
	#pragma omp parallel for
	for (int i = 0; i < (int) this->angles.size(); ++i) {
		for (int j = 0; j < (int) (this->xgrid.size() - 1); ++j) {
			if ( (this->angles.at(i) > this->xgrid.at(j)) && (this->angles.at(i) < this->xgrid.at(j + 1)) ) {
				#pragma omp critical
				this->histogram[j] += histogram_normalizer;
				break;
			}
		}
	}
}
Gce::Gce(py::list &hist, py::list &grid, int frames) : resolution(py::len(hist)), t_star(0) {
	this->histogram = new double[py::len(hist)] ();
	for (int i = 0; i < py::len(hist); ++i) {
		this->histogram[i] = py::extract<double>(hist[i]);
	}
	for (int i = 0; i < py::len(grid); ++i) {
		this->xgrid.push_back(py::extract<double>(grid[i]));
	}
	this->n_frames = frames;
}

Gce::Gce(py::list &l) : Gce(l, -1) {}

double Gce::integrate_p(py::str &intName, double min, double max) {
	std::string integralName = std::string(py::extract<char *>(intName));
	return this->integrate_c(integralName, min, max);
}

#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

/********************************************************************************
 * 										*
 * The python modules								*
 * 										*
 ********************************************************************************/
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
