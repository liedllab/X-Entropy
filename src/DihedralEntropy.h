#pragma once


#include <vector>

#include <boost/python.hpp>

#include "Integrators.h"
#include "OMPExceptionHandler.h"
#include "kde.h"


class DihedralEntropy {
private:
	double m_entropy;
	int m_res;
	int m_length;
	std::vector<double> m_angles;
	void integrate(const std::unique_ptr<IIntegration>& t);
	void integrate();
public:
	DihedralEntropy(void) {}
	DihedralEntropy(boost::python::list &, int);
	DihedralEntropy(boost::python::list &, boost::python::str &);
	DihedralEntropy(boost::python::list &, int, boost::python::str &);
	DihedralEntropy(boost::python::list &);
	double getEntropy(void);
};