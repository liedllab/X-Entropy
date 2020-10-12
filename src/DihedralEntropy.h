#pragma once


#include <vector>

#include "Integrators.h"
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
	DihedralEntropy(const std::vector<double> &, int);
	DihedralEntropy(const std::vector<double> &, const std::string &);
	DihedralEntropy(const std::vector<double> &, int, const std::string &);
	DihedralEntropy(const std::vector<double> &);
	double getEntropy(void);
};