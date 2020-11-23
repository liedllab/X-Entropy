#pragma once

#include <vector>
#include <string>
#include <iostream>
#include <limits>

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
    double fixedpoint(const std::vector<double> &, const std::vector<double> &);
    std::pair<double, double> extrema(const std::vector<double> &) const noexcept;
    double csiszar(double, int, const std::vector<double> &, const std::vector<double> &) const noexcept;
    double calcIntPow(double, int) const noexcept;
    double calcHistogramNormalizer(int) const;
    double calcHistogramNormalizer(const std::vector<double> &) const;
    void calculateHistogram(double histogramNormalizer, const std::vector<double> &weights);

public:
    Gce(double *, int, int);
    Gce(const std::vector<double> &, int);
    Gce(const std::vector<double> &histogram, const std::vector<double> &xGrid, int res);

    void dct(double *, double *);
    void idct(double *, double *);
    void calculate(void);
    double integrate(const std::string &, double, double);
    double entropy(const std::string &type, double min, double max);
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
