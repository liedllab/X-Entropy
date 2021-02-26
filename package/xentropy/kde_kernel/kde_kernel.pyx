import cython
from libcpp.vector cimport vector
from libcpp.string cimport string

import numpy as np

cdef extern from "../../src/kde.cpp":
    pass

cdef extern from "../../src/kde.h":
    cdef cppclass KDE:
        KDE(vector[double]&, int) except +
        KDE(vector[double]&, vector[double] &, int) except +

        void calculate() except +
        double integrate(string &, double, double) except +
        double entropy(string &, double, double) except +

        vector[double] getDensityEstimation()
        vector[double] getCenters()
        double getBandwidth()

@cython.embedsignature(True)
cdef class _kde_kernel:
    """
    This is a Docstring.
    """
    def __init__(self, data, resolution, weights=None):
        pass

    cdef KDE* m_kde_kernel
    
    def __cinit__(self, data, resolution, weights=None):
        if not (weights is None):
            self.m_kde_kernel = new KDE(data, weights, resolution)
        else:
            self.m_kde_kernel = new KDE(data, resolution)
    
    def __dealloc__(self):
        del self.m_kde_kernel
    
    def calculate(self):
        self.m_kde_kernel.calculate()

    def integrate(self, start, end, method="Simpson"):
        return self.m_kde_kernel.integrate(method.encode("utf-8"), start, end)
    
    def calculate_entropy(self, start, end, method="Simpson"):
        return self.m_kde_kernel.entropy(method.encode("utf-8"), start, end)

    def get_grid(self):
        return np.array(self.m_kde_kernel.getCenters())

    def get_pdf(self):
        return np.array(self.m_kde_kernel.getDensityEstimation())
    
    def get_bandwidth(self):
        return self.m_kde_kernel.getBandwidth()