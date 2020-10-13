import cython
from libcpp.vector cimport vector
from libcpp.string cimport string

import numpy as np

cdef extern from "../../src/kde.cpp":
    pass

cdef extern from "../../src/kde.h":
    cdef cppclass Gce:
        Gce(vector[double]&, int) except +

        void calculate() except +
        double integrate(string &, double, double) except +

        vector[double] getDensityEstimation()
        vector[double] getGrid()
        double getBandwidth()

@cython.embedsignature(True)
cdef class __kde_kernel:
    """
    This is a Docstring.
    """
    def __init__(self, data, resolution):
        pass

    cdef Gce* m_gce_kernel

    def __cinit__(self, vector[double] data, int resolution):
        self.m_gce_kernel = new Gce(data, resolution)
    
    def __dealloc__(self):
        del self.m_gce_kernel
    
    def calculate(self):
        self.m_gce_kernel.calculate()

    def integrate(self, start, end, method="Simpson"):
        return self.m_gce_kernel.integrate(method.encode("utf-8"), start, end)

    def getGrid(self):
        return np.array(self.m_gce_kernel.getGrid())

    def getPDF(self):
        return np.array(self.m_gce_kernel.getDensityEstimation())
    
    def getBandwidth(self):
        return self.m_gce_kernel.getBandwidth()