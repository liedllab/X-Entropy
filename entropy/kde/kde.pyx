# distutils: language = c++
from libcpp.vector cimport vector
from libcpp.string cimport string

cdef extern from "../../src/kde.cpp":
    pass

cdef extern from "../../src/kde.h":
    cdef cppclass Gce:
        Gce(vector[double]&, int) except +

        void calculate()
        double integrate(string &, double, double)

cdef class GCE:
    cdef Gce*m_gce_kernel

    def __cinit__(self, vector[double] input, int resolution):
        self.m_gce_kernel = new Gce(input, resolution)
    
    def __dealloc__(self):
        del self.m_gce_kernel
    
    def calculate(self):
        self.m_gce_kernel.calculate()