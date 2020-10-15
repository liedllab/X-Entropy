#from distutils.core import setup, Extension
#import sys
#from glob import glob
from setuptools import Extension, setup
from Cython.Build import cythonize

setup(
        name = 'DihedralEntropy',
        version = '1.3',
        description = "Calculate the entropy of a given set of data, using the fftw3 library.",
        author = "Johannes Kraml",
        author_email="johannes.kraml@uibk.ac.at",

        py_modules=['entropy.dihedrals', 'entropy.kde', 'entropy.internal.resolution', 'entropy.internal.pre_post_processing'],

        ext_modules = cythonize([
            Extension(
                'entropy.kde_kernel',
                sources = [
                  'entropy/kde_kernel/kde_kernel.pyx',
                  'src/Integrators.cpp'
                ],
                libraries=['fftw3', 'm', 'gomp'],
                extra_compile_args=['-O3', '-fopenmp'],
                language = "c++",
            ),
        ])
)
