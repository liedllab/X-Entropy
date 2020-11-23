#from distutils.core import setup, Extension
#import sys
#from glob import glob
from setuptools import Extension, setup
from Cython.Build import cythonize
#from datetime import datetime

win_flags = ['/O2', '/openmp']
win_libs = ['fftw3']

linux_flags = ['-O3', '-fopenmp']
linux_libs = ['fftw3', 'm', 'gomp']

setup(
        name = 'DihedralEntropy',
        version = '1.8',
        description = "Calculate the entropy of a given set of data, using the fftw3 library.",
        author = "Johannes Kraml",
        author_email="johannes.kraml@uibk.ac.at",

        py_modules=['entropy.dihedrals', 'entropy.kde', 'entropy.constants', 'entropy.entropy',
            'entropy.reweighting','entropy.internal.resolution', 
            'entropy.internal.pre_post_processing'],

        ext_modules = cythonize([
            Extension(
                'entropy.kde_kernel',
                sources = [
                  'entropy/kde_kernel/kde_kernel.pyx',
                  'src/Integrators.cpp'
                ],
                libraries=win_libs,
                extra_compile_args=win_flags,
                language = "c++",
            ),
        ])
)
#print("Installed at:\n{}".format(datetime.now()))

