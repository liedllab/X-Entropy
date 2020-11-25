#from distutils.core import setup, Extension
#import sys
#from glob import glob
from setuptools import Extension, setup
from Cython.Build import cythonize
import platform
import os
#from datetime import datetime

if platform.system() == 'Windows':
    flags = ['/O2', '/openmp']
    libs = ['fftw3']
elif platform.system() == 'Linux':
    flags = ['-fopenmp']
    libs = ['fftw3', 'm', 'gomp']
elif platform.system() == 'Darwin':
    flags = ['-fopenmp']
    libs = ['fftw3', 'm', 'gomp']

setup(
        name = 'XEntropy',
        version = '1.8',
        description = "Calculate the entropy of a given set of data, using the fftw3 library.",
        author = "Johannes Kraml",
        author_email="johannes.kraml@uibk.ac.at",

        py_modules=['xentropy.dihedrals', 'xentropy.kde', 'xentropy.constants', 'xentropy.xentropy',
            'xentropy.reweighting','xentropy.internal.resolution', 
            'xentropy.internal.pre_post_processing'],

        ext_modules = cythonize([
            Extension(
                'xentropy.kde_kernel',
                sources = [
                  'xentropy/kde_kernel/kde_kernel.pyx',
                  'src/Integrators.cpp'
                ],
                libraries=libs,
                extra_compile_args=flags,
                language = "c++",
            ),
        ])
)
#print("Installed at:\n{}".format(datetime.now()))

