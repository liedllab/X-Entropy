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
    inc_dirs=[]
    lib_dirs=[]
elif platform.system() == 'Linux':
    flags = ['-O3', '-fopenmp']
    libs = ['fftw3', 'm', 'gomp']
    inc_dirs=[]
    lib_dirs=[]
    os.environ['CC'] = 'gcc'
    os.environ['CXX'] = 'g++'
elif platform.system() == 'Darwin':
    flags = ['-O3', '-fopenmp']
    libs = ['fftw3', 'm', 'gomp']
    inc_dirs=[]
    lib_dirs=[]
    os.environ['CC'] = 'gcc'
    os.environ['CXX'] = 'g++'

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
                include_dirs=[],
                library_dirs=[],
                libraries=libs,
                extra_compile_args=flags,
                language = "c++",
            ),
        ])
)
#print("Installed at:\n{}".format(datetime.now()))

