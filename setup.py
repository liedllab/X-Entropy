#from distutils.core import setup, Extension
#import sys
#from glob import glob
from setuptools import Extension, setup, find_packages
from Cython.Build import cythonize

setup(
        name = 'DihedralEntropy',
        version = '1.3',
        description = "Calculate the entropy of a given set of data, using the fftw3 library.",
        author = "Johannes Kraml",
        author_email="johannes.kraml@uibk.ac.at",

        #py_modules=['entropy.dihedrals', 'entropy.kde', 'entropy.resolution'],
        packages = find_packages(),

        ext_modules = [
            Extension(
                'entropy.kde',
                [
                  'src/kde.cpp',
                  'src/Integrators.cpp',
                  'src/DihedralEntropy.cpp',
                  'src/Module.cpp'
                ],
                libraries=['fftw3', 'm', 'boost_python-py{}'.format(version), 'gomp'],
                extra_compile_args=['-O3', '-fopenmp'],
            ),
        ]
)
