from distutils.core import setup, Extension
import sys


if sys.version_info.major == 3:
    version = (3,8)
else:
    version = (2,7)


setup(
        name = 'DihedralEntropy',
        version = '1.3',
        description = "Calculate the entropy of a given set of data, using the fftw3 library.",
        author = "Johannes Kraml",
        author_email="johannes.kraml@uibk.ac.at",

        py_modules=['entropy.dihedrals'],

        ext_modules = [
            Extension(
                'entropy.kde',
                [
                  'src/kde.cpp',
                  'src/Integrators.cpp',
                  'src/DihedralEntropy.cpp',
                  'src/Module.cpp'
                ],
                libraries=['fftw3', 'm', 'boost_python-py{}{}'.format(version[0], version[1]), 'gomp'],
                extra_compile_args=['-O3', '-fopenmp'],
            ),
        ]
)
