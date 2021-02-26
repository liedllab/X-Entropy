# X-Entropy

![C++:CMake](https://github.com/jokr91/dihedral_entropy/workflows/C++:CMake/badge.svg) ![Python:Pip](https://github.com/jokr91/dihedral_entropy/workflows/Python:Pip/badge.svg)

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This library is primarily meant to calculate dihedral entropies from MD
simulation data. For this, we use a KDE with automatic bandwidth selection as
suggested by K. Botev et al. We tried to keep the package as generalized as
possible, therefore, the package can be used to calculate the entropy of any
data, or also to simply calculate the KDE.

## Installation

### Using Conda

This is probably the easiest way to install this package. After cloning this
repository (or downloading the conda directory directly), you can easily install
the package using conda. Be aware that in this case you need Python 3.8 installed
for Linux and MacOS and Python 3.7 for Windows.
```bash
conda install -c ${ABSOLUTE_PATH_TO_XENTROPY_REPOSITORY}/conda xentropy
```
This should work for Windows 10, Linux, MacOS. Be aware, that for Windows you
need to specify the full path to the conda directory. If it for some reason does
not work, you might have to build the conda installer yourself. If you want to
do that, you need to install conda-build and conda-verify and then just build the
installer in the xentropy root diretory with
```bash
conda-build package
```
On MacOS you will need to follow the instructions on
https://docs.conda.io/projects/conda-build/en/latest/resources/compiler-tools.html#macos-sdk
before building the package.

### Building with pip

Be aware that you are responsible for the dependencies, if you decide to build 
this way.

Please ensure that you have a python version higher than 3.3 for this package to
work properly. First change to the directory containing the entropy module. The
package can then be installed via

```bash
pip install .
```

### Building with setup.py

Be aware that you are responsible for the dependencies, if you decide to build 
this way.

```bash
python setup.py build
python setup.py install
```

### Import

After this, the package is callable by using

```python
import xentropy
```

### Dependencies

- openmp
- fftw3
- numpy
- Cython

If installed using conda, these are added automatically.

## Usage

For usage examples, please check the examples folder.

## Contributors

 - Johannes Kraml
 - Patrick Quoika
 - Florian Hofer

## Bugs

As of now, there are no known bugs, if you find any, please contact me.

johannes@kraml.dev

## License Info

This module was parallelized using the openmp toolkit (license:
<https://github.com/llvm-mirror/openmp/blob/master/LICENSE.txt>). Parts of the
module use the fftw3 library (<http://www.fftw.org/>).

The KDE in this package is based on the publication by Z. Botev: 
[Kernel density estimation via diffusion
Z. I. Botev, J. F. Grotowski, and D. P. Kroese (2010)
Annals of Statistics, Volume 38, Number 5, pages 2916-2957
doi:10.1214/10-AOS799]
(https://projecteuclid.org/journals/annals-of-statistics/volume-38/issue-5/Kernel-density-estimation-via-diffusion/10.1214/10-AOS799.full).
