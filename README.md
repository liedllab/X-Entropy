# Dihedral Entropy

![C++:CMake](https://github.com/jokr91/dihedral_entropy/workflows/C++:CMake/badge.svg) ![Python:Pip](https://github.com/jokr91/dihedral_entropy/workflows/Python:Pip/badge.svg)

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This is an adaption of the kernel density estimator (kde) of Z. Botev, based on a matlab implementation. It was ported to C++ by Johannes Kraml with the help of Florian Hofer (mostly debugging and stuff), both based on the matlab implementation (the names of the variables were kept for better comparison) and the C implementation of Roland Huber.

Features were added to ensure better parallelism for the kde, as well as for the subsequent integration. Before the kernel density estimation is performed, the dihedral angles are mirrored, in order to circumvent boundary conditions. If this behaviour is undesired, a call to the kde class itself is necessary to calculate the density estimation alone.

## Installation

### Building with pip

Please ensure that you are having a python version higher than 3.3 for this package to work properly. First change to the directory containing the entropy module. The package can then be installed via

```bash
pip install .
```

### Building with setup.py

```bash
python setup.py build
python setup.py install
```

### Using Conda

This is probably the easiest way to install this package. After cloning this repository 
(or downloading the conda directory directly), you can easily install the package using
conda.
```bash
conda install -c /PATH/TO/XENTROPY/CONDA/DIRECTORY xentropy
```
This should work for windows 10, linux, MacOS. If it for some reason does not work, you
might have to build the conda installer yourself. If you want to do that, you need to
install conda-build and then just build the installer when you are in the xentropy root 
diretory with
```bash
conda-build .
```
On MacOS you will need to follow the instructions on 
https://docs.conda.io/projects/conda-build/en/latest/resources/compiler-tools.html#macos-sdk

### Test if installation worked

After this, the package is callable by using

```python
import entropy
```

### Dependencies

- openmp
- fftw3

If installed using conda, these are added automatically.

## Bugs

As of now, there are no known bugs, if you find any, please contact me.

johannes@kraml.dev

## License Info

This module was parallelized using the openmp toolkit (license: <https://github.com/llvm-mirror/openmp/blob/master/LICENSE.txt>). Parts of the module use the fftw3 library (<http://www.fftw.org/>).
