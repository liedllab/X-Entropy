# Dihedral Entropy

![C++:CMake](https://github.com/jokr91/dihedral_entropy/workflows/C++:CMake/badge.svg) ![Python:Pip](https://github.com/jokr91/dihedral_entropy/workflows/Python:Pip/badge.svg)

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This is a library is primarily meant to calculate the dihedral entropy. However, we are using a KDE with automatic 
bandwidth selection as suggested by K. Botev et al. We tried to keep the package as general as possible, therefore,
the package can be used to calculate the entropy of any data, or also to simply calculate the KDE.

## Installation

### Building with pip

Please ensure that you are having a python version higher than 3.3 for this package to work properly. First change 
to the directory containing the entropy module. The package can then be installed via

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
import xentropy
```

### Dependencies

- openmp
- fftw3

If installed using conda, these are added automatically.

## Contributors

 - Johannes Kraml
 - Patrick Quoika
 - Florian Hofer

## Bugs

As of now, there are no known bugs, if you find any, please contact me.

johannes@kraml.dev

## License Info

This module was parallelized using the openmp toolkit (license: <https://github.com/llvm-mirror/openmp/blob/master/LICENSE.txt>). Parts of the module use the fftw3 library (<http://www.fftw.org/>).
