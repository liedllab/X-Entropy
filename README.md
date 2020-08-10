# Dihedral Entropy

![C++:CMake](https://github.com/jokr91/dihedral_entropy/workflows/C++:CMake/badge.svg) ![Python:Pip](https://github.com/jokr91/dihedral_entropy/workflows/Python:Pip/badge.svg)

[![License: GPL v2](https://img.shields.io/badge/License-GPL%20v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)

This is an adaption of the kernel density estimator (kde) of Z. Botev, based on a matlab implementation. It was ported to C++ by Johannes Kraml with the help of Florian Hofer (mostly debugging and stuff), both based on the matlab implementation (the names of the variables were kept for better comparison) and the C implementation of Roland Huber.

Features were added to ensure better parallelism for the kde, as well as for the subsequent integration. Before the kernel density estimation is performed, the dihedral angles are mirrored, in order to circumvent boundary conditions. If this behaviour is undesired, a call to the kde class itself is necessary to calculate the density estimation alone.

## Installation

Please ensure that you are having a python version higher than 3.3 for this package to work properly. First change to the directory containing the entropy module. The package can then be installed via

```bash
pip install .
```

If pip is not working correctly, the command

```bash
python setup.py build
python setup.py install
```

After this, the package is callable by using

```python
import entropy
```

Dependencies:
The dependencies are not yet handled correctly, thus the user is asked to ensure he fulfills the necessary dependencies on his system.

- openmp
- fftw3

As of now, there are no known bugs, if you find any, please contact me.

johannes@kraml.dev

## License Info

This module was parallelized using the openmp toolkit (license: <https://github.com/llvm-mirror/openmp/blob/master/LICENSE.txt>). Parts of the module use the fftw3 library (<http://www.fftw.org/>).

The module entropy is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version. (<https://www.gnu.org/licenses/gpl-2.0.html>)
