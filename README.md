# Dihedral Entropy

This is an adaption of the kernel density estimator (kde) of Z. Botev, based on a matlab implementation. It was ported to C++ by Johannes Kraml with the help of Florian Hofer (mostly debugging and stuff), both based on the matlab implementation (the names of the variables were kept for better comparison) and the C implementation of Roland Huber.

Features were added to ensure better parallelism for the kde, as well as for the subsequent integration. Before the kernel density estimation is performed, the dihedral angles are mirrored, in order to circumvent boundary conditions. If this behaviour is undesired, a call to the kde class itself is necessary to calculate the density estimation alone.

## Installation

The package can be installed via

```bash
pip install /PATH/TO/ENTROPY
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

Because of some stupidity between Python versions 2.x and 3.3+, if your Python version is either 2.x or < 3.3, the user has to add a __init__.py to the dihedrals.py directory (this file can be empty). On Python 3.3+ the user should not do that. A fix is in the making.

As of now, there are no known bugs, if you find any, please contact me.

johannes@kraml.dev

## License Info

This module was parallelized using the openmp toolkit (license: <https://github.com/llvm-mirror/openmp/blob/master/LICENSE.txt>). Parts of the module use the fftw3 library (<http://www.fftw.org/>).

The module entropy is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version. (<https://www.gnu.org/licenses/gpl-2.0.html>)
