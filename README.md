# fits2x
[![pipeline status](https://gitlab.peulen.xyz/tpeulen/fit2x/badges/master/pipeline.svg)](https://gitlab.peulen.xyz/tpeulen/fit2x/badges/master/pipeline.svg)
[![Anaconda-Server Badge](https://anaconda.org/tpeulen/fit2x/badges/installer/conda.svg)](https://anaconda.org/tpeulen/fit2x)

## Warning
fit2x is currently in early development.

## General description
fit2x is a collection of models that use maximum likelihood
estimators for polarization and time-resolved fluorescence decays.

fit2x implement burst integrated fluorescence lifetime fits (BIFL)
with scatter. The library can be used in conjuncture with tttrlib to 
analyze and  process single-molecule FRET (smFRET) experiments and 
confocal laser  scanning microscopy (CLSM) data. The fit2x shared libary 
can be used from LabView and is wrapped by SWIG (Simplified Wrapper and 
Interface Generator) for common scripting languages as Python as main 
target language. 

![fit2x smfit23][1]


## Design goals
*   Low memory footprint (keep objective large datasets, e.g.  FLIM in memory).
*   Platform independent C/C++ library with interfaces for scripting libraries

## Capabilities
*   Polarization and time-resolved analysis 
*   Stable analysis results with minimum photon counts 
*   Robust and thoroughly tested maximum likelihood estimators

## Examples

```python
import fit2x
import numpy as np

irf = np.array(
    [0, 0, 0, 260, 1582, 155, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 22, 1074, 830, 10, 0, 0,
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
     0, 0], dtype=np.float64
)

data = np.array(
    [0, 0, 0, 1, 9, 7, 5, 5, 5, 2, 2, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 2, 2, 2, 2, 3, 0, 1, 0,
     1, 1, 1, 2, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
)

settings = {
    'dt': 0.5079365079365079,
    'g_factor': 1.0,
    'l1': 0.1,
    'l2': 0.2,
    'convolution_stop': 31,
    'irf': irf,
    'period': 16.0,
    'background': np.zeros_like(irf)
}
fit23 = fit2x.Fit23(**settings)

# initial values
tau, gamma, r0, rho = 2.2, 0.01, 0.38, 1.22
x0 = np.array([tau, gamma, r0, rho])
fixed = np.array([0, 1, 1, 0])
# pass data and initial values to fit. 
# The return value contains the fit result 
r = fit23(
    data=data,
    initial_values=x0,
    fixed=fixed
)
```

## Implementation
Pure pure C/C++ high performance algorithms for analysis of low photon count 
data.

## Building and Installation

### C++ shared library

The C++ shared library can be installed from source with [cmake](https://cmake.org/):

```console
git clone --recursive https://github.com/fluorescence-tools/fit2x.git
mkdir fit2x/build; cd fit2x/build
cmake ..
sudo make install
```

### Python bindings
The Python bindings can be either be installed by downloading and compiling the 
source code or by using a precompiled distribution for Python anaconda environment.

The following commands can be used to download and compile the source code:

```console
git clone --recursive https://github.com/fluorescence-tools/fit2x.git
cd fit2x
sudo python setup.py install
```

In an [anaconda](https://www.anaconda.com/) environment the library can 
be installed by the following command: 
```console
conda install -c tpeulen fit2x
```

For most users the latter approach is recommended. Currently, pre-compiled 
packages for the anaconda distribution system are available for 64bit Windows, macOS,
and Linux.

Legacy 32-bit platforms and versions of programming languages, e.g, Python 2.7 
are not supported.

## Documentation
The API of fit2x as well as some use cases are documented 
on its [web page](https://fluorescence-tools.github.io/fit2x) 

In case you notice unusual behaviour do not hesitate to contact the authors. 
    
## License
fit2x was developed at the [Seidel Lab](<https://www.mpc.uni-duesseldorf.de/>) (Heinrich Heine University). 
and is maintained by Thomas Peulen. fit2x is released under the open source [MIT license](<https://opensource.org/licenses/MIT>).

## Citation
If you have used fit2x in a scientific publication, we would appreciate citations to the following paper: 
[![DOI for citing fit2x](https://img.shields.io/badge/10.1021/ac000877g-blue.svg)](https://doi.org/10.1021/ac000877g)
> Michael Maus, Mircea Cotlet, Johan Hofkens, Thomas Gensch, Frans C. De Schryver, J. Schaffer, and C. A. M. Seidel, 2001. An Experimental Comparison of the Maximum Likelihood Estimation and Nonlinear Least-Squares Fluorescence Lifetime Analysis of Single Molecules. Anal. Chem., 73, 9, pp2078–2086.

 [1]: https://raw.githubusercontent.com/Fluorescence-Tools/fit2x/gh-pages/_images/plot_fit23_1.png "Fit23 single molecule MLE"
 