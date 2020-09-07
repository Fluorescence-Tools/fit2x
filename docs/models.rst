Fits
====
General
-------
The fits covered by fit2x (e.g. fit23, fit24, and fit25) all have their own
model function, target (objective) function, fit function, and parameter
correction function. Briefly, the purpose of these functions is as follows.

+--------------------------+--------------------------------------------+
| Function                 | Description                                |
+==========================+============================================+
| Model function           | Computes the model for a set of parameters |
+--------------------------+--------------------------------------------+
| Target/objective function| Computes ob and returns a score for data   |
+--------------------------+--------------------------------------------+
| Fit function             | Optimizes the parameters to the data       |
+--------------------------+--------------------------------------------+
| Correction function      | Assures reasonable range of parameters     |
+--------------------------+--------------------------------------------+

The model function computes for a set of input parameters the corresponding
realization of the model (the model fluorescence decay curve(s)). The target
function (more often called objective function) computes the disagreement of the
data and the model function for a given set of model parameters. The fitting
function optimizes a set of selected input parameters to the data by minimizing
the disagreement of the data and the model. Finally, the estimates of the model
parameters are corrected by a correction function. These functions have names
that relate to the name of the model, e.g., `target23` corresponds to `fit23`.
For computationally efficiency, the functions are hardcoded in C.

The models of the different fits differ in their parameters. The parameters of
the fits are described in detail below. The target functions explicitly consider
the counting noise and provide an estimate of the likelihood of the parameters.
The initial values of the parameters are specified by the user. Next, the likelihood
of the target is maximized by the Limited-Memory Broyden-Fletcher-Goldfarb-Shanno
Algorithm (L-BFGS). L-BFGS is an iterative method for solving unconstrained nonlinear
optimization problems. Gradients for the L-BFGS algorithm are approximated by linear
interpolation to handle the differentiable components, i.e., the model parameters
that are free. The parameters of a model ca be either specified as free or as fixed.
Only free parameters are optimized.

Python
------
For a use in python the ``fit2x`` module exposes a set of C functions that can be
used to (1) compute model and (2) target/objective function and (3) optimize
model parameters to experimental data. Besides the exposed functions the fit models
are accessible via a simplified object-based interface that reduces the number of
lines of code that need to be written for analyzing fluorescence decay histograms.
The code blocks that are used below to illustrate the fit2x functionality are
extracts from the tests located in the test folder of the fit2x repository. The
test can be used as a more detailed refernce on how to use fit2x.

Model functions can be computed using the ``modelf2x`` function of the ``fit2x``
module. Here, the ``x`` represents a particular model function. The use of a
``fit2x`` model function is for the model function 23 (``modelf23``) below.
model function

.. literalinclude:: ../examples/fluorescence_decay/plot_fit23_1.py
   :language: python
   :lines: 36-44
   :linenos:

To compute a model, first a set of model parameters and a set of corrections need
to be specified. All input parameters for ``modelf23`` are numpy arrays. In addition
to the model parameters and the corrections ``modelf23`` requires an instrument
response function (irf) and a background pattern. The model functions operate on
numpy arrays and modify the numpy array for a model in-place. This means, that the
output of the function is written to the input numpy-array of the model. In the
example above the output is written to the array ``model``.


To compute the value of a target for a realization of model parameters ``fit2x``
provides the functions ``targetf2x``. The inputs of a target function (here
``tarferf23``) are an numpy array containing the model parameter values and a
structure that contains the corrections and all other necessary data to compute
the value of the objective function (for fit23 i.e. data, irf, backgound, time
resolution).

.. literalinclude:: ../test/test_fit23.py
   :language: python
   :lines: 171-179
   :linenos:

The data needed in addition to the model parameters are passed to the target function
using ``fit2x.MParam`` objects that can be created by the ``fit2x.CreateMParam``
function from numpy arrays. The return value of a target function is the score
of the model parameters

Model parameters can be optimized to the data by fit functions for fit 23 the
fit function ``fit2x.fit23`` is used.

.. literalinclude:: ../test/test_fit23.py
   :language: python
   :lines: 250-253
   :linenos:

The fit functions takes like the target function an object of the type
``fit2x.MParam`` in addition to the initial values, and a list of fixed model
parameters as an input. The array containing the initial values of the model
parameters are modified in-place buy the fit function.

Alternatively, a simplified python interface can be used to optimize a set of
model as illustrated by the source code embedded in the plot below. The simplified
interface handles the creation of auxiliary data structures such as ``fit2x.MParam``.

.. plot:: ../examples/fluorescence_decay/plot_fit23_mini_example.py

    Analysis result of fit23 by the simplified python interface profided
    by ``fit2x``

In the example shown above, first a fit object of the type ``fit2x.Fit23`` is
created. All necessary data except for the experimental data for a fit is passed
to the fit object when it is created. To perform a fit on experimental data for
a set for a set of initial values, the fit object is called using the inital values
and the data as parameters.

Description of fits
===================
fit23
-----
fit23 optimizes a single rotational correlation time :math:`\rho` and a
fluorescence lifetime :math:`\tau` to a polarization resolved fluorescence decay
considering the fraction of scattered light and instrument response function in
the two detection channels for the parallel and perpendicular fluorescence. Fit23
operates on fluorescence decays in the :term:`Jordi-format`.

.. plot:: ../examples/fluorescence_decay/plot_fit23_1.py

    Simulation and analysis of low photon count data by ``fit2x.fit23``.


fit23 is intended to be used for data with very few photons, e.g. for pixel analysis
in fluorescence lifetime image microscopy (FLIM) or for single-molecule spectroscopy.
The fit implements a maximum likelihood estimimator as previously described :cite:`maus_experimental_2001`.
Briefly, the MLE fit quality parameter 2I* = :math:`-2\ln L(n,g)` (where :math:`L`
is the likelihood function, :math:`n` are the experimental counts, and :math:`g`
is the model function) is minimized. The model function :math:`g` under magic-angle
is given by:

:math:`g_i=N_g \left[ (1-\gamma) \frac{irf_i \ast \exp(iT/k\tau) + c}{\sum_{i=0}^{k}irf_i \ast \exp(iT/k\tau) + c} + \gamma \frac{bg_i}{\sum_i^{k} bg_i} \right]`

:math:`N_e` is not a fitting parameter but set to the experimental number of
photons :math:`N`, :math:`\ast` is the convolution operation, :math:`\tau` is the
fluorescence lifetime, :math:`irf` is the instrument response function, :math:`i`
is the channel number, :math:`bg_i` is the background count in the channel :math:`i`,

The convolution by fit23 is computed recusively and accounts for high repetition
rates:

:math:`irf_i \ast \exp(iT/k\tau) = \sum_{j=1}^{min(i,l)}irf_j\exp(-(i-j)T/k\tau) + \sum_{j=i+1}^{k}irf_j\exp(-(i+k-j)T/k\tau) `

The anisotropy treated as previously described :cite:`schaffer_identification_1999`.
The correction factors needed for a correct anisotropy used by ``fit2x`` are
defined in the glossary (:term:`Anisotropy`).


fit24
-----
Fit24 optimizes is a bi-exponential model function with two fluorescence lifetimes
:math:`\tau_1`, :math:`\tau_2`, and amplitude of the second lifetime :math:`a_2`,
the fraction scattered light :math:`\gamma`, and a constant offset to experimental
data. Fit24 does not describe anisotropy. The decays passed to the fit in the
:term:`Jordi` format. The two decays in the Jordi array are both treated by the
same model function and convoluted with the corresponding
:term:`instrument response function<IRF>`. The model function is

.. :math:

    M_i =(1-a_2) \exp(i\Delta t/\tau_1) + a_2 \exp(i\Delta t/\tau_2)

where :math:`\Delta t` is the time per micro time channel, :math:`i` is the micro time
channel number, :math:`a_2` is the fraction of the second species. The model function
is corrected for the fraction of background and a constant offset.

.. :math:

    g_i =(1 - \gamma) \frac{M_i}{\sum_i M_i} + \gamma \frac{B_i}{\sum_i B_i} + c

Where,  :math:`c` is a constant offset, :math:`B` the background pattern, :math:`M`
the model function and :math:`\gamma` the fraction of scattered light.

fit25
-----
Selects the lifetime out of a set of 4 fixed lifetimes that best describes the data.
Works with polarization resolved Jordi stacks, computes rotational correlation time
by the anisotropy. This function selects out of a set of 4 lifetimes tau the lifetime
that fits best the data.

fit26
-----
Pattern fit. Determines the fraction :math:`f` of two mixed patterns.

.. :math:

    g_i = f \cdot pattern_{1,i} + (1-f) \cdot pattern_{2,i}

(No convolution of patterns, area of pattern is normalized by fit)