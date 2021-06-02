#####
fit2x
#####

General description
===================
fits2x is a collection of models and fits that use maximum likelihood
estimators for polarization and time-resolved fluorescence decays
(see :cite:`maus_experimental_2001`).

fit2x implement burst integrated fluorescence lifetime fits (BIFL)
with scatter. The library can be used in conjuncture with tttrlib to
analyze and  process single-molecule FRET (smFRET) experiments and
confocal laser  scanning microscopy (CLSM) data. The fit2x shared libary
can be used from LabView and is wrapped by SWIG (Simplified Wrapper and
Interface Generator) for common scripting languages as Python as main
target language.

## Design goals
*   Low memory footprint (keep objective large datasets, e.g.  FLIM in memory).
*   Platform independent C/C++ library with interfaces for scripting libraries

## Capabilities
*   Polarization and time-resolved analysis
*   Stable analysis results with minimum photon counts
*   Robust thoroughly tested maximum likelihood estimators


Documentation
=============
Documentation for the latest releases is available `here <https://fluorescence-tools.github.io/fit2x/>`_,.
Previous releases are available through the read-the-docs page
`page <https://fit2x.readthedocs.io/>`_,.


Contents
========

.. toctree::
   :maxdepth: 1

   installation
   introduction
   auto_examples/index
   cpp_api
   glossary
   zreferences


Indices and tables
==================

* :ref:`genindex`
* :ref:`search`


License
=======
fit2x was developed at the Seidel Lab (Heinrich Heine University) and is maintained
by Thomas Peulen. fit2x is released under the open source `MIT license <https://opensource.org/licenses/MIT>`_.
