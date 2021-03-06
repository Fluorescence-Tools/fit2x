************
Installation
************

``fit2x`` can either be installed from prebuilt binaries or from the source code.
For most users it is recommended to install ``fit2x`` using the prebuilt binaries.
Below the installation using prebuilt binaries and the prerequisites to compile
``fit2x`` are briefly outlined.


Prebuilt binaries
=================
It is recommended to install ``fit2x`` for Python environments via ``conda`` by.

.. code-block:: bash

    conda install -c tpeulen fit2x


Alternatively, ``fit2x`` can be installed via ``pip``.

.. code-block:: bash

    pip install fit2x


Compliation
===========

``fit2x`` can be compiled and installed using the source code provided in the
git repository.

.. code-block:: bash

    git clone --recursive https://github.com/tpeulen/fit2x.git
    cd fit2x
    sudo python setup.py install

To compile ``fit2x`` a set of prerequisites need to be fulfilled:

    1. An installed compiler.
    2. The `HDF5 <https://www.hdfgroup.org/>`_ library with C/C++ include files.
    3. A recent 64bit `Python <https://www.python.org/>`_ installation with include files.
    4. `cmake <https://cmake.org/>`_
    5. `SWIG <http://www.swig.org/>`_

To debug and analyze the ouput of ``fit2x`` in more detail it can be useful to
build the library in Debug mode that outputs additional information to the standard
standard logging stream. fit2x can be build in debugging mode by.

.. code-block:: bash

    python setup.py build_ext --debug


Windows
-------
On windows ``fit2x`` is best compiled with the `Visual Studio 2017 <https://visualstudio.microsoft.com/>`_. For
compilation the Visual Studion Community edition is sufficient. In addition to
Visual Studio the libraries and the include files as listed above need to be
installed. The prebuilt binaries are compiled on Windows 10 with using 64bit anaconda
Python environments `miniconda <https://docs.conda.io/en/latest/miniconda.html>`_
using the conda build recipe that is provided with the source code in the ``conda-recipe``
folder.

macOS
-----
For MacOS the prebuilt binaries are compiled on MacOS 10.13 with the Apple clang
compiler using a anaconda distribution and the provided ``conda-recipe``.

Linux
-----
The Linux prebuilt binaries are compiled on Ubuntu 18.04 in an anaconda distribution
and the provided ``conda-recipe``.

Conda
-----

A conda recipe is provided in the folder 'conda-recipe' to build the fit2x library with the
`conda build <https://docs.conda.io/projects/conda-build/en/latest/>`_ environment.

To build the library download the recipe, install the conda build package and use
the provided recipe to build the library.

.. code-block:: bash

    conda install conda-build
    conda build conda-recipe

