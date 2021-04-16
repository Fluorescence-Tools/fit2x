#!/usr/bin/env bash
doxygen ../include/Doxyfile
python doxy2swig.py ../docs/_build/api/xml/index.xml ../ext/python/documentation.i
