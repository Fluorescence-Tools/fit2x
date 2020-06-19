%module(directors="1", package="fit2x") fit2x
%feature("kwargs", 1);
%{
// This fixes numpy int casting to std::vector,int>
// (see: https://github.com/swig/swig/issues/888)
#define SWIG_PYTHON_CAST_MODE

// This is needed for numpy as you need SWIG_FILE_WITH_INIT
#define SWIG_FILE_WITH_INIT
#include "../include/fits2x.h"
%}

%include "documentation.i"
%include "exception.i"
%include "std_vector.i";
%include "numpy.i"

%init %{
import_array();
%}

// Templates
%template(VectorDouble) std::vector<double>;
%template(VectorInt32) std::vector<int>;

%apply (int DIM1, double* INPLACE_ARRAY1) {(int len5, double* mfunction)}
%apply (int DIM1, double* IN_ARRAY1)
{
    (int len1, double* param),
    (int len2, double* irf),
    (int len3, double* bg),
    (int len4, double* corrections)
}
%apply (int DIM1, double* INPLACE_ARRAY1) {(int len1, double* x)}
%apply (int DIM1, short* IN_ARRAY1) {(int len2, short* fixed)}

%include "../include/fits2x.h"
%pythoncode "../ext/python/fit2x.py"
%include "lvarray.i"

%include "fit23.i"
%include "fit24.i"
%include "fit25.i"
%include "fit26.i"

