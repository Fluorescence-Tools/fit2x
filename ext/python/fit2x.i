%module(directors="1", package="fit2x") fit2x
%feature("kwargs", 1);
%{
// This fixes numpy int casting to std::vector,int>
// (see: https://github.com/swig/swig/issues/888)
#define SWIG_PYTHON_CAST_MODE

// This is needed for numpy as you need SWIG_FILE_WITH_INIT
#define SWIG_FILE_WITH_INIT
#include "../include/fits2x.h"
#include "tttrlib/tttr.h"
%}


%include "documentation.i"
%include "exception.i"
%include "std_vector.i";
%include "std_string.i";
%include "std_wstring.i";
%include "cpointer.i"
%include "numpy.i"
%include "attribute.i"
%include "std_shared_ptr.i";
%include <typemaps.i>

// Use shared_prt for TTTR to pass TTTR around
%shared_ptr(TTTR)

%init %{
import_array();
%}

// Generic input arrays
// floating numbers
%apply(double* IN_ARRAY1, int DIM1) {(double *input, int n_input)}
%apply(double* IN_ARRAY2, int DIM1, DIM2) {(double *input, int n_input1, int n_input2)}
// integers
%apply(char* IN_ARRAY1, int DIM1) {(char *input, int n_input)}
%apply(short* IN_ARRAY1, int DIM1) {(short* input, int n_input)}
%apply(unsigned short* IN_ARRAY1, int DIM1) {(unsigned short* input, int n_input)}
%apply(int* IN_ARRAY1, int DIM1) {(int* input, int n_input)}
%apply(long long* IN_ARRAY1, int DIM1) {(long long *input, int n_input)}
%apply(unsigned long long* IN_ARRAY1, int DIM1) {(unsigned long long *input, int n_input)}

// Generic output arrays views
// floating points
%apply(double** ARGOUTVIEW_ARRAY1, int* DIM1) {(double** output_view, int* n_output)}

// Generic output memory managed arrays
// floating points
%apply(double** ARGOUTVIEWM_ARRAY1, int* DIM1) {(double** output, int* n_output)}
%apply(double** ARGOUTVIEWM_ARRAY2, int* DIM1, int* DIM2) {(double** output, int* n_output1, int* n_output2)}
%apply(double** ARGOUTVIEWM_ARRAY1, int* DIM1) {(double** output, int* n_output)}
%apply (double** ARGOUTVIEWM_ARRAY3, int* DIM1, int* DIM2, int* DIM3) {(double** output, int* dim1, int* dim2, int* dim3)}
%apply (float** ARGOUTVIEWM_ARRAY4, int* DIM1, int* DIM2, int* DIM3, int* DIM4) {(float** output, int* dim1, int* dim2, int* dim3, int* dim4)}
// integers
%apply(long long** ARGOUTVIEWM_ARRAY1, int* DIM1) {(long long **output, int *n_output)}
%apply(unsigned long long** ARGOUTVIEWM_ARRAY1, int* DIM1) {(unsigned long long** output, int* n_output)}
%apply(int** ARGOUTVIEWM_ARRAY1, int* DIM1) {(int** output, int* n_output)}
%apply(unsigned int** ARGOUTVIEWM_ARRAY1, int* DIM1) {(unsigned int** output, int* n_output)}
%apply(short** ARGOUTVIEWM_ARRAY1, int* DIM1) {(short** output, int* n_output)}
%apply(unsigned short** ARGOUTVIEWM_ARRAY1, int* DIM1) {(unsigned short** output, int* n_output)}
%apply(char** ARGOUTVIEWM_ARRAY1, int* DIM1) {(char** output, int* n_output)}
%apply(signed char** ARGOUTVIEW_ARRAY1, int* DIM1) {(signed char** output, int* n_output)}
%apply (unsigned int** ARGOUTVIEWM_ARRAY2, int* DIM1, int* DIM2) {(unsigned int** output, int* dim1, int* dim2)}
%apply (unsigned int** ARGOUTVIEWM_ARRAY3, int* DIM1, int* DIM2, int* DIM3) {(unsigned int** output, int* dim1, int* dim2, int* dim3)}
%apply (unsigned char** ARGOUTVIEWM_ARRAY4, int* DIM1, int* DIM2, int* DIM3, int* DIM4) {(unsigned char** output, int* dim1, int* dim2, int* dim3, int* dim4)}

// Generic inplace arrays
%apply(double* INPLACE_ARRAY1, int DIM1) {(double* inplace_output, int n_output)}


// Templates
%template(VectorDouble) std::vector<double>;
%template(VectorInt32) std::vector<int>;

%apply (int DIM1, double* INPLACE_ARRAY1) {(int len1, double* x)}
%apply (int DIM1, short* IN_ARRAY1) {(int len2, short* fixed)}

%include "../include/fits2x.h"
%pythoncode "../ext/python/fit2x.py"

%include "lvarray.i"
%include "fsconv.i"
%include "fit23.i"
%include "fit24.i"
%include "fit25.i"
%include "fit26.i"
%include "decay.i"
%include "phasor.i"
