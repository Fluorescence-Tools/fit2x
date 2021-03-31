%module(directors="1", package="fit2x") fit2x
%feature("kwargs", 1);
%{
// This fixes numpy int casting to std::vector,int>
// (see: https://github.com/swig/swig/issues/888)
#define SWIG_PYTHON_CAST_MODE
// This is needed for numpy as you need SWIG_FILE_WITH_INIT
#define SWIG_FILE_WITH_INIT
#include <assert.h>
#include "../include/fits2x.h"
%}

%include <typemaps.i>
%include <std_shared_ptr.i>
%include <cpointer.i>
%include <std_vector.i>
%include <attribute.i>
%include "exception.i"
%include "numpy.i"
%include "documentation.i"


%init %{
import_array();
%}

// Python code that should be included at the begining (import, Base class, etc)
%pythoncode "../ext/python/fit2x.py"

// Templates
%template(VectorDouble) std::vector<double>;
%template(VectorInt32) std::vector<int>;

// shared_prt
%shared_ptr(TTTR) // to pass TTTR around
%shared_ptr(std::string)

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
%apply(float** ARGOUTVIEW_ARRAY1, int* DIM1, int* DIM2, int* DIM3, int* DIM4) {(float **output, int *dim1, int *dim2, int *dim3, int *dim4)}

// Generic output memory managed arrays
// float and double
%apply(double** ARGOUTVIEWM_ARRAY1, int* DIM1) {(double** output, int* n_output)}
%apply(double** ARGOUTVIEWM_ARRAY2, int* DIM1, int* DIM2) {(double** output, int* n_output1, int* n_output2)}
%apply(double** ARGOUTVIEWM_ARRAY1, int* DIM1) {(double** output, int* n_output)}
%apply (double** ARGOUTVIEWM_ARRAY3, int* DIM1, int* DIM2, int* DIM3) {(double** output, int* dim1, int* dim2, int* dim3)}
%apply (float** ARGOUTVIEWM_ARRAY4, int* DIM1, int* DIM2, int* DIM3, int* DIM4) {(float **output, int *dim1, int *dim2, int *dim3, int *dim4)}

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

// for fit23, fit24, etc.
%apply (int DIM1, double* INPLACE_ARRAY1) {(int len1, double* x)}
%apply (int DIM1, short* IN_ARRAY1) {(int len2, short* fixed)}


/* Convolution and LabView interface*/
%include "lvarray.i"
%include "fsconv.i"
%include "phasor.i"

/* Fits and Decay*/
%include "fit23.i"
%include "fit24.i"
%include "fit25.i"
%include "fit26.i"
%include "decay.i"

%include "../include/fits2x.h"
