%{
#include "DecayCurve.h"
%}

%apply(double* IN_ARRAY1, int DIM1) {(double* curve1, int n_curve1), (double* curve2, int n_curve2)}
%attribute(DecayCurve, double, acquisition_time, get_acquisition_time, set_acquisition_time);
%attribute(DecayCurve, double, shift, get_shift, set_shift);

%include "DecayCurve.h"
%extend DecayCurve{%pythoncode "../ext/python/DecayCurve.py"}

