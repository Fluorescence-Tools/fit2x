%{
#include "../include/DecayCurve.h"
%}

//// https://github.com/swig/swig/issues/562
//// Fix double free with inplace operators
//%feature("del","") *::operator+=;
//%feature("new","") *::operator+=;
//%feature("del","") *::operator*=;
//%feature("new","") *::operator*=;

// Input arrays
%apply(double* IN_ARRAY1, int DIM1) {
    (double* curve1, int n_curve1),
    (double* curve2, int n_curve2)
}

%include "../include/DecayCurve.h"
%attribute(
        DecayCurve,
        double, acquisition_time,
        get_acquisition_time, set_acquisition_time
);
%extend DecayCurve{
        %pythoncode "../ext/python/DecayCurve.py"
}
