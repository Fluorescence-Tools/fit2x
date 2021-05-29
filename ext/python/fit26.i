%{
#include "fit26.h"
%}

%rename (targetf26) my_targetf26;
%exception my_targetf26{$action if (PyErr_Occurred()) SWIG_fail;}

%inline %{
double my_targetf26(int len1, double* x, MParam* p){
    if (len1 != 8) {
        PyErr_Format(PyExc_ValueError,
                     "The length of the parameter vector must of length 8. "
                     "Arrays of length (%d) given",
                     len1);
        return 0.0;
    }
    return targetf26(x, p);
}
%}

%rename (fit26) my_fit26;
%inline %{
double my_fit26(int len1, double* x, int len2, short* fixed, MParam* p){
    if (len1 != 1) {
        PyErr_Format(PyExc_ValueError,
                     "The length of the parameter vector must of length 1. "
                     "Arrays of length (%d) given",
                     len1);
        return 0.0;
    }
    if (len2 < 1) {
        PyErr_Format(PyExc_ValueError,
                     "The length of the vector fixed be at least of length 1. "
                     "Array of lengths (%d) given",
                     len2);
        return 0.0;
    }
    return fit26(x, fixed, p);
}
%}

%pythoncode "../ext/python/fit26.py"
