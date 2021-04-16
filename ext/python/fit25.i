%{
#include "../include/fit25.h"
%}


%rename (targetf25) my_targetf25;
%exception my_targetf25{$action if (PyErr_Occurred()) SWIG_fail;}

%inline %{
double my_targetf25(int len1, double* x, MParam* p){
    if (len1 != 8) {
        PyErr_Format(PyExc_ValueError,
                     "The length of the parameter vector must of length 8. "
                     "Arrays of length (%d) given",
                     len1);
        return 0.0;
    }
    return targetf25(x, p);
}
%}

%rename (fit25) my_fit25;
%inline %{
double my_fit25(int len1, double* x, int len2, short* fixed, MParam* p){
    if (len1 != 9) {
        PyErr_Format(PyExc_ValueError,
                     "The length of the parameter vector must of length 9. "
                     "Arrays of length (%d) given",
                     len1);
        return 0.0;
    }
    if (len2 < 5) {
        PyErr_Format(PyExc_ValueError,
                     "The length of the vector fixed be at least of length 5. "
                     "Array of lengths (%d) given",
                     len2);
        return 0.0;
    }
    return fit25(x, fixed, p);
}
%}

%pythoncode "../ext/python/fit25.py"

