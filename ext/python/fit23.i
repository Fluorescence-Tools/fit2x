
%exception my_modelf23{$action if (PyErr_Occurred()) SWIG_fail;}
%rename (modelf23) my_modelf23;
%exception my_modelf23{
        $action
        if (PyErr_Occurred()) SWIG_fail;
}


%inline %{
int my_modelf23(
        int len1, double* param,
        int len2, double* irf,
        int len3, double* bg,
        double dt,
        int len4, double* corrections,
        int len5, double* mfunction
){
    if (len2 != len3) {
        PyErr_Format(PyExc_ValueError,
                     "IRF and Bg array should have same length. "
                     "Arrays of lengths (%d,%d) given",
                     len2, len3);
        return 0.0;
    }
    if (len5 != len3) {
        PyErr_Format(PyExc_ValueError,
                     "Output array should be of length inputs. "
                     "Arrays of lengths (%d,%d) given",
                     len5, len3);
        return 0.0;
    }
    if (len1 != 4) {
        PyErr_Format(PyExc_ValueError,
                     "Parameter array should be of length 4. "
                     "Arrays of length (%d) given",
                     len1);
        return 0.0;
    }
    if (len4 != 5) {
        PyErr_Format(PyExc_ValueError,
                     "Corrections array should be of length 4. "
                     "Arrays of length (%d) given",
                     len1);
        return 0.0;
    }
    return modelf23(param, irf, bg, len5 / 2, dt, corrections, mfunction);
}
%}

%rename (targetf23) my_targetf23;
%exception my_targetf23{$action if (PyErr_Occurred()) SWIG_fail;}
%inline %{
double my_targetf23(int len1, double* x, MParam* p){
    if (len1 != 8) {
        PyErr_Format(PyExc_ValueError,
                     "The length of the parameter vector must of length 8. "
                     "Arrays of length (%d) given",
                     len1);
        return 0.0;
    }
    return targetf23(x, p);
}
%}

%rename (fit23) my_fit23;
%inline %{
double my_fit23(int len1, double* x, int len2, short* fixed, MParam* p){
    if (len1 != 8) {
        PyErr_Format(PyExc_ValueError,
                     "The length of the parameter vector must of length 8"
                     "Arrays of length (%d) given",
                     len1);
        return 0.0;
    }
    if (len2 < 4) {
        PyErr_Format(PyExc_ValueError,
                     "The length of the vector fixed be at least of length 6"
                     "Arrays of lengths (%d) given",
                     len2);
        return 0.0;
    }
    return fit23(x, fixed, p);
}
%}


%pythoncode "../ext/python/fit23.py"
