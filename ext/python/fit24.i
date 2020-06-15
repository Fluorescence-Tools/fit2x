%rename (modelf24) my_modelf24;
%exception my_modelf24{
        $action
        if (PyErr_Occurred()) SWIG_fail;
}
%inline %{
int my_modelf24(int len1, double* param,
                int len2, double* irf,
                int len3, double* bg,
                double dt,
                int len4, double* corrections,
                int len5, double* mfunction)
{
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
    if (len1 != 5) {
        PyErr_Format(PyExc_ValueError,
                     "Parameter array should be of length 5. "
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
    return modelf24(param, irf, bg, len5 / 2, dt, corrections, mfunction);
}
%}


