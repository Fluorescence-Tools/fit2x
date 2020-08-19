%include "typemaps.i"

// manually added instead of including header file as all other functions
// are also manually added
%apply (double* INPLACE_ARRAY1, int DIM1) {
    (double* fit, int len1),
    (double* model, int len3),
    (double *model, int n_model),
    (double *time_axis, int n_time_axis),
    (double *instrument_response_function, int n_instrument_response_function),
    (double *lifetime_spectrum, int n_lifetime_spectrum),
    (double *data, int n_data),
    (double* w_sq, int len3),
    (double* x, int len3),
    (double* decay, int len2),
    (double* irf, int len2),
    (double* lamp, int len1), // used by shift_lamp
    (double* lampsh, int len2) // used by shift_lamp
}

void fconv_per_cs_time_axis(
        double *model, int n_model,
        double *time_axis, int n_time_axis,
        double *instrument_response_function, int n_instrument_response_function,
        double *lifetime_spectrum, int n_lifetime_spectrum,
        int convolution_start = 0,
        int convolution_stop = -1,
        double period = 100.0
);

void add_pile_up_to_model(
        double* model, int n_model,
        double* data, int n_data,
        double repetition_rate,
        double dead_time,
        double measurement_time,
        std::string pile_up_model
);

//// rescale
//////////////

%rename (rescale) my_rescale;
%exception my_rescale{$action if (PyErr_Occurred()) SWIG_fail;}
%inline %{
double my_rescale(
        double* fit, int len1,
        double* decay, int len2,
        int start = 0,
        int stop = -1
){
    if (len1 != len2) {
        PyErr_Format(PyExc_ValueError,
                     "Model and decay array should have same length. "
                     "Arrays of lengths (%d,%d) given",
                     len1, len2);
    }
    if(start < 0){
        PyErr_Format(PyExc_ValueError,
                     "Start index needs to be larger or equal to zero."
        );
    }
    if(stop < 0){
        stop = len1;
    }
    if (start >= len1) {
        PyErr_Format(PyExc_ValueError,
                     "Start index (%d) too large for array of lengths (%d).",
                     start, len1);
    }
    if (stop > len1) {
        PyErr_Format(PyExc_ValueError,
                     "Stop index (%d) too large for array of lengths (%d).",
                     stop, len1);
    }
    double scale = 0.0;
    rescale(fit, decay, &scale, start, stop);
    return scale;
}
%}

//// rescale_w
////////////////

%rename (rescale_w) my_rescale_w;
%exception my_rescale{$action if (PyErr_Occurred()) SWIG_fail;}
%inline %{
double my_rescale_w(
        double* fit, int len1,
        double* decay, int len2,
        double* w_sq, int len3,
        int start = 0,
        int stop = -1
){
    if (len1 != len2) {
        PyErr_Format(PyExc_ValueError,
                     "Model and decay array should have same length. "
                     "Arrays of lengths (%d,%d) given",
                     len1, len2);
    }
    if (len3 != len2) {
        PyErr_Format(PyExc_ValueError,
                     "Weight and decay array should have same length. "
                     "Arrays of lengths (%d,%d) given",
                     len2, len3);
    }
    if(start < 0){
        PyErr_Format(PyExc_ValueError,
                     "Start index needs to be larger or equal to zero."
        );
    }
    if(stop < 0){
        stop = len1;
    }
    if (start > len1) {
        PyErr_Format(PyExc_ValueError,
                     "Start index (%d) too large for array of lengths (%d).",
                     start, len1);
    }
    if (stop > len1) {
        PyErr_Format(PyExc_ValueError,
                     "Stop index (%d) too large for array of lengths (%d).",
                     stop, len1);
    }
    double scale = 0.0;
    rescale_w(fit, decay, w_sq, &scale, start, stop);
    return scale;
}
%}

/// rescale_w_bg
///////////////////
%rename (rescale_w_bg) my_rescale_w_bg;
%exception my_rescale_w_bg{$action if (PyErr_Occurred()) SWIG_fail;}
%inline %{
double my_rescale_w_bg(
        double* fit, int len1,
        double* decay, int len2,
        double* w_sq, int len3,
        double bg = 0.0,
        int start = 0,
        int stop = -1
){
    if (len1 != len2) {
        PyErr_Format(PyExc_ValueError,
                     "Model and decay array should have same length. "
                     "Arrays of lengths (%d,%d) given",
                     len1, len2);
    }
    if (len3 != len2) {
        PyErr_Format(PyExc_ValueError,
                     "Weight and decay array should have same length. "
                     "Arrays of lengths (%d,%d) given",
                     len2, len3);
    }
    if(start < 0){
        PyErr_Format(PyExc_ValueError,
                     "Start index needs to be larger or equal to zero."
        );
    }
    if(stop < 0){
        stop = len1;
    }
    if (start > len1) {
        PyErr_Format(PyExc_ValueError,
                     "Start index (%d) too large for array of lengths (%d).",
                     start, len1);
    }
    if (stop > len1) {
        PyErr_Format(PyExc_ValueError,
                     "Stop index (%d) too large for array of lengths (%d).",
                     stop, len1);
    }
    double scale = 0.0;
    rescale_w_bg(fit, decay, w_sq, bg, &scale, start, stop);
    return scale;
}
%}

//// fconv
///////////////////
%rename (fconv) my_fconv;
%exception my_fconv{$action if (PyErr_Occurred()) SWIG_fail;}
%inline %{
void my_fconv(
        double* fit, int len1,
        double* irf, int len2,
        double* x, int len3,
        int start = 0,
        int stop = -1,
        double dt = 1.0
){
    if (len1 != len2) {
        PyErr_Format(PyExc_ValueError,
                     "Model and decay array should have same length. "
                     "Arrays of lengths (%d,%d) given",
                     len1, len2);
    }
    if(start < 0){
        PyErr_Format(PyExc_ValueError,
                     "Start index needs to be larger or equal to zero."
        );
    }
    if(stop < 0){
        stop = len1;
    }
    if (start > len1) {
        PyErr_Format(PyExc_ValueError,
                     "Start index (%d) too large for array of lengths (%d).",
                     start, len1);
    }
    if (stop > len1) {
        PyErr_Format(PyExc_ValueError,
                     "Stop index (%d) too large for array of lengths (%d).",
                     stop, len1);
    }
    fconv(fit, x, irf, len3 / 2, start, stop, dt);
}
%}


//// fconv_avx
///////////////////
%rename (fconv_avx) my_fconv_avx;
%exception my_fconv_avx{$action if (PyErr_Occurred()) SWIG_fail;}
%inline %{
void my_fconv_avx(
        double* fit, int len1,
        double* irf, int len2,
        double* x, int len3,
        int start,
        int stop,
        double dt
){
    if (len1 != len2) {
        PyErr_Format(PyExc_ValueError,
                     "Model and decay array should have same length. "
                     "Arrays of lengths (%d,%d) given",
                     len1, len2);
    }
    if(start < 0){
        PyErr_Format(PyExc_ValueError,
                     "Start index needs to be larger or equal to zero."
        );
    }
    if(stop < 0){
        stop = len1;
    }
    if (start > len1) {
        PyErr_Format(PyExc_ValueError,
                     "Start index (%d) too large for array of lengths (%d).",
                     start, len1);
    }
    if (stop > len1) {
        PyErr_Format(PyExc_ValueError,
                     "Stop index (%d) too large for array of lengths (%d).",
                     stop, len1);
    }
    fconv_avx(fit, x, irf, len3 / 2, start, stop, dt);
}
%}


//// fconv_per
///////////////////
%rename (fconv_per) my_fconv_per;
%exception my_fconv_per{$action if (PyErr_Occurred()) SWIG_fail;}
%inline %{
void my_fconv_per(
        double* fit, int len1,
        double* irf, int len2,
        double* x, int len3,
        double period,
        int start = 0,
        int stop = -1,
        double dt = 1.0
){
    if (len1 != len2) {
        PyErr_Format(PyExc_ValueError,
                     "Model and decay array should have same length. "
                     "Arrays of lengths (%d,%d) given",
                     len1, len2);
    }
    if(start < 0){
        PyErr_Format(PyExc_ValueError,
                     "Start index needs to be larger or equal to zero."
        );
    }
    if(stop < 0){
        stop = len1;
    }
    if (start > len1) {
        PyErr_Format(PyExc_ValueError,
                     "Start index (%d) too large for array of lengths (%d).",
                     start, len1);
    }
    if (stop > len1) {
        PyErr_Format(PyExc_ValueError,
                     "Stop index (%d) too large for array of lengths (%d).",
                     stop, len1);
    }
    fconv_per(fit, x, irf, len3 / 2, start, stop, len1, period, dt);
}
%}

//// fconv_per_avx
///////////////////
%rename (fconv_per_avx) my_fconv_per_avx;
%exception my_fconv_per_avx{$action if (PyErr_Occurred()) SWIG_fail;}
%inline %{
void my_fconv_per_avx(
        double* fit, int len1,
        double* irf, int len2,
        double* x, int len3,
        double period,
        int start = 0,
        int stop = -1,
        double dt = 1.0
){
    if (len1 != len2) {
        PyErr_Format(PyExc_ValueError,
                     "Model and decay array should have same length. "
                     "Arrays of lengths (%d,%d) given",
                     len1, len2);
    }
    if(start < 0){
        PyErr_Format(PyExc_ValueError,
                     "Start index needs to be larger or equal to zero."
        );
    }
    if(stop < 0){
        stop = len1;
    }
    if (start > len1) {
        PyErr_Format(PyExc_ValueError,
                     "Start index (%d) too large for array of lengths (%d).",
                     start, len1);
    }
    if (stop > len1) {
        PyErr_Format(PyExc_ValueError,
                     "Stop index (%d) too large for array of lengths (%d).",
                     stop, len1);
    }
    fconv_per_avx(fit, x, irf, len3 / 2, start, stop, len1, period, dt);
}
%}


//// fconv_per_cs
///////////////////
%rename (fconv_per_cs) my_fconv_per_cs;
%exception my_fconv_per_cs{$action if (PyErr_Occurred()) SWIG_fail;}
%inline %{
void my_fconv_per_cs(
        double* fit, int len1,
        double* irf, int len2,
        double* x, int len3,
        double period,
        int conv_stop = -1,
        int stop = -1,
        double dt = 1.0
){
    if (len1 != len2) {
        PyErr_Format(PyExc_ValueError,
                     "Model and decay array should have same length. "
                     "Arrays of lengths (%d,%d) given",
                     len1, len2);
    }
    if(stop < 0){
        stop = len1;
    }
    if(conv_stop < 0){
        conv_stop = len1;
    }
    if (stop > len1) {
        PyErr_Format(PyExc_ValueError,
                     "Stop index (%d) too large for array of lengths (%d).",
                     stop, len1);
    }
    fconv_per_cs(fit, x, irf, len3 / 2, stop, len1, period, conv_stop);
}
%}

//// fconv_ref
///////////////////
%rename (fconv_ref) my_fconv_ref;
%exception my_fconv_ref{$action if (PyErr_Occurred()) SWIG_fail;}
%inline %{
void my_fconv_ref(
        double* fit, int len1,
        double* irf, int len2,
        double* x, int len3,
        double tauref,
        int start = 0,
        int stop = -1,
        double dt = 1.0
){
    if (len1 != len2) {
        PyErr_Format(PyExc_ValueError,
                     "Model and decay array should have same length. "
                     "Arrays of lengths (%d,%d) given",
                     len1, len2);
    }
    if(start < 0){
        PyErr_Format(PyExc_ValueError,
                     "Start index needs to be larger or equal to zero."
        );
    }
    if(stop < 0){
        stop = len1;
    }
    if (start > len1) {
        PyErr_Format(PyExc_ValueError,
                     "Start index (%d) too large for array of lengths (%d).",
                     start, len1);
    }
    if (stop > len1) {
        PyErr_Format(PyExc_ValueError,
                     "Stop index (%d) too large for array of lengths (%d).",
                     stop, len1);
    }
    fconv_ref(fit, x, irf, len3 / 2, start, stop, tauref);
}
%}

//// sconv
///////////////////
%rename (sconv) my_sconv;
%exception my_sconv{$action if (PyErr_Occurred()) SWIG_fail;}
%inline %{
void my_sconv(
        double* fit, int len1,
        double* irf, int len2,
        double* model, int len3,
        int start = 0,
        int stop = -1
){
    if (len1 != len2) {
        PyErr_Format(PyExc_ValueError,
                     "Model and decay array should have same length. "
                     "Arrays of lengths (%d,%d) given",
                     len1, len2);
    }
    if (len3 != len2) {
        PyErr_Format(PyExc_ValueError,
                     "Model and fit array should have same length. "
                     "Arrays of lengths (%d,%d) given",
                     len3, len2);
    }
    if(start < 0){
        PyErr_Format(PyExc_ValueError,
                     "Start index needs to be larger or equal to zero."
        );
    }
    if(stop < 0){
        stop = len1;
    }
    if (start > len1) {
        PyErr_Format(PyExc_ValueError,
                     "Start index (%d) too large for array of lengths (%d).",
                     start, len1);
    }
    if (stop > len1) {
        PyErr_Format(PyExc_ValueError,
                     "Stop index (%d) too large for array of lengths (%d).",
                     stop, len1);
    }
    sconv(fit, model, irf, start, stop);
}
%}

//// shift_lamp
///////////////////
%rename (shift_lamp) my_shift_lamp;
%exception my_shift_lamp{$action if (PyErr_Occurred()) SWIG_fail;}
%inline %{
void my_shift_lamp(
        double* lamp, int len1,
        double* lampsh, int len2,
        double ts
){
    if (len1 != len2) {
        PyErr_Format(PyExc_ValueError,
                     "IRF and shifted IRF array should have same length. "
                     "Arrays of lengths (%d,%d) given",
                     len1, len2);
    }
    shift_lamp(lampsh, lamp, ts, len1);
}
%}

