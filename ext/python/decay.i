%{
#include "../include/decay.h"
#include "tttrlib/tttr.h"
%}
// Input arrays
%apply(double* IN_ARRAY1, int DIM1) {
    (double* time_axis, int n_time_axis),
    (double* squared_weights, int n_weights),
    (double* instrument_response_function, int n_instrument_response_function),
    (double* lifetime_spectrum, int n_lifetime_spectrum),
    (double* data, int n_data),
    (double* curve1, int n_curve1),
    (double* curve2, int n_curve2),
    (double *lifetime_spectrum, int n_lifetime_spectrum)
}

%apply(short* IN_ARRAY1, int DIM1){
    (short* fixed, int n_fixed)
}

// Input output arrays
%apply(double* INPLACE_ARRAY1, int DIM1) {
    (double* model_function, int n_model_function),
    (double* x, int n_x)
}

%attribute(Decay, double, background, get_constant_background, set_constant_background);
%attribute(Decay, double, irf_shift, get_irf_shift_channels, set_irf_shift_channels);
%attribute(Decay, double, irf_background, get_irf_background_counts, set_irf_background_counts);
%attribute(Decay, int, convolution_stop, get_convolution_stop, set_convolution_stop);
%attribute(Decay, int, convolution_start, get_convolution_start, set_convolution_start);
%attribute(Decay, bool, correct_pile_up, get_correct_pile_up, set_correct_pile_up);
%attribute(Decay, double, scatter_fraction, get_areal_scatter_fraction, set_areal_scatter_fraction);



%include "../include/decay.h"

%extend Decay{
    %pythoncode "../ext/python/decay/decay_extension.py"
}

