%{
#include "../include/decay.h"
%}

// Input arrays
%apply(double* IN_ARRAY1, int DIM1) {
    (double* time_axis, int n_time_axis),
    (double* squared_data_weights, int n_squared_data_weights),
    (double* scatter, int n_scatter),
    (double* weights, int n_weights),
    (double* data_weights, int n_data_weights),
    (double* irf_histogram, int n_irf_histogram),
    (double* linearization, int n_linearization),
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

%include "../include/decay.h"

%attribute(Decay, double, constant_offset, get_constant_offset, set_constant_offset);
%attribute(Decay, double, irf_shift_channels, get_irf_shift_channels, set_irf_shift_channels);
%attribute(Decay, double, irf_background_counts, get_irf_background_counts, set_irf_background_counts);
%attribute(Decay, int, convolution_stop, get_convolution_stop, set_convolution_stop);
%attribute(Decay, int, convolution_start, get_convolution_start, set_convolution_start);
%attribute(Decay, bool, use_pile_up_correction, get_use_pile_up_correction, set_use_pile_up_correction);
%attribute(Decay, bool, abs_lifetime_spectrum, get_abs_lifetime_spectrum, set_abs_lifetime_spectrum);
%attribute(Decay, bool, use_amplitude_threshold, get_use_amplitude_threshold, set_use_amplitude_threshold);
%attribute(Decay, double, amplitude_threshold, get_amplitude_threshold, set_amplitude_threshold);
%attribute(Decay, double, scatter_fraction, get_scatter_fraction, set_scatter_fraction);
%attribute(Decay, bool, scale_model_to_data, get_scale_model_to_data, set_scale_model_to_data);
%attribute(Decay, bool, use_corrected_irf_as_scatter, get_use_corrected_irf_as_scatter, set_use_corrected_irf_as_scatter);
%attribute(Decay, bool, use_linearization, get_use_linearization, set_use_linearization);
%attribute(Decay, bool, is_valid, get_is_valid);
%attribute(Decay, int, convolution_method, get_convolution_method, set_convolution_method);

%attribute(Decay, double, acquisition_time, get_acquisition_time, set_acquisition_time);
%attribute(Decay, double, instrument_dead_time, get_instrument_dead_time, set_instrument_dead_time);
%attribute(Decay, double, number_of_photons, get_number_of_photons, set_number_of_photons);
%attribute(Decay, double, excitation_period, get_excitation_period, set_excitation_period);

%extend Decay{
    %pythoncode "../ext/python/decay_extension.py"
}

