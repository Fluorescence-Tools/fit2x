%{
#include "../include/Decay.h"
%}

//// Input arrays
%apply(double* IN_ARRAY1, int DIM1) {
    (double *lifetime_spectrum, int n_lifetime_spectrum),
    (double *lifetime_spectrum = nullptr, int n_lifetime_spectrum = -1)
}

%include "../include/Decay.h"

%attribute(Decay, bool, is_valid, get_is_valid);

// Convolution
%attribute(Decay, int, convolution_method, get_convolution_method, set_convolution_method);
%attribute(Decay, double, excitation_period, get_excitation_period, set_excitation_period);
%attribute(Decay, double, irf_shift_channels, get_irf_shift_channels, set_irf_shift_channels);
%attribute(Decay, double, scatter_fraction, get_scatter_fraction, set_scatter_fraction);
%attribute(Decay, int, convolution_start, get_convolution_start, set_convolution_start);
%attribute(Decay, int, convolution_stop, get_convolution_stop, set_convolution_stop);
%attribute(Decay, double, irf_background_counts, get_irf_background_counts, set_irf_background_counts);
%attribute(Decay, bool, use_corrected_irf_as_scatter, get_use_corrected_irf_as_scatter, set_use_corrected_irf_as_scatter);

// Lifetime spectrum
%attribute(Decay, bool, abs_lifetime_spectrum, get_abs_lifetime_spectrum, set_abs_lifetime_spectrum);
%attribute(Decay, bool, use_amplitude_threshold, get_use_amplitude_threshold, set_use_amplitude_threshold);
%attribute(Decay, double, amplitude_threshold, get_amplitude_threshold, set_amplitude_threshold);

// Score
// TODO: score_type

// Background
%attribute(Decay, double, constant_offset, get_constant_offset, set_constant_offset);

// Pileup
%attribute(Decay, bool, use_pile_up_correction, get_use_pile_up_correction, set_use_pile_up_correction);
%attribute(Decay, double, instrument_dead_time, get_instrument_dead_time, set_instrument_dead_time);
%attribute(Decay, double, repetition_rate, get_repetition_rate, set_repetition_rate);
%attributestring(Decay, std::string, pile_up_model, get_pile_up_model, set_pile_up_model);

// Experimental data
%attribute(Decay, double, acquisition_time, get_acquisition_time, set_acquisition_time);

// Scale
%attribute(Decay, bool, scale_model_to_data, get_scale_model_to_data, set_scale_model_to_data);
%attribute(Decay, double, number_of_photons, get_number_of_photons, set_number_of_photons);
// scale_start/stop set by convolution_start/stop

// Linearization
%attribute(Decay, bool, use_linearization, get_use_linearization, set_use_linearization);

// Scoring
//%attributeval(Decay, std::vector<int>, score_range, get_score_range, set_score_range);
%attribute(Decay, double, score, get_score);
%attributestring(Decay, std::string, score_type, get_score_type, set_score_type);

%extend Decay{
    %pythoncode "../ext/python/Decay.py"
}

