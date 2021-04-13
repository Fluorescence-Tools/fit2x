%feature("kwargs", 1);
%{
#include "../include/DecayConvolution.h"
%}

%include "../include/DecayConvolution.h"
%attribute(DecayConvolution, bool, use_corrected_irf_as_scatter, get_use_corrected_irf_as_scatter, set_use_corrected_irf_as_scatter);
%attribute(DecayConvolution, double, scatter_fraction, get_scatter_fraction, set_scatter_fraction);
%attribute(DecayConvolution, int, convolution_stop, get_convolution_stop, set_convolution_stop);
%attribute(DecayConvolution, int, convolution_start, get_convolution_start, set_convolution_start);
%attribute(DecayConvolution, int, convolution_method, get_convolution_method, set_convolution_method);
%attribute(DecayConvolution, double, excitation_period, get_excitation_period, set_excitation_period);
%attribute(DecayConvolution, double, irf_shift_channels, get_irf_shift_channels, set_irf_shift_channels);
%attribute(DecayConvolution, double, irf_background_counts, get_irf_background_counts, set_irf_background_counts);

