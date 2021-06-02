%{
#include "DecayConvolution.h"
%}

%attribute(DecayConvolution, DecayCurve, corrected_irf, get_corrected_irf);
%attribute(DecayConvolution, int, convolution_method, get_convolution_method, set_convolution_method);
%attribute(DecayConvolution, double, excitation_period, get_excitation_period, set_excitation_period);
%attribute(DecayConvolution, double, irf_shift_channels, get_irf_shift_channels, set_irf_shift_channels);
%attribute(DecayConvolution, double, irf_background_counts, get_irf_background_counts, set_irf_background_counts);

%include "DecayConvolution.h"
