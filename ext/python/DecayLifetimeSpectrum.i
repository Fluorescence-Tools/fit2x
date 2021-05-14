%{
#include "../include/DecayLifetimeSpectrum.h"
%}

//// Input arrays
%apply(double* IN_ARRAY1, int DIM1) {
    (double *lifetime_spectrum, int n_lifetime_spectrum)
}
%include "../include/DecayLifetimeSpectrum.h"
%shared_ptr(DecayLifetimeSpectrum);

%attribute(DecayLifetimeSpectrum, bool, abs_lifetime_spectrum, get_abs_lifetime_spectrum, set_abs_lifetime_spectrum);
%attribute(DecayLifetimeSpectrum, bool, use_amplitude_threshold, get_use_amplitude_threshold, set_use_amplitude_threshold);
%attribute(DecayLifetimeSpectrum, double, amplitude_threshold, get_amplitude_threshold, set_amplitude_threshold);

%extend DecayLifetimeSpectrum{
        %pythoncode "../ext/python/DecayLifetimeSpectrum.py"
}
