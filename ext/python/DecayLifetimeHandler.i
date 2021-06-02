%{
#include "DecayLifetimeHandler.h"
%}

%rename(get_lifetime_spectrum_vec) get_lifetime_spectrum();

//// Input arrays
%attribute(DecayLifetimeHandler, bool, abs_lifetime_spectrum, get_abs_lifetime_spectrum, set_abs_lifetime_spectrum);
%attribute(DecayLifetimeHandler, bool, use_amplitude_threshold, get_use_amplitude_threshold, set_use_amplitude_threshold);
%attribute(DecayLifetimeHandler, double, amplitude_threshold, get_amplitude_threshold, set_amplitude_threshold);

%include "DecayLifetimeHandler.h"
%extend DecayLifetimeHandler{%pythoncode "../ext/python/DecayLifetimeHandler.py"}
