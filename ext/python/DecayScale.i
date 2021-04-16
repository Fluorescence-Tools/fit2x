%feature("kwargs", 1);
%{
#include "../include/DecayScale.h"
%}

%include "../include/DecayScale.h"
%attribute(DecayScale, double, number_of_photons, get_number_of_photons, set_number_of_photons);

%attribute(DecayScale, int, scale_start, get_scale_start, set_scale_start);
%attribute(DecayScale, int, scale_stop, get_scale_stop, set_scale_stop);
%attribute(DecayScale, bool, scale_model_to_data, get_scale_model_to_data, set_scale_model_to_data);
