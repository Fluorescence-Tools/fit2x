%{
#include "DecayScale.h"
%}

%attribute(DecayScale, double, number_of_photons, get_number_of_photons);
%attribute(DecayScale, double, constant_background, get_constant_background, set_constant_background);

%include "DecayScale.h"
